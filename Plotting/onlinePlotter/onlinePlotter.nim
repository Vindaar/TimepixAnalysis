import plotly, fsmonitor2, nimhdf5, shell
import ingrid / [raw_data_manipulation, reconstruction, tos_helpers, ingrid_types]
from .. / karaPlot / plotData import readEvent, read
import helpers/utils
import sequtils, algorithm, strutils, strformat, os, tables, heapqueue, options
import json, strscans, asyncdispatch, posix
import threadpool_simple
import zero_functional
import arraymancer
import ws
import asynchttpserver, asyncnet
import protocol

# - read all files already in folder
# - keep monitoring directory for more files to be added

# - store table (?) of files as a global state w/ bool whehter already in H5 file
# - async processing new events coming from monitoring.
# - async create new plot send via websocket to plotter
# - plotter very simple:
#   - read data packet w/ JSON from websocket, redraw that plot. loop
type
  InGridFile = object
    fname: string
    evNumber: int
  Event = object
    data: Tensor[float]
    idx: int
  Histo[T] = object
    # store as 1D sequence and recalculate the histogram each time
    # or store bins?
    # First go with recalculation and see how performance is
    data: seq[T]
    idx: Slice[int]

  OnlinePlotter = object
    files: HeapQueue[InGridFile]  ## seq containing filenames
    processed: int                ## number of files processed
    numReco: int                  ## size of reconstructed datasets
    processedPlot: int            ## number of files processed in plot
    chip: int                     ## chip to show
    runNumber: int
    event: Event
    pixHisto: Histo[int]
    plt: Grid
    hasNew: bool                  ## flag indicating new data is avilable

import docopt

const docStr = """
Simple online event display + # pix histogram viewer

Usage:
  onlinePlotter <runFolder> [--out <h5file>] [--chip <chipNum>][options]

Options:
  --out <h5file>         Optional name for H5 file in which to store
                         reconstruced data.
  --chip <chipNum>       Chip to look at. Default is 3.
  -h, --help             Show this help
  -v, --version          Show the version number

Hand the path to the run folder for which to do online monitoring.
"""
const doc = withDocopt(docStr)


var channel: Channel[string]
var eventChannel: Channel[InGridFile]
channel.open()
eventChannel.open()



proc `<`(a, b: InGridFile): bool = a.evNumber < b.evNumber

proc processFiles(h5f: var H5FileObj,
                  files: seq[string],
                  runNumber: int,
                  rfKind: RunFolderKind,
                  attrsWritten: bool): bool =
  ## receives a slice of `OnlinePlotters` files and processes them
  ## returns a bool indicating whether we actually wrote something
  ## to the H5 file (not a single broken file encountered)
  const batchsize = 50000
  let r = readAndProcessInGrid(files, runNumber, rfKind)
  if r.events.len > 0:
    # only continue if any events survived
    let nChips = r.nChips

    if attrsWritten == false:
      writeInGridAttrs(h5f, r, rfKind, rtNone)
      # create datasets in H5 file
      initInGridInH5(h5f, runNumber, nChips, batchsize)

    writeProcessedRunToH5(h5f, r)
    #
    result = true
  else:
    result = false

proc getEvNum(fname: string): Option[int] =
  ## extracts event number from InGrid filename
  var evNumStr: string
  var dummy: string
  let fnameNoPath = fname.extractFilename
  # Not needed since we extract the filename beforehand
  #if scanf(fnameNoPath, "$*/" & NewVirtexEventScanfNoPath, dummy, evNumStr):
  #  result = some(evNumStr.parseInt)
  if scanf(fnameNoPath, NewVirtexEventScanfNoPath, evNumStr):
    result = some(evNumStr.parseInt)
  else:
    echo "Fname ", fname
    result = none[int]()

proc initialRead(h5file: string,
                 path: string,
                 runNumber: int,
                 rfKind: RunFolderKind): seq[InGridFile] =
  var h5f = H5open(h5file, "rw")
  # first get list of files; snapshot at this point in time
  var files = getSortedListOfFiles(path,
                                   EventSortType.fname,
                                   EventType.InGridType,
                                   rfKind)
  let good = processFiles(h5f, files, runNumber, rfKind, false)
  # there should be no reason why the initial read should fail
  doAssert good == true, "Initial processing failed!"
  # convert files to seq[InGridFile] and return
  result = files.mapIt(InGridFile(fname: it, evNumber: it.getEvNum.get))
  discard h5f.close()

proc initialReco(h5file: string, runNumber: int, chip: int): int =
  ## intial reconstruction of all events already written to disk
  var h5f = H5open(h5file, "rw")
  h5f.visit_file()
  var reco_run: seq[FlowVar[ref RecoEvent[Pix]]]
  for ch, pixdata in h5f.readDataFromH5(runNumber):
    let reco = reconstructSingleChip(pixdata, runNumber, ch)
    reco_run.add reco
    # set inital number of processed events for online plotter
    if ch == chip:
      result = reco.len
  h5f.writeRecoRunToH5(h5f, reco_run, runNumber)
  discard h5f.close()

proc push[T: seq[InGridFile] | InGridFile](onPlt: var OnlinePlotter, file: T) =
  when T is InGridFile:
    onPlt.files.push file
  else:
    for f in file:
      onPlt.files.push f

proc `[]`(onPlt: OnlinePlotter, idx: Slice[int]): seq[string] =
  result = newSeqOfCap[string](idx.len)
  for i in idx.a .. idx.b:
    result.add onPlt.files[i].fname

proc len(onPlt: OnlinePlotter): int = result = onPlt.files.len
proc high(onPlt: OnlinePlotter): int = result = onPlt.files.len - 1

proc createPlots(onPlt: var OnlinePlotter, h5f: var H5FileObj) =
  ## creates the initial plots based on the already written data in the
  ## H5 file object
  let layout = Layout(title: "Septempower is over 9000!!!",
                      width: 1800,
                      height: 850)
  var grid = createGrid(numPlots = 2, layout = layout)
  # read initial data for pixel histogram
  let hits = h5f.read(onPlt.runNumber, "hits", onPlt.chip, dtype = int)
  # get the last event. Just use number of hits to know how many events
  let event = readEvent(h5f, onPlt.runNumber, onPlt.chip, @[hits.high]).squeeze

  onPlt.event = Event(data: event,
                      idx: hits.high)
  onPlt.pixHisto = Histo[int](data: hits,
                              idx: 0 .. hits.high)

  let histp = histPlot(onPlt.pixHisto.data)
    .binRange(0, 300)
    .binSize(2.0)
  let ev = heatmap(event.toRawSeq.reshape([256, 256]))
  grid[0] = histp
  grid[1] = ev
  onPlt.plt = grid
  onPlt.processedPlot = hits.high

proc updatePlots(onPlt: var OnlinePlotter, h5f: var H5FileObj) =
  if onPlt.processedPlot < onPlt.numReco:
    let idxToRead = onPlt.pixHisto.idx.b + 1 ..< onPlt.numReco
    let hits = h5f.read(onPlt.runNumber, "hits", onPlt.chip, dtype = int,
                        idx = toSeq(idxToRead))
    # get the last event. Just use number of hits to know how many events
    let event = readEvent(h5f, onPlt.runNumber, onPlt.chip, @[onPlt.numReco - 1]).squeeze

    onPlt.event = Event(data: event,
                        idx: onPlt.numReco - 1)
    let
      idxLow = onPlt.pixHisto.idx.a
      idxHigh = idxToRead.b
      combinedHits = concat(onPlt.pixHisto.data, hits)
    onPlt.pixHisto = Histo[int](data: combinedHits,
                                 idx: idxLow .. idxHigh)

    let histp = histPlot(onPlt.pixHisto.data)
      .binRange(0, 300)
      .binSize(2.0)
    let ev = heatmap(event.toRawSeq.reshape([256, 256]))
    onPlt.plt[0] = histp
    onPlt.plt[1] = ev

    onPlt.processedPlot = onPlt.numReco

proc writeNewEvents(onPlt: var OnlinePlotter,
                    h5f: var H5FileObj,
                    reco: seq[RecoEvent]) =
  ## adds the needed information for the newly added events
  let
    runNumber = onPlt.runNumber
    chip = onPlt.chip
    grpBase = recoDataChipBase(runNumber) & $chip
  var
    # reco len is a *lower* bound on the sequence size, since there may be
    # several clusters in a single frame!
    x  = newSeqOfCap[seq[uint8]](reco.len)
    y  = newSeqOfCap[seq[uint8]](reco.len)
    ch = newSeqOfCap[seq[uint16]](reco.len)
    ev = newSeqOfCap[int](reco.len)
    hits = newSeqOfCap[int](reco.len)
    xD = h5f[(grpBase / "x").dset_str]
    yD = h5f[(grpBase / "y").dset_str]
    chD = h5f[(grpBase / "ToT").dset_str]
    evD = h5f[(grpBase / "eventNumber").dset_str]
    hitsD = h5f[(grpBase / "hits").dset_str]
  for event in reco:
    let
      num = event.event_number
      chip = event.chip_number
    for i, cl in event.cluster:
      x.add(newSeq[uint8](cl.data.len))
      y.add(newSeq[uint8](cl.data.len))
      ch.add(newSeq[uint16](cl.data.len))
      ev.add num
      hits.add cl.data.len
      for j in 0..cl.data.high:
        x[^1][j]  = cl.data[j].x
        y[^1][j]  = cl.data[j].y
        ch[^1][j] = cl.data[j].ch
  let oldsize = xD.shape[0]
  # use event sequence length to get number of clusters
  let numClusters = ev.len
  let newsize = oldSize + numClusters
  echo "Old size ", oldsize, " newsize ", newsize, " last ev ", ev[^1]
  xD.resize((newsize, 1))
  yD.resize((newsize, 1))
  chD.resize((newsize, 1))
  evD.resize((newsize, 1))
  hitsD.resize((newsize, 1))
  xD.write_hyperslab(x, offset = @[oldsize, 0], count = @[numClusters, 1])
  yD.write_hyperslab(y, offset = @[oldsize, 0], count = @[numClusters, 1])
  chD.write_hyperslab(ch, offset = @[oldsize, 0], count = @[numClusters, 1])
  evD.write_hyperslab(ev, offset = @[oldsize, 0], count = @[numClusters, 1])
  hitsD.write_hyperslab(hits, offset = @[oldsize, 0], count = @[numClusters, 1])

  # set max index
  onPlt.numReco = newsize

proc updateData(onPlt: var OnlinePlotter, h5f: var H5FileObj,
                rfKind: RunFolderKind) =
  ## update the data based on new files in watched directory
  if onPlt.processed < onPlt.len:
    echo "New file was added apparently ", onPlt.processed, " and ", onPlt.len
    # new data was added
    let files = onPlt[onPlt.processed ..< onPlt.len]
    let numNew = onPlt.len - onPlt.processed
    echo "Files :", files, " is ", numNew
    # read files and reconstruct individual events
    let good = processFiles(h5f, files, onPlt.runNumber, rfKind, true)
    if good:
      # read new events into a `seq[(Pixels, int)]` type
      let grpName = rawDataChipBase(onPlt.runNumber) & $(onPlt.chip)
      let
        special_u8 = special_type(uint8)
        special_u16 = special_type(uint16)
      let
        xDset = h5f[(grpName / "raw_x").dset_str]
        yDset = h5f[(grpName / "raw_y").dset_str]
        chDset = h5f[(grpName / "raw_ch").dset_str]
        evNumDset =  h5f[(grpName.parentDir / "eventNumber").dset_str]
      let toRead = toSeq(onPlt.processed ..< onPlt.len)
      let
        xs = xDset[special_u8, uint8, toRead]
        ys = yDset[special_u8, uint8, toRead]
        chs = chDset[special_u16, uint16, toRead]
        evNum = evNumDset[toRead, int]
      var pixData: seq[(Pixels, int)]
      for i, ev in evNum:
        let
          xi = xs[i]
          yi = ys[i]
          chi = chs[i]
        let pix: Pixels = zip(xi, yi, chi) --> map((x: it[0], y: it[1], ch: it[2]))
        pixData.add (pix, ev)
      let reco = reconstructSingleChip(pixdata,
                                       onPlt.runNumber,
                                       onPlt.chip)
        .mapIt((^it)[])
      # still need to write back to H5 file
      writeNewEvents(onPlt, h5f, reco)
      # finally set processed files to new value
      onPlt.processed = onPlt.len
    else:
      echo "No events could be read to reconstruct!"
      # sleep for a short time to make sure we don't encounter a problem again,
      # if the reason we end up here is that a single event was read, and that
      # failed parsing (i.e. broken InGrid file error)
      sleep(50)

proc sendPacket[T](channel: var Channel[T], onPlt: OnlinePlotter) =
  ## send current plot via channel
  var data = ""
  toUgly(data, (% onPlt.plt))
  let dp = initDataPacket(kind = Messages.Data,
                          payload = data)
  channel.send(dp.asData.string)

proc processClient(req: Request) {.async.} =
  ## handle a single client
  let ws = await ws.newWebSocket(req)
  if ws.isNil:
    echo "WS negotiation failed: ", ws.readyState
    await req.respond(Http400, "Websocket negotiation failed: " & $ws.readyState)
    req.client.close()
    return
  else:
    # receive connection successful package
    let (opcodeConnect, dataConnect) = await ws.receivePacket()
    let dp = parseDataPacket(dataConnect)
    case dp.kind
    of Messages.Connected:
      echo "New websocket customer arrived! ", dataConnect
    else:
      echo "Received wrong packet, quit early"
      return

  var dp: DataPacket
  var
    hasData = false
    dpData = ""
  while true:
    let (opcode, data) = await ws.receivePacket()
    echo "(opcode: ", opcode, ", data length: ", data

    # reset has data bool
    hasData = false

    case opcode
    of Opcode.Text, Opcode.Cont:
      # parse the `DataPacket` we just received
      dp = parseDataPacket(data)
      # now just send this DataPacket and await response from worker
      doAssert dp.kind == Messages.Request

      while not hasData:
        (hasData, dpData) = channel.tryRecv()
        sleep(50)
      # send data as packet
      waitFor ws.send(dpData)
    of Opcode.Close:
      ws.close()
      echo "Socket went away, close code: ", ws.readyState
      req.client.close()
      return
    else:
      echo "Unkown error: ", opcode

  echo "This client dies now!"
  ws.close()
  req.client.close()

proc serveClient(server: AsyncHttpServer) {.async.} =
  var
    stopAvailable = false
    stop = false
  var clientFut = server.serve(Port(8080), processClient)
  while not stopAvailable:
    if not clientFut.finished:
      # client still connected, continue
      poll(500)
    else:
      # this client disconnected early, so accept another one
      clientFut = server.serve(Port(8080), processClient)

proc serve() =
  var server = newAsyncHttpServer()
  shell:
    "../../Tools/karaRun/karaRun" "-d:release -r client.nim"
  waitFor server.serveClient()

proc watcher(path: string) {.gcsafe.} =
  var
    monitor = newMonitor()
  monitor.add(path, {MonitorCreate})
  #monitor.register(
  proc handle(ev: MonitorEvent) =
    case ev.kind
    of MonitorCreate:
      # file created, check valid event file
      let evNum = ev.name.getEvNum

      if evNum.isSome:
        eventChannel.send InGridFile(fname: path / ev.name, evNumber: evNum.get)
        echo "New data available!"
    else: discard

  while true:
    let fut = monitor.read()
    while not fut.finished():
      sleep(50)
    #for cb in monitor.handleEvents:
    for action in fut.read():
      handle(action)
  #monitor.watch()
  #runForever()

proc worker(h5file, path: string,
            runNumber, chip: int,
            rfKind: RunFolderKind) =
  var onPlt = OnlinePlotter(chip: chip,
                            runNumber: runNumber)
  # perform initial read and process
  onPlt.push(initialRead(h5file, path, runNumber, rfKind))
  sleep(1000)
  onPlt.processed = initialReco(h5file, runNumber, onPlt.chip)

  # open file and send first plot
  var h5f = H5open(h5file, "rw")
  createPlots(onPlt, h5f)
  channel.sendPacket(onPlt)
  # now watch channel for new data
  var
    newFiles: seq[string]
    hasData = true
    data: InGridFile
  while true:
    # reset has data
    hasData = true
    sleep(50)
    while hasData:
      (hasData, data) = eventChannel.tryRecv()
      if hasData:
        onPlt.push data
        onPlt.hasNew = true
    if onPlt.hasNew:
      onPlt.hasNew = false
      updateData(onPlt, h5f, rfKind)
      updatePlots(onPlt, h5f)
      # take current plot and send via channel
      # TODO: check if new data, only send then
      channel.sendPacket(onPlt)
    else:
      echo "No new data available!"

  discard h5f.close()

proc main =
  let args = docopt(doc)
  echo args
  let path = $args["<runFolder>"]
  var h5file = $args["--out"]

  var chip = 3
  if $args["--chip"] != "nil":
    chip = ($args["--chip"]).parseInt

  let (isRunFolder, runNumber, rfKind, contains_run_folder) = isTosRunFolder(path)
  if h5file == "nil":
    h5file = "run_" & $runNumber & ".h5"

  if isRunFolder:
    # start the server thread
    var thr: Thread[void]
    thr.createThread(serve)

    var watcherThr: Thread[string]
    watcherThr.createThread(watcher, path)

    # now perform main work of reading, reconstruction and storing data
    worker(h5file, path, runNumber, chip, rfKind)

    # TODO: properly clean up on Ctrl+C
    joinThread(thr)
    joinThread(watcherThr)



when isMainModule:
  main()
