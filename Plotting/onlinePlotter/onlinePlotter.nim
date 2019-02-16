import plotly, fsmonitor2, nimhdf5
import ingrid / [raw_data_manipulation, reconstruction, tos_helpers, ingrid_types]
import helpers/utils
import zero_functional
import sequtils, algorithm, strutils, strformat, os, tables, heapqueue, options
import strscans
import asyncdispatch
import threadpool_simple
import posix

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
  OnlinePlotter = object
    files: HeapQueue[InGridFile]  ## seq containing filenames
    processed: int                ## number of files processed
    chip: int                     ## chip to show
    runNumber: int

import docopt

const docStr = """
Simple online event display + # pix histogram viewer

Usage:
  debugDiff <runFolder> [--out <h5file>] [--chip <chipNum>][options]

Options:
  --out <h5file>         Optional name for H5 file in which to store
                         reconstruced data.
  --chip <chipNum>       Chip to look at. Default is 3.
  -h, --help             Show this help
  -v, --version          Show the version number

Hand the path to the run folder for which to do online monitoring.
"""
const doc = withDocopt(docStr)

proc `<`(a, b: InGridFile): bool = a.evNumber < b.evNumber

proc processFiles(h5f: var H5FileObj,
                  files: seq[string],
                  runNumber: int,
                  rfKind: RunFolderKind,
                  attrsWritten: bool) =
  ## receives a slice of `OnlinePlotters` files and processes them
  const batchsize = 50000
  let r = readAndProcessInGrid(files, runNumber, rfKind)
  let nChips = r.nChips

  if attrsWritten == false:
    writeInGridAttrs(h5f, r, rfKind, rtNone)
    # create datasets in H5 file
    initInGridInH5(h5f, runNumber, nChips, batchsize)

  writeProcessedRunToH5(h5f, r)

proc getEvNum(fname: string): Option[int] =
  ## extracts event number from InGrid filename
  var evNumStr: string
  var dummy: string
  if scanf(fname, NewVirtexEventScanfNoPath, evNumStr):
    result = some(evNumStr.parseInt)
  elif scanf(fname, "$*" & NewVirtexEventScanfNoPath, dummy, evNumStr):
    result = some(evNumStr.parseInt)
  else:
    echo "Fname ", fname
    result = none[int]()

proc initialRead(h5file: string,
                 path: string,
                 runNumber: int,
                 rfKind: RunFolderKind): seq[InGridFile] =
  var h5f = H5file(h5file, "rw")
  # first get list of files; snapshot at this point in time
  var files = getSortedListOfFiles(path,
                                   EventSortType.fname,
                                   EventType.InGridType,
                                   rfKind)
  processFiles(h5f, files, runNumber, rfKind, false)

  # convert files to seq[InGridFile] and return
  result = files.mapIt(InGridFile(fname: it, evNumber: it.getEvNum.get))
  discard h5f.close()

proc initialReco(h5file: string, runNumber: int, chip: int): int =
  ## intial reconstruction of all events already written to disk
  var h5f = H5file(h5file, "rw")
  h5f.visit_file()
  var reco_run: seq[FlowVar[ref RecoEvent]]
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

proc writeNewEvents(h5f: var H5FileObj,
                    runNumber, chip: int,
                    reco: seq[RecoEvent]) =
  ## adds the needed information for the newly added events
  let grpBase = recoDataChipBase(runNumber) & $chip
  var
    x  = newSeqOfCap[seq[uint8]](reco.len)
    y  = newSeqOfCap[seq[uint8]](reco.len)
    ch = newSeqOfCap[seq[uint16]](reco.len)
    ev = newSeqOfCap[int](reco.len)
    xD = h5f[(grpBase / "x").dset_str]
    yD = h5f[(grpBase / "y").dset_str]
    chD = h5f[(grpBase / "ToT").dset_str]
    evD = h5f[(grpBase / "eventNumber").dset_str]
  for ev in reco:
    let
      num = ev.event_number
      chip = ev.chip_number
    for i, cl in ev.cluster:
      x.add(newSeq[uint8](cl.data.len))
      y.add(newSeq[uint8](cl.data.len))
      ch.add(newSeq[uint16](cl.data.len))
      ev.add(num)
      for j in 0..cl.data.high:
        x[^1][j]  = cl.data[j].x
        y[^1][j]  = cl.data[j].y
        ch[^1][j] = cl.data[j].ch
  let oldsize = xD.shape[0]
  let newsize = oldSize + reco.len
  xD.resize((newsize, 1))
  yD.resize((newsize, 1))
  chD.resize((newsize, 1))
  xD.write_hyperslab(x, offset = @[oldsize, 0], count = @[reco.len, 1])
  yD.write_hyperslab(y, offset = @[oldsize, 0], count = @[reco.len, 1])
  chD.write_hyperslab(ch, offset = @[oldsize, 0], count = @[reco.len, 1])
  evD.write_hyperslab(ev, offset = @[oldsize, 0], count = @[reco.len, 1])

#proc updatePlots(h5f: var H5FileObj, runNumber: int, rfKind: RunFolderKind) =
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
    processFiles(h5f, files, onPlt.runNumber, rfKind, true)

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
    h5f.writeNewEvents(onPlt.runNumber, onPlt.chip, reco)
    # finally set processed files to new value
    onPlt.processed = onPlt.len

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
    var
      monitor = newMonitor()
      onPlt = OnlinePlotter(chip: chip,
                            runNumber: runNumber)
    monitor.add(path)

    monitor.register(
      proc (ev: MonitorEvent) =
        case ev.kind
        of MonitorCreate:
          # file created, check valid event file
          let evNum = ev.fullName.getEvNum
          if evNum.isSome:
            onPlt.files.push InGridFile(fname: ev.fullName, evNumber: evNum.get)
        else: discard
    )
    monitor.watch()
    poll()
    onPlt.push(initialRead(h5file, path, runNumber, rfKind))
    onPlt.processed = initialReco(h5file, runNumber, onPlt.chip)

    var h5f = H5file(h5file, "rw")
    while true:
      poll()
      updateData(onPlt, h5f, rfKind)
      #updatePlots(onPlt, h5f)
    discard h5f.close()


when isMainModule:
  main()
