## this script performs the raw data manipulation of a given run or list of runs
## and outputs the resulting data to a HDF5 file
## steps which are performed
## - reading all data*.txt and data*.txt-fadc files and writing them to a
##   HDF5 file, one group for each run
## - calculating the occupancy of each run
## - calculating the num_pix / event histogram
## - caluclating the FADC signal depth / event histogram
## - calculate the ToT per pixel histogram
## standard lib
import std / [os, osproc, logging, sequtils, sugar, algorithm, tables, times,
              strutils, strformat, rdstdin]

# InGrid-module
import fadc_helpers
import helpers/utils
import tos_helpers
import ingrid_types

# other modules
import seqmath
import nimhdf5
import arraymancer
import parsetoml
import docopt

type
  RawFlagKind = enum
    rfIgnoreRunList, rfOverwrite, rfNoFadc, rfTpx3

const FILE_BUFSIZE = 75_000

##############################
# create globals for 2014/15 run list
##############################
# projectDefs contains `OldTosRunListPath` among others
import projectDefs
## TODO: this check is broken! Oh no, it's not broken, but if the binary is moved after compilation
## the CT check is invalid!
var
  oldTosRunListFound = false
  oldTosCalibRuns: set[uint16] = {}
  oldTosBackRuns: set[uint16]  = {}
  oldTosXrayRuns: set[uint16]  = {}
if fileExists(OldTosRunListPath):
  oldTosCalibRuns = parseOldTosRunlist(OldTosRunListPath, rtCalibration)
  oldTosBackRuns  = parseOldTosRunlist(OldTosRunListPath, rtBackground)
  oldTosXrayRuns  = parseOldTosRunlist(OldTosRunListPath, rtXrayFinger)
  oldTosRunListFound = true

const docStr = """
InGrid raw data manipulation.

Usage:
  raw_data_manipulation <folder> [options]
  raw_data_manipulation <folder> --runType <type> [options]
  raw_data_manipulation <folder> --out=<name> [--nofadc] [--runType=<type>] [--ignoreRunList] [options]
  raw_data_manipulation <folder> --nofadc [options]
  raw_data_manipulation --tpx3 <H5File> [options]
  raw_data_manipulation --tpx3 <H5File> --runType <type> [options]
  raw_data_manipulation --tpx3 <H5File> --runType <type> --out=<name> [options]
  raw_data_manipulation -h | --help
  raw_data_manipulation --version

Options:
  --tpx3 <H5File>     Convert data from a Timepix3 H5 file to TPA format
  --runType=<type>    Select run type (Calib | Back | Xray)
                      The following are parsed case insensetive:
                      Calib = {"calib", "calibration", "c"}
                      Back = {"back", "background", "b"}
                      Xray = {"xray", "xrayfinger", "x"}
  --out=<name>        Filename of output file
  --nofadc            Do not read FADC files
  --ignoreRunList     If set ignores the run list 2014/15 to indicate
                      using any rfOldTos run
  --overwrite         If set will overwrite runs already existing in the
                      file. By default runs found in the file will be skipped.
                      HOWEVER: overwriting is assumed, if you only hand a
                      run folder!
  -h --help           Show this help
  --version           Show version.

"""
const doc = withDocopt(docStr)

# define the compression filter we use
#let filter = H5Filter(kind: fkZlib, zlibLevel: 6) # on run 146 took: 0.889728331565857 min
let filter = H5Filter(kind: fkZlib, zlibLevel: 4) # on run 146 took: 0.5987015843391419 min
#let filter = H5Filter(kind: fkZlib, zlibLevel: 9) # on run 146 took:  2.170824865500132 min
# let filter = H5Filter(kind: fkBlosc, bloscLevel: 4,
#                       doShuffle: false,
#                       compressor: BloscCompressor.LZ4) # on run 146 took: 0.7583944002787272 min
# let filter = H5Filter(kind: fkBlosc, bloscLevel: 9,
#                       doShuffle: false,
#                       compressor: BloscCompressor.LZ4) # on run 146 took: 0.6799025972684224 min
# let filter = H5Filter(kind: fkBlosc, bloscLevel: 4,
#                       doShuffle: true,
#                       compressor: BloscCompressor.LZ4) # on run 146 took: 0.6656568646430969 min
# let filter = H5Filter(kind: fkBlosc, bloscLevel: 9,
#                       doShuffle: true,
#                       compressor: BloscCompressor.LZ4) # on run 146 took: 0.4899360696474711 min
# let filter = H5Filter(kind: fkNone)

################################################################################

# set up the logger
var L = newConsoleLogger()
if not dirExists("logs"):
  createDir("logs")
var fL = newFileLogger("logs/raw_data_manipulation.log", fmtStr = verboseFmtStr)
when isMainModule:
  addHandler(L)
  addHandler(fL)

################################################################################

template ch_len(): int = 2560
template all_ch_len(): int = ch_len() * 4

proc parseTomlConfig(): TomlValueRef =
  const sourceDir = currentSourcePath().parentDir
  result = parseToml.parseFile(sourceDir / "config.toml")

proc specialTypesAndEvKeys(): (hid_t, hid_t, array[7, string]) =
  let
    # create datatypes for variable length data
    ev_type_xy = special_type(uint8)
    ev_type_ch = special_type(uint16)
    # dataset names corresponding to event header keys
    eventHeaderKeys = ["eventNumber", "useHvFadc", "fadcReadout", "timestamp",
                       "szint1ClockInt", "szint2ClockInt", "fadcTriggerClock"]
  result[0] = ev_type_xy
  result[1] = ev_type_ch
  result[2] = eventHeaderKeys

proc getTotHitOccDsets(h5f: var H5FileObj, chipGroups: seq[H5Group]):
                      (seq[H5DataSet], seq[H5DataSet], seq[H5DataSet]) =
  proc fromChipGroups(h5f: var H5FileObj, chipGroups: seq[H5Group], name: string): seq[H5DataSet] =
    chipGroups.mapIt(h5f[(it.name & name).dset_str])
  var
    totDset = h5f.fromChipGroups(chipGroups, "/ToT")
    hitDset = h5f.fromChipGroups(chipGroups, "/Hits")
    occDset = h5f.fromChipGroups(chipGroups, "/Occupancy")
  result = (totDset, hitDset, occDset)

template batchFiles(files: var seq[string], bufsize, actions: untyped): untyped =
  ## this is a template, which can be used to batch a set of files into chunks
  ## of size `bufsize`. This is done by deleting elements from `files` until it
  ## is empty. A block of code `actions` is given to the template, which will be
  ## performed.
  ## The variable ind_high is injected into the calling space, to allow the
  ## `actions` block to slice the current files
  ##
  ## Example:
  ##
  ## .. code-block:
  ##   let
  ##     fname_base = "data$#.txt"
  ##     files = mapIt(toSeq(0..<1000), fname_base & $it)
  ##   batch_files(files, 100):
  ##     # make use of batching for calc with high memory consumption, making use of
  ##     # injected `ind_high` variable
  ##     echo memoryConsumptuousCalc(files[0..ind_high])
  while len(files) > 0:
    # variable to set last index to read to
    var ind_high {.inject.} = min(files.high, bufsize - 1)
    # perform actions as desired
    actions

    info "... removing read elements from list"
    # sequtils.delete removes the element with ind_high as well!
    files.delete(0, ind_high)

proc batchFileReading[T](files: var seq[string],
                         rfKind: RunFolderKind = rfNewTos,
                         bufsize: int = FILE_BUFSIZE):
                          seq[T] {.inline.} =
  # and removes all elements in the file list until all events have been read and the seq
  # is empty
  let
    t0 = epochTime()
    n_files = len(files)
  var count = 0
  # initialize sequence
  result = @[]

  var buf_seq: type(result)
  batchFiles(files, bufsize):
    # read files into buffer sequence, `ind_high` is an injected variable of the template
    # after each iteration the `files` variable is modified. Read files are deleted.
    when T is Event:
      case rfKind
      of rfOldTos, rfNewTos, rfSrsTos:
        buf_seq = readListOfInGridFiles(files[0 .. ind_high], rfKind)
      else:
        raise newException(IOError, "Unknown run folder kind. Cannot read " &
          "event files!")
    elif T is FadcFile:
      buf_seq = readListOfFadcFiles(files[0 .. ind_high])

    info "... and concating buffered sequence to result"
    result = concat(result, buf_seq)
    count += bufsize
  info "All files read. Number = " & $len(result)
  info "Reading took $# seconds" % $(epochTime() - t0)
  info "Compared with starting files " & $n_files


proc readRawInGridData*(listOfFiles: seq[string],
                        rfKind: RunFolderKind):
                          seq[Event] =
  ## given a run_folder it reads all event files (data<number>.txt) and returns
  ## a sequence of Events, which store the raw event data
  ## Intermediately we receive FlowVars to ref Events after reading. We read via
  ## inodes, which may be scramled, so we sort the data and get the FlowVar values.
  ## NOTE: this procedure does the reading of the data in parallel, thanks to
  ## using spawn
  # get a sorted list of files, sorted by filename first
  var files: seq[string] = sortByInode(listOfFiles)
  # split the sorted files into batches, and sort each batch by inode
  result = batchFileReading[Event](files, rfKind)

proc sortReadInGridData(rawInGrid: seq[Event],
                        rfKind: RunFolderKind): seq[Event] =
  ## sorts the seq of FlowVars according to event numbers again, otherwise
  ## h5 file is all mangled
  info "Sorting data..."
  # case on the old TOS' data storage and new TOS' version
  let t0 = epochTime()
  case rfKind
  of rfNewTos, rfOldTos, rfSrsTos:
    # in this case there may be missing events, so we simply sort by the indices themselves
    # sorting is done the following way:
    # - extract list of eventNumbers from `Events`
    # - create list of tuples of (`EventNumber`, `AtIndex`)
    # - sort tuples by `EventNumber`
    # - add elements to result taken from `AtIndex` in `Events`
    let
      numEvents = raw_ingrid.len
      # get event numbers
      numList = mapIt(raw_ingrid, it.evHeader["eventNumber"].parseInt)
      # zip event numbers with indices in raw_ingrid (unsorted!)
      zipped = zip(numList, toSeq(0 ..< numEvents))
      # sort tuples by event numbers (indices thus mangled, but in "correct" order
      # for insertion)
      sortedNums = zipped.sortedByIt(it[0])
    info &"Min event number {min(numList)} and max number {max(numList)}"
    # insert elements into result
    result = newSeqOfCap[Event](raw_ingrid.len)
    for i in sortedNums:
      result.add raw_ingrid[i[1]]
  else:
    # we'll never end up here with rfUnknown, unless something bad happens
    logging.fatal("Unkown error. Ended up with unknown run folder kind " &
      "in `sortReadInGridData`. Stopping program")
    quit(1)

  let t1 = epochTime()
  info &"...Sorting done, took {$(t1 - t0)} seconds"


proc processRawInGridData(run: Run): ProcessedRun =
  ## procedure to process the raw data read from the event files by readRawInGridData
  ## inputs:
  ##    ch: seq[Event]] = seq of Event objects, which each store raw data of a single event.
  ##        We read normal events, perform calculations
  ##        to obtain ToT per pixel, number of hits and occupancies of that data
  ##    runNumber: int = the run number of the current run
  ##    runHeader = The header valid for the whole run.
  ## outputs:
  ##   ProcessedRun containing:
  ##    events:    seq[Event] = the raw data from the seq of FlowVars saved in seq of Events
  ##    tuple of:
  ##      tot:  seq[seq[int]] = ToT values for each chip of Septemboard for run
  ##      hits: seq[seq[int]] = number of hits for each chip of Septemboard fo run
  ##      occ:    Tensor[int] = (nChips, 256, 256) tensor containing occupancies of all chips for
  ##        this data.

  let
    ch = run.events
    runHeader = run.runHeader
  # get number of chips from header
  let nChips = parseInt(runHeader["numChips"])

  # variable to count number of processed files
  var
    count = 0
    # store ToT data of all events for each chip
    # Note: still contains a single seq for each event! Need to concat
    # these seqs at the end
    tot_run: seq[seq[seq[uint16]]] = newSeq[seq[seq[uint16]]](nChips)
    # store occupancy frames for each chip
    # TODO: allow for other values than 7 chips!
    occ = zeros[int64](nChips, 256, 256)
    # store number of hits for each chip
    hits = newSeq[seq[uint16]](nChips)
    # initialize the events sequence of result, since we add to this sequence
    # instead of copying ch to it!
    events = newSeq[Event](len(ch))
  let
    # get the run specific time and shutter mode
    time = parseFloat(runHeader["shutterTime"])
    mode = float(parseShutterMode(runHeader["shutterMode"]))

  # set the run number
  result.runNumber = run.runNumber

  # initialize empty sequences. Input to anonymous function is var
  # as we change each inner sequence in place with newSeq
  apply(tot_run, (x: var seq[seq[uint16]]) => newSeq[seq[uint16]](x, len(ch)))
  apply(hits, (x: var seq[uint16]) => newSeq[uint16](x, len(ch)))

  info "starting to process events..."
  count = 0
  for i in 0 .. ch.high:
  # assign the length field of the ref object
    # for the rest, get a copy of the event
    var a: Event = ch[i]
    # TODO: think about parallelizing here by having proc, which
    # works processes the single event?
    a.length = calcLength(a, time, mode)
    events[i] = a
    let chips = a.chips
    for c in chips:
      let
        num = c.chip.number
        pixels = c.pixels
      addPixelsToOccupancySeptem(occ, pixels, num)
      let tot_event = pixelsToTOT(pixels)
      tot_run[num][i] = tot_event
      let n_pix = len(tot_event).uint16
      if n_pix > 0'u16:
        # if the compiler flag (-d:CUT_ON_CENTER) is set, we cut all events, which are
        # in the center 4.5mm^2 square of the chip
        when defined(CUT_ON_CENTER):
          if isNearCenterOfChip(pixels) == true:
            hits[num][i] = n_pix
        else:
          hits[num][i] = n_pix
    echoCount(count, msg = " files processed.")

  # use first event of run to fill event header. Fine, because event
  # header is contained in every file
  result.runHeader = runHeader
  result.chips = run.chips
  result.nChips = nChips
  result.events = events
  result.tots = newSeq[seq[uint16]](nChips)
  for i, tot in tot_run:
    result.tots[i] = concat tot
  result.hits = hits
  result.occupancies = occ

proc processFadcData(fadcFiles: seq[FadcFile]): ProcessedFadcData {.inline.} =
  ## proc which performs all processing needed to be done on the raw FADC
  ## data. Starting from conversion of FadcFiles -> FadcData, but includes
  ## calculation of minimum and check for noisy events
  # sequence to store the indices needed to extract the 0 channel
  let
    #fadc_ch0_indices = getCh0Indices()
    #ch_len = ch_len()
    #pedestal_run = getPedestalRun()
    nEvents = fadcFiles.len
    # we demand at least 4 dips, before we can consider an event as noisy
    #n_dips = 4
    # the percentile considered for the calculation of the minimum
    #min_percentile = 0.95

  # convert FlowVars of FadcFiles to sequence of FadcFiles
  result.raw_fadc_data = newSeq[seq[uint16]](nEvents)
  #result.fadc_data = zeros[float]([nEvents, ch_len])
  result.trigRecs = newSeq[int](nEvents)
  #result.noisy = newSeq[int](nEvents)
  #result.minVals = newSeq[float](nEvents)
  result.eventNumber = newSeq[int](nEvents)
  let t0 = epochTime()
  # TODO: parallelize this somehow so that it's faster!
  for i, ev in fadcFiles:
    result.raw_fadc_data[i] = ev.data
    result.trigRecs[i]      = ev.trigRec
    result.eventNumber[i]   = ev.eventNumber
    #let fadc_dat = ev.fadcFileToFadcData(pedestal_run, fadc_ch0_indices).data
    #result.fadc_data[i, _]  = fadc_dat.reshape([1, ch_len])
    #result.noisy[i]         = fadc_dat.isFadcFileNoisy(n_dips)
    #result.minVals[i]       = fadc_dat.calcMinOfPulse(min_percentile)

  # this parallel solution seems to be slower, instead of faster ?! well, probably
  # because we only have these two spawns and one of these functions is much slower
  # than the other
  # parallel:
  #   for i, event in fadcFiles:
  #     let ev = (^event)[]
  #     if i < result.raw_fadc_data.len:
  #       result.raw_fadc_data[i] = ev.data
  #     let fadc_dat = ev.fadcFileToFadcData(pedestal_run, fadc_ch0_indices).data
  #     if i < result.trigRecs.len:
  #       result.trigRecs[i]      = ev.trigrec
  #       result.fadc_data[i, _]  = fadc_dat.reshape([1, ch_len])
  #     if i < result.noisy.len:
  #       result.noisy[i]         = spawn fadc_dat.isFadcFileNoisy(n_dips)
  #     if i < result.minVals.len:
  #       result.minVals[i]       = spawn fadc_dat.calcMinOfPulse(min_percentile)
  #     sync()
  info "Calculation of $# events took $# seconds" % [$nEvents, $(epochTime() - t0)]

proc initFadcInH5(h5f: var H5FileObj, runNumber, batchsize: int, filename: string) =
  # proc to initialize the datasets etc in the HDF5 file for the FADC. Useful
  # since we don't want to do this every time we call the write function
  const
    ch_len = ch_len()
    all_ch_len = all_ch_len()

  let groupName = fadcRawPath(runNumber)
  template datasetCreation(h5f, name, shape, `type`: untyped): untyped =
    ## inserts the correct data set creation parameters
    h5f.create_dataset(name,
                       shape,
                       dtype = `type`,
                       chunksize = @[batchsize, shape[1]],
                       maxshape = @[int.high, shape[1]],
                       filter = filter)
  var
    runGroup = h5f.create_group(groupName)
    # create the datasets for raw data etc
    # NOTE: we initialize all datasets with a size of 0. This means we need to extend
    # it immediately. However, this allows us to always (!) simply extend and write
    # the data to dset.len onwards!
    raw_fadc_dset    = h5f.datasetCreation(rawFadcBasename(runNumber), (0, all_ch_len), uint16)
    #fadc_dset        = h5f.datasetCreation(fadcDataBasename(runNumber), (0, ch_len), float)
    trigrec_dset     = h5f.datasetCreation(trigrecBasename(runNumber), (0, 1), int)
    # dataset of eventNumber
    eventNumber_dset = h5f.datasetCreation(eventNumberBasenameRaw(runNumber), (0, 1), int)
    # dataset stores flag whether FADC event was a noisyo one (using our algorithm)
    #noisy_dset       = h5f.datasetCreation(noiseBasename(runNumber), (0, 1), int)
    # dataset stores minima of each FADC event, dip voltage
    #minVals_dset     = h5f.datasetCreation(minValsBasename(runNumber), (0, 1), float)

  # write attributes to FADC groups
  # read the given FADC file and extract that information from it
  let fadc_for_attrs = readFadcFile(filename)
  # helper sequence to loop over both groups to write attrs
  runGroup.attrs["posttrig"] = fadc_for_attrs.posttrig
  runGroup.attrs["pretrig"] = fadc_for_attrs.pretrig
  runGroup.attrs["n_channels"] = fadc_for_attrs.n_channels
  runGroup.attrs["channel_mask"] = fadc_for_attrs.channel_mask
  runGroup.attrs["frequency"] = fadc_for_attrs.frequency
  runGroup.attrs["sampling_mode"] = fadc_for_attrs.sampling_mode
  runGroup.attrs["pedestal_run"] = if fadc_for_attrs.pedestal_run == true: 1 else: 0

proc writeFadcDataToH5(h5f: var H5FileObj, runNumber: int, f_proc: ProcessedFadcData) =
  # proc to write the current FADC data to the H5 file
  # now write the data
  let
    raw_name = rawFadcBasename(runNumber)
    trigRec_name = trigRecBasename(runNumber)
    eventNumber_name = eventNumberBasenameRaw(runNumber)
    ch_len = ch_len()
    all_ch_len = all_ch_len()
    nEvents = f_proc.raw_fadc_data.len
  var
    raw_fadc_dset = h5f[raw_name.dset_str]
    trigRec_dset = h5f[trigRec_name.dset_str]
    eventNumber_dset = h5f[eventNumber_name.dset_str]

  info raw_fadc_dset.shape
  info raw_fadc_dset.maxshape
  info trigRec_dset.shape
  info trigRec_dset.maxshape
  info eventNumber_dset.shape
  info eventNumber_dset.maxshape
  # first need to extend the dataset, as we start with a size of 0.
  let oldsize = raw_fadc_dset.shape[0]
  let newsize = oldsize + nEvents
  # TODO: currently this is somewhat problematic. We simply resize always. In a way this is
  # fine, because we need to resize to append. But in case we start this program twice in
  # a row, without deleting the file, we simply extend the dataset further, because we read
  # the current (final!) shape from the file
  # NOTE: one way to mitigate t his, would be to set oldsize as a {.global.} variable
  # in which case we simply set it to 0 on the first call and afterwards extend it by
  # the size we add
  info "Adding to FADC datasets, from $# to $#" % [$oldsize, $newsize]
  raw_fadc_dset.add f_proc.raw_fadc_data
  trigRec_dset.add f_proc.trigRecs
  eventNumber_dset.add f_proc.eventNumber

  # now write the data
  let t0 = epochTime()
  # write using hyperslab
  #fadc_dset.write_hyperslab(f_proc.fadc_data.toRawSeq,
  #                          offset = @[oldsize, 0],
  #                          count = @[nEvents, ch_len])
  #noisy_dset.write_hyperslab(f_proc.noisy,
  #                           offset = @[oldsize, 0],
  #                           count = @[nEvents, 1])
  #minVals_dset.write_hyperslab(f_proc.minVals,
  #                             offset = @[oldsize, 0],
  #                             count = @[nEvents, 1])
  info "Writing of FADC data took $# seconds" % $(epochTime() - t0)

proc readWriteFadcData(run_folder: string, runNumber: int, h5f: var H5FileObj) =
  ## given a run_folder it reads all fadc files (data<number>.txt-fadc),
  ## processes it (FadcFile -> FadcRaw) and writes it to the HDF5 file
  # get a sorted list of files, sorted by inode
  var
    files: seq[string] = getSortedListOfFiles(run_folder,
                                              EventSortType.fname,
                                              EventType.FadcType,
                                              RunFolderKind.rfUnknown)
    raw_fadc_data: seq[FadcFile]
    # variable to store the processed FADC data
    f_proc: ProcessedFadcData
    # in case of FADC data we cannot afford to read all files into memory before
    # writing some to HDF5, because the memory overhead from storing all files
    # in seq[string] is too large (17000 FADC files -> 10GB needed!)
    # thus already perform batching here
    files_read: seq[string]
  if files.len == 0:
    # in case there are no FADC files, return from this proc
    return
  const batchsize = 2500
  # before we start iterating over the files, initialize the H5 file
  h5f.initFadcInH5(runNumber, batchsize, files[0])
  batchFiles(files, batchsize):
    # batch in 1000 file pieces
    var mfiles = files[0 .. ind_high]
    info "Starting with file $# and ending with file $#" % [$mfiles[0], $mfiles[^1]]
    files_read = files_read.concat(mfiles)
    raw_fadc_data = batchFileReading[FadcFile](mfiles)

    # TODO: read FADC files also by inode and then sort the fadc
    # we just read here. NOTE: for that have to change the writeFadcDataToH5
    # proc to accomodate that!

    # given read files, we now need to append this data to the HDF5 file, before
    # we can process more data, otherwise we might run out of RAM
    f_proc = raw_fadc_data.processFadcData
    info "Number of FADC files in this batch ", raw_fadc_data.len

    h5f.writeFadcDataToH5(runNumber, f_proc)

  info "Number of files read: ", files_read.toSet.len
  # finally finish writing to the HDF5 file
  # finishFadcWriteToH5(h5f, runNumber)

proc createChipGroups(h5f: var H5FileObj,
                      runNumber: int,
                      nChips: int = 0): seq[H5Group] =
  let chipGroupName = getGroupNameRaw(runNumber) & "/chip_$#"
  result = mapIt(toSeq(0 ..< nChips), h5f.create_group(chipGroupName % $it))

proc initInGridInH5*(h5f: var H5FileObj, runNumber, nChips, batchsize: int) =
  ## This proc creates the groups and dataset for the InGrid data in the H5 file
  ## inputs:
  ##   h5f: H5file = the H5 file object of the writeable HDF5 file
  ##   ?
  # create variables for group names (NOTE: dirty template!)
  let groupName = getGroupNameRaw(runNumber)
  let chipGroups = createChipGroups(h5f, runNumber, nChips)
  let (ev_type_xy, ev_type_ch, eventHeaderKeys) = specialTypesAndEvKeys()

  template datasetCreation(h5f: untyped, name: untyped, `type`: untyped): untyped =
    ## inserts the correct data set creation parameters
    h5f.create_dataset(name,
                       (0, 1),
                       dtype = `type`,
                       chunksize = @[batchsize, 1],
                       maxshape = @[int.high, 1],
                       filter = filter)

  var
    # datasets are chunked in the batchsize we read. Size originally 0
    x_dsets  = mapIt(chip_groups, h5f.datasetCreation(it.name & "/raw_x", ev_type_xy))
    y_dsets  = mapIt(chip_groups, h5f.datasetCreation(it.name & "/raw_y", ev_type_xy))
    ch_dsets = mapIt(chip_groups, h5f.datasetCreation(it.name & "/raw_ch", ev_type_ch))

    # datasets to store the header information for each event
    evHeadersDsetTab = eventHeaderKeys.mapIt(
      (it,
       h5f.datasetCreation(groupName & "/" & it, int))
    ).toTable
    # TODO: add string of datetime as well
    #dateTimeDset = h5f.create_dataset(joinPath(group_name, "dateTime"), nEvents, string)

    # other single column data
    durationDset = h5f.datasetCreation(joinPath(groupName, "eventDuration"), float)
  let names = chipGroups.mapIt(it.name)
  var
    totDset = mapIt(names, h5f.datasetCreation(it & "/ToT", uint16))
    hitDset = mapIt(names, h5f.datasetCreation(it & "/Hits", uint16))
    # use normal dataset creation proc, due to static size of occupancies
    occDset = mapIt(names, h5f.create_dataset(it & "/Occupancy", (256, 256), int))

proc getCenterChipAndName(run: ProcessedRun): (int, string) =
  ## returns the chip number and the name of the center chip
  # TODO: Find nicer solution!
  var centerChip = 0
  case run.nChips
  of 1:
    centerChip = 0
  of 7:
    centerChip = 3
  else:
    warn &"This number of chips ({run.nChips}) currently unsupported for" &
      " `centerChip` determination. Will be set to 0."
  let centerName = run.chips[centerChip].name
  result = (centerChip, centerName)

proc writeRawAttrs*(h5f: var H5FileObj,
                    run: ProcessedRun,
                    rfKind: RunFolderKind,
                    runType: RunTypeKind) =
  # finally write run type to base runs group
  var rawG = h5f["runs".grp_str]
  rawG.attrs["runType"] = $runType
  rawG.attrs["runFolderKind"] = $rfKind
  let (centerChip, centerName) = getCenterChipAndName(run)
  rawG.attrs["centerChip"] = centerChip
  rawG.attrs["centerChipName"] = centerName
  let start = run.events[0].evHeader["dateTime"]
  let stop = run.events[^1].evHeader["dateTime"]
  rawG.attrs["runStart"] = start
  rawG.attrs["runStop"] = stop
  # NOTE: the run length will be wrong by the duration of the last event!
  rawG.attrs["totalRunDuration"] = (parseTOSDateString(stop) -
                                    parseTOSDateString(start)).inSeconds

proc writeRunGrpAttrs*(h5f: var H5FileObj, group: var H5Group,
                       runType: RunTypeKind,
                       run: ProcessedRun) =
  ## writes all attributes to given `group` that can be extracted from
  ## the `ProcessedRun`, `rfKind` and `runType`.
  # now write attribute data (containing the event run header, for a start
  # NOTE: unfortunately we cannot write all of it simply using applyIt,
  # because we need to parse some numbers as ints, leave some as strings
  let asInt = ["runNumber", "runTime", "runTimeFrames", "numChips", "shutterTime",
               "runMode", "fastClock", "externalTrigger"]
  let asString = ["pathName", "dateTime", "shutterMode"]

  # write run header
  # Note: need to check for existence of the keys, because for old TOS data,
  # not all keys are set!
  for it in asInt:
    if it in run.runHeader:
      let att = parseInt(run.runHeader[it])
      group.attrs[it]  = att
  for it in asString:
    if it in run.runHeader:
      let att = run.runHeader[it]
      group.attrs[it] = att
  let (centerChip, centerName) = getCenterChipAndName(run)
  group.attrs["centerChipName"] = centerName
  group.attrs["centerChip"] = centerChip
  # initialize the attribute for the current number of stored events to 0
  group.attrs["numEventsStored"] = 0
  group.attrs["runType"] = $runType

proc writeChipAttrs*(h5f: var H5FileObj,
                     chipGroups: var seq[H5Group],
                     run: ProcessedRun) =
  # write attributes for each chip
  for i, grp in mpairs(chip_groups):
    grp.attrs["chipNumber"] = run.chips[i].number
    grp.attrs["chipName"]   = run.chips[i].name
    # initialize the attribute for the current number of stored events to 0
    grp.attrs["numEventsStored"] = 0

proc writeInGridAttrs*(h5f: var H5FileObj, run: ProcessedRun,
                       rfKind: RunFolderKind, runType: RunTypeKind) =
  # writes all attributes into the output file. This includes
  # - "runs" group attributes
  # - individual run group attributes
  # - chip group attributes
  # "runs" group
  writeRawAttrs(h5f, run, rfKind, runType)
  # individual run group
  let groupName = getGroupNameRaw(run.runNumber)
  var group = h5f[groupName.grp_str]
  writeRunGrpAttrs(h5f, group, runType, run)
  # chip groups
  var chipGroups = createChipGroups(h5f, run.runNumber, run.nChips)
  writeChipAttrs(h5f, chipGroups, run)

proc fillDataForH5(x, y: var seq[seq[seq[uint8]]],
                   ch: var seq[seq[seq[uint16]]],
                   evHeaders: var Table[string, seq[int]],
                   duration: var seq[float],
                   events: seq[Event],
                   startEvent: int) =
  let
    nEvents = events.len
    # take 0 event to get number of chips, since same for whole run
    nChips = events[0].nChips
  for i in 0 ..< nChips:
    x[i]  = newSeq[seq[uint8]](nEvents)
    y[i]  = newSeq[seq[uint8]](nEvents)
    ch[i] = newSeq[seq[uint16]](nEvents)
  for i, event in events:
    duration[i] = event.length
    # add event header information
    for key in keys(evHeaders):
      try:
        evHeaders[key][i] = parseInt(event.evHeader[key])
      except KeyError:
        #echo "Event $# with evHeaders does not contain key $#" % [$event, $key]
        discard
      except IndexError:
        logging.error "Event $# with evHeaders does not contain key $#" % [$event, $key]
    # add raw chip pixel information
    for chp in event.chips:
      let
        num = chp.chip.number
      let hits = chp.pixels.len
      var
        xl: seq[uint8] = newSeq[uint8](hits)
        yl: seq[uint8] = newSeq[uint8](hits)
        chl: seq[uint16] = newSeq[uint16](hits)
      for i, p in chp.pixels:
        xl[i] = uint8(p[0])
        yl[i] = uint8(p[1])
        chl[i] = uint16(p[2])
      x[num][i] = xl
      y[num][i] = yl
      ch[num][i] = chl

proc writeProcessedRunToH5*(h5f: var H5FileObj, run: ProcessedRun) =
  ## this procedure writes the data from the processed run to a HDF5
  ## (opened already) given by h5file_id
  ## inputs:
  ##   h5f: H5file = the H5 file object of the writeable HDF5 file
  ##   run: ProcessedRun = a tuple of the processed run data
  let
    nEvents = run.events.len
    runNumber = run.runNumber
    nChips = run.nChips

  # TODO: write the run information into the meta data of the group
  info "Create data to write to HDF5 file"
  let t0 = epochTime()
  # first write the raw data
  # get the names of the groups
  let groupName = getGroupNameRaw(runNumber)
  var runGroup = h5f[groupName.grp_str]
  var chipGroups = createChipGroups(h5f, runNumber, nChips)
  let (ev_type_xy, ev_type_ch, eventHeaderKeys) = specialTypesAndEvKeys()

  var
    # get the datasets from the file in `chipGroups`
    x_dsets  = mapIt(chipGroups, h5f[(it.name & "/raw_x" ).dset_str])
    y_dsets  = mapIt(chipGroups, h5f[(it.name & "/raw_y" ).dset_str])
    ch_dsets = mapIt(chipGroups, h5f[(it.name & "/raw_ch").dset_str])

    # datasets to store the header information for each event
    evHeadersDsetTab = eventHeaderKeys.mapIt(
        (it, h5f[(groupName & "/" & it).dset_str])
      ).toTable
    # TODO: add string of datetime as well
    #dateTimeDset = h5f.create_dataset(joinPath(groupName, "dateTime"), nEvents, string)

    # other single column data
    durationDset = h5f[(joinPath(groupName, "eventDuration")).dset_str]
    duration = newSeq[float](nEvents)

    evHeaders = initTable[string, seq[int]]()
    x  = newSeq[seq[seq[uint8]]](nChips)
    y  = newSeq[seq[seq[uint8]]](nChips)
    ch = newSeq[seq[seq[uint16]]](nChips)

  # prepare event header keys and value (init seqs)
  for key in eventHeaderKeys:
    evHeaders[key] = newSeq[int](nEvents)

  ##############################
  ##### Fill the data seqs #####
  ##############################

  # use
  let oldsize = runGroup.attrs["numEventsStored", int]
  let newsize = oldsize + nEvents
  # set new size as attribute
  runGroup.attrs["numEventsStored"] = newsize

  # call proc to write the data from the events to the seqs, tables
  fillDataForH5(x, y, ch, evHeaders, duration, run.events, oldsize)

  ##############################
  ###### Write the data ########
  ##############################

  info "Writing all dset x data "
  for i in 0 ..< nChips:
    withDebug:
      info "Writing dsets ", i, " size x ", x_dsets.len
    x_dsets[i].add x[i]
    y_dsets[i].add y[i]
    ch_dsets[i].add ch[i]
    withDebug:
      info "Shape of x ", x[i].len, " ", x[i].shape
      info "Shape of dset ", x_dsets[i].shape

  for key, dset in mpairs(evHeadersDsetTab):
    withDebug:
      info "Writing $# in $#" % [$key, $dset]
    dset.add evHeaders[key]

  # write other single column datasets
  durationDset.add duration
  info "took a total of $# seconds" % $(epochTime() - t0)

  ####################
  #  Reconstruction  #
  ####################

  # TODO: maybe this needs to be done in a pass after everything has been done?
  # at least for occupancy?

  info "ToTs shape is ", run.tots.shape
  info "hits shape is ", run.hits.shape
  # into the reco group name we now write the ToT and Hits information
  var (totDsets, hitDsets, occDsets) = getTotHitOccDsets(h5f, chipGroups)

  for chip in 0 ..< run.nChips:
    # since not every chip has hits on each event, we need to create one group
    # for each chip and store each chip's data in these
    # in this chip group store:
    # - occupancy
    # - ToTs
    # - Hits
    # ... more later
    let
      tot = run.tots[chip]
      hit = run.hits[chip]
      occ = run.occupancies[chip, _, _].squeeze.clone
    var
      totDset = totDsets[chip]
      hitDset = hitDsets[chip]
      occDset = occDsets[chip]
    totDset.add tot
    hitDset.add hit
    # before writing the occupancy dataset, we need to read the old, stack the current
    # occupancy on it and finally write the result
    # TODO: this does not seem to make sense to me. We're iterating over all chips in a whole run.
    # Why would there be data in the occupancy dataset for us to read?
    let stackOcc = occDset[int64].toTensor.reshape([256, 256]) .+ occ
    occDset.unsafeWrite(stackOcc.get_data_ptr, stackOcc.size)

#proc linkRawToReco(h5f: var H5FileObj, runNumber, nChips: int) =
#  ## perform linking from raw group to reco group
#  let (groupName,
#       recoGroupName,
#       chipGroupName,
#       combineGroupName) = inGridRawGroupNames(runNumber)
#  let (_, _, eventHeaderKeys) = specialTypesAndEvKeys()
#  let (totDsetNames,
#       hitDsetNames,
#       occDsetNames) = getTotHitOccDsetNames(chipGroupName, nChips)
#  let
#    durationDsetName = joinPath(groupName, "eventDuration")
#
#  # link over to reconstruction group
#  h5f.create_hardlink(durationDsetName, recoGroupName / extractFilename(durationDsetName))
#  # create hard links of header data to reco group
#  for key in eventHeaderKeys:
#    h5f.create_hardlink(joinPath(groupName, key), joinPath(recoGroupName, key))

proc createRun(runHeader: Table[string, string],
               runNumber: int,
               events: seq[Event],
               rfKind: RunFolderKind): Run =
  ## performs the conversion from a sequence `Event` plus meta
  ## information to a `Run`, which combines everything and deals
  ## with differences of data storage for SRS, old V6 and new V6
  result.runHeader = runHeader
  result.runNumber = runNumber
  result.events = events
  # now extract the chips correclty depending on data type
  case rfKind
  of rfOldTos, rfNewTos:
    # just get it from any event
    let ev = events[0]
    # extract each `Chip` from each `ChipEvent`
    result.chips = ev.chips.mapIt(it.chip)
  of rfSrsTos:
    # in this case need to extract information from the
    # run header
    if not result.runHeader.hasKey(SrsRunIncomplete) and
       not result.runHeader.hasKey(SrsNoChipId):
      let nChips = result.runHeader["numChips"].parseInt
      result.chips = newSeq[Chip](nChips)
      for i in 0 ..< nChips:
        let name = result.runHeader[&"chip_{i}"]
        result.chips[i] = (name, i)
    else:
      # in this case take bad information from chips too
      let ev = events[0]
      result.chips = ev.chips.mapIt(it.chip)
  of rfUnknown:
    raise newException(Exception, "Creation of a `Run` is impossible for an " &
      "unknown run folder kind at the moment!")

proc calcTimeOfEvent(evDuration: float,
                     eventNumber: int,
                     timestamp: string,
                     runStart: DateTime): DateTime =
  let
    tstamp = timestamp
      .align(count = 9, padding = '0')
      .parse("HHmmssfff")
    tdiff = (evDuration * eventNumber.float).round(places = 3)
  result = runStart + initDuration(seconds = tdiff.round.int,
                                   milliseconds = (tdiff - tdiff.trunc).round.int)
  # replace time of evDate
  result.second = tstamp.second
  result.minute = tstamp.minute
  result.hour = tstamp.hour

proc calcTimeOfEvent(runTime, totalEvents, eventNumber: int,
                     timestamp: string,
                     runStart: DateTime): DateTime =
  let evDuration = runTime.float / totalEvents.float
  result = calcTimeOfEvent(evDuration, eventNumber, timestamp, runStart)

proc fixOldTosTimestamps(runHeader: Table[string, string],
                         events: var seq[Event]) =
  ## applies the correct time stamp to the events in `events` for old TOS
  ## data, since each event only has the 24h time associated to it.
  let
    runTime = runHeader["runTime"].parseInt
    totalEvents = runHeader["totalEvents"].parseInt
    evDuration = runTime.float / totalEvents.float
    runStart = runHeader["dateTime"].parse(TosDateString)
  for ev in mitems(events):
    # walk over all events and replace the `timestamp` field with a corrected value
    let
      evNumber = ev.evHeader["eventNumber"].parseInt
      evDate = calcTimeOfEvent(evDuration, evNumber, ev.evHeader["timestamp"],
                               runStart)
    ev.evHeader["timestamp"] = $evDate.toTime.toUnix

proc readAndProcessInGrid*(listOfFiles: seq[string],
                           runNumber: int,
                           rfKind: RunFolderKind): ProcessedRun =
  ## Calls the procs to read InGrid data and hands it to the processing proc
  ## inputs:
  ##   listOfFiles: all the files to be read
  ##   runNumber: run we're reading from
  ##   rfKind: old or new TOS data
  # read the raw event data into a seq of FlowVars
  info "list of files ", listOfFiles.len
  let ingrid = readRawInGridData(listOfFiles, rfKind)
  var sortedIngrid = sortReadInGridData(ingrid, rfKind)

  # only continue, if any fails in input
  # This may not be the case if proc is used in a dynamic environment!
  if sortedIngrid.len > 0:
    # to extract the run header, we only need any element of
    # the data. For `rfNewTos` and `rfOldTos` the run header is
    # equivalent to the event header. For `rfSrsTos` it's different
    let runHeader = getRunHeader(sortedIngrid[0], runNumber, rfKind)
    case rfKind
    of rfOldTos: fixOldTosTimestamps(runHeader, sortedIngrid)
    else: discard
    # process the data read into seq of FlowVars, save as result
    let run = createRun(runHeader, runNumber, sortedIngrid, rfKind)
    result = processRawInGridData(run)

proc processAndWriteFadc(run_folder: string, runNumber: int, h5f: var H5FileObj) =
  # for the FADC we call a single function here, which works on
  # the FADC files in a buffered way, always reading 1000 FADC
  # filsa at a time.
  let mem1 = getOccupiedMem()
  info "occupied memory before fadc $# \n\n" % [$mem1]
  readWriteFadcData(run_folder, runNumber, h5f)
  info "FADC took $# data" % $(getOccupiedMem() - mem1)

proc createProcessedTpx3Run(data: seq[Tpx3Data]): ProcessedRun =
  # add new cluster if diff in time larger than 50 clock cycles
  ## TODO: make config.toml adjustable?
  const cutoff = 50
  result = computeTpx3RunParameters(data, clusterTimeCutoff = cutoff)
  result.nChips = 1 ## TODO: allow multiple chips, find out where to best read from input file
  result.chips = @[(name: "W15 E5", number: 0)]
  result.runNumber = 0
  #result.runHeader =

proc processAndWriteSingleRun(h5f: var H5FileObj, run_folder: string,
                              flags: set[RawFlagKind], runType: RunTypeKind = rtNone) =
  ## proc to process and write a single run
  ## inputs:
  ##     h5f: var H5FileObj = mutable copy of the H5 file object to which we will write
  ##         the data
  ##     flags: set[RawFlagKind] = flags indicating different settings, e.g. `nofadc`
  const batchsize = FILE_BUFSIZE * 2
  var attrsWritten = false
  var nChips: int
  # parse config toml file
  let cfgTable = parseTomlConfig()
  let plotOutPath = cfgTable["RawData"]["plotDirectory"].getStr


  let (_, runNumber, rfKind, _) = isTosRunFolder(runFolder)
  var files = getSortedListOfFiles(run_folder,
                                   EventSortType.fname,
                                   EventType.InGridType,
                                   rfKind)
  let plotDirPrefix = h5f.genPlotDirname(plotOutPath, PlotDirRawPrefixAttr)
  batchFiles(files, batchsize):
    let r = readAndProcessInGrid(files[0 .. ind_high], runNumber, rfKind)
    if r.events.len > 0:
      nChips = r.nChips

      if not attrsWritten:
        # create datasets in H5 file
        initInGridInH5(h5f, runNumber, nChips, batchsize)
        # now init attributes
        writeInGridAttrs(h5f, r, rfKind, runType)
        attrsWritten = true

      for chip in 0 ..< nChips:
        plotOccupancy(squeeze(r.occupancies[chip,_,_]),
                      plotDirPrefix, r.runNumber, chip)
      writeProcessedRunToH5(h5f, r)
      info "Size of total ProcessedRun object = ", sizeof(r)
    else:
      warn "Skipped writing to file, since ProcessedRun contains no events!"

  ####################
  # Create Hardlinks #
  ####################
  # linkRawToReco(h5f, runNumber, nChips)

  # dump sequences to file
  #dumpToTandHits(folder, runType, r.tots, r.hits)

  if rfNoFadc notin flags:
    processAndWriteFadc(runFolder, runNumber, h5f)

  # finally once we're done, add `rawDataFinished` attribute
  runFinished(h5f, runNumber)
  # TODO: write all other settings to file too? e.g. `nofadc`,
  # `ignoreRunList` etc?

proc askNoRunListContinue(): bool =
  var buf: string
  while true:
    buf = readLineFromStdin("Do you want to continue? (Y/n)")
    case buf.normalize
    of "": return true
    of "y", "yes": return true
    of "n", "no": return false
    else: continue

proc handleTimepix1(folder: string, runType: RunTypeKind, outfile: string,
                    flags: set[RawFlagKind], configFile: string) =
  # first check whether given folder is valid run folder
  let (is_run_folder, runNumber, rfKind, contains_run_folder) = isTosRunFolder(folder)
  info "Is run folder       : ", is_run_folder
  info "Contains run folder : ", contains_run_folder

  var flags = flags # mutable local copy
  if rfKind == rfOldTos:
    # in case of old TOS runs, there never was a detector with an FADC
    # so force `nofadc`
    info "runKind is " & $rfOldTos & ", hence `nofadc` -> true"
    flags.incl rfNoFadc

  if is_run_folder == true and contains_run_folder == false:
    # hand H5FileObj to processSingleRun, because we need to write intermediate
    # steps to the H5 file for the FADC, otherwise we use too much RAM
    # in order to write the processed run and FADC data to file, open the HDF5 file
    var h5f = H5open(outfile, "rw")
    case rfKind
    of rfOldTos:
      if not oldTosRunListFound:
        warn "Old TOS run list was not found! Expected in " & OldTosRunListPath & "."
        if not askNoRunListContinue():
          # stopping without runlist
          info "Stopping on user desire, due to missing old TOS runlist."
          return
      case runType
      of rtCalibration:
        if rfIgnoreRunList in flags or runNumber.uint16 in oldTosCalibRuns:
          processAndWriteSingleRun(h5f, folder, flags, runType)
        else:
          info &"Run {runNumber} with path {folder} is an invalid run for type {runType}!"
      of rtBackground:
        if rfIgnoreRunList in flags or runNumber.uint16 in oldTosBackRuns:
          processAndWriteSingleRun(h5f, folder, flags, runType)
        else:
          info &"Run {runNumber} with path {folder} is an invalid run for type {runType}!"
      of rtXrayFinger:
        if rfIgnoreRunList in flags or runNumber.uint16 in oldTosXrayRuns:
          processAndWriteSingleRun(h5f, folder, flags, runType)
        else:
          info &"Run {runNumber} with path {folder} is an invalid run for type {runType}!"
      else:
        info &"Run {runNumber} with path {folder} is invalid for type {runType}"
    of rfNewTos, rfSrsTos:
      processAndWriteSingleRun(h5f, folder, flags, runType)
    else:
      raise newException(IOError, "Unknown run folder kind. Cannot read " &
        "events!")
    info "free memory ", getFreeMem()
    info "occupied memory so far $# \n\n" % [$getOccupiedMem()]
    info "Closing h5file with code ", h5f.close()

  elif is_run_folder == false and contains_run_folder == true:
    # in this case loop over all folder again and call processSingleRun() for each
    # run folder
    # open H5 output file to check if run already exists
    var h5f: H5FileObj
    if fileExists(outfile):
      h5f = H5open(outfile, "r")
    else:
      h5f = H5open(outfile, "rw")
    for kind, path in walkDir(folder):
      case kind
      of pcDir, pcLinkToDir:
        # only support run folders, not nested run folders
        echo "occupied memory before run $# \n\n" % [$getOccupiedMem()]
        let (is_rf, runNumber, _, contains_rf) = isTosRunFolder(path)
        if is_rf == true and contains_rf == false:
          if rfOverwrite notin flags and hasRawRun(h5f, runNumber):
            # skip this run if no overwrite or run not in file
            info &"Run number {runNumber} already exists in file, skipping."
            continue
          else:
            info &"Starting raw data manipulation for run number {runNumber}"
            # close h5 file so that subprocess can access it
            discard h5f.close()

          let program = getAppFilename()
          var command = program & " " & path
          for i in 2 .. paramCount():
            let c = paramStr(i)
            command = command & " " & c
          info "Calling command ", command
          let errC = execCmd(command)
          if errC != 0:
            quit("Subprocess failed with " & $errC)

          # reopen the hdf5 file
          h5f = H5open(outfile, "r")

          # TODO: the following is the normal code. However, it leaks memory. That's
          # why we currently just call this script on the subfolder
          # var h5f = H5file(outfile, "rw")
          # processAndWriteSingleRun(h5f, path, nofadc)
          # dumpNumberOfInstances()
          # GC_fullCollect()
          # dumpNumberOfInstances()
        info "occupied memory after gc $#" % [$getOccupiedMem()]
      else:
        # other types will be skipped
        discard
  elif is_run_folder == true and contains_run_folder == true:
    logging.error "Currently not implemented to run over nested run folders."
    quit()
  else:
    logging.error "No run folder found in given path."
    quit()

proc handleTimepix3(h5file: string, runType: RunTypeKind, outfile: string, flags: set[RawFlagKind]) =
  ## handles converting a Timepix3 input file from tpx3-daq / basil format to required TPA format
  info "Converting Tpx3 data from " & $h5file & " and storing it in " & $outfile
  var h5fout = H5File(outfile, "rw")

  const nChips = 1 ## NOTE: so far we just use 1 chip for simplicity
  # parse config toml file
  let cfgTable = parseTomlConfig()
  let plotOutPath = cfgTable["RawData"]["plotDirectory"].getStr

  let plotDirPrefix = h5fout.genPlotDirname(plotOutPath, PlotDirRawPrefixAttr)

  var h5f = H5File(h5file, "rw")
  const tpx3Buf = 50_000_000
  let dset = h5f["interpreted/hit_data_0".dset_str]
  let pixNum = dset.shape[0]
  let batches = ceil(pixNum.float / tpx3Buf.float).int
  var oldIdx = 0
  var countIdx = min(tpx3Buf, pixNum)
  var attrsWritten = false
  var runNumber = 0
  for idx in 0 ..< batches:
    var data = dset.read_hyperslab(Tpx3Data, @[oldIdx, 0],
                                   count = @[countIdx, 1])
    oldIdx += countIdx
    countIdx = if oldIdx + tpx3Buf < pixNum: tpx3Buf
               else: pixNum - oldIdx
    let r = createProcessedTpx3Run(data)
    if r.events.len > 0:
      if not attrsWritten:
        runNumber = r.runNumber
        #nChips = processedRun.nChips
        # create datasets in H5 file
        initInGridInH5(h5fout, runNumber, nChips, batchsize = FILE_BUFSIZE)
        # now init attributes
        writeInGridAttrs(h5fout, r, rfUnknown, runType)
        for chip in 0 ..< nChips:
          plotOccupancy(squeeze(r.occupancies[chip,_,_]),
                        plotDirPrefix, runNumber, chip)
        attrsWritten = true
      writeProcessedRunToH5(h5fout, r)
      info "Size of total ProcessedRun object = ", sizeof(r)
    else:
      warn "Skipped writing to file, since ProcessedRun contains no events!"

  discard h5f.close()
  # finally once we're done, add `rawDataFinished` attribute
  runFinished(h5fout, runNumber)
  discard h5fout.close()

proc main() =

  # use the usage docstring to generate an CL argument table
  let args = docopt(doc)
  echo args

  #echo oldTosBackRuns

  let folder = $args["<folder>"]
  var runTypeStr = $args["--runType"]
  var runType: RunTypeKind
  var flags: set[RawFlagKind]
  var outfile = $args["--out"]
  if runTypeStr != "nil":
    runType = parseRunType(runTypeStr)
  if outfile == "nil":
    outfile = "run_file.h5"
  if $args["--nofadc"] == "true":
    flags.incl rfNoFadc
  if $args["--ignoreRunList"] != "false":
    flags.incl rfIgnoreRunList
  if $args["--overwrite"] == "true":
    flags.incl rfOverwrite
  let tpx3File = if $args["--tpx3"] == "nil": "" else: $args["--tpx3"]
  if tpx3File.len > 0:
    flags.incl rfTpx3

  echo &"Flags are is {flags}"
  let t0 = epochTime()

  if rfTpx3 notin flags:
    handleTimepix1(folder, runType, outfile, flags)
  else:
    handleTimepix3(tpx3File, runType, outfile, flags)

  info "Processing all given runs took $# minutes" % $( (epochTime() - t0) / 60'f )

when isMainModule:
  main()
