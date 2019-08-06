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
import os, osproc, logging
import re
import sequtils, sugar
import algorithm
import tables
import times
#import threadpool
import threadpool_simple
import memfiles
import strutils, strformat, parseutils
import docopt
import typetraits
import sets
import macros
#import nimprof

# InGrid-module
import fadc_helpers
import helpers/utils
import tos_helpers
import ingrid_types

# other modules
import seqmath
import nimhdf5
import arraymancer
import zero_functional

# global experimental pragma to use parallel: statement in readRawInGridData()
{.experimental.}

## TODO:
## change the we write to H5 to appropriate types, e.g.
## x / y pixel coordinates may as well be written as uint8
## FADC raw as uint16
## all flags as uint8
## ToT values as uint16
## hits as uint16

type
  RawFlagKind = enum
    rfIgnoreRunList, rfOverwrite, rfNoFadc

const FILE_BUFSIZE = 25000

##############################
# create globals for 2014/15 run list
##############################
# TOD: should this actually be in here?
const OldTosRunlist = "Runlist-CAST-D03-W0063.csv"
const AppDir = getProjectPath()
var TpxDir {.compileTime.} = ""
static:
  discard parseUntil(AppDir, TpxDir, "TimepixAnalysis")
  TpxDir = TpxDir / "TimepixAnalysis/resources/" / OldTosRunList
const OldTosRunListPath = TpxDir
when fileExists(OldTosRunListPath):
  let oldTosCalibRuns = parseOldTosRunlist(OldTosRunListPath, rtCalibration)
  let oldTosBackRuns  = parseOldTosRunlist(OldTosRunListPath, rtBackground)
  let oldTosXrayRuns  = parseOldTosRunlist(OldTosRunListPath, rtXrayFinger)
else:
  static:
    hint("Compiling without 2014/15 run list")
  const oldTosCalibRuns = ""
  const oldTosBackRuns  = ""
  const oldTosXrayRuns  = ""

const docStr = """
InGrid raw data manipulation.

Usage:
  raw_data_manipulation <folder> [options]
  raw_data_manipulation <folder> --runType <type> [options]
  raw_data_manipulation <folder> --out=<name> [--nofadc] [--runType=<type>] [--ignoreRunList] [options]
  raw_data_manipulation <folder> --nofadc [options]
  raw_data_manipulation -h | --help
  raw_data_manipulation --version

Options:
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

template inGridGroupNames(runNumber: int): (string, string, string, string) =
  ## returns the names of all groups in a H5 file related to InGrid
  let
    # group name for raw data
    groupName = getGroupNameForRun(runNumber)
    # group name for reconstructed data
    recoGroupName = getRecoNameForRun(runNumber)
    # create datatypes for variable length data
    chipGroupName = group_name & "/chip_$#" #% $chp
    combineGroupName = getRawCombineName()
  var result = (groupName, recoGroupName, chipGroupName, combineGroupName)
  result

template inGridGroups(h5f: var H5FileObj,
                      nChips: int = 0,
                      forFadc: static[bool] = false):
         (H5Group, H5Group, H5Group, seq[H5Group]) =
  ## template to get the H5 groups for
  ## Note: some variables need to be defined in the calling scope! This
  ## is why this is a template!
  var
    # create the groups for the run and reconstruction data
    runGroup = h5f.create_group(groupName)
    recoGroup = h5f.create_group(recoGroupName)
  when forFadc == false:
    var
      # combined data group
      combineGroup = h5f.create_group(combineGroupName)
      # create group for each chip
      chipGroups = mapIt(toSeq(0 ..< nChips), h5f.create_group(chipGroupName % $it))
    var result = (runGroup, recoGroup, combineGroup, chipGroups)
  else:
    var
      fadcCombine = h5f.create_group(combineRecoBasenameFadc())
    var result = (runGroup, recoGroup, fadcCombine, newSeq[H5Group](0))
  result

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

proc getTotHitOccDsetNames(chipGroupName: string,
                           nChips: int):
                          (seq[string], seq[string], seq[string]) =
  let
    totDsetNames = toSeq(0 ..< nChips).mapIt((chipGroupName % $it) & "/ToT")
    hitDsetNames = toSeq(0 ..< nChips).mapIt((chipGroupName % $it) & "/Hits")
    occDsetNames = toSeq(0 ..< nChips).mapIt((chipGroupName % $it) & "/Occupancy")
  result = (totDsetNames, hitDsetNames, occDsetNames)

proc getTotHitOccDsets(h5f: var H5FileObj, chipGroupName: string, nChips: int):
                      (seq[H5DataSet], seq[H5DataSet], seq[H5DataSet]) =
  let (totDsetNames,
       hitDsetNames,
       occDsetNames) = getTotHitOccDsetNames(chipGroupName, nChips)

  var
    totDset = totDsetNames.mapIt(h5f[it.dset_str])
    hitDset = hitDsetNames.mapIt(h5f[it.dset_str])
    occDset = occDsetNames.mapIt(h5f[it.dset_str])
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
    var ind_high {.inject.} = bufsize
    if files.high < bufsize:
      ind_high = files.high

    # perform actions as desired
    actions

    info "... removing read elements from list"
    # sequtils.delete removes the element with ind_high as well!
    files.delete(0, ind_high)

proc batchFileReading[T](files: var seq[string],
                         rfKind: RunFolderKind = rfNewTos,
                         bufsize: int = FILE_BUFSIZE):
                          seq[FlowVar[ref T]] {.inline.} =
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
        buf_seq = readListOfInGridFiles(files[0..ind_high], rfKind)
      else:
        raise newException(IOError, "Unknown run folder kind. Cannot read " &
          "event files!")
    elif T is FadcFile:
      buf_seq = readListOfFadcFiles(files[0..ind_high])

    info "... and concating buffered sequence to result"
    result = concat(result, buf_seq)
    count += bufsize
  info "All files read. Number = " & $len(result)
  info "Reading took $# seconds" % $(epochTime() - t0)
  info "Compared with starting files " & $n_files


proc readRawInGridData*(listOfFiles: seq[string],
                        rfKind: RunFolderKind):
                          seq[FlowVar[ref Event]] =
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

proc sortReadInGridData(rawIngridNil: seq[FlowVar[ref Event]],
                        rfKind: RunFolderKind): seq[Event] =
  ## sorts the seq of FlowVars according to event numbers again, otherwise
  ## h5 file is all mangled
  info "Sorting data..."
  # case on the old TOS' data storage and new TOS' version
  let t0 = epochTime()

  # first filter out all nil references (which can happen due to
  # broken files)
  var rawIngrid = newSeqOfCap[Event](rawInGridNil.len)
  for evFut in rawInGridNil:
    if not evFut.isNil:
      let ev = ^evFut
      if not ev.isNil:
        rawIngrid.add ev[]
  if rawInGridNil.len != rawInGrid.len:
    warn &"Removed {rawInGridNil.len - rawInGrid.len} broken InGrid files!"

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
    echoFilesCounted(count, msg = " files processed.")

  # use first event of run to fill event header. Fine, because event
  # header is contained in every file
  result.runHeader = runHeader
  result.chips = run.chips
  result.nChips = nChips
  result.events = events
  result.tots = tot_run -->> map(it -->
                                 flatten().
                                 to(seq[uint16])) --> to(seq[seq[uint16]])# flatten 1 level
  result.hits = hits
  result.occupancies = occ

{.experimental.}
proc processFadcData(fadcFilesNil: seq[FlowVar[ref FadcFile]]): ProcessedFadcData {.inline.} =
  ## proc which performs all processing needed to be done on the raw FADC
  ## data. Starting from conversion of FadcFiles -> FadcData, but includes
  ## calculation of minimum and check for noisy events
  # sequence to store the indices needed to extract the 0 channel

  # filter out possible nil refs, due to broken files
  var fadcFiles = newSeqOfCap[FadcFile](fadcFilesNil.len)
  for evFut in fadcFilesNil:
    if not evFut.isNil:
      let ev = ^evFut
      if not ev.isNil:
        fadcFiles.add ev[]

  if fadcFilesNil.len != fadcFiles.len:
    warn &"Removed {fadcFilesNil.len - fadcFiles.len} broken FADC files!"

  let
    fadc_ch0_indices = getCh0Indices()
    ch_len = ch_len()
    pedestal_run = getPedestalRun()
    nEvents = fadcFiles.len
    # we demand at least 4 dips, before we can consider an event as noisy
    n_dips = 4
    # the percentile considered for the calculation of the minimum
    min_percentile = 0.95

  # convert FlowVars of FadcFiles to sequence of FadcFiles
  result.raw_fadc_data = newSeq[seq[uint16]](nEvents)
  result.fadc_data = zeros[float]([nEvents, ch_len])
  result.trigRecs = newSeq[int](nEvents)
  result.noisy = newSeq[int](nEvents)
  result.minVals = newSeq[float](nEvents)
  result.eventNumber = newSeq[int](nEvents)
  let t0 = epochTime()
  # TODO: parallelize this somehow so that it's faster!
  for i, ev in fadcFiles:
    result.raw_fadc_data[i] = ev.data
    result.trigRecs[i]      = ev.trigRec
    result.eventNumber[i]   = ev.eventNumber
    let fadc_dat = ev.fadcFileToFadcData(pedestal_run, fadc_ch0_indices).data
    result.fadc_data[i, _]  = fadc_dat.reshape([1, ch_len])
    result.noisy[i]         = fadc_dat.isFadcFileNoisy(n_dips)
    result.minVals[i]       = fadc_dat.calcMinOfPulse(min_percentile)

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
  let
    ch_len = ch_len()
    all_ch_len = all_ch_len()

  let
    group_name = getGroupNameForRun(runNumber) & "/fadc"
    reco_group_name = getRecoNameForRun(runNumber) & "/fadc"

  # use dirty template to get groups for the group names
  var (runGroup, recoGroup, fadcCombine, _) = inGridGroups(h5f, forFadc = true)
  var
    # create the datasets for raw data etc
    # NOTE: we initialize all datasets with a size of 0. This means we need to extend
    # it immediately. However, this allows us to always (!) simply extend and write
    # the data to dset.len onwards!
    raw_fadc_dset    = h5f.create_dataset(rawFadcBasename(runNumber), (0, all_ch_len),
                                          uint16,
                                          chunksize = @[batchsize, all_ch_len],
                                          maxshape = @[int.high, all_ch_len],
                                          filter = filter)
    fadc_dset        = h5f.create_dataset(fadcDataBasename(runNumber), (0, ch_len),
                                          float,
                                          chunksize = @[batchsize, ch_len],
                                          maxshape = @[int.high, ch_len],
                                          filter = filter)
    trigrec_dset     = h5f.create_dataset(trigrecBasename(runNumber), (0, 1),
                                          int,
                                          chunksize = @[batchsize, 1],
                                          maxshape = @[int.high, 1],
                                          filter = filter)
    # dataset of eventNumber
    eventNumber_dset = h5f.create_dataset(eventNumberBasename(runNumber), (0, 1),
                                          int,
                                          chunksize = @[batchsize, 1],
                                          maxshape = @[int.high, 1],
                                          filter = filter)

    # dataset stores flag whether FADC event was a noisy one (using our algorithm)
    noisy_dset       = h5f.create_dataset(noiseBasename(runNumber), (0, 1),
                                          int,
                                          chunksize = @[batchsize, 1],
                                          maxshape = @[int.high, 1],
                                          filter = filter)
    # dataset stores minima of each FADC event, dip voltage
    minVals_dset     = h5f.create_dataset(minValsBasename(runNumber), (0, 1),
                                          float,
                                          chunksize = @[batchsize, 1],
                                          maxshape = @[int.high, 1],
                                          filter = filter)

  # write attributes to FADC groups
  # read the given FADC file and extract that information from it
  let fadc_for_attrs = readFadcFile(filename)
  # helper sequence to loop over both groups to write attrs
  var group_seq = @[run_group, reco_group]
  for group in mitems(group_seq):
    group.attrs["posttrig"] = fadc_for_attrs.posttrig
    group.attrs["pretrig"] = fadc_for_attrs.pretrig
    group.attrs["n_channels"] = fadc_for_attrs.n_channels
    group.attrs["channel_mask"] = fadc_for_attrs.channel_mask
    group.attrs["frequency"] = fadc_for_attrs.frequency
    group.attrs["sampling_mode"] = fadc_for_attrs.sampling_mode
    group.attrs["pedestal_run"] = if fadc_for_attrs.pedestal_run == true: 1 else: 0

proc writeFadcDataToH5(h5f: var H5FileObj, runNumber: int, f_proc: ProcessedFadcData) =
  # proc to write the current FADC data to the H5 file
  # now write the data
  let
    reco_group_name = getRecoNameForRun(runNumber)
    raw_name = rawFadcBasename(runNumber)
    reco_name = fadcDataBasename(runNumber)
    trigRec_name = trigRecBasename(runNumber)
    eventNumber_name = eventNumberBasename(runNumber)
    ch_len = ch_len()
    all_ch_len = all_ch_len()
    nEvents = f_proc.raw_fadc_data.len
  var
    raw_fadc_dset = h5f[raw_name.dset_str]
    fadc_dset = h5f[reco_name.dset_str]
    trigRec_dset = h5f[trigRec_name.dset_str]
    eventNumber_dset = h5f[eventNumber_name.dset_str]
    noisy_dset = h5f[noiseBasename(runNumber).dset_str]
    minVals_dset = h5f[minValsBasename(runNumber).dset_str]

  info raw_fadc_dset.shape
  info raw_fadc_dset.maxshape
  info fadc_dset.shape
  info fadc_dset.maxshape
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
  raw_fadc_dset.resize((newsize, all_ch_len))
  fadc_dset.resize((newsize, ch_len))
  trigRec_dset.resize((newsize, 1))
  eventNumber_dset.resize((newsize, 1))
  noisy_dset.resize((newsize, 1))
  minVals_dset.resize((newsize, 1))

  # now write the data
  let t0 = epochTime()
  # TODO: speed this up
  # write using hyperslab
  info "Trying to write using hyperslab! from $# to $#" % [$oldsize, $newsize]
  raw_fadc_dset.write_hyperslab(f_proc.raw_fadc_data,
                                offset = @[oldsize, 0],
                                count = @[nEvents, all_ch_len])
  fadc_dset.write_hyperslab(f_proc.fadc_data.toRawSeq,
                            offset = @[oldsize, 0],
                            count = @[nEvents, ch_len])
  trigRec_dset.write_hyperslab(f_proc.trigRecs,
                               offset = @[oldsize, 0],
                               count = @[nEvents, 1])
  eventNumber_dset.write_hyperslab(f_proc.eventNumber,
                               offset = @[oldsize, 0],
                               count = @[nEvents, 1])
  noisy_dset.write_hyperslab(f_proc.noisy,
                             offset = @[oldsize, 0],
                             count = @[nEvents, 1])
  minVals_dset.write_hyperslab(f_proc.minVals,
                               offset = @[oldsize, 0],
                               count = @[nEvents, 1])
  info "Writing of FADC data took $# seconds" % $(epochTime() - t0)

proc finishFadcWriteToH5(h5f: var H5FileObj, runNumber: int) =
  # proc to finalize the last things we need to write for the FADC data
  # for now only hardlinking to combine group
  let
    noisy_target = noiseBasename(runNumber)
    minVals_target = minValsBasename(runNumber)
    noisy_link = combineRecoBasenameNoisy(runNumber)
    minVals_link = combineRecoBasenameMinVals(runNumber)

  h5f.create_hardlink(noisy_target, noisy_link)
  h5f.create_hardlink(minVals_target, minVals_link)

proc readProcessWriteFadcData(run_folder: string, runNumber: int, h5f: var H5FileObj) =
  ## given a run_folder it reads all fadc files (data<number>.txt-fadc),
  ## processes it (FadcFile -> FadcData) and writes it to the HDF5 file

  # get a sorted list of files, sorted by inode
  var
    files: seq[string] = getSortedListOfFiles(run_folder,
                                              EventSortType.fname,
                                              EventType.FadcType,
                                              RunFolderKind.rfUnknown)
    raw_fadc_data: seq[FlowVar[ref FadcFile]] = @[]
    # variable to store the processed FADC data
    f_proc: ProcessedFadcData
    # in case of FADC data we cannot afford to read all files into memory before
    # writing some to HDF5, because the memory overhead from storing all files
    # in seq[string] is too large (17000 FADC files -> 10GB needed!)
    # thus already perform batching here
    files_read: seq[string] = @[]
  # use batchFiles template to work on 1000 files per batch

  if files.len == 0:
    # in case there are no FADC files, return from this proc
    return

  const batchsize = 1000

  # before we start iterating over the files, initialize the H5 file
  h5f.initFadcInH5(runNumber, batchsize, files[0])
  batchFiles(files, batchsize - 1):
    # batch in 1000 file pieces
    var mfiles = files[0..ind_high]
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
  finishFadcWriteToH5(h5f, runNumber)

proc initInGridInH5*(h5f: var H5FileObj, runNumber, nChips, batchsize: int) =
  ## This proc creates the groups and dataset for the InGrid data in the H5 file
  ## inputs:
  ##   h5f: H5file = the H5 file object of the writeable HDF5 file
  ##   ?
  # create variables for group names (NOTE: dirty template!)
  let (groupName, recoGroupName, chipGroupName, combineGroupName) = inGridGroupNames(runNumber)
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
    # group for raw data
    run_group = h5f.create_group(group_name)
    # group for reconstructed data
    reco_group = h5f.create_group(reco_group_name)
    # combined data group
    combine_group = h5f.create_group(combine_group_name)

    # create group for each chip
    chip_groups = mapIt(toSeq(0 ..< nChips), h5f.create_group(chipGroupName % $it))
    # datasets are chunked in the batchsize we read. Size originally 0
    x_dsets  = mapIt(chip_groups, h5f.datasetCreation(it.name & "/raw_x", ev_type_xy))
    y_dsets  = mapIt(chip_groups, h5f.datasetCreation(it.name & "/raw_y", ev_type_xy))
    ch_dsets = mapIt(chip_groups, h5f.datasetCreation(it.name & "/raw_ch", ev_type_ch))

    # datasets to store the header information for each event
    evHeadersDsetTab = eventHeaderKeys.mapIt(
      (it,
       h5f.datasetCreation(group_name & "/" & it, int))
    ).toTable
    # TODO: add string of datetime as well
    #dateTimeDset = h5f.create_dataset(joinPath(group_name, "dateTime"), nEvents, string)

    # other single column data
    durationDset = h5f.datasetCreation(joinPath(group_name, "eventDuration"), float)
  let names = mapIt(toSeq(0 ..< nChips), chipGroupName % $it)
  var
    totDset = mapIt(names, h5f.datasetCreation(it & "/ToT", uint16))
    hitDset = mapIt(names, h5f.datasetCreation(it & "/Hits", uint16))
    # use normal dataset creation proc, due to static size of occupancies
    occDset = mapIt(names, h5f.create_dataset(it & "/Occupancy", (256, 256), int))

proc writeInGridAttrs*(h5f: var H5FileObj, run: ProcessedRun,
                       rfKind: RunFolderKind,
                       runType: RunTypeKind) =

  # use dirty template to define variables of group names
  let (groupName,
       recoGroupName,
       chipGroupName,
       combineGroupName) = inGridGroupNames(run.runNumber)
  # use another dirty template to get the groups for the names
  var (runGroup, recoGroup, combineGroup, chipGroups) = inGridGroups(h5f, run.nChips)
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
      run_group.attrs[it]  = att
      reco_group.attrs[it] = att
  for it in asString:
    if it in run.runHeader:
      let att = run.runHeader[it]
      run_group.attrs[it] = att
      reco_group.attrs[it] = att

  # initialize the attribute for the current number of stored events to 0
  run_group.attrs["numEventsStored"] = 0
  # write attributes for each chip
  var i = 0
  for grp in mitems(chip_groups):
    grp.attrs["chipNumber"] = run.chips[i].number
    grp.attrs["chipName"]   = run.chips[i].name
    # initialize the attribute for the current number of stored events to 0
    grp.attrs["numEventsStored"] = 0
    inc i

  # finally write run type to base runs / reconstruction groups
  var
    rawG = h5f["runs".grp_str]
    recoG = h5f["reconstruction".grp_str]
  rawG.attrs["runType"] = $runType
  recoG.attrs["runType"] = $runType
  rawG.attrs["runFolderKind"] = $rfKind
  recoG.attrs["runFolderKind"] = $rfKind
  # Currently hardcode the number of chips we use.
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
  rawG.attrs["centerChip"] = centerChip
  recoG.attrs["centerChip"] = centerChip
  run_group.attrs["centerChip"] = centerChip
  reco_group.attrs["centerChip"] = centerChip
  rawG.attrs["centerChipName"] = centerName
  recoG.attrs["centerChipName"] = centerName
  run_group.attrs["centerChipName"] = centerName
  reco_group.attrs["centerChipName"] = centerName

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
    # given chip number, we can write the data of this event to the correct dataset
    # need to subtract number of events already in file (min(evNumber) > 0!)
    # NOTE: we now simply enumerate the number of events we read. There should be no reason
    # why we should have to consider the real event numbers
    # TODO: if you read this in the future and everything works, remove the
    # next two lines :)
    let evNumberRaw = parseInt(event.evHeader["eventNumber"])
    let evNumber = parseInt(event.evHeader["eventNumber"]) - startEvent

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
  let (groupName,
       recoGroupName,
       chipGroupName,
       combineGroupName) = inGridGroupNames(run.runNumber)
  # another dirty template to get the groups for the names
  var (runGroup, recoGroup, combineGroup, chipGroups) = inGridGroups(h5f, nChips)
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
    #dateTimeDset = h5f.create_dataset(joinPath(group_name, "dateTime"), nEvents, string)

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
  #let all = x_dsets[0].all

  template writeHyper(dset: untyped, data: untyped): untyped =
    dset.write_hyperslab(data, offset = @[oldsize, 0], count = @[nEvents, 1])

  for i in 0 ..< nChips:
    withDebug:
      info "Writing dsets ", i, " size x ", x_dsets.len
    x_dsets[i].resize((newsize, 1))
    y_dsets[i].resize((newsize, 1))
    ch_dsets[i].resize((newsize, 1))
    withDebug:
      info "Shape of x ", x[i].len, " ", x[i].shape
      info "Shape of dset ", x_dsets[i].shape
    x_dsets[i].writeHyper(x[i])
    y_dsets[i].writeHyper(y[i])
    ch_dsets[i].writeHyper(ch[i])

  for key, dset in mpairs(evHeadersDsetTab):
    withDebug:
      info "Writing $# in $#" % [$key, $dset]
    dset.resize((newsize, 1))
    dset.writeHyper(evHeaders[key])

  # write other single column datasets
  durationDset.resize((newsize, 1))
  durationDset.writeHyper(duration)
  info "took a total of $# seconds" % $(epochTime() - t0)

  ####################
  #  Reconstruction  #
  ####################

  # TODO: maybe this needs to be done in a pass after everything has been done?
  # at least for occupancy?

  info "ToTs shape is ", run.tots.shape
  info "hits shape is ", run.hits.shape
  # into the reco group name we now write the ToT and Hits information
  var (totDsets, hitDsets, occDsets) = getTotHitOccDsets(h5f, chipGroupName, nChips)

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

    let newTotSize = totDset.shape[0] + tot.len
    let newHitSize = hitDset.shape[0] + hit.len

    let totOldSize = totDset.shape[0]
    totDset.resize((newTotSize, 1))
    hitDset.resize((newHitSize, 1))
    # need to handle ToT dataset differently
    totDset.write_hyperslab(tot.reshape([tot.len, 1]), offset = @[totOldSize, 0], count = @[tot.len, 1])
    hitDset.writeHyper(hit.reshape([hit.len, 1]))
    # before writing the occupancy dataset, we need to read the old, stack the current
    # occupancy on it and finally write the result
    let stackOcc = occDset[int64].toTensor.reshape([256, 256]) .+ occ
    occDset.unsafeWrite(stackOcc.get_data_ptr, stackOcc.size)


proc linkRawToReco(h5f: var H5FileObj, runNumber, nChips: int) =
  ## perform linking from raw group to reco group
  let (groupName,
       recoGroupName,
       chipGroupName,
       combineGroupName) = inGridGroupNames(runNumber)
  let (_, _, eventHeaderKeys) = specialTypesAndEvKeys()
  let (totDsetNames,
       hitDsetNames,
       occDsetNames) = getTotHitOccDsetNames(chipGroupName, nChips)
  let
    durationDsetName = joinPath(groupName, "eventDuration")

  # link over to reconstruction group
  h5f.create_hardlink(durationDsetName, recoGroupName / extractFilename(durationDsetName))
  # create hard links of header data to reco group
  for key in eventHeaderKeys:
    h5f.create_hardlink(joinPath(groupName, key), joinPath(recoGroupName, key))
  for chip in 0 ..< nChips:
    # create hardlinks for ToT and Hits
    h5f.create_hardlink(totDsetNames[chip], combineRawBasenameToT(chip, runNumber))
    h5f.create_hardlink(hitDsetNames[chip], combineRawBasenameHits(chip, runNumber))

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
  readProcessWriteFadcData(run_folder, runNumber, h5f)
  info "FADC took $# data" % $(getOccupiedMem() - mem1)

proc processAndWriteSingleRun(h5f: var H5FileObj, run_folder: string,
                              flags: set[RawFlagKind], runType: RunTypeKind = rtNone) =
  ## proc to process and write a single run
  ## inputs:
  ##     h5f: var H5FileObj = mutable copy of the H5 file object to which we will write
  ##         the data
  ##     flags: set[RawFlagKind] = flags indicating different settings, e.g. `nofadc`
  const batchsize = 50000
  var attrsWritten = false
  var nChips: int

  let (_, runNumber, rfKind, _) = isTosRunFolder(runFolder)
  var files = getSortedListOfFiles(run_folder,
                                   EventSortType.fname,
                                   EventType.InGridType,
                                   rfKind)


  batchFiles(files, batchsize - 1):
    let r = readAndProcessInGrid(files[0 .. ind_high], runNumber, rfKind)
    if r.events.len > 0:
      nChips = r.nChips

      if attrsWritten == false:
        writeInGridAttrs(h5f, r, rfKind, runType)
        # create datasets in H5 file
        initInGridInH5(h5f, runNumber, nChips, batchsize)
        attrsWritten = true

      let a = squeeze(r.occupancies[0,_,_])
      dumpFrameToFile("tmp/frame.txt", a)
      writeProcessedRunToH5(h5f, r)
      info "Size of total ProcessedRun object = ", sizeof(r)
    else:
      warn "Skipped writing to file, since ProcessedRun contains no events!"

  ####################
  # Create Hardlinks #
  ####################
  linkRawToReco(h5f, runNumber, nChips)

  # dump sequences to file
  #dumpToTandHits(folder, runType, r.tots, r.hits)

  if rfNoFadc notin flags:
    processAndWriteFadc(runFolder, runNumber, h5f)

  # finally once we're done, add `rawDataFinished` attribute
  runFinished(h5f, runNumber)
  # TODO: write all other settings to file too? e.g. `nofadc`,
  # `ignoreRunList` etc?

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

  echo &"Flags are is {flags}"

  # first check whether given folder is valid run folder
  let (is_run_folder, runNumber, rfKind, contains_run_folder) = isTosRunFolder(folder)
  info "Is run folder       : ", is_run_folder
  info "Contains run folder : ", contains_run_folder

  if rfKind == rfOldTos:
    # in case of old TOS runs, there never was a detector with an FADC
    # so force `nofadc`
    info "runKind is " & $rfOldTos & ", hence `nofadc` -> true"
    flags.incl rfNoFadc

  let t0 = epochTime()
  if is_run_folder == true and contains_run_folder == false:
    # hand H5FileObj to processSingleRun, because we need to write intermediate
    # steps to the H5 file for the FADC, otherwise we use too much RAM
    # in order to write the processed run and FADC data to file, open the HDF5 file
    var h5f = H5file(outfile, "rw")
    case rfKind
    of rfOldTos:
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
      h5f = H5file(outfile, "r")
    else:
      h5f = H5file(outfile, "rw")
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
          h5f = H5file(outfile, "r")

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

  info "Processing all given runs took $# minutes" % $( (epochTime() - t0) / 60'f )

when isMainModule:
  main()
