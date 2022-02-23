# this file contains the necessary code for reconstruction of X-rays in the raw data
# consists of
#   - finding clusters in data
#   - calculating properties of Events

# nim stdlib
import os, sequtils, sugar, math, tables, strutils, strformat, macros, options
#import threadpool
when not defined(gcDestructors):
  import threadpool_simple
elif true:
  import std / threadpool
else:
  import weave
import logging
import docopt
import times
import sets

# external modules
import nimhdf5
import parsetoml

# custom modules
import tos_helpers
import helpers/utils
import ingrid_types
import calibration
import fadc_analysis

type
  RecoFlagKind = enum
    rfNone, rfCreateFe, rfOnlyEnergy, rfOnlyEnergyElectrons,
    rfOnlyCharge, rfOnlyFeSpec, rfOnlyFadc, rfOnlyGasGain, rfOnlyGainFit,
    rfReadAllRuns

  ConfigFlagKind = enum
    cfNone, cfShowPlots

  DocoptTab = Table[string, docopt.Value]

when defined(linux):
  const commitHash = staticExec("git rev-parse --short HEAD")
else:
  const commitHash = ""
# get date using `CompileDate` magic
const currentDate = CompileDate & " at " & CompileTime

const docTmpl = """
Version: $# built on: $#
InGrid reconstruction and energy calibration.

NOTE: When calling `reconstruction` without any of the `--only_*` flags, the input file
has to be a H5 file resulting from `raw_data_manipulation`. In the other cases the input is
simply a file resulting from a prior `reconstruction` call!

The `--only_*` flags below are given roughly in the order in which the full analysis chain
requires them to be run. If unsure on the order, check the runAnalysisChain.nim file.

Usage:
  reconstruction <HDF5file> --out <name> [options]
  reconstruction <HDF5file> --out <name> [--runNumber <number>] [--create_fe_spec] [options]
  reconstruction <HDF5file> [--runNumber <number>] --only_fadc [options]
  reconstruction <HDF5file> [--runNumber <number>] --only_fe_spec [options]
  reconstruction <HDF5file> [--runNumber <number>] --only_charge [options]
  reconstruction <HDF5file> [--runNumber <number>] --only_gas_gain [options]
  reconstruction <HDF5file> [--runNumber <number>] --only_gain_fit [options]
  reconstruction <HDF5file> [--runNumber <number>] --only_energy_from_e [options]
  reconstruction <HDF5file> [--runNumber <number>] --only_energy <factor> [options]
  reconstruction -h | --help
  reconstruction --version


Options:
  --out <name>            Filename and path of output file
  --runNumber <number>    Only work on this run
  --create_fe_spec        Toggle to create Fe calibration spectrum based on cuts
                          Takes precedence over --calib_energy if set!
  --only_energy <factor>  Toggle to /only/ perform energy calibration using the given factor.
                          Takes precedence over --create_fe_spec if set.
                          If no runNumber is given, performs energy calibration on all runs
                          in the HDF5 file.
  --only_energy_from_e    Toggle to /only/ calculate the energy for each cluster based on
                          the Fe charge spectrum vs gas gain calibration
  --only_fe_spec          Toggle to /only/ create the Fe spectrum for this run and perform the
                          fit of it. Will try to perform a charge calibration, if possible.
  --only_charge           Toggle to /only/ calculate the charge for each TOT value based on
                          the TOT calibration. The `ingridDatabase.h5` needs to be present.
  --only_fadc             If this flag is set, the reconstructed FADC data is used to calculate
                          FADC values such as rise and fall times among others, which are written
                          to the H5 file.
  --config <path>         Path to the configuration file to use.
  -h --help               Show this help
  --version               Show version.
"""
const doc = docTmpl % [commitHash, currentDate]

# the filter we use globally in this file
let filter = H5Filter(kind: fkZlib, zlibLevel: 4)
#let filter = H5Filter(kind: fkNone)
const Chunksize = 10000


################################################################################

# set up the logger
var L = newConsoleLogger()
if not dirExists("logs"):
  createDir("logs")
var fL = newFileLogger("logs/reconstruction.log", fmtStr = verboseFmtStr)
when isMainModule:
  addHandler(L)
  addHandler(fL)

################################################################################

template ch_len(): int = 2560
template all_ch_len(): int = ch_len() * 4

template benchmark(num: int, actions: untyped) {.dirty.} =
  for i in 0 ..< num:
    actions

macro setTabFields(tab: Table[string, seq[seq[typed]]],
                   names: typed, # some array[N, string]
                   chip: int,
                   obj: untyped): untyped =
  ## taking some table `tab` sets all fields of the table taken from an array `names` to the
  ## fields of the object `obj` of the same names as the names in the table
  result = newStmtList()
  let namesImpl = names.getImpl
  for name in namesImpl:
    let field = parseExpr(name.strVal)
    result.add quote do:
      `tab`[`name`][`chip`].add `obj`.`field`

proc initDataTab[T: (float | int), N: int](tab: var Table[string, seq[seq[T]]],
                                           nchips: int,
                                           names: array[N, string]) =
  ## convenienve proc to initialize the seqs inside a table `tab` storing the
  ## data for all `nchips` of the `names`
  for dset in names:
    tab[dset] = newSeq[seq[T]](nchips)
    # initialize the sequences
    for s in mitems(tab[dset]):
      s = @[]

proc createDatasets[N: int](dset_tab: var Table[string, seq[H5DataSet]],
                            h5f: H5File,
                            names: array[N, string],
                            nchips: int,
                            lengths: seq[int],
                            groups: seq[H5Group],
                            dtype: typedesc) =
  ## creates the actual H5 datasets from the names and datatypes and stores them
  ## in the given `dset_tab`. One dataset for each element in `names` of type `dtype`
  ## is created of length `lengths` for each chip.
  for dset in names:
    dset_tab[dset] = newSeq[H5DataSet](nchips)
    for chip in 0 ..< nchips:
      dset_tab[dset][chip] = h5f.create_dataset(groups[chip].name / dset, lengths[chip],
                                                dtype = dtype,
                                                chunksize = @[Chunksize, 1],
                                                maxshape = @[int.high, 1],
                                                filter = filter)

proc writeRecoRunToH5*[T: SomePix](h5f: H5File,
                                   h5fraw: H5File,
                                   reco_run: seq[RecoEvent[T]],
                                   runNumber: int,
                                   timepixVersion: TimepixVersion) =
  ## proc which writes the reconstructed event data from a single run into
  ## the given H5 file. Called after every processed run
  ## inputs:
  ##     `h5f`: H5File = the H5 file object in which to store data
  ##     `h5fraw`: H5File = the H5 file objcect containing the raw data,
  ##       may be the same file as `h5f`.
  ##     `reco_run`: seq[RecoEvent] = sequence of reconstructed events
  ##       to write to file.
  ## outputs: -
  ## throws:
  ##     potentially throws HDF5LibraryError, if a call to the H5 library fails
  # now we need to write into reco group for each chip
  var rawGroup = h5fraw[getGroupNameRaw(runNumber).grp_str]
  let nChips = rawGroup.attrs["numChips", int]

  # for now hardcode the number of chips. easy to change by getting the number
  # simply from a file or whatever
  info "Number of events in total ", reco_run.len

  let
    # start time for timing the write
    t0 = epochTime()
    # group name for reconstructed data
    reco_group_name = getGroupNameReco(runNumber)
    chip_group_name = reco_group_name / "chip_$#"

  const
    # define the names for the datasets which we want to write
    int_cluster_names = getIntClusterNames()
    int_dset_names = getIntDsetNames()
    # name of float datasets, part of geometry, cluster object
    float_geometry_names = getFloatGeometryNames()
    float_cluster_names = getFloatClusterNames()
    float_dset_names = getFloatDsetNames()
  # now parsing all the data is really fucking ugly, thanks to the tons of
  # different variables, which we want to write :( Unfortunately, we cannot
  # simply make that a compound datatype or something. Well
  # create group for each chip
  var chip_groups = mapIt(toSeq(0..<nChips), h5f.create_group(chip_group_name % $it))

  # create a table containing the sequences for int datasets and corresponding names
  var int_data_tab = initTable[string, seq[seq[int]]]()
  var float_data_tab = initTable[string, seq[seq[float]]]()
  int_data_tab.initDataTab(nChips, int_dset_names)
  float_data_tab.initDataTab(nChips, float_dset_names)

  var
    # before we create the datasets, first parse the data we will write
    x  = newSeq[seq[seq[uint8]]](nChips)
    y  = newSeq[seq[seq[uint8]]](nChips)
    ch = newSeq[seq[seq[uint16]]](nChips)
    toa = newSeq[seq[seq[uint16]]](nChips)
    toaCombined = newSeq[seq[seq[uint64]]](nChips)

  for chip in 0 ..< nChips:
    x[chip]  = @[]
    y[chip]  = @[]
    ch[chip] = @[]
  var count = 0
  var tUnpack = epochTime()
  for event in reco_run:
    let
      num = event.event_number
      chip = event.chip_number

    for i, cl in event.cluster:
      x[chip].add(newSeq[uint8](cl.data.len))
      y[chip].add(newSeq[uint8](cl.data.len))
      ch[chip].add(newSeq[uint16](cl.data.len))
      for j in 0..cl.data.high:
        x[chip][^1][j] = cl.data[j].x
        y[chip][^1][j] = cl.data[j].y
        ch[chip][^1][j] = cl.data[j].ch

      int_data_tab.setTabFields(int_cluster_names, chip, cl)
      float_data_tab.setTabFields(float_cluster_names, chip, cl)
      float_data_tab.setTabFields(float_geometry_names, chip, cl.geometry)

      # add event number individually, since it's not part of some object we can
      # use our macro for
      int_data_tab["eventNumber"][chip].add num

      # add the found clusters in ToA
      if timepixVersion == Timepix3:
        toa[chip].add cl.toa
        toaCombined[chip].add cl.toaCombined

  # now that we have the data and now how many elements each type has
  # we can create the datasets
  info "Now creating datasets"
  # define type for variable length pixel data
  let ev_type_xy = special_type(uint8)
  let ev_type_ch = special_type(uint16)

  template datasetCreation(h5f: untyped, name, dlen, `type`: untyped): untyped =
    ## inserts the correct data set creation parameters
    info "Creating dataset ", name
    h5f.create_dataset(name,
                       dlen,
                       dtype = `type`,
                       chunksize = @[Chunksize, 1],
                       maxshape = @[int.high, 1],
                       filter = filter)

  var
    # datasets for x, y and charge
    int_dsets = initTable[string, seq[H5DataSet]]()
    float_dsets = initTable[string, seq[H5DataSet]]()
    # x, y and charge datasets
    x_dsets = mapIt(toSeq(0..<nChips),
                    h5f.datasetCreation(chip_groups[it].name & "/x",
                                        x[it].len,
                                        ev_type_xy))
    y_dsets = mapIt(toSeq(0..<nChips),
                    h5f.datasetCreation(chip_groups[it].name & "/y",
                                        y[it].len,
                                        ev_type_xy))
    ch_dsets = mapIt(toSeq(0..<nChips),
                    h5f.datasetCreation(chip_groups[it].name & "/ToT",
                                        ch[it].len,
                                        ev_type_ch))
    toa_dsets: seq[H5DataSet]
    toa_combined_dsets: seq[H5DataSet]
  if timepixVersion == Timepix3:
    toa_dsets = mapIt(toSeq(0..<nChips),
                      h5f.datasetCreation(chip_groups[it].name & "/ToA",
                                          ch[it].len,
                                          ev_type_ch))
    toa_combined_dsets = mapIt(toSeq(0..<nChips),
                      h5f.datasetCreation(chip_groups[it].name & "/ToACombined",
                                          ch[it].len,
                                          ev_type_ch))

  # variable to store number of events for each chip
  let eventsPerChip = mapIt(x, it.len)
  # now create all datasets and store them in the dataset tables
  int_dsets.createDatasets(h5f, int_dset_names, nChips, eventsPerChip, chip_groups, int)
  float_dsets.createDatasets(h5f, float_dset_names, nChips, eventsPerChip, chip_groups, float)

  # now that we have the datasets, write everything...
  let all = x_dsets[0].all
  # get locations of raw data groups, so that we can copy
  # the attributes
  info "Writing data to datasets"
  let raw_groups = rawDataChipBase(runNumber)
  for chip in 0 ..< nChips:
    for dset in int_dset_names:
      int_dsets[dset][chip][all] = int_data_tab[dset][chip]
    for dset in float_cluster_names:
      float_dsets[dset][chip][all] = float_data_tab[dset][chip]
    for dset in float_geometry_names:
      float_dsets[dset][chip][all] = float_data_tab[dset][chip]
    # pixel specific seqs
    x_dsets[chip][all] = x[chip]
    y_dsets[chip][all] = y[chip]
    ch_dsets[chip][all] = ch[chip]

    if timepixVersion == Timepix3:
      toa_dsets[chip][all] = toa[chip]
      toa_combined_dsets[chip][all] = toaCombined[chip]

    # anyways, write the chip dataset attributes
    let raw_chip_group_name = raw_groups & $chip
    var raw_group = h5fraw[raw_chip_group_name.grp_str]
    let ch_numb = raw_group.attrs["chipNumber", int]
    let ch_name = raw_group.attrs["chipName", string]
    # and write these to the current group
    chip_groups[ch_numb].attrs["chipNumber"] = ch_numb
    chip_groups[ch_numb].attrs["chipName"] = ch_name

iterator readDataFromH5*(h5f: H5File, runNumber: int,
                         timepixVersion: TimepixVersion):
    tuple[chip: int,
          eventData: RecoInputData[Pix]] =
  ## proc to read data from the HDF5 file from `group`
  ## returns the chip number and a sequence containing the pixel data for this
  ## event and its event number
  var chip_base = rawDataChipBase(runNumber)
  for grp in keys(h5f.groups):
    # get the event numbers for this run
    if chip_base in grp:
      let evNumbers = h5f[grp.parentDir / "eventNumber", int]
      # now can start reading, get the group containing the data for this chip
      var group = h5f[grp.grp_str]
      # get the chip number from the attributes of the group
      let chip_number = group.attrs["chipNumber", int]
      # given group and chip number, we can read vlen data
      let vlen_xy = special_type(uint8)
      let vlen_ch = special_type(uint16)
      let
        raw_x  = h5f[grp / "raw_x", vlen_xy, uint8]
        raw_y  = h5f[grp / "raw_y", vlen_xy, uint8]
        raw_ch = h5f[grp / "raw_ch", vlen_ch, uint16]
      var
        raw_toa: seq[seq[uint16]]
        raw_toa_combined: seq[seq[uint64]]
      if timepixVersion == Timepix3:
        raw_toa = h5f[grp / "raw_toa", vlen_ch, uint16]
        raw_toa_combined = h5f[grp / "raw_toa_combined", special_type(uint64), uint64]
      var runPix = newSeq[RecoInputEvent[Pix]](raw_x.len)
      for i in 0 ..< runPix.len:
        let rpix = zipEm(raw_x[i], raw_y[i], raw_ch[i])
        case timepixVersion
        of Timepix1:
          runPix[i] = (pixels: rpix, eventNumber: evNumbers[i],
                       toa: newSeq[uint16](), toaCombined: newSeq[uint64]())
        of Timepix3:
          runPix[i] = (pixels: rpix, eventNumber: evNumbers[i],
                       toa: raw_toa[i], toaCombined: raw_toa_combined[i])

      # and yield them
      yield (chip: chip_number, eventData: run_pix)

{.experimental.}
proc reconstructSingleChip*(data: RecoInputData[Pix],
                            run, chip, searchRadius: int,
                            dbscanEpsilon: float,
                            clusterAlgo: ClusteringAlgorithm,
                            timepixVersion: TimepixVersion
                           ): seq[RecoEvent[Pix]] =
  ## procedure which receives pixel data for a given chip and run number
  ## and performs the reconstruction on it
  ## inputs:
  ##    data: seq[Pixels] = data of pixels for this chip containing pixels for each event
  ##    run: int = run number of run
  ##    chip: int = chip number working on
  info &"Working on chip {chip} in run {run}"
  info &"We have {data.len} events to reconstruct"
  var count = 0
  let numElems = data.len
  result = newSeq[RecoEvent[Pix]](numElems)
  when not defined(gcDestructors):
    let p = newThreadPool()
    var res = newSeq[FlowVar[RecoEvent[Pix]]](numElems)
    for event in 0 ..< numElems:
      if event < result.len:
        res[event] = p.spawn recoEvent(data[event], chip, run, searchRadius,
                                       dbscanEpsilon = dbscanEpsilon,
                                       clusterAlgo = clusterAlgo,
                                       timepixVersion = timepixVersion)
      echoCount(count, 5000, msg = " clusters reconstructed")
    p.sync()
    count = 0
    for event in 0 ..< numElems:
      result[event] = ^(res[event])
      echoCount(count, 5000, msg = " cluster FlowVars unpacked")
  elif true:
    {.error: "Compiling with `--gc:arc` is currently unsupported as we don't have a working threadpool " &
      "that can return objects containing seq / string like data that works with arc...".}
    var res = newSeq[FlowVar[RecoEvent[Pix]]](numElems)
    for event in 0 ..< numElems:
      if event < result.len:
        res[event] = spawn recoEvent(data[event], chip, run, searchRadius,
                                     dbscanEpsilon = dbscanEpsilon,
                                     clusterAlgo = clusterAlgo,
                                     timepixVersion = timepixVersion)
      echoCount(count, 5000, msg = " clusters reconstructed")
    sync()
    count = 0
    for event in 0 ..< numElems:
      result[event] = ^(res[event])
      echoCount(count, 5000, msg = " cluster FlowVars unpacked")
  else:
    {.error: "Compiling with `--gc:arc` is currently unsupported as we don't have a working threadpool " &
      "that can return objects containing seq / string like data that works with arc...".}
    init(Weave)
    var resBuf = cast[ptr UncheckedArray[RecoEvent[Pix]]](result[0].addr)
    parallelFor event in 0 ..< numElems:
      captures: {resBuf, data, chip, run, searchRadius, dbscanEpsilon, clusterAlgo}
      resBuf[event] = recoEvent(data[event], chip, run, searchRadius,
                                dbscanEpsilon = dbscanEpsilon,
                                clusterAlgo = clusterAlgo,
                                timepixVersion = timepixVersion)
      echoCounted(event, 5000, msg = " clusters reconstructed")

proc createAndFitFeSpec(h5f: H5File,
                        runNumber: int,
                        fittingOnly: bool) =
  ## create the Fe spectrum for the run, apply the charge calibration if possible
  ## and then fit to the Fe spectrum
  var centerChip = h5f.getCenterChip(runNumber)
  h5f.createFeSpectrum(runNumber, centerChip)
  try:
    h5f.applyChargeCalibration(runNumber)
  except KeyError as e:
    warn "No charge calibration possible for current one or " &
         "more chips. Exception message:\n" & e.msg
  h5f.fitToFeSpectrum(runNumber, centerChip, fittingOnly)

proc initRecoFadcInH5(h5f, h5fout: H5File, runNumber, batchsize: int) =
  # proc to initialize the datasets etc in the HDF5 file for the FADC. Useful
  # since we don't want to do this every time we call the write function

  if fadcRawPath(runNumber) notin h5f:
    # means `raw_data_manipulation` was run with `--nofadc` or does not have any
    # FADC data (e.g. 2014 data)
    return

  const
    ch_len = ch_len()
    all_ch_len = all_ch_len()
  let groupName = fadcRecoPath(runNumber)
  template datasetCreation(h5f, name, shape, `type`: untyped): untyped =
    ## inserts the correct data set creation parameters
    h5f.create_dataset(name,
                       shape,
                       dtype = `type`,
                       chunksize = @[batchsize, shape[1]],
                       maxshape = @[int.high, shape[1]],
                       filter = filter)
  var
    # NOTE: we initialize all datasets with a size of 0. This means we need to extend
    # it immediately. However, this allows us to always (!) simply extend and write
    # the data to dset.len onwards!
    recoGroup = h5fout.create_group(groupName)
    fadc_dset        = h5fout.datasetCreation(fadcDataBasename(runNumber), (0, ch_len), float)
    # dataset of eventNumber
    eventNumber_dset = h5fout.datasetCreation(eventNumberBasenameReco(runNumber), (0, 1), int)
    # dataset stores flag whether FADC event was a noisyo one (using our algorithm)
    noisy_dset       = h5fout.datasetCreation(noiseBasename(runNumber), (0, 1), int)
    # dataset stores minima of each FADC event, dip voltage
    minVals_dset     = h5fout.datasetCreation(minValsBasename(runNumber), (0, 1), float)

  # write attributes to FADC groups
  let fadcRaw = h5f[fadcRawPath(runNumber).grp_str]
  recoGroup.attrs["posttrig"] = fadcRaw.attrs["posttrig", int64]
  recoGroup.attrs["pretrig"] = fadcRaw.attrs["pretrig", int64]
  recoGroup.attrs["n_channels"] = fadcRaw.attrs["n_channels", int64]
  recoGroup.attrs["channel_mask"] = fadcRaw.attrs["channel_mask", int64]
  recoGroup.attrs["frequency"] = fadcRaw.attrs["frequency", int64]
  recoGroup.attrs["sampling_mode"] = fadcRaw.attrs["sampling_mode", int64]
  recoGroup.attrs["pedestal_run"] = fadcRaw.attrs["pedestal_run", int64]

proc copyOverDataAttrs(h5f, h5fout: H5File, runNumber: int) =
  template copyAttrs(pout, pin: untyped): untyped =
    let recoGrp = h5fout.create_group(pout) # [pout.grp_str]
    let rawGrp = h5f[pin.grp_str]
    recoGrp.copy_attributes(rawGrp.attrs)
  # copy attributes of `runs/run_XYZ -> reconstruction/run_XYZ`
  copyAttrs("/reconstruction", "/runs")
  copyAttrs((recoBase() & $runNumber), (rawDataBase() & $runNumber))
  template copyOver(path, dtype: untyped): untyped =
    # don't care for the returned group here
    discard h5fout.write_dataset(recoBase() & path, h5f[rawDataBase() & path, dtype])
  h5f.visit_file() # is this required?
  let rawGrp = h5f[(rawDataBase() & $runNumber).grp_str]
  info "Copying over data from input ", h5f.name, " to ", h5fout.name
  for dset in items(rawGrp): # copy one level (nested is not copied)
    info "Copying dataset ", dset
    case dset.dtypeAnyKind
    of dkInt64:
      copyOver($runNumber / dset.name.extractFilename, int64)
    of dkFloat64:
      copyOver($runNumber / dset.name.extractFilename, float64)
    else:
      doAssert false, "Unexpected dataset " & $dset & " with base type " & $dset.dtypeAnyKind

proc calcTriggerFractions(h5f: H5File, runNumber: int) =
  ## calculates the fraction of events within a given run of events with
  ## - FADC triggers
  ## - non trivial scinti 1 and 2 triggers
  ## and writes it as attributes to the reconstruction run group
  const dsetNames = ["fadcReadout", "szint1ClockInt", "szint2ClockInt"]

  let fadcData = h5f[recoBase & $runNumber / "fadcReadout", int64]
  let szint1Data = h5f[recoBase & $runNumber / "szint1ClockInt", int64]
  let szint2Data = h5f[recoBase & $runNumber / "szint2ClockInt", int64]

  let fadcRatio = fadcData.filterIt(it == 1).len.float / fadcData.len.float
  let szint1Ratio = szint1Data.filterIt(it > 0).len.float / szint1Data.len.float
  let szint2Ratio = szint2Data.filterIt(it > 0).len.float / szint2Data.len.float
  let s1NonTrivial = szint1Data.filterIt(it > 0 and it < 4095).len.float / szint1Data.len.float
  let s2NonTrivial = szint2Data.filterIt(it > 0 and it < 4095).len.float / szint2Data.len.float

  let grp = h5f[(recoBase & $runNumber).grp_str]
  grp.attrs["fadcTriggerRatio"] = fadcRatio
  grp.attrs["szint1TriggerRatio"] = szint1Ratio
  grp.attrs["szint2TriggerRatio"] = szint2Ratio
  grp.attrs["nonTrivial_szint1TriggerRatio"] = s1NonTrivial
  grp.attrs["nonTrivial_szint2TriggerRatio"] = s2NonTrivial

template recordIterRuns*(base: string, body: untyped): untyped =
  ## Helper template, which iterates all runs in the file, records the runs iterated
  ## over and injects the current `runNumber` and current `grp` into the calling scope
  let t0 = epochTime()
  # iterate over all raw data groups
  var runNumbersIterated: set[uint16]
  var runNumbersDone: set[uint16]

  # check if run type is stored in group, read it
  var runType = rtNone
  for num, curGrp in runs(h5f, base):
    # now read some data. Return value will be added later
    let grp {.inject.} = curGrp
    let runNumber {.inject.} = num
    runNumbersIterated.incl runNumber.uint16
    body
    runNumbersDone.incl runNumber.uint16

  info "Reconstruction of all runs in $# with flags: $# took $# seconds" % [$h5f.name,
                                                                            $flags,
                                                                            $(epochTime() - t0)]
  info "Performed reconstruction of the following runs:"
  info $runNumbersDone
  info "while iterating over the following:"
  info $runNumbersIterated
  # add all flags that were processed

proc reconstructRunsInFile(h5f: H5File,
                           h5fout: H5File,
                           flags: set[RecoFlagKind],
                           cfgFlags: set[ConfigFlagKind],
                           searchRadius: int,
                           dbscanEpsilon: float,
                           clusterAlgo: ClusteringAlgorithm,
                           runNumberArg: Option[int] = none[int]()) =
  ## proc which performs reconstruction of runs in a given file (all by default). It only takes
  ## care of the general purpose conversion from raw data to reconstructed ingrid clusters
  ## plus FADC data. More complicated calibrations are handled by `applyCalibrationSteps` below.
  ## inputs:
  ##   `h5f`: the file from which we read the raw data
  ##   `h5fout`: the file to which we write the reconstructed data. May be the same file
  ##   `flags`: stores the command line arguments
  ##   `runNumberArg`: optional run number, if given only this run is reconstructed
  ##   `h5fout`: output file to which the reconstructed data is written
  const batchsize = 5000
  var reco_run: seq[RecoEvent[Pix]] = @[]
  let showPlots = if cfShowPlots in cfgFlags: true else: false
  # read the timepix version from the input file
  let timepixVersion = h5f.timepixVersion()
  recordIterRuns(rawDataBase()):
    if (recoBase() & $runNumber) in h5fout:
      # check attributes whether this run was actually finished
      let h5grp = h5fout[(recoBase() & $runNumber).grp_str]
      if "RawTransferFinished" in h5grp.attrs and
         h5grp.attrs["RawTransferFinished", string] == "true":
        continue
      # else we redo it
    # check whether all runs are read, if not if this run is correct run number
    let runType = parseEnum[RunTypeKind](
      h5f[(rawDataBase() & $runNumber).grp_str].attrs["runType", string]
    )
    if (runNumberArg.isSome and runNumber == runNumberArg.get) or
       rfReadAllRuns in flags:
      # intersection of `flags` with all `"only flags"` has to be empty
      doAssert (flags * {rfOnlyEnergy .. rfOnlyGainFit}).card == 0
      # initialize groups in `h5fout`
      initRecoFadcInH5(h5f, h5fout, runNumber, batchsize)
      copyOverDataAttrs(h5f, h5fout, runNumber)

      var runGroupForAttrs = h5f[grp.grp_str]
      let nChips = runGroupForAttrs.attrs["numChips", int]
      # TODO: we can in principle perform energy calibration in one go
      # together with creation of spectrum, if we work as follows:
      # 1. calibration runs:
      #    - need to interface with Python code, i.e. call fitting procedure,
      #      which returns the value to the Nim program as its return value
      let t1 = epochTime()
      for chip, pixdata in h5f.readDataFromH5(runNumber, timepixVersion):
        # given single runs pixel data, call reconstruct run proc
        # NOTE: the data returned from the iterator contains all
        # events in ascending order of event number, i.e.
        # [0] -> eventNumber == 0 and so on
        reco_run.add reconstructSingleChip(pixdata, runNumber, chip,
                                           searchRadius, dbscanEpsilon, clusterAlgo,
                                           timepixVersion)
        info &"Reco run now contains {reco_run.len} elements"
      info "Reconstruction of run $# took $# seconds" % [$runNumber, $(epochTime() - t1)]
      # finished run, so write run to H5 file
      h5fout.writeRecoRunToH5(h5f, reco_run, runNumber, timepixVersion)
      # now flush both files
      h5fout.flush
      h5f.flush
      # set reco run length back to 0
      reco_run.setLen(0)
      # calculate fractions of FADC / Scinti triggers per run
      if fadcRawPath(runNumber) in h5f:
        h5fout.calcTriggerFractions(runNumber)
      # now check whether create iron spectrum flag is set
      # or this is a calibration run, then always create it
      if rfCreateFe in flags or runType == rtCalibration:
        createAndFitFeSpec(h5fout, runNumber, not showPlots)

      # add flag that this run is finished
      let h5grp = h5fout[(recoBase() & $runNumber).grp_str]
      h5grp.attrs["RawTransferFinished"] = "true"

proc applyCalibrationSteps(h5f: H5File,
                           flags: set[RecoFlagKind],
                           cfgFlags: set[ConfigFlagKind],
                           cfgTable: TomlValueRef,
                           runNumberArg = none[int](),
                           calib_factor = none[float]()) =
  ## inputs:
  ##   `h5f`: the file from which we read the raw data
  ##   `h5fout`: the file to which we write the reconstructed data. May be the same file
  ##   `flags`: stores the command line arguments
  ##   `runNumberArg`: optional run number, if given only this run is reconstructed
  ##   `calib_factor`: optional factor to use to calculate energy of clusters based on number
  ##     of pixels in the cluster.
  ##   `h5fout`: optional file to which the reconstructed data is written instead
  ## NOTE: While the flags are called "--only*" the ``only`` portion simply refers to
  ## not converting raw data to reconstructed data. Different flags do not exclude one another!
  let showPlots = if cfShowPlots in cfgFlags: true else: false
  if rfOnlyEnergyElectrons in flags:
    #h5fout.calcEnergyFromPixels(runNumber, calib_factor)
    let interval = cfgTable["Calibration"]["gasGainInterval"].getFloat
    let gcKind = parseEnum[GasGainVsChargeCalibKind](
      cfgTable["Calibration"]["gasGainEnergyKind"].getStr
    )
    h5f.calcEnergyFromCharge(interval, gcKind)
  if rfOnlyGainFit in flags:
    let interval = cfgTable["Calibration"]["gasGainInterval"].getFloat
    let gcKind = parseEnum[GasGainVsChargeCalibKind](
      cfgTable["Calibration"]["gasGainEnergyKind"].getStr
    )
    h5f.performChargeCalibGasGainFit(interval, gcKind)
  recordIterRuns(recoBase()):
    if (runNumberArg.isSome and runNumber == runNumberArg.get) or
       rfReadAllRuns in flags:
      # intersection of `flags` with all `"only flags"` must not be empty
      doAssert (flags * {rfOnlyEnergy .. rfOnlyGainFit}).card > 0
      # only perform energy calibration of the reconstructed runs in file
      # check if reconstructed run exists
      if hasKey(h5f.groups, (recoBase & $runNumber)) == true:
        if rfOnlyEnergy in flags:
          # TODO: take this out
          doAssert calib_factor.isSome, "Need calibration factor to calculate " &
            "energy from number of pixels!"
          h5f.calcEnergyFromPixels(runNumber, calib_factor.get)
        if rfOnlyCharge in flags:
          let toDelete = cfgTable["Calibration"]["deleteChargeDset"].getBool
          h5f.applyChargeCalibration(runNumber, toDelete = toDelete)
        if rfOnlyGasGain in flags:
          let interval = cfgTable["Calibration"]["gasGainInterval"].getFloat
          let minInterval = cfgTable["Calibration"]["minimumGasGainInterval"].getFloat
          h5f.calcGasGain(runNumber, interval, minInterval)
        if rfOnlyFadc in flags:
          h5f.calcRiseAndFallTimes(runNumber)
        if rfOnlyFeSpec in flags:
          createAndFitFeSpec(h5f, runNumber, not showPlots)
      else:
        warn "No reconstructed run found for $#" % $grp

proc parseTomlConfig(configFile: string): (TomlValueRef, set[ConfigFlagKind]) =
  ## parses our config.toml file and returns a set of flags
  ## corresponding to different settings and the full toml table
  # TODO: concat together from `TpxDir`
  var config: TomlValueRef
  if configFile.len == 0:
    const sourceDir = currentSourcePath().parentDir
    config = parseToml.parseFile(sourceDir / "config.toml")
  else:
    config = parseToml.parseFile(configFile)
  var flags: set[ConfigFlagKind]
  if config["Calibration"]["showPlots"].getBool:
    flags.incl cfShowPlots
  result = (config, flags)

proc parseOnlyFlags(args: DocoptTab): set[RecoFlagKind] =
  if $args["--only_charge"] == "true":
    result.incl rfOnlyCharge
  if $args["--only_fadc"] == "true":
    result.incl rfOnlyFadc
  if $args["--only_fe_spec"] == "true":
    result.incl rfOnlyFeSpec
  if $args["--only_gas_gain"] == "true":
    result.incl rfOnlyGasGain
  if $args["--only_gain_fit"] == "true":
    result.incl rfOnlyGainFit
  if $args["--only_energy_from_e"] == "true":
    result.incl rfOnlyEnergyElectrons

proc flagsValid(h5f: H5File, flags: set[RecoFlagKind]): bool =
  ## Checks whether the flags are actually valid for the given file
  let grp = h5f[recoGroupGrpStr()]
  result = true
  if "runType" in grp.attrs:
    let runType = parseEnum[RunTypeKind](grp.attrs["runType", string])
    if rfOnlyGainFit in flags and runType != rtCalibration:
      warn "Fit to charge calibration / gas gain only possible for " &
        "calibration runs!"
      return false

proc main() =

  # create command line arguments using docopt
  let args = docopt(doc)
  echo args

  let
    h5f_name = $args["<HDF5file>"]
    create_fe_arg = $args["--create_fe_spec"]
    configFile = if $args["--config"] != "nil": $args["--config"]
                 else: ""
  var
    flags: set[RecoFlagKind]
    runNumberArg = $args["--runNumber"]
    runNumber: Option[int]
    calibFactorStr = $args["--only_energy"]
    calibFactor: Option[float]
    outfile = $args["--out"]
  if runNumberArg == "nil":
    flags.incl rfReadAllRuns
  else:
    runNumber = some(parseInt(runNumberArg))
  if outfile != "nil":
    outfile = $args["--out"]
    info &"Set outfile to {outfile}"
  if create_fe_arg == "true":
    flags.incl rfCreateFe
  if calibFactorStr != "nil":
    flags.incl rfOnlyEnergy
    calibFactor = some(parseFloat(calibFactorStr))

  # parse the explicit `--only_...` flags
  flags = flags + parseOnlyFlags(args)

  # parse config toml file
  let (cfgTable, cfgFlags) = parseTomlConfig(configFile)

  let plotOutPath = cfgTable["Calibration"]["plotDirectory"].getStr
  var h5f = H5open(h5f_name, "rw")
  let plotDirPrefix = h5f.genPlotDirname(plotOutPath, PlotDirPrefixAttr)
  let searchRadius = cfgTable["Reconstruction"]["searchRadius"].getInt
  let dbscanEpsilon = cfgTable["Reconstruction"]["epsilon"].getFloat
  let clusterAlgo = parseEnum[ClusteringAlgorithm](cfgTable["Reconstruction"]["clusterAlgo"].getStr)

  if (flags * {rfOnlyEnergy .. rfOnlyGainFit}).card == 0:
    # `reconstruction` call w/o `--only-*` flag
    # visit the whole file to read which groups exist
    h5f.visitFile
    var h5fout: H5File
    if outfile != "None":
      h5fout = H5open(outfile, "rw")
      h5fout.visitFile
      # copy over `PlotDirPrefixAttr` to h5fout
      h5fout.attrs[PlotDirPrefixAttr] = plotDirPrefix
      reconstructRunsInFile(h5f, h5fout, flags, cfgFlags,
                            searchRadius = searchRadius,
                            dbscanEpsilon = dbscanEpsilon,
                            clusterAlgo = clusterAlgo,
                            runNumberArg = runNumber)

    var err = h5f.close()
    if err != 0:
      logging.error &"Failed to close H5 file {h5f.name}"
    err = h5fout.close()
    if err != 0:
      logging.error &"Failed to close H5 file {h5fout.name}"
  else:
    var h5f = H5open(h5f_name, "rw")
    if flagsValid(h5f, flags):
      applyCalibrationSteps(h5f, flags, cfgFlags, cfgTable,
                            runNumberArg = runNumber,
                            calibFactor = calibFactor)
    else:
      logging.warn &"Invalid flags given for file {h5f_name}: {flags}"


when isMainModule:
  main()
