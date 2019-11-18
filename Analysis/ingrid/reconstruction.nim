# this file contains the necessary code for reconstruction of X-rays in the raw data
# consists of
#   - finding clusters in data
#   - calculating properties of Events

# nim stdlib
import os
import sequtils, future
#import threadpool
import threadpool_simple
import logging
import math
import tables
import docopt
import strutils
import strformat
import typeinfo
import re
import times
import stats
import sets
import macros

# external modules
import nlopt
import seqmath
import nimhdf5
import parsetoml

# custom modules
import tos_helpers
import helpers/utils
import ingrid_types
import calibration
import fadc_analysis

type
  # fit object, which is handed to the NLopt library in the
  # `VarStruct` -> to the eccentricity function
  FitObject[T: SomePix] = object
    cluster: Cluster[T]
    xy: tuple[x, y: float64]

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

Usage:
  reconstruction <HDF5file> --out <name> [options]
  reconstruction <HDF5file> --out <name> [--runNumber <number>] [--create_fe_spec] [options]
  reconstruction <HDF5file> [--runNumber <number>] --only_energy <factor> [options]
  reconstruction <HDF5file> [--runNumber <number>] --only_energy_from_e [options]
  reconstruction <HDF5file> [--runNumber <number>] --only_fe_spec [options]
  reconstruction <HDF5file> [--runNumber <number>] --only_charge [options]
  reconstruction <HDF5file> [--runNumber <number>] --only_fadc [options]
  reconstruction <HDF5file> [--runNumber <number>] --only_gas_gain [options]
  reconstruction <HDF5file> [--runNumber <number>] --only_gain_fit [options]
  reconstruction -h | --help
  reconstruction --version


Options:
  --out <name>            Filename and path of output file
  --runNumber <number>    Only work on this run
  --create_fe_spec        Toggle to create Fe calibration spectrum based on cuts
                          Takes precedence over --calib_energy if set!
  --only_energy <factor>  Toggle to /only/ perform energy calibration using the given factor.
                          Takes precedence over --create_fe_spec and --calib_energy if set.
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
                            h5f: var H5FileObj,
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

proc writeRecoRunToH5*[T: SomePix](h5f: var H5FileObj,
                                   h5fraw: var H5FileObj,
                                   reco_run: seq[FlowVar[ref RecoEvent[T]]],
                                   runNumber: int) =
  ## proc which writes the reconstructed event data from a single run into
  ## the given H5 file. Called after every processed run
  ## inputs:
  ##     `h5f`: var H5FileObj = the H5 file object in which to store data
  ##     `h5fraw`: var H5FileObj = the H5 file objcect containing the raw data,
  ##       may be the same file as `h5f`.
  ##     `reco_run`: seq[FlowVar[ref RecoEvent]] = sequence of reconstructed events
  ##       to write to file. FlowVar[ref Event], due to using threadpool / spawn
  ##       to reconstruct the runs in parallel
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

  for chip in 0 ..< nChips:
    x[chip]  = @[]
    y[chip]  = @[]
    ch[chip] = @[]

  for event_f in reco_run:
    let
      # get the RecoEvent from the FlowVar ref
      event = (^event_f)[]
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

      expandMacros:
        int_data_tab.setTabFields(int_cluster_names, chip, cl)
        float_data_tab.setTabFields(float_cluster_names, chip, cl)
        float_data_tab.setTabFields(float_geometry_names, chip, cl.geometry)

      # add event number individually, since it's not part of some object we can
      # use our macro for
      int_data_tab["eventNumber"][chip].add num

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

    # anyways, write the chip dataset attributes
    let raw_chip_group_name = raw_groups & $chip
    var raw_group = h5fraw[raw_chip_group_name.grp_str]
    let ch_numb = raw_group.attrs["chipNumber", int]
    let ch_name = raw_group.attrs["chipName", string]
    # and write these to the current group
    chip_groups[ch_numb].attrs["chipNumber"] = ch_numb
    chip_groups[ch_numb].attrs["chipName"] = ch_name

iterator readDataFromH5*(h5f: var H5FileObj, runNumber: int):
         (int, seq[(Pixels, int)]) =
  ## proc to read data from the HDF5 file from `group`
  ## returns the chip number and a sequence containing the pixel data for this
  ## event and its event number
  var chip_base = rawDataChipBase(runNumber)
  let raw_data_basename = rawDataBase()
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
      var
        raw_x_dset  = h5f[(grp / "raw_x").dset_str]
        raw_y_dset  = h5f[(grp / "raw_y").dset_str]
        raw_ch_dset = h5f[(grp / "raw_ch").dset_str]
      let
        raw_x  = raw_x_dset[vlen_xy, uint8]
        raw_y  = raw_y_dset[vlen_xy, uint8]
        raw_ch = raw_ch_dset[vlen_ch, uint16]

      var runPix = newSeq[(Pixels, int)](raw_x.len)
      for i in 0 ..< raw_x.len:
        let rpix = zipEm(raw_x[i], raw_y[i], raw_ch[i])
        runPix[i] = (rpix, evNumbers[i])
      # and yield them
      yield (chip_number, run_pix)

proc newClusterGeometry(): ClusterGeometry =
  result = ClusterGeometry(rmsLongitudinal: Inf,
                           rmsTransverse: Inf,
                           eccentricity: Inf,
                           rotationAngle: Inf,
                           skewnessLongitudinal: Inf,
                           skewnessTransverse: Inf,
                           kurtosisLongitudinal:Inf,
                           kurtosisTransverse: Inf,
                           length: Inf,
                           width: Inf,
                           fractionInTransverseRms: Inf)


proc newClusterObject[T: SomePix](c: Cluster[T]): ClusterObject[T] =
  ## initialize variables with Inf for now
  # TODO: should we initialize geometry values by Inf as well?
  let geometry = ClusterGeometry()
  result = ClusterObject[T](data: c,
                            centerX: Inf,
                            centerY: Inf,
                            energy: Inf,
                            geometry: geometry)

proc eccentricity[T: SomePix](p: seq[float], func_data: FitObject[T]): float =
  ## this function calculates the eccentricity of a found pixel cluster using nimnlopt.
  ## Since no proper high level library is yet available, we need to pass a var pointer
  ## of func_data, which contains the x and y arrays in which the data is stored, in
  ## order to calculate the RMS variables
  # first recover the data from the pointer to func_data, by casting the
  # raw pointer to a Cluster object
  let fit = func_data
  let c = fit.cluster
  let (centerX, centerY) = fit.xy

  var
    sum_x: float = 0
    sum_y: float = 0
    sum_x2: float = 0
    sum_y2: float = 0

  for i in 0..<len(c):
    let
      new_x = cos(p[0]) * (c[i].x.float - centerX) * PITCH - sin(p[0]) * (c[i].y.float - centerY) * PITCH
      new_y = sin(p[0]) * (c[i].x.float - centerX) * PITCH + cos(p[0]) * (c[i].y.float - centerY) * PITCH
    sum_x += new_x
    sum_y += new_y
    sum_x2 += (new_x * new_x)
    sum_y2 += (new_y * new_y)

  let
    n_elements: float = len(c).float
    rms_x: float = sqrt( (sum_x2 / n_elements) - (sum_x * sum_x / n_elements / n_elements))
    rms_y: float = sqrt( (sum_y2 / n_elements) - (sum_y * sum_y / n_elements / n_elements))

  # calc eccentricity from RMS
  let exc = rms_x / rms_y
  result = -exc

#proc eccentricityInt(p: seq[float], func_data: FitObjectInt): float =
#  eccentricityImpl()

proc calcGeometry[T: SomePix](cluster: Cluster[T],
                              pos_x, pos_y, rot_angle: float): ClusterGeometry =
  ## given a cluster and the rotation angle of it, calculate the different
  ## statistical moments, i.e. RMS, skewness and kurtosis in longitudinal and
  ## transverse direction
  ## done by rotating the cluster by the angle to define the two directions
  let npix = cluster.len

  var
    xRot = newSeq[float](npix)
    yRot = newSeq[float](npix)
    radius: float
    x_max, x_min: float
    y_max, y_min: float
    i = 0
  for p in cluster:
    let (x, y) = applyPitchConversion(p.x, p.y)
    xRot[i] = cos(-rot_angle) * (x - pos_x) - sin(-rot_angle) * (y - pos_y)
    yRot[i] = sin(-rot_angle) * (x - pos_x) + cos(-rot_angle) * (y - pos_y)

    # calculate distance from center
    let dist = distance(xRot[i], yRot[i])
    if dist > radius:
      radius = dist
    inc i
  # define statistics objects
  var
    stat_x: RunningStat
    stat_y: RunningStat
  # and push our new vectors to them
  stat_x.push(xRot)
  stat_y.push(yRot)

  # now we have all data to calculate the geometric properties
  result.length               = max(xRot) - min(xRot)
  result.width                = max(yRot) - min(yRot)
  result.rmsTransverse        = stat_y.standardDeviation()
  result.rmsLongitudinal      = stat_x.standardDeviation()
  result.skewnessTransverse   = stat_y.skewness()
  result.skewnessLongitudinal = stat_x.skewness()
  result.kurtosisTransverse   = stat_y.kurtosis()
  result.kurtosisLongitudinal = stat_x.kurtosis()
  result.rotationAngle        = rot_angle
  result.eccentricity         = result.rmsLongitudinal / result.rmsTransverse
  # get fraction of all pixels within the transverse RMS, by filtering all elements
  # within the transverse RMS radius and dividing by total pix
  # when not defined(release):
  #   # DEBUG
  #   echo "rms trans is ", result.rmsTransverse
  #   echo "std is ", stat_y.variance()
  #   echo "thus filter is ", filterIt(zip(xRot, yRot), distance(it.a, it.b) <= result.rmsTransverse)
  result.lengthDivRmsTrans = result.length / result.rmsTransverse
  result.fractionInTransverseRms = (
    filterIt(zip(xRot, yRot),
             distance(it[0], it[1]) <= result.rmsTransverse).len
  ).float / float(npix)

proc isPixInSearchRadius[T: SomeInteger](p1, p2: Coord[T], search_r: int): bool =
  ## given two pixels, p1 and p2, we check whether p2 is within one square search
  ## of p1
  ## inputs:
  ##   p1: Pix = pixel from which to start search
  ##   p2: Pix = pixel for which to check
  ##   search_r: int = search radius (square) in which to check for p2 in (p1 V search_r)
  ## outpuits:
  ##   bool = true if within search_r
  ##          false if not
  let
    # determine boundary of search space
    right = p1.x.int + search_r
    left  = p1.x.int - search_r
    up    = p1.y.int + search_r
    down  = p1.y.int - search_r
  # NOTE: for performance we may want to use the fact that we may know that
  # p1 is either to the left (if in the same row) or below (if not in same row)
  var
    in_x: bool = false
    in_y: bool = false

  if p2.x.int < right and p2.x.int > left:
    in_x = true
  if p2.y.int < up and p2.y.int > down:
    in_y = true
  result = if in_x == true and in_y == true: true else: false

proc findSimpleCluster*[T: SomePix](pixels: seq[T]): seq[Cluster[T]] =
  ## this procedure searches for clusters based on a fixed search radius, whether
  ## a pixel is found within that boundary, e.g. searchRadius = 50:
  ## if within 50 pixels another pixel is found, add pixel to cluster, continue
  ## inputs:
  ##   -
  ## ouputs:
  ##   -

  # - iterate over all hit pixels in event
  # - check next pixel, is it within search bound? yes,
  #   add pixel to hit
  var
    # sequence in which to store the hits, which are part of a cluster
    c: Cluster[T] = @[]
    # create copy of pixels so that we can remove elements from it
    raw_event = pixels
    # counter
    i = 0
  result = @[] #new seq[Cluster]

  when defined(onlySingleCluster):
    let
      search_r = 128
      cutoff_size = 1
  else:
    let
      search_r = 50
      cutoff_size = 5

  # add the first pixel of the given sequence to have a starting pixel, from which we
  # look for other pixels in the cluster
  c.add(pixels[0])
  raw_event.deleteIntersection(@[pixels[0]])
  while raw_event.len > 0 and i < c.len:
    let p1: Coord[type(c[i].x)] = (x: c[i].x, y: c[i].y)
    # alternatively:
    let t = raw_event.filterIt(isPixInSearchRadius(p1, (it.x, it.y), search_r))
    #let t = filter(raw_event, (p: tuple[x, y: uint8, ch: uint16]) ->
    #               bool => isPixInSearchRadius(p1, (p.x, p.y), search_r))

    # add all found pixels to current cluster
    c = concat(c, t)
    # remove elements from t in raw_event
    deleteIntersection(raw_event, t)

    if i == c.len - 1 and raw_event.len > 0:
      # if we are at the last hit pixel in the event, but raw_events is not empty,
      # this means there is more than 1 cluster in the frame. Thus, add the current
      # cluster 'c' to the seq of clusters and call this function recursively with
      # raw_events as the starting parameter
      if len(c) > cutoff_size:
        result.add(c)
      if len(raw_event) > cutoff_size:
        result.add(findSimpleCluster(raw_event))
    elif raw_event.len == 0 and len(c) > cutoff_size:
      result.add(c)
    inc i

proc eccentricityNloptOptimizer[T: SomePix](fitObject: FitObject[T]):
  NloptOpt[FitObject[T]] =
  ## returns the already configured Nlopt optimizer to fit the rotation angle /
  ## eccentricity
  var
    # set the boundary values corresponding to range of 360 deg
    lb = (-4.0 * arctan(1.0), 4.0 * arctan(1.0))
  type tFitObj = type(fitObject)
  result = newNloptOpt[tFitObj](LN_BOBYQA, 1, @[lb])
  # hand the function to fit as well as the data object we need in it
  # NOTE: workaround for https://github.com/nim-lang/Nim/issues/11778
  var varStruct = VarStruct[tFitObj](userFunc: eccentricity, data: fitObject,
                                     kind: FuncKind.NoGrad)
  result.setFunction(varStruct)
  # set relative precisions of x and y, as well as limit max time the algorithm
  # should take to 1 second
  # these default values have proven to be working
  result.xtol_rel = 1e-8
  result.ftol_rel = 1e-8
  result.maxtime  = 1.0
  result.initial_step = 0.02

proc fitRotAngle[T: SomePix](cl_obj: ClusterObject[T],
                             rotAngleEstimate: float): (float, float) =
  ## Performs the fitting of the rotation angle on the given ClusterObject
  ## `cl_obj` and returns the final parameters as well as the minimum
  ## value at those parameters.
  # set the fit object with which we hand the necessary data to the
  # eccentricity function
  var fit_object = FitObject[T](cluster: cl_obj.data,
                                xy: (x: cl_obj.centerX, y: cl_obj.centerY))
  var opt = eccentricityNloptOptimizer(fit_object)
  # start minimization
  var p = @[rotAngleEstimate]
  let (params, min_val) = opt.optimize(p)
  if opt.status < NLOPT_SUCCESS:
    info opt.status
    warn "nlopt failed!"
  # clean up optimizer
  destroy(opt)
  # now return the optimized parameters and the corresponding min value
  result = (params[0], min_val)

proc recoCluster[T: SomePix](c: Cluster[T]): ClusterObject[T] {.gcsafe.} =
  result = newClusterObject(c)
  let
    clustersize: int = len(c)
    (sum_x, sum_y, sumTotInCluster) = sum(c)
    # using map and sum[Pix] we can calculate sum of x^2, y^2 and x*y in one line
    (sum_x2, sum_y2, sum_xy) = sum(c.mapIt((it.x.int * it.x.int,
                                            it.y.int * it.y.int,
                                            it.x.int * it.y.int)))
    #(sum_x2, sum_y2, sum_xy) = sum(map(c, (p: Pix) ->
    #                               (int, int, int) => (p.x.int * p.x.int,
    #                                                   p.y.int * p.y.int,
    #                                                   p.x.int * p.y.int)))
    pos_x = float64(sum_x) / float64(clustersize)
    pos_y = float64(sum_y) / float64(clustersize)
  var
    rms_x = sqrt(float64(sum_x2) / float64(clustersize) - pos_x * pos_x)
    rms_y = sqrt(float64(sum_y2) / float64(clustersize) - pos_y * pos_y)
    rotAngleEstimate = arctan( (float64(sum_xy) / float64(clustersize)) -
                               pos_x * pos_y / (rms_x * rms_x))

  # set the total "charge" in the cluster (sum of ToT values), can be
  # converted to electrons with ToT calibration
  result.sum_tot = sumTotInCluster.int
  # set number of hits in cluster
  result.hits = clustersize
  # set the position
  (result.centerX, result.centerY) = applyPitchConversion(pos_x, pos_y)
  #(float(NPIX) - float(pos_x) + 0.5) * PITCH
  #result.pos_y = (float(pos_y) + 0.5) * PITCH
  # prepare rot angle fit
  if rotAngleEstimate < 0:
    #echo "correcting 1"
    rotAngleEstimate += 8 * arctan(1.0)
  if rotAngleEstimate > 4 * arctan(1.0):
    #echo "correcting 2"
    rotAngleEstimate -= 4 * arctan(1.0)
  elif classify(rotAngleEstimate) != fcNormal:
    warn "Rot angle estimate is NaN, vals are ", $rms_x, " ", $rms_y
    # what do we do in this case with the geometry?!
    #raise newException(ValueError, "Rotation angle estimate returned bad value")
    warn "Fit will probably fail!"

  # else we can minimize the rotation angle and calc the eccentricity
  let (rot_angle, eccentricity) = fitRotAngle(result, rotAngleEstimate)

  # now we still need to use the rotation angle to calculate the different geometric
  # properties, i.e. RMS, skewness and kurtosis along the long axis of the cluster
  result.geometry = calcGeometry(c, result.centerX, result.centerY, rot_angle)

proc recoEvent*[T: SomePix](data: (seq[T], int), chip: int): ref RecoEvent[T] {.gcsafe.} =
  result = new RecoEvent[T]
  result.event_number = data[1]
  result.chip_number = chip
  if data[0].len > 0:
    let cluster = findSimpleCluster(data[0])
    result.cluster = newSeq[ClusterObject[T]](cluster.len)
    for i, cl in cluster:
      result.cluster[i] = recoCluster(cl)

{.experimental.}
#proc reconstructSingleChip(data: seq[Pixels], run, chip: int): seq[ref RecoEvent] =
proc reconstructSingleChip*(data: seq[(Pixels, int)], run, chip: int): seq[FlowVar[ref RecoEvent[Pix]]] =
                           #evNumbers: seq[int]): seq[FlowVar[ref RecoEvent]] =
  ## procedure which receives pixel data for a given chip and run number
  ## and performs the reconstruction on it
  ## inputs:
  ##    data: seq[Pixels] = data of pixels for this chip containing pixels for each event
  ##    run: int = run number of run
  ##    chip: int = chip number working on
  info &"Working on chip {chip} in run {run}"
  info &"We have {data.len} events to reconstruct"
  var count = 0
  result = newSeq[FlowVar[ref RecoEvent[Pix]]](data.len)
  #parallel:
  let p = newThreadPool()
  for event in 0..data.high:
    if event < result.len:
      result[event] = p.spawn recoEvent(data[event], chip)
    echoFilesCounted(count, 2500)
  p.sync()

#iterator matchingGroup(h5f: var H5FileObj,

proc createAndFitFeSpec(h5f: var H5FileObj,
                        h5fout: var H5FileObj,
                        runNumber: int,
                        fittingOnly: bool) =
  ## create the Fe spectrum for the run, apply the charge calibration if possible
  ## and then fit to the Fe spectrum, writing results to `h5fout`
  var centerChip = h5f.getCenterChip(runNumber)
  h5fout.createFeSpectrum(runNumber, centerChip)
  try:
    h5fout.applyChargeCalibration(runNumber)
  except KeyError as e:
    warn "No charge calibration possible for current one or " &
         "more chips. Exception message:\n" & e.msg
  h5fout.fitToFeSpectrum(runNumber, centerChip, fittingOnly)

proc initRecoFadcInH5(h5f, h5fout: var H5FileObj, runNumber, batchsize: int) =
  # proc to initialize the datasets etc in the HDF5 file for the FADC. Useful
  # since we don't want to do this every time we call the write function
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

proc copyOverDataAttrs(h5f, h5fout: var H5FileObj, runNumber: int) =
  let recoGrp = h5fout[(recoBase() & $runNumber).grp_str]
  let rawGrp = h5f[(rawDataBase() & $runNumber).grp_str]
  # copy attributes of `runs/run_XYZ -> reconstruction/run_XYZ`
  recoGrp.copy_attributes(rawGrp.attrs)
  template copyOver(path, dtype: untyped): untyped =
    # don't care for the returned group here
    discard h5fout.write_dataset(recoBase() & path, h5f[rawDataBase() & path, dtype])

  for dset in items(rawGrp, start_path = rawDataBase()):
    echo "Copying: ", dset
    case dset.dtypeAnyKind
    of akInt64:
      copyOver($runNumber / dset.name.extractFilename, int64)
    of akFloat64:
      copyOver($runNumber / dset.name.extractFilename, float64)
    else:
      doAssert false, "Unexpected dataset " & $dset & " with base type " & $dset.dtypeAnyKind

proc calcTriggerFractions(h5f: var H5FileObj, runNumber: int) =
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

proc reconstructRunsInFile(h5f: var H5FileObj,
                           h5fout: var H5FileObj,
                           flags: set[RecoFlagKind],
                           cfgFlags: set[ConfigFlagKind],
                           runNumberArg: int = -1,
                           calib_factor: float = 1.0) =
  ## proc which performs reconstruction of runs in a given file (all by default)
  ## if the --only_energy command line argument is set, we skip the reconstruction
  ## and only perform an energy calibration of the existing (!) reconstructed
  ## runs using the calibration factor given
  ## inputs:
  ##   `h5f`: the file from which we read the raw data
  ##   `h5fout`: the file to which we write the reconstructed data. May be the same file
  ##   `flags`: stores the command line arguments
  ##   `runNumberArg`: optional run number, if given only this run is reconstructed
  ##   `calib_factor`: factor to use to calculate energy of clusters
  ##   `h5fout`: optional file to which the reconstructed data is written instead

  let
    raw_data_basename = rawDataBase()
    t0 = epochTime()

  const batchsize = 5000
  var reco_run: seq[FlowVar[ref RecoEvent[Pix]]] = @[]
  let showPlots = if cfShowPlots in cfgFlags: true else: false

  # iterate over all raw data groups
  var runNumbersIterated: set[uint16]
  var runNumbersDone: set[uint16]

  # check if run type is stored in group, read it
  var runType = rtNone
  var rawGroup = h5f[rawGroupGrpStr]
  if "runType" in rawGroup.attrs:
    runType = parseEnum[RunTypeKind](rawGroup.attrs["runType", string])
    if rfOnlyGainFit in flags and runType != rtCalibration:
      warn "Fit to charge calibration / gas gain only possible for " &
        "calibration runs!"
      return
    elif rfOnlyGainFit in flags:
      h5f.performChargeCalibGasGainFit()
      # return early
      # TODO: move this whole stuff somewhere else! Does not belong into
      # reconstruction like this!
      return

  for num, grp in runs(h5f, rawDataBase()):
    # now read some data. Return value will be added later
    let runNumber = parseInt(num)
    runNumbersIterated.incl runNumber.uint16
    # check whether all runs are read, if not if this run is correct run number
    if rfReadAllRuns in flags or runNumber == runNumberArg:
      # check if intersection of `flags` with all `"only flags"` is empty
      if (flags * {rfOnlyEnergy .. rfOnlyGainFit}).card == 0:
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
        for chip, pixdata in h5f.readDataFromH5(runNumber):
          # given single runs pixel data, call reconstruct run proc
          # NOTE: the data returned from the iterator contains all
          # events in ascending order of event number, i.e.
          # [0] -> eventNumber == 0 and so on
          reco_run.add reconstructSingleChip(pixdata, runNumber, chip)
          info &"Reco run now contains {reco_run.len} elements"

        info "Reconstruction of run $# took $# seconds" % [$runNumber, $(epochTime() - t1)]
        # finished run, so write run to H5 file
        h5fout.writeRecoRunToH5(h5f, reco_run, runNumber)

        # now flush both files
        h5fout.flush
        h5f.flush
        # set reco run length back to 0
        reco_run.setLen(0)

        # calculate fractions of FADC / Scinti triggers per run
        h5fout.calcTriggerFractions(runNumber)

        runNumbersDone.incl runNumber.uint16

        # now check whether create iron spectrum flag is set
        # or this is a calibration run, then always create it
        if rfCreateFe in flags or runType == rtCalibration:
          createAndFitFeSpec(h5f, h5fout, runNumber, not showPlots)
      else:
        # only perform energy calibration of the reconstructed runs in file
        # check if reconstructed run exists
        if hasKey(h5fout.groups, (recoBase & $runNumber)) == true:
          if rfOnlyEnergy in flags:
            # TODO: take this out
            h5fout.calcEnergyFromPixels(runNumber, calib_factor)
          if rfOnlyCharge in flags:
            h5fout.applyChargeCalibration(runNumber)
          if rfOnlyGasGain in flags:
            h5fout.calcGasGain(runNumber, showPlots)
          if rfOnlyFadc in flags:
            h5fout.calcRiseAndFallTimes(runNumber)
          if rfOnlyFeSpec in flags:
            createAndFitFeSpec(h5f, h5fout, runNumber, not showPlots)
        else:
          warn "No reconstructed run found for $#" % $grp

  if rfOnlyEnergyElectrons in flags:
    #h5fout.calcEnergyFromPixels(runNumber, calib_factor)
    h5fout.calcEnergyFromCharge()

  info "Reconstruction of all runs in $# took $# seconds" % [$h5f.name, $(epochTime() - t0)]

  info "Performed reconstruction of the following runs:"
  info $runNumbersDone
  info "while iterating over the following:"
  info $runNumbersIterated

proc reconstructRunsInFile(h5f: var H5FileObj,
                           flags: set[RecoFlagKind],
                           cfgFlags: set[ConfigFlagKind],
                           runNumberArg: int = -1,
                           calib_factor: float = 1.0) =
  ## this proc is a wrapper around the one above, which is called in case
  ## no output H5 file is given. In that case we simply hand the input file
  ## as the output file to the proc
  # in case we want to write the reconstructed data to the same file,
  # simply set h5fout to h5f. This creates a copy of h5f, but we don't care
  # since it points to the same file and the important group / dset tables
  # are stored as references anyways
  reconstructRunsInFile(h5f, h5f, flags, cfgFlags, runNumberArg, calib_factor)

proc reconstructSingleRunFolder(folder: string)
    {.deprecated: "This proc is deprecated! Run raw_data_manipulation before!".} =
  ## procedure which receives path to a run folder and reconstructs the objects
  ## in that folder
  ## inputs:
  ##    folder: string = the run folder from which to reconstruct events

  # TODO: Either remove this proc soon or revive it, but with proper run support
  # and calling the actual reconstruction procs from the data!
  # Latter might be a good idea actually for quick convenience!

  info "Starting to read list of files"
  let
    files = getSortedListOfFiles(folder, EventSortType.inode, EventType.InGridType,
                                 rfUnknown)
    f_to_read = if files.high < 30000: files.high else: 30000
    data = readListOfInGridFiles(files[0..f_to_read], rfUnknown)
  var
    min_val = 10.0
    min_seq = newSeq[seq[float64]](7)
  apply(min_seq, (x: seq[float64]) -> seq[float64] => newSeqOfCap[float64](files.high))
  for e in data:
    let
      a: Event = (^e)[]
      chips = a.chips
    for c in chips:
      let
        num: int = c.chip.number
        cluster = findSimpleCluster(c.pixels)
      for cl in cluster:
        #echo "Starting reco of ", a.evHeader["eventNumber"]
        let ob = recoCluster(cl)
        if ob.geometry.rotationAngle < min_val:
          min_val = ob.geometry.rotationAngle
        min_seq[num].add ob.geometry.rotationAngle
  dumpRotAngle(min_seq)

proc parseTomlConfig(): set[ConfigFlagKind] =
  ## parses our config.toml file and (for now) just returns a set of flags
  ## corresponding to different settings
  # TODO: concat together from `TpxDir`
  const sourceDir = currentSourcePath().parentDir
  let config = parseToml.parseFile(sourceDir / "config.toml")
  let showPlot = config["Calibration"]["showPlots"].getBool
  if showPlot:
    result.incl cfShowPlots

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

proc main() =

  # create command line arguments using docopt
  let args = docopt(doc)
  echo args

  let
    h5f_name = $args["<HDF5file>"]
    create_fe_arg = $args["--create_fe_spec"]
    calib_energy_arg = $args["--calib_energy"]

  var
    flags: set[RecoFlagKind]
    runNumber = $args["--runNumber"]
    outfile = $args["--out"]
    calib_factor_str = $args["--only_energy"]
    calib_factor: float = Inf
  if runNumber == "nil":
    flags.incl rfReadAllRuns
  if outfile == "nil":
    outfile = "None"
  else:
    outfile = $args["--out"]
    info &"Set outfile to {outfile}"
    warn "WARNING: writing to a different H5 file is not quite finished yet."
    warn "\t The resulting file will miss the attributes of the reco run groups"
    warn "\t as well as the common datasets like timestamps!"
  if create_fe_arg == "true":
    flags.incl rfCreateFe
  if calib_factor_str != "nil":
    flags.incl rfOnlyEnergy
    calib_factor = parseFloat(calib_factor_str)

  # parse the explicit `--only_...` flags
  flags = flags + parseOnlyFlags(args)

  # parse config toml file
  let cfgFlags = parseTomlConfig()

  var h5f = H5file(h5f_name, "r")
  # visit the whole file to read which groups exist
  h5f.visitFile
  var h5fout: H5FileObj
  if outfile != "None":
    h5fout = H5file(outfile, "rw")
    h5fout.visitFile

  let raw_data_basename = rawDataBase()
  if rfReadAllRuns in flags and outfile == "None":
    reconstructRunsInFile(h5f, flags, cfgFlags, calib_factor = calib_factor)
  elif rfReadAllRuns notin flags and outfile == "None":
    reconstructRunsInFile(h5f, flags, cfgFlags, parseInt(runNumber), calib_factor)
  elif rfReadAllRuns in flags and outfile != "None":
    reconstructRunsInFile(h5f, h5fout, flags, cfgFlags, calib_factor = calib_factor)
  elif rfReadAllRuns notin flags and outfile != "None":
    reconstructRunsInFile(h5f, h5fout, flags, cfgFlags, parseInt(runNumber), calib_factor = calib_factor)

  var err: herr_t
  err = h5f.close()
  if err != 0:
    logging.error &"Failed to close H5 file {h5f.name}"
  if outfile != "None":
    err = h5fout.close()
    if err != 0:
      logging.error &"Failed to close H5 file {h5fout.name}"

  # NOTE: there's no point for this anymore, at least not at the moment
  # # first check whether given folder is valid run folder
  # let (is_run_folder, contains_run_folder) = isTosRunFolder(folder)
  # echo "Is run folder       : ", is_run_folder
  # echo "Contains run folder : ", contains_run_folder
  # reconstructSingleRun(folder)

when isMainModule:
  main()
