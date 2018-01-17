# this file contains the necessary code for reconstruction of X-rays in the raw data
# consists of
#   - finding clusters in data
#   - calculating properties of Events

# nim stdlib
import os
import sequtils, future
import threadpool
import math
import tables
import docopt
import strutils
import typetraits
import re
import times
import stats
import sets

# external modules
import nlopt
import nlopt_wrapper
import nimhdf5

# custom modules
import tos_helper_functions
import helper_functions
import ingrid_types
import energy_calibration
import fadc_analysis

type
  FitObject = object
    cluster: Cluster
    xy: tuple[x, y: float64]


let doc = """
InGrid reconstruction and energy calibration.

Usage:
  reconstruction <HDF5file> [options]
  reconstruction <HDF5file> --run_number <number> [options]
  reconstruction <HDF5file> --run_number <number> --only_energy <factor> [options]
  reconstruction <HDF5file> --run_number <number> --only_fadc [options]
  reconstruction <HDF5file> (--create_fe_spec | --calib_energy) [options]
  reconstruction <HDF5file> --only_energy <factor> [options]
  reconstruction <HDF5file> --only_fadc [options]
  reconstruction <HDF5file> --out <name> [options]
  reconstruction -h | --help
  reconstruction --version


Options:
  --run_number <number>   Only work on this run
  --create_fe_spec        Toggle to create Fe calibration spectrum based on cuts
                          Takes precedence over --calib_energy if set!
  --calib_energy          Toggle to perform energy calibration. Conversion factors needed!
  --only_energy <factor>  Toggle to /only/ perform energy calibration using the given factor. 
                          Takes precedence over --create_fe_spec and --calib_energy if set.
                          If no run_number is given, performs energy calibration on all runs
                          in the HDF5 file.
  --only_fadc             If this flag is set, the reconstructed FADC data is used to calculate
                          FADC values such as rise and fall times among others, which are written
                          to the H5 file.
  --out <name>            Filename and path of output file
  -h --help               Show this help
  --version               Show version.
"""

template benchmark(num: int, actions: untyped) {.dirty.} =
  for i in 0 ..< num:
    actions

proc writeRecoRunToH5(h5f: var H5FileObj, reco_run: seq[FlowVar[ref RecoEvent]], run_number: int) =
  ## proc which writes the reconstructed event data from a single run into
  ## the given H5 file. Called after every processed run
  ## inputs:
  ##     h5f: var H5FileObj = the H5 file object in which to store data
  ##     reco_run: seq[FlowVar[ref RecoEvent]] = sequence of reconstructed events
  ##       to write to file. FlowVar[ref Event], due to using threadpool / spawn
  ##       to reconstruct the runs in parallel
  ## outputs: -
  ## throws:
  ##     potentially throws HDF5LibraryError, if a call to the H5 library fails

  # TODO!!!!!!! holy fuck, this is the worst piece of code in a long time
  # rewrite this!!!!!

  
  # for now hardcode the number of chips. easy to change by getting the number
  # simply from a file or whatever
  const nchips = 7
  echo "Number of events in total ", reco_run.len
  # now we need to write into reco group for each chip
  # as well as combined
  # that is it, right?
  let t0 = epochTime()
  # first write the raw data
  let
    # group name for reconstructed data
    reco_group_name = getRecoNameForRun(run_number)
    ev_type = special_type(int)
    chip_group_name = reco_group_name / "chip_$#"
    combine_group_name = getRecoCombineName()

    # define the names for the datasets which we want to write
    int_dset_names = ["hits", "sumToT", "eventNumber"]
    float_dset_names = ["centerX", "centerY", "energyFromPixel", "rmsLongitudinal", "rmsTransverse", 
                        "skewnessLongitudinal", "skewnessTransverse", "kurtosisLongitudinal", "kurtosisTransverse",
                        "eccentricity", "rotationAngle", "length", "width", "fractionInTransveseRms"]
  # now parsing all the data is really fucking ugly, thanks to the tons of
  # different variables, which we want to write :( Unfortunately, we cannot
  # simply make that a compound datatype or something. Well
  var
    # create group for each chip
    chip_groups = mapIt(toSeq(0..<nchips), h5f.create_group(chip_group_name % $it))
    combine_group = h5f.create_group(combine_group_name)

    # NOTE: can we not handle this with two tables instead?
  var
    # before we create the datasets, first parse the data we will write
    x  = newSeq[seq[seq[int]]](nchips)
    y  = newSeq[seq[seq[int]]](nchips)
    ch = newSeq[seq[seq[int]]](nchips)
    ev_numbers = newSeq[seq[int]](nchips)
    hits = newSeq[seq[int]](nchips)
    sum_tot = newSeq[seq[int]](nchips)
    # float seqs
    pos_x = newSeq[seq[float]](nchips)
    pos_y = newSeq[seq[float]](nchips)
    rms_long = newSeq[seq[float]](nchips)
    rms_trans = newSeq[seq[float]](nchips)
    eccentricity = newSeq[seq[float]](nchips)
    rot_angle = newSeq[seq[float]](nchips)
    skew_long = newSeq[seq[float]](nchips)
    skew_trans = newSeq[seq[float]](nchips)
    kurt_long = newSeq[seq[float]](nchips)
    kurt_trans = newSeq[seq[float]](nchips)
    length = newSeq[seq[float]](nchips)
    width = newSeq[seq[float]](nchips)
    fraction_transverse_rms = newSeq[seq[float]](nchips)

  for chip in 0 ..< nchips:
    pos_x[chip] = @[]
    pos_y[chip] = @[]
    hits[chip] = @[]
    sum_tot[chip] = @[]
    ev_numbers[chip] = @[]
    # float seqs
    rms_long[chip] = @[]
    rms_trans[chip] = @[]
    eccentricity[chip] = @[]
    rot_angle[chip] = @[]
    skew_long[chip] = @[]
    skew_trans[chip] = @[]
    kurt_long[chip] = @[]
    kurt_trans[chip] = @[]
    length[chip] = @[]
    width[chip] = @[]
    fraction_transverse_rms[chip] = @[]
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
      x[chip].add(newSeq[int](cl.data.len))
      y[chip].add(newSeq[int](cl.data.len))
      ch[chip].add(newSeq[int](cl.data.len))
      for j in 0..cl.data.high:
        x[chip][^1][j] = cl.data[j].x
        y[chip][^1][j] = cl.data[j].y
        ch[chip][^1][j] = cl.data[j].ch
      # int seqs
      pos_x[chip].add cl.pos_x
      pos_y[chip].add cl.pos_y
      hits[chip].add cl.hits
      sum_tot[chip].add cl.sum_tot
      ev_numbers[chip].add num
      # float seqs
      rms_long[chip].add cl.geometry.rms_long
      rms_trans[chip].add cl.geometry.rms_trans 
      eccentricity[chip].add cl.geometry.eccentricity
      rot_angle[chip].add cl.geometry.rot_angle 
      skew_long[chip].add cl.geometry.skew_long 
      skew_trans[chip].add cl.geometry.skew_trans
      kurt_long[chip].add cl.geometry.kurt_long 
      kurt_trans[chip].add cl.geometry.kurt_trans
      length[chip].add cl.geometry.length 
      width[chip].add cl.geometry.width
      fraction_transverse_rms[chip].add cl.geometry.fraction_transverse_rms

  # now that we have the data and now how many elements each type has
  # we can create the datasets
  echo "Now creating dset stuff"      
  var
    # datasets for x, y and charge
    int_dsets = initTable[string, seq[H5DataSet]]()
    float_dsets = initTable[string, seq[H5DataSet]]()    
    # x, y and charge datasets
    x_dsets  = mapIt(toSeq(0..<nchips),
                     h5f.create_dataset(chip_groups[it].name & "/x", x[it].len, dtype = ev_type))
    y_dsets  = mapIt(toSeq(0..<nchips),
                     h5f.create_dataset(chip_groups[it].name & "/y", y[it].len, dtype = ev_type))
    ch_dsets  = mapIt(toSeq(0..<nchips),
                     h5f.create_dataset(chip_groups[it].name & "/ToT", ch[it].len, dtype = ev_type))
  
  for dset in int_dset_names:
    int_dsets[dset] = newSeq[H5DataSet](nchips)
    for chip in 0 ..< nchips:
      int_dsets[dset][chip] = h5f.create_dataset(chip_groups[chip].name / dset, x[chip].len, dtype = int)
  for dset in float_dset_names:
    float_dsets[dset] = newSeq[H5DataSet](nchips)
    for chip in 0 ..< nchips:
      float_dsets[dset][chip] = h5f.create_dataset(chip_groups[chip].name / dset, x[chip].len, dtype = float)  
    
  # now that we have the datasets, write everything...
  let all = x_dsets[0].all
  # get locations of raw data groups, so that we can copy
  # the attributes
  let raw_groups = rawDataChipBase(run_number)    
  for chip in 0 ..< nchips:
    int_dsets["hits"][chip][all] = hits[chip]
    int_dsets["sumToT"][chip][all] = sum_tot[chip]
    int_dsets["eventNumber"][chip][all] = ev_numbers[chip]
    # float seqs
    float_dsets["centerX"][chip][all] = pos_x[chip]
    float_dsets["centerY"][chip][all] = pos_y[chip]
    float_dsets["rmsLongitudinal"][chip][all] = rms_long[chip]
    float_dsets["rmsTransverse"][chip][all] = rms_trans[chip]
    float_dsets["eccentricity"][chip][all] = eccentricity[chip]
    float_dsets["rotationAngle"][chip][all] = rot_angle[chip]
    float_dsets["skewnessLongitudinal"][chip][all] = skew_long[chip]
    float_dsets["skewnessTransverse"][chip][all] = skew_trans[chip]
    float_dsets["kurtosisLongitudinal"][chip][all] = kurt_long[chip]
    float_dsets["kurtosisTransverse"][chip][all] = kurt_trans[chip]
    float_dsets["length"][chip][all] = length[chip]
    float_dsets["width"][chip][all] = width[chip]
    float_dsets["fractionInTransveseRms"][chip][all] = fraction_transverse_rms[chip]
    x_dsets[chip][all] = x[chip]
    y_dsets[chip][all] = y[chip]
    ch_dsets[chip][all] = ch[chip]
  
    # anyways, write the chip dataset attributes
    let raw_chip_group_name = raw_groups & $chip
    var raw_group = h5f[raw_chip_group_name.grp_str]
    let ch_numb = raw_group.attrs["chipNumber", int]
    let ch_name = raw_group.attrs["chipName", string]
    # and write these to the current group
    chip_groups[ch_numb].attrs["chipNumber"] = ch_numb
    chip_groups[ch_numb].attrs["chipName"] = ch_name
  
  # what in the actual fuck


iterator readDataFromH5(h5f: var H5FileObj, group: string, run_number: int): (int, seq[Pixels]) =
  # proc to read data from the HDF5 file from `group`
  var chip_base = rawDataChipBase(run_number)
  let raw_data_basename = rawDataBase()
  for grp in keys(h5f.groups):
    if chip_base in grp:
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
      # combine the raw data into a sequence of pixels
      let run_pix = map(toSeq(0..raw_x.high), (i: int) -> seq[(int, int, int)] =>
                        mapIt(toSeq(0..raw_x[i].high), (int(raw_x[i][it]), int(raw_y[i][it]), int(raw_ch[i][it]))))
      # and yield them
      yield (chip_number, run_pix)

proc newClusterGeometry(): ClusterGeometry =
  result = ClusterGeometry(rms_long: Inf,
                           rms_trans: Inf,
                           eccentricity: Inf,
                           rot_angle: Inf,
                           skew_long: Inf,
                           skew_trans: Inf,
                           kurt_long: Inf,
                           kurt_trans: Inf,
                           length: Inf,
                           width: Inf,
                           fraction_transverse_rms: Inf)


proc newClusterObject(c: Cluster): ClusterObject =
  # initialize variables with Inf for now
  # TODO: should we initialize geometry values by Inf as well?
  let geometry = ClusterGeometry()
  result = ClusterObject(data: c,
                         pos_x: Inf,
                         pos_y: Inf,
                         energy: Inf,
                         geometry: geometry)


proc excentricity(n: cuint, p: array[1, cdouble], grad: var array[1, cdouble], func_data: var pointer): cdouble {.cdecl.} =
  # this function calculates the excentricity of a found pixel cluster using nimnlopt.
  # Since no proper high level library is yet available, we need to pass a var pointer
  # of func_data, which contains the x and y arrays in which the data is stored, in
  # order to calculate the RMS variables

  # first recover the data from the pointer to func_data, by casting the
  # raw pointer to a Cluster object
  let fit = cast[FitObject](func_data)
  let c = fit.cluster
  let (x, y) = fit.xy

  var
    sum_x: cdouble = 0
    sum_y: cdouble = 0
    sum_x2: cdouble = 0
    sum_y2: cdouble = 0

  # echo "Starting new calc with param ", p[0]
  for i in 0..<len(c):
    let
      new_x = cos(p[0]) * (cdouble(c[i].x) - cdouble(x)) * PITCH - sin(p[0]) * (cdouble(c[i].y) - cdouble(y)) * PITCH
      new_y = sin(p[0]) * (cdouble(c[i].x) - cdouble(x)) * PITCH + cos(p[0]) * (cdouble(c[i].y) - cdouble(y)) * PITCH
    # echo "new_x is ", new_x, " and new_y ", new_y
    sum_x += new_x
    sum_y += new_y
    sum_x2 += (new_x * new_x)
    sum_y2 += (new_y * new_y)
  
  let
    n_elements: cdouble = cdouble(len(c))
    rms_x: cdouble = sqrt( (sum_x2 / n_elements) - (sum_x * sum_x / n_elements / n_elements))
    rms_y: cdouble = sqrt( (sum_y2 / n_elements) - (sum_y * sum_y / n_elements / n_elements))    

  #echo "sum_x2 / elements ", sum_x2 / n_elements, "sum_x * sum_x / elements / elements ", sum_x * sum_x / n_elements / n_elements
  let exc = rms_x / rms_y

  # need to check whether grad is nil. Only used for some algorithms, otherwise a
  # NULL pointer is handed in C
  if unlikely(addr(grad) != nil):
    # normally we'd calculate the gradient for the current parameters, but we're
    # not going to use it. Can also remove this whole if statement
    discard

  #echo "parameters ", sum_x, " ", sum_y, " ", sum_x2, " ", sum_y2, " ", rms_x, " ", rms_y, " ", n_elements, " exc ", exc
  result = -exc

proc calcGeomtry(cluster: Cluster, pos_x, pos_y, rot_angle: float): ClusterGeometry =
  # given a cluster and the rotation angle of it, calculate the different
  # statistical moments, i.e. RMS, skewness and kurtosis in longitudinal and
  # transverse direction
  # done by rotating the cluster by the angle to define the two directions
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
  result.length       = max(xRot) - min(xRot)
  result.width        = max(yRot) - min(yRot)
  result.rms_long     = stat_x.standardDeviation()
  result.rms_trans    = stat_y.standardDeviation()
  result.skew_long    = stat_x.skewness()
  result.skew_trans   = stat_y.skewness()
  result.kurt_long    = stat_x.kurtosis()
  result.kurt_trans   = stat_y.kurtosis()
  result.rot_angle    = rot_angle
  result.eccentricity = result.rms_long / result.rms_trans
  # get fraction of all pixels within the transverse RMS, by filtering all elements
  # within the transverse RMS radius and dividing by total pix
  # DEBUG
  # echo "rms trans is ", result.rms_trans
  # echo "std is ", stat_y.variance()
  # echo "thus filter is ", filterIt(zip(xRot, yRot), distance(it.a, it.b) <= result.rms_trans)
  result.fraction_transverse_rms = float(filterIt(zip(xRot, yRot), distance(it.a, it.b) <= result.rms_trans).len) / float(npix)

proc isPixInSearchRadius(p1, p2: Coord, search_r: int): bool =
  # given two pixels, p1 and p2, we check whether p2 is within one square search
  # of p1
  # inputs:
  #   p1: Pix = pixel from which to start search
  #   p2: Pix = pixel for which to check
  #   search_r: int = search radius (square) in which to check for p2 in (p1 V search_r)
  # outpuits:
  #   bool = true if within search_r
  #          false if not
  let
    # determine boundary of search space
    right = p1.x + search_r
    left  = p1.x - search_r
    up    = p1.y + search_r
    down  = p1.y - search_r
  # NOTE: for performance we may want to use the fact that we may know that
  # p1 is either to the left (if in the same row) or below (if not in same row)
  var
    in_x: bool = false
    in_y: bool = false

  if p2.x < right and p2.x > left:
    in_x = true
  if p2.y < up and p2.y > down:
    in_y = true
  result = if in_x == true and in_y == true: true else: false

proc findSimpleCluster*(pixels: Pixels): seq[Cluster] =
  # this procedure searches for clusters based on a fixed search radius, whether
  # a pixel is found within that boundary, e.g. searchRadius = 50:
  # if within 50 pixels another pixel is found, add pixel to cluster, continue
  # inputs:
  #   -
  # ouputs:
  #   - 

  # - iterate over all hit pixels in event
  # - check next pixel, is it within search bound? yes,
  #   add pixel to hit
  var
    # sequence in which to store the hits, which are part of a cluster
    c: Cluster = @[]
    # create copy of pixels so that we can remove elements from it
    raw_event = pixels
    # counter
    i = 0 
  result = @[] #new seq[Cluster]
  
  let
    search_r = 50
    cutoff_size = 5

  # add the first pixel of the given sequence to have a starting pixel, from which we
  # look for other pixels in the cluster
  c.add(pixels[0])
  raw_event.deleteIntersection(@[pixels[0]])
  
  while raw_event.len > 0 and i < c.len:
    let p1: Coord = (x: c[i].x, y: c[i].y)
    # alternatively:
    let t = filter(raw_event, (p: tuple[x, y, ch: int]) -> bool => isPixInSearchRadius(p1, (p.x, p.y), search_r))

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

template eccentricityNloptOptimizer(fit_object: FitObject): NloptOpt =
  # returns the already configured Nlopt optimizer to fit the rotation angle /
  # eccentricity
  # set  the values of the fit objectn
  var
    # set the boundary values corresponding to range of 360 deg
    lb = (-4.0 * arctan(1.0), 4.0 * arctan(1.0))
  var opt = newNloptOpt("LN_COBYLA", lb)
  # hand the function to fit as well as the data object we need in it
  opt.setFunction(excentricity, fit_object)
  # set relative precisions of x and y, as well as limit max time the algorithm
  # should take to 1 second
  # these default values have proven to be working
  opt.xtol_rel = 1e-8
  opt.ftol_rel = 1e-8
  opt.maxtime  = 1.0
  opt.initial_step = 0.02
  opt

template fitRotAngle(cl_obj: ClusterObject, rotAngleEstimate: float): (float, float) =
  # simple template which wraps the optimization of the rotation angle /
  # eccentricity
  # cluster object is handed as var to avoid any copying

  # TODO: think about what to do with 11810 pixels, throw them out
    
  var
    # set the fit object with which we hand the necessary data to the
    # excentricity function
    fit_object: FitObject
    # the resulting fit parameter  
    p = @[rotAngleEstimate]
    
  fit_object.cluster = cl_obj.data
  fit_object.xy = (x: cl_obj.pos_x, y: cl_obj.pos_y)
  var opt = eccentricityNloptOptimizer(fit_object)
  # start minimization
  let (params, min_val) = opt.optimize(p)
  if opt.status < NLOPT_SUCCESS:
    echo opt.status
    echo "nlopt failed!\n"
  # clean up optimizer
  nlopt_destroy(opt.optimizer)
  # now return the rotation angle and eccentricity
  (params[0], min_val)
  
proc recoCluster(c: Cluster): ClusterObject =
  result = newClusterObject(c)
  let
    clustersize: int = len(c)
    (sum_x, sum_y, sumTotInCluster) = sum(c)
    # using map and sum[Pix] we can calculate sum of x^2, y^2 and x*y in one line
    (sum_x2, sum_y2, sum_xy) = sum(map(c, (p: Pix) -> Pix => (p.x * p.x,
                                                              p.y * p.y,
                                                              p.x * p.y)))
    pos_x = float64(sum_x) / float64(clustersize)
    pos_y = float64(sum_y) / float64(clustersize)
  var
    rms_x = sqrt(float64(sum_x2) / float64(clustersize) - pos_x * pos_x)
    rms_y = sqrt(float64(sum_y2) / float64(clustersize) - pos_y * pos_y)
    rotAngleEstimate = arctan( (float64(sum_xy) / float64(clustersize)) -
                               pos_x * pos_y / (rms_x * rms_x))

  # set the total "charge" in the cluster (sum of ToT values), can be
  # converted to electrons with ToT calibration
  result.sum_tot = sumTotInCluster
  # set number of hits in cluster
  result.hits = clustersize
  # set the position
  (result.pos_x, result.pos_y) = applyPitchConversion(pos_x, pos_y)
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
    echo "Rot angle estimate is NaN, vals are ", rms_x, " ", rms_y
    # what do we do in this case with the geometry?!
    #raise newException(ValueError, "Rotation angle estimate returned bad value")
    echo "Fit will probably fail!"
    
  # else we can minimize the rotation angle and calc the eccentricity
  let (rot_angle, eccentricity) = fitRotAngle(result, rotAngleEstimate)

  # now we still need to use the rotation angle to calculate the different geometric
  # properties, i.e. RMS, skewness and kurtosis along the long axis of the cluster
  result.geometry = calcGeomtry(c, result.pos_x, result.pos_y, rot_angle)

proc recoEvent(data: Pixels, event, chip: int): ref RecoEvent =
  result = new RecoEvent
  result.event_number = event
  result.chip_number = chip
  if data.len > 0:
    let cluster = findSimpleCluster(data)
    result.cluster = newSeq[ClusterObject](cluster.len)
    for i, cl in cluster:
      result.cluster[i] = recoCluster(cl)

{.experimental.}
#proc reconstructSingleChip(data: seq[Pixels], run, chip: int): seq[ref RecoEvent] =
proc reconstructSingleChip(data: seq[Pixels], run, chip: int): seq[FlowVar[ref RecoEvent]] =
  # procedure which receives pixel data for a given chip and run number
  # and performs the reconstruction on it
  # inputs:
  #    data: seq[Pixels] = data of pixels for this chip containing pixels for each event
  #    run: int = run number of run
  #    chip: int = chip number working on
  echo "Working on chip $# in run $# " % [$chip, $run]
  echo "We have $# events to reconstruct" % $data.len
  var count = 0
  result = newSeq[FlowVar[ref RecoEvent]](data.len)
  #result = newSeq[ref RecoEvent](data.len)  
  parallel:
    for event in 0..data.high:
      if event < result.len:
        result[event] = spawn recoEvent(data[event], event, chip)
      echoFilesCounted(count, 2500)
  sync()

#iterator matchingGroup(h5f: var H5FileObj,   

proc reconstructAllRunsInFile(h5f: var H5FileObj, flags_tab: Table[string, bool], calib_factor: float = 1.0) =
  ## proc which performs reconstruction of all runs in a given file
  ## if the --only_energy command line argument is set, we skip the reconstruction
  ## and only perform an energy calibration of the existing (!) reconstructed
  ## runs using the calibration factor given
  
  let
    raw_data_basename = rawDataBase()
    run_regex = re(raw_data_basename & r"(\d+)$")
    t0 = epochTime()
  var run: array[1, string]
  var reco_run: seq[FlowVar[ref RecoEvent]] = @[]
  for grp in keys(h5f.groups):
    if grp.match(run_regex, run) == true:
      # now read some data. Return value will be added later
      let run_number = parseInt(run[0])

      if flags_tab["only_energy"] == false and flags_tab["only_fadc"] == false:
        # TODO: we can in principle perform energy calibration in one go
        # together with creation of spectrum, if we work as follows:
        # 1. calibration runs:
        #    - need to interface with Python code, i.e. call fitting procedure,
        #      which returns the value to the Nim program as its return value
        
        let t1 = epochTime()      
        for chip, pixdata in h5f.readDataFromH5(grp, run_number):
          # given single runs pixel data, call reconstruct run proc
          # NOTE: the data returned from the iterator contains all
          # events in ascending order of event number, i.e.
          # [0] -> eventNumber == 0 and so on
          reco_run.add reconstructSingleChip(pixdata, run_number, chip)
          
        echo "Reconstruction of run $# took $# seconds" % [$run_number, $(epochTime() - t1)]
        # finished run, so write run to H5 file
        h5f.writeRecoRunToH5(reco_run, run_number)
        # set reco run length back to 0
        reco_run.setLen(0)

        # now check whether create iron spectrum flag is set
        if flags_tab["create_fe"] == true:
          createFeSpectrum(h5f, run_number)
        elif flags_tab["calib_energy"] == true:
          applyEnergyCalibration(h5f, run_number, 1.1)
      else:
        # only perform energy calibration of the reconstructed runs in file
        # check if reconstructed run exists
        if hasKey(h5f.groups, (recoBase & $run_number)) == true:
          if flags_tab["only_energy"] == true:
            applyEnergyCalibration(h5f, run_number, calib_factor)
          if flags_tab["only_fadc"] == true:
            calcRiseAndFallTimes(h5f, run_number)
        else:
          echo "No reconstructed run found for $#" % $grp
    else:
      # this is the case in which group was not matched to regex
      discard
  echo "Reconstruction of all runs in $# took $# seconds" % [$h5f.name, $(epochTime() - t0)]

proc reconstructSingleRunInFile(h5f: var H5FileObj,
                                run_number: int,
                                flags_tab: Table[string, bool],
                                calib_factor: float = 1.0) =
  ## proc which performs reconstruction of a single run in a given file
  ## if the --only_energy command line argument is set, we skip the reconstruction
  ## and only perform an energy calibration of the existing (!) reconstructed
  ## runs using the calibration factor given

  ## TODO: combine this with proc above
  let
    raw_data_basename = rawDataBase()
    run_regex = re(raw_data_basename & $run_number & r"$")
    t0 = epochTime()
  echo "Reading group name $#" % $(raw_data_basename & $run_number & r"$")
  var reco_run: seq[FlowVar[ref RecoEvent]] = @[]
  for grp in keys(h5f.groups):
    if grp.match(run_regex) == true:
      # now read some data. Return value will be added later
      if flags_tab["only_energy"] == false and flags_tab["only_fadc"] == false:
        # TODO: we can in principle perform energy calibration in one go
        # together with creation of spectrum, if we work as follows:
        # 1. calibration runs:
        #    - need to interface with Python code, i.e. call fitting procedure,
        #      which returns the value to the Nim program as its return value
        
        let t1 = epochTime()      
        for chip, pixdata in h5f.readDataFromH5(grp, run_number):
          # given single runs pixel data, call reconstruct run proc
          # NOTE: the data returned from the iterator contains all
          # events in ascending order of event number, i.e.
          # [0] -> eventNumber == 0 and so on
          reco_run.add reconstructSingleChip(pixdata, run_number, chip)
          
        echo "Reconstruction of run $# took $# seconds" % [$run_number, $(epochTime() - t1)]
        # finished run, so write run to H5 file
        h5f.writeRecoRunToH5(reco_run, run_number)
        # set reco run length back to 0
        reco_run.setLen(0)

        # now check whether create iron spectrum flag is set
        if flags_tab["create_fe"] == true:
          createFeSpectrum(h5f, run_number)
        elif flags_tab["calib_energy"] == true:
          applyEnergyCalibration(h5f, run_number, 1.1)
      else:
        # only perform energy calibration of the reconstructed runs in file
        # check if reconstructed run exists
        if hasKey(h5f.groups, (recoBase & $run_number)) == true:
          if flags_tab["only_energy"] == true:
            applyEnergyCalibration(h5f, run_number, calib_factor)
          if flags_tab["only_fadc"] == true:
            calcRiseAndFallTimes(h5f, run_number)
        else:
          echo "No reconstructed run found for $#" % $grp

      # since we found the run, we can break from the loop
      break
    else:
      # this is the case in which group was not matched to regex
      discard
  echo "Reconstruction of run $# in $# took $# seconds" % [$run_number, $h5f.name, $(epochTime() - t0)]

proc reconstructSingleRunFolder(folder: string) =
  # procedure which receives path to a run folder and reconstructs the objects
  # in that folder
  # inputs:
  #    folder: string = the run folder from which to reconstruct events
  echo "Starting to read list of files"
  let    
    files = getSortedListOfFiles(folder, EventSortType.inode, EventType.InGridType)
    regex_tup = getRegexForEvents()
    f_to_read = if files.high < 30000: files.high else: 30000
    data = readListOfInGridFiles(files[0..f_to_read], regex_tup)
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
        if ob.geometry.rot_angle < min_val:
          min_val = ob.geometry.rot_angle
        min_seq[num].add ob.geometry.rot_angle
  dumpRotAngle(min_seq)

proc main() =

  # create command line arguments using docopt
  let args = docopt(doc)
  echo args
  
  let
    h5f_name = $args["<HDF5file>"]
    create_fe_arg = $args["--create_fe_spec"]
    calib_energy_arg = $args["--calib_energy"]

  var
    run_number = $args["--run_number"]
    outfile = $args["--out"]
    create_fe_flag: bool
    calib_energy_flag: bool
    calib_factor_str = $args["--only_energy"]
    calib_factor: float = Inf
    only_energy_flag: bool
    only_fadc: bool
  if run_number == "nil":
    run_number = ""
  if outfile == "nil":
    echo "Currently not implemented. What even for?"
    outfile = "None"
  if create_fe_arg == "nil":
    create_fe_flag = false
  else:
    create_fe_flag = true
  if calib_energy_arg == "nil":
    calib_energy_flag = false
  else:
    calib_energy_flag = true
  if calib_factor_str == "nil":
    only_energy_flag = false
  else:
    only_energy_flag = true
    calib_factor = parseFloat(calib_factor_str)
  if $args["--only_fadc"] != "nil":
    only_fadc = true
  else:
    only_fadc = false
      

  let flags_tab = { "create_fe": create_fe_flag,
                    "calib_energy": calib_energy_flag,
                    "only_energy": only_energy_flag,
                    "only_fadc": only_fadc}.toTable                  
    

  var h5f = H5file(h5f_name, "rw")
  # visit the whole file to read which groups exist
  h5f.visitFile
  let raw_data_basename = rawDataBase()
  if run_number == "":
    reconstructAllRunsInFile(h5f, flags_tab, calib_factor)
  else:
    reconstructSingleRunInFile(h5f, parseInt(run_number), flags_tab, calib_factor)
  
  # NOTE: there's no point for this anymore, at least not at the moment
  # # first check whether given folder is valid run folder
  # let (is_run_folder, contains_run_folder) = isTosRunFolder(folder)
  # echo "Is run folder       : ", is_run_folder
  # echo "Contains run folder : ", contains_run_folder
  # reconstructSingleRun(folder)

when isMainModule:
  main()
