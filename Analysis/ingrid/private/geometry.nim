import math, stats
import ingrid / ingrid_types
import helpers / utils
import logging
import ingrid / private / [pure, cdl_cuts]

import seqmath, sequtils, nlopt

#[
This module contains (among others) the actual ``reconstruction`` calculations
done by `reconstruction.nim`, that is those calculation, which are not part of
the calibration steps (via `--only_*` flags), but instead concern the conversion
from raw data to reconstructed data. This is basically cluster finding and
rotation angle / geometrical property calcuations.
]#

type
  # fit object, which is handed to the NLopt library in the
  # `VarStruct` -> to the eccentricity function
  FitObject[T: SomePix] = object
    cluster: Cluster[T]
    xy: tuple[x, y: float64]

macro hijackMe(procImpl: untyped): untyped =
  when defined(activateHijack) and declared(replaceBody):
    replaceBody(procImpl)
  else:
    procImpl

################################################################################
############# Geometry calculation related procs ###############################
################################################################################

proc newClusterGeometry*(): ClusterGeometry =
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


proc newClusterObject*[T: SomePix](c: Cluster[T]): ClusterObject[T] =
  ## initialize variables with Inf for now
  # TODO: should we initialize geometry values by Inf as well?
  let geometry = ClusterGeometry()
  result = ClusterObject[T](data: c,
                            centerX: Inf,
                            centerY: Inf,
                            energy: Inf,
                            geometry: geometry)

template withSeptemXY*(chipNumber: int, actions: untyped): untyped =
  ## injects the x0, y0 coordinates of the given chip number embedded into
  ## the septem frame
  var
    x0 {.inject.}: int
    y0 {.inject.}: int
  case chipNumber
  of 0:
    # chip bottom left of board
    y0 = 0 # top end of bottom row
    x0 = 128 # shifted by half a chip to the right
  of 1:
    # chip bottom right
    y0 = 0
    x0 = 128 + 256
  of 2:
    # middle left
    y0 = 256
    x0 = 0
  of 3:
    # middle middle
    y0 = 256
    x0 = 256
  of 4:
    # middle right
    y0 = 256
    x0 = 2 * 256
  of 5:
    # top right chip
    y0 = 3 * 256
    x0 = 2 * 256 + 128
  of 6:
    # top left chip
    y0 = 3 * 256
    x0 = 128 + 256
  actions

func determineChip*[T:  SomePix](p: T): int =
  ## determines which chip the given septem pixel coordinate corresponds to
  if p.y in 0 .. 255 and p.x in 128 .. 128 + 255:
    # bottom left
    result = 0
  elif p.y in 0 .. 255:
    # bottom right
    result = 1
  elif p.y in 256 .. 511 and p.x in 0 .. 255:
    # middle left
    result = 2
  elif p.y in 256 .. 511 and p.x in 256 .. 511:
    # center
    result = 3
  elif p.y in 256 .. 511:
    # middle right
    result = 4
  elif p.x in 128 + 256 .. 512 + 127:
    # top right
    result = 5
  elif p.x in 128 .. 128 + 255:
    # top left
    result = 6
  else:
    raise newException(Exception, "This chip should not exist! " & $p)

func septemPixToChpPix*[T: SomePix](p: T, chipNumber: range[0 .. 6]): T =
  ## inverse of chpPixToSeptemPix
  result = p
  case chipNumber
  of 0:
    # chip bottom left of board
    result.x = result.x - 128 # shifted by half a chip to the right
  of 1:
    # chip bottom right
    result.x = result.x - (128 + 256)
  of 2:
    # middle left
    result.y = result.y - 256
  of 3:
    # middle middle
    result.y = result.y - 256
    result.x = result.x - 256
  of 4:
    # middle right
    result.y = result.y - 256
    result.x = result.x - 512
  of 5:
    # top right chip
    result.y = -(result.y - 3 * 256)
    result.x = -(result.x - (2 * 256 + 128))
  of 6:
    # top left chip
    result.y = -(result.y - 3 * 256)
    result.x = -(result.x - (128 + 256))

func chpPixToSeptemPix*(p: Pix, chipNumber: range[0 .. 6]): PixInt =
  ## converts the given local chip pixel to the full septem frame coordinate system
  withSeptemXY(chipNumber):
    var xIdx, yIdx: int
    case chipNumber
    of 0, 1, 2, 3, 4:
      xIdx = x0 + p.x.int
      yIdx = y0 + p.y.int
    of 5, 6:
      xIdx = x0 - p.x.int
      yIdx = y0 - p.y.int
    result = (x: xIdx, y: yIdx, ch: p.ch.int)

func chpPixToSeptemPix*(pix: Pixels, chipNumber: range[0 .. 6]): PixelsInt =
  ## converts the given local chip pixels to the full septem frame coordinate system
  result.setLen(pix.len)
  for i, p in pix:
    let pp = chpPixToSeptemPix(p, chipNumber)
    result[i] = pp

# proc sum*[T: tuple](s: seq[T]): T {.inline.} =
#   # this procedure sums the given array along the given axis
#   # if T is itself e.g. a tuple, we will return a tuple, one
#   # element for each field in the tuple
#   assert s.len > 0, "Can't sum empty sequences"
#   var sum_t: T
#   for p in s:
#     for n, f in fieldPairs(p):
#       sum_t[f] += p[n]

proc calcCentroidOfEvent*(pix: Pixels): tuple[x, y: float] =
  ## proc to calc centroid of the given pixels
  ## inputs:
  ##    pixels object (seq[tuple[x, y, ch: int]]) containing raw event
  ## outputs:
  ##    tuple[x, y: int]: tuple containing centroid x and y position
  ## let x = map(pix, (p: tuple[x, y, ch: int]) -> int => p.x)
  ## let y = map(pix, (p: tuple[x, y, ch: int]) -> int => p.y)
  ## let sum_x = foldl(x, a + b)
  ## let sum_y = foldl(y, a + b)
  var
    sum_x: int = 0
    sum_y: int = 0
  for p in pix:
    sum_x += p.x.int
    sum_y += p.y.int
  #let (sum_x, sum_y, sum_ch) = sum(pix)
  result.x = float(sum_x) / float(len(pix))
  result.y = float(sum_y) / float(len(pix))


proc isNearCenterOfChip*(pix: Pixels): bool =
  ## proc to check whether event is located around center of chip
  ## inputs:
  ##    pixels object (seq[tuple[x, y, ch: int]]) containing raw event
  ## outputs:
  ##    true if within 4.5mm center square, false otherwise
  let (center_x, center_y) = calcCentroidOfEvent(pix)
  # pitch in um
  let pitch = 0.05
  let n_pix_to_bound = 2.25 / pitch
  # center pixel is (127, 127)
  let center_pix = 127'f
  var
    in_x = false
    in_y = false
  if center_x > (center_pix - n_pix_to_bound) and center_x < (center_pix + n_pix_to_bound):
    in_x = true
  if center_y > (center_pix - n_pix_to_bound) and center_y < (center_pix + n_pix_to_bound):
    in_y = true
  if in_x == true and in_y == true:
    result = true
  else:
    result = false

# template which calculates euclidean distance between 2 points
template distance*(x, y): float = sqrt(x * x + y * y)

# template which returns pitch converted positions on chip pixel values
# to mm from center of chip
# constants are:
# const NPIX = 256
# const PITCH = 0.0055 (see ingrid_types)
template applyPitchConversion*[T: (float | SomeInteger)](x, y: T): (float, float) =
  ## template which returns the converted positions on a Timepix
  ## pixel position --> position from center in mm
  ((float(NPIX) - float(x) + 0.5) * PITCH, (float(y) + 0.5) * PITCH)


func inRegion*(centerX, centerY: float, region: ChipRegion): bool {.inline.} =
  ## returns the result of a cut on a certain chip `region`. Inputs the
  ## `centerX` and `centerY` position of a cluster and returns true if
  ## the cluster is within the region
  const centerChip = 7.0
  case region
  of crGold:
    # make sure this is only initialized once somehow...
    let regCut = getRegionCut(region)
    result = if centerX >= regCut.xMin and
                centerX <= regCut.xMax and
                centerY >= regCut.yMin and
                centerY <= regCut.yMax:
               true
             else:
               false
  of crAll:
    # simply always return good
    result = true
  else:
    # make sure this is only initialized once somehow...
    let regCut = getRegionCut(region)
    # silver and bronze region only different by radius
    let
      xdiff = (centerX - centerChip)
      ydiff = (centerY - centerChip)
      radius = distance(xdiff, ydiff)
    # TODO: gold cut is NOT part of the silver region (see C. Krieger PhD p. 133)
    result = if radius <= regCut.radius: true else : false

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

proc calcGeometry*[T: SomePix](cluster: Cluster[T],
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
                                     kind: nlopt.FuncKind.NoGrad)
  result.setFunction(varStruct)
  # set relative precisions of x and y, as well as limit max time the algorithm
  # should take to 1 second
  # these default values have proven to be working
  result.xtol_rel = 1e-8
  result.ftol_rel = 1e-8
  result.maxtime  = 1.0
  result.initial_step = 0.02

proc fitRotAngle*[T: SomePix](cl_obj: ClusterObject[T],
                              rotAngleEstimate: float): (float, float) = #{.hijackMe.} =
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

proc recoCluster*[T: SomePix](c: Cluster[T]): ClusterObject[T] {.gcsafe.} =
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

proc recoEvent*[T: SomePix](dat: tuple[pixels: seq[T], eventNumber: int],
                            chip, run: int): ref RecoEvent[T] {.gcsafe, hijackMe.} =
  result = new RecoEvent[T]
  result.event_number = data.eventNumber
  result.chip_number = chip
  if data[0].len > 0:
    let cluster = findSimpleCluster(data.pixels)
    result.cluster = newSeq[ClusterObject[T]](cluster.len)
    for i, cl in cluster:
      result.cluster[i] = recoCluster(cl)
