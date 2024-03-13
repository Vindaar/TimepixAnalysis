import std / [math, stats, strutils, envvars]
import .. / ingrid_types
import helpers / utils
import logging
import pure, cdl_cuts, clustering

import sequtils, nlopt
import arraymancer # for tensor and dbscan

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


proc newClusterObject*[T: SomePix](timepix: TimepixVersion): ClusterObject[T] =
  ## initialize variables with Inf for now
  # TODO: should we initialize geometry values by Inf as well?
  let geometry = ClusterGeometry()
  result = ClusterObject[T](centerX: Inf,
                            centerY: Inf,
                            energy: Inf,
                            geometry: geometry,
                            version: timepix)

proc to*[T: SomePix; U: SomePix](c: Cluster[T], _: typedesc[U]): Cluster[U] =
  ## Converts the input pix type to the output
  ## May throw away information
  when T is U: result = c
  elif T is Pix and U is PixTpx3:
    # return with empty `toa`, `toaCombined`
    warn("Conversion from `Pix` to `PixTpx3` adds empty ToA data!")
    result = newSeq[U](c.len)
    for i in 0 ..< result.len:
      result[i] = (x: c[i].x, y: c[i].y, ch: c[i].ch, toa: 0'u16, toaCombined: 0'u64)
  elif T is PixTpx3 and U is Pix:
    warn("Conversion from `PixTpx3` to `Pix` throws away ToA information!")
    result = newSeq[U](c.len)
    for i in 0 ..< result.len:
      result[i] = (x: c[i].x, y: c[i].y, ch: c[i].ch)
  elif T is PixInt or U is PixInt:
    error("Currently unsupported for `PixInt` type! Need to make sure we perform " &
      "coordinate transformation correctly!")

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
    # top right chip (- 1 as we start at top/right, which would be out of bounds if x/y == 0
    #                 for other chips x/y == 0 leads to first idx on next chip)
    y0 = 3 * 256 - 1
    x0 = 2 * 256 + 127
  of 6:
    # top left chip
    y0 = 3 * 256 - 1
    x0 = 256 + 127
  else: doAssert false, "Invalid chip number encountered in `withSeptemXY`"
  actions

const
  ## All the sizes are given in MilliMeter
  Width* = 14.1
  Height* = 14.1
  BondHeight* = 2.0
  FullHeight* = Height + BondHeight
  # helper distances
  YRow1Row2* = 0.38 # distance between row 1 and 2
  YRow2Row3* = 3.1 # distance between row 2 and 3
  # distances in X between chips on this row
  Row1XDist* = 0.85
  Row2XDist* = 0.35
  Row3XDist* = 0.35
  XSize* = 3 * Width + (3 - 1) * Row2XDist # 43mm
  YSize* = FullHeight + 2 * Height + YRow1Row2 + YRow2Row3 #3 * FullHeight + YRow1Row2 + YRow2Row3
  # offsets of the x positions of each row
  Row1XOffset* = 6.95
  Row2XOffset* = 0.0
  Row3XOffset* = 7.22
  # sizes in pixel of the full layout
  #           v--- size of tight layout in pixel
  #                       v--- size of chips in mm ('real size' of tight layout)
  #                                   v--- size in mm of the real layout
  XSizePix* = ((256 * 3) / (Width * 3) * XSize).ceil.int
  YSizePix* = ((256 * 3) / (Height * 3) * YSize).ceil.int

# subtract from 768 as we do the same in `applyPitchConversion` (why?)
# normalize by full width (14 mm * 3 chips) and scale to all pixels
proc toRealXPix*(x: float): int =
  clamp(((x / XSize) * XSizePix.float).round.int, 0, XSizePix-1)

proc toRealYPix*(y: float): int =
  clamp(((y / YSize) * YSizePix.float).round.int, 0, YSizePix-1)

proc toRealXPos*(x: int|float): float =
  (float(x) + 0.5) * (XSize / XSizePix.float)

proc toRealYPos*(y: int|float): float =
  (float(y) + 0.5) * (YSize / YSizePix.float)

# subtract from 768 as we do the same in `applyPitchConversion` to
# invert the view from 'top down at' the detector to 'looking through'
# the detector (like a camera)
# normalize by full width (14 mm * 3 chips) and scale to all pixels
## XXX: this is not actually used! I think because inversion here is not needed. The
proc toXPixAlt*(x: float): int = ## <- yeah, we don't use it.
  clamp((768 - (x / (TimepixSize * 3.0)) * 768.0).int, 0, 767)

proc toXPix*(x: float): int =
  clamp(((x / (TimepixSize * 3.0)) * 768.0).int, 0, 767)

proc toYPix*(y: float): int =
  clamp(((y / (TimepixSize * 3.0)) * 768.0).int, 0, 767)

template withRealSeptemXY*(chipNumber: int, actions: untyped): untyped =
  ## injects the x0, y0 coordinates of the given chip number embedded into
  ## the septem frame using the real distances between the chips (converted into
  ## single pixels though).
  var
    x0 {.inject.}: float
    y0 {.inject.}: float
  case chipNumber
  of 0:
    # chip bottom left of board
    y0 = 0.0 # top end of bottom row
    x0 = Row3XOffset # shifted by half a chip to the right
  of 1:
    # chip bottom right
    y0 = 0.0
    x0 = (Row3XOffset + Width + Row3XDist)
  of 2:
    # middle left
    y0 = (FullHeight + YRow2Row3)
    x0 = Row2XOffset
  of 3:
    # middle middle
    y0 = (FullHeight + YRow2Row3)
    x0 = (Row2XOffset + Width + Row2XDist)
  of 4:
    # middle right
    y0 = (FullHeight + YRow2Row3)
    x0 = (Row2XOffset + 2 * Width + 2 * Row2XDist)
  of 5:
    # top right chip (- 1 as we start at top/right, which would be out of bounds if x/y == 0
    #                 for other chips x/y == 0 leads to first idx on next chip)
    y0 = (FullHeight + 2 * Height + YRow1Row2 + YRow2Row3)
    x0 = (Row1XOffset + 2 * Width + Row1XDist) # why one more chip, i.e. starting from right?
  of 6:
    # top left chip
    y0 = (FullHeight + 2 * Height + YRow1Row2 + YRow2Row3)
    x0 = (Row1XOffset + Width) # why one more chip, i.e. starting from right?
  else: doAssert false, "Invalid chip number encountered in `withRealSeptemXY`"
  actions

func determineChip*[T:  SomePix](p: T, allowOutsideChip = false, useRealLayout = false): int =
  ## determines which chip the given septem pixel coordinate corresponds to
  ##
  ## If `allowOutsideChip` is true, we return `-1` for any position that is not in the
  ## valid chip pixels. This can happen for e.g. cluster centers on SeptemEvents
  ##
  ## Note: this function assumes the layout between the chips is 'tight', i.e. no spacing
  ## at all!
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
    if allowOutsideChip:
      result = -1
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
    result.y = -(result.y - (3 * 256 - 1))
    result.x = -(result.x - (2 * 256 + 127))
  of 6:
    # top left chip
    result.y = -(result.y - (3 * 256 - 1))
    result.x = -(result.x - (127 + 256))

func septemToChp*(p: tuple[x, y: float], chipNumber: range[0 .. 6]): tuple[x, y: float] =
  ## Converts global coordinates of the regular septemboard layout (tight, i.e. no spacing)
  ## into a chip local coordinate.
  result = p
  case chipNumber
  of 0:
    # chip bottom left of board
    result.x = result.x - 7.0 # shifted by half a chip to the right
  of 1:
    # chip bottom right
    result.x = result.x - (7.0 + 14.0)
  of 2:
    # middle left
    result.y = result.y - 14.0
  of 3:
    # middle middle
    result.y = result.y - 14.0
    result.x = result.x - 14.0
  of 4:
    # middle right
    result.y = result.y - 14.0
    result.x = result.x - 14.0*2.0
  of 5:
    # top right chip
    result.y = -(result.y - 3 * 14.0)
    result.x = -(result.x - (2 * 14.0 + 7.0))
  of 6:
    # top left chip
    result.y = -(result.y - 3 * 14.0)
    result.x = -(result.x - (7.0 * 14.0))


func chpPixToSeptemPix*[T: SomePix](p: T, chipNumber: range[0 .. 6], realLayout = false): PixInt =
  ## converts the given local chip pixel to the full septem frame coordinate system
  if realLayout:
    withRealSeptemXY(chipNumber):
      var xIdx, yIdx: int
      case chipNumber
      of 0, 1, 2, 3, 4:
        xIdx = x0.toRealXPix() + p.x.int
        yIdx = y0.toRealYPix() + p.y.int
      of 5, 6:
        xIdx = x0.toRealXPix() - p.x.int
        yIdx = y0.toRealYPix() - p.y.int
      result = (x: min(xIdx, XSizePix), y: min(yIdx, YSizePix), ch: p.ch.int)
  else:
    withSeptemXY(chipNumber):
      var xIdx, yIdx: int
      case chipNumber
      of 0, 1, 2, 3, 4:
        xIdx = x0 + p.x.int
        yIdx = y0 + p.y.int
      of 5, 6:
        xIdx = x0 - p.x.int
        yIdx = y0 - p.y.int
      result = (x: min(xIdx, 767), y: min(yIdx, 767), ch: p.ch.int)

proc septemPixToRealPix*(p: PixInt): PixInt =
  let chip = p.determineChip(allowOutsideChip = false)
  result = septemPixToChpPix(p, chip).chpPixToSeptemPix(chip, realLayout = true)

proc septemPixToRealPix*(p: PixelsInt): PixelsInt =
  ## Converts all given septemboard pixels (of the tight layout) to the real layout
  result = newSeq[PixInt](p.len)
  for i in 0 ..< p.len:
    result[i] = septemPixToRealPix(p[i])

proc tightToReal*(c: tuple[x, y: float]): tuple[x, y: float] =
  ## Converts a given set of x, y coordinates for the given chip number into
  ## coordinates on the real septemboard layout.
  let chipNumber = determineChip((x: c.x.toXPix(), y: c.y.toYPix(), ch: 0))
  let p = septemToChp(c, chipNumber)
  withRealSeptemXY(chipNumber):
    case chipNumber
    of 0, 1, 2, 3, 4:
      result.x = (XSize - (x0 + p.x)) # invert position here
      result.y = y0 + p.y
    of 5, 6:
      result.x = (XSize - (x0 - p.x))
      result.y = y0 - p.y
    else:
      doAssert false, "Invalid chip: " & $chipNumber

proc chpPixToRealPix*(a: int, isX: bool, chipNumber: range[0 .. 6]): int =
  ## Note: This is not particularly efficient... But doing it better requires
  ## a rewrite of the templates above.
  let pix = if isX: (x: a, y: 0, ch: 0)
            else: (x: 129, y: a, ch: 0) # 129 is a dummy that corresponds to a point that is
                                        # for sure on a chip, independent of the septem row (indep. of y)
  let sepCoord = pix.chpPixToSeptemPix(chipNumber, realLayout = true)
  result = if isX: sepCoord.x
           else: sepCoord.y

proc chpPixToSeptemPix*(a: int, isX: bool, chipNumber: range[0 .. 6]): int =
  ## Note: This is not particularly efficient... But doing it better requires
  ## a rewrite of the templates above.
  let pix = if isX: (x: a, y: 0, ch: 0)
            else: (x: 129, y: a, ch: 0) # 129 is a dummy that corresponds to a point that is
                                        # for sure on a chip, independent of the septem row (indep. of y)
  let sepCoord = pix.chpPixToSeptemPix(chipNumber, realLayout = false)
  result = if isX: sepCoord.x
           else: sepCoord.y

proc determineRealChip*(p: PixInt, allowOutsideChip = false): int =
  ## Determines the chip for an input pixel on the real layout.
  ## XXX: why cannot be `let`?
  ## XXX2: can merge this with `determineChip`, no?
  let chipEdges = block:
    var res = newSeq[tuple[left, bottom, right, top: int]]()
    for chip in 0 ..< 7:
      let # replae by `chpPixToSeptemPix` to do regular layout!
        left = chpPixToRealPix(0, isX = true, chip)
        bottom = chpPixToRealPix(0, isX = false, chip)
        right = chpPixToRealPix(255, isX = true, chip)
        top =  chpPixToRealPix(255, isX = false, chip)
      if chip in [5, 6]: # top two chips are inverted in X
        res.add (right, top, left, bottom)
      else:
        res.add (left, bottom, right, top)
    res
  for chip, bounds in chipEdges:
    if p.x >= bounds.left and p.x <= bounds.right and
       p.y >= bounds.bottom and p.y <= bounds.top:
      return chip
  if not allowOutsideChip:
    doAssert false, "The pixel : " & $p & " does not belong to a chip!"
  else:
    result = -1

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
  if true: quit("`isNearCenterOfChip` is broken!")
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
template distance*(x, y: float): float = sqrt(x * x + y * y)

# template which returns pitch converted positions on chip pixel values
# to mm from center of chip
# constants are:
# const NPIX = 256
# const PITCH = 0.055 (see ingrid_types)
func applyPitchConversion*[T: (float | SomeInteger)](x, y: T, npix: int): (float, float) =
  ## template which returns the converted positions on a Timepix
  ## pixel position --> absolute position from pixel center in mm.
  ## Note that the x axis is 'inverted'! This follows Christoph's code.
  ## It essentially switches the data from a 'top down view' onto the detector
  ## to a "camera-like" view 'through' the detector. I prefer the former, but
  ## I only much thought about this in the context of the limit calculation.
  ##
  ## Note: this should really use `npix - 1`. `x` can be maximum 255, so we have
  ## a shift off by 1 here. But this would require rerunning *everything* for a
  ## absolute tiny difference. Hence I'll leave it for now.
  ## Also: `PITCH` is currently `0.055`, but more correct would be something that
  ## really uses the Timepix size.
  ((float(npix) - float(x) + 0.5) * PITCH, (float(y) + 0.5) * PITCH)

func inRegion*(centerX, centerY: float, region: ChipRegion): bool {.inline.} =
  ## returns the result of a cut on a certain chip `region`. Inputs the
  ## `centerX` and `centerY` position of a cluster and returns true if
  ## the cluster is within the region
  const centerChip = TimepixSize / 2.0
  case region
  of crGold:
    # make sure this is only initialized once somehow...
    let regCut = getRegionCut(region)
    result = centerX >= regCut.xMin and
             centerX <= regCut.xMax and
             centerY >= regCut.yMin and
             centerY <= regCut.yMax
  of crAll:
    # return true `iff` the coordinate is actually a valid coordinate on the chip!
    result = centerX >= 0.0 and centerX <= TimepixSize and
             centerY >= 0.0 and centerY <= TimepixSize
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

import unchained
defUnit(cm², toExport = true)
proc areaOf*(region: ChipRegion): cm² =
  case region
  of crGold: result = pow(0.95 - 0.45, 2).cm² # area of gold region!
  of crSilver: result = (π * (4.5.mm * 4.5.mm)).to(cm²)
  of crBronze: result = (π * (5.5.mm * 5.5.mm)).to(cm²)
  of crAll: result = (14.1.mm * 14.1.mm).to(cm²)

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
                               pos_x, pos_y, rot_angle: float,
                               useRealLayout = false): ClusterGeometry =
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
  let rotAngle = if not useRealLayout: -rot_angle # for single chip we invert `x ↦ -x`. Thus need to invert angle too.
                 else: rotAngle
  for p in cluster:
    when T is Pix or T is PixTpx3:
      let (x, y) = applyPitchConversion(p.x, p.y, NPIX) # `useRealLayout` has no point here
    elif T is PixInt or T is PixIntTpx3:
      var x, y: float
      if useRealLayout: ## Assumes input is already in real septemboard (pixel) coordinates
        (x, y) = (toRealXPos(p.x), toRealYPos(p.y))
      else:
        (x, y) = applyPitchConversion(p.x, p.y, NPIX * 3)
    else:
      error("Invalid type: " & $T)
    xRot[i] = cos(rotAngle) * (x - pos_x) - sin(rotAngle) * (y - pos_y)
    yRot[i] = sin(rotAngle) * (x - pos_x) + cos(rotAngle) * (y - pos_y)

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
  let xRms  = stat_x.standardDeviation() # this way we simply take the longer axis, as a mirroring
  let yRms  = stat_y.standardDeviation() # does not matter
  let takeX = if xRms >= yRms: true # x is the long axis
              else: false          # y is the long axis
  # doAssert takeX, "Xr > Yr? " & $(xRms > yRms) & " ? " & $xRms & " vs " & $yRms & " pixels ? " & $xRot & " vs " & $yRot & " for rot angle\n " & $rotAngle & " \n cluster " & $cluster
  template takeIt(x, y, takeX: untyped): untyped =
    if takeX: x else: y
  let xLong = max(xRot) - min(xRot)     # helpers are there because in septem board data
  let yLong = max(yRot) - min(yRot)     # using real layout, the axes are sometimes switched.
  let xSkew = stat_x.skewness()
  let ySkew = stat_y.skewness()
  let xKurt = stat_x.kurtosis()
  let yKurt = stat_y.kurtosis()
  # for each of the variables now take the correct x / y variable depending
  # on which axis was the longer one
  result.length               = takeIt(xLong, yLong, takeX)
  result.width                = takeIt(xLong, yLong, not takeX)
  result.rmsLongitudinal      = takeIt(xRms ,  yRms, takeX)
  result.rmsTransverse        = takeIt(xRms ,  yRms, not takeX)
  result.skewnessLongitudinal = takeIt(xSkew, ySkew, takeX)
  result.skewnessTransverse   = takeIt(xSkew, ySkew, not takeX)
  result.kurtosisLongitudinal = takeIt(xKurt, yKurt, takeX)
  result.kurtosisTransverse   = takeIt(xKurt, yKurt, not takeX)
  # We don't care about the sign of the angle. 0-180° covers all angle information
  # we deduce based on nlopt fit anyway! This might be changed if we ever try to
  # deduce the directionality of a track!
  # This is important because our `rotAngle` modification above makes positive angles negative
  result.rotationAngle        = abs(rotAngle)
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

proc calcToAGeometry*[T: SomePix](cluster: var ClusterObject[T]): ToAGeometry =
  ## Given a cluster, computes different ToA based "geometric" properties (i.e.
  ## the length in ToA etc.) and also (hence `cluster` is `var`) modifies the
  ## `toa` field such that each cluster starts at 0.

  ## XXX: `toaLength` is already computed in `raw_data_manipulation` as `length` field of
  ## the `ProcessedRun`!
  var minToA = uint16.high
  var maxToA = 0'u16
  for i, toa in cluster.toa:
    minToA = min(minToA, toa)
    maxToA = max(maxToA, toa)
  ## use min ToA knowledge to push subtracted values to stat and modify `toa`
  var
    stat: RunningStat
  for i, toa in mpairs(cluster.toa):
    let toaZ = toa.int - minToA.int
    stat.push(toaZ.float)
    doAssert toaZ >= 0
    toa = toaZ.uint16
  ## Cannot safely treat `toaLength` as uint16 due to underflow danger
  result.toaLength = (maxToA.float - minToA.float)
  result.toaMean = stat.mean()
  result.toaRms = stat.standardDeviation()
  result.toaSkewness = stat.skewness()
  result.toaKurtosis = stat.kurtosis()
  result.toaMin = minToA

proc wrapDbscan(p: Tensor[float], eps: float, minSamples: int): seq[int] =
  ## This is a wrapper around `dbscan`. Without it for some reason `seqmath's` `arange`
  ## is bound in the context of the `kdtree` code for some reason (binding manually is
  ## no help, neither here nor in `reconstruction.nim`)
  dbscan(p, eps, minSamples)

proc findClusterDBSCAN*[T: SomePix](pixels: seq[T], eps: float = 65.0,
                                    minSamples: int = 3): seq[Cluster[T]] =
  var pT = newTensorUninit[float]([pixels.len, 2])
  for i, tup in pixels:
    pT[i, _] = [tup.x.float, tup.y.float].toTensor.unsqueeze(axis = 0)
  if pixels.len == 0: return

  let clusterIdxs = wrapDbscan(pT, eps, minSamples)
  for i, clIdx in clusterIdxs:
    if clIdx == -1: continue
    if clIdx >= result.len:
      result.setLen(clIdx+1)
    if result[clIdx].len == 0:
      result[clIdx] = newSeqOfCap[T](pixels.len)
    result[clIdx].add pixels[i]

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


#proc recoCluster*(c: Cluster[Pix]): ClusterObject[Pix] {.gcsafe.} =
proc recoCluster*[T: SomePix; U: SomePix](c: Cluster[T],
                                          timepix: TimepixVersion = Timepix1,
                                          _: typedesc[U],
                                          useRealLayout = false
                                         ): ClusterObject[U] {.gcsafe, hijackMe.} =
  result = newClusterObject[U](timepix)

  let clustersize: int = len(c)
  ##
  const NeedConvert = T is PixTpx3 and U is Pix or
                      T is PixIntTpx3 and U is PixInt

  var cl = newSeq[U](clustersize)
  var
    sum_x, sum_x2: int
    sum_y, sum_y2, sum_xy: int
    sum_ToT: int
    minToA: uint16
  when T is PixTpx3:
    result.toa = newSeq[uint16](clustersize)
    result.toaCombined = newSeq[uint64](clustersize)
  for i in 0 ..< clustersize:
    let ci = c[i]
    sum_x  += ci.x.int
    sum_y  += ci.y.int
    sumToT += ci.ch.int
    sum_x2 += ci.x.int * ci.x.int
    sum_y2 += ci.y.int * ci.y.int
    sum_xy += ci.x.int * ci.y.int
    when NeedConvert:
      cl[i] = (x: ci.x, y: ci.y, ch: ci.ch)
      result.toa[i] = ci.toa
      result.toaCombined[i] = ci.toaCombined
  when NeedConvert:
    result.data = cl
  else:
    result.data = c
  let
    pos_x = float64(sum_x) / float64(clustersize)
    pos_y = float64(sum_y) / float64(clustersize)
  var
    rms_x = sqrt(float64(sum_x2) / float64(clustersize) - pos_x * pos_x)
    rms_y = sqrt(float64(sum_y2) / float64(clustersize) - pos_y * pos_y)
    rotAngleEstimate = arctan( (float64(sum_xy) / float64(clustersize)) -
                               pos_x * pos_y / (rms_x * rms_x))

  # set the total "charge" in the cluster (sum of ToT values), can be
  # converted to electrons with ToT calibration
  result.sum_tot = sumTot
  # set number of hits in cluster
  result.hits = clustersize
  # set the position
  when T is Pix or T is PixTpx3:
    (result.centerX, result.centerY) = applyPitchConversion(pos_x, pos_y, NPIX)
  elif T is PixInt or T is PixIntTpx3:
    if useRealLayout:
      (result.centerX, result.centerY) = (toRealXPos(pos_x), toRealYPos(pos_y))
    else:
      (result.centerX, result.centerY) = applyPitchConversion(pos_x, pos_y, NPIX * 3)
  else:
    error("Invalid type: " & $T)
  # prepare rot angle fit
  if rotAngleEstimate < 0:
    #echo "correcting 1"
    rotAngleEstimate += 8 * arctan(1.0)
  if rotAngleEstimate > 4 * arctan(1.0):
    #echo "correcting 2"
    rotAngleEstimate -= 4 * arctan(1.0)
  elif classify(rotAngleEstimate) != fcNormal:
    warn "Rot angle estimate is NaN, vals are ", $rms_x, " ", $rms_y, " from data: ", $c
    # what do we do in this case with the geometry?!
    #raise newException(ValueError, "Rotation angle estimate returned bad value")
    warn "Fit will probably fail!"

  # else we can minimize the rotation angle and calc the eccentricity
  let (rot_angle, eccentricity) = fitRotAngle(result, rotAngleEstimate)

  # now we still need to use the rotation angle to calculate the different geometric
  # properties, i.e. RMS, skewness and kurtosis along the long axis of the cluster
  result.geometry = calcGeometry(c, result.centerX, result.centerY, rot_angle, useRealLayout)
  when T is PixTpx3:
    result.toaGeometry = calcToAGeometry(result)


proc getPixels[T; U](dat: RecoInputEvent[U], _: typedesc[T],
                     useRealLayout = false): seq[T] =
  ## Returns the pixel data of the input event taking into account a
  ## possible conversion from Tpx1 -> Tpx3 data or transformation to the
  ## real Septemboard layout including spacing.
  when T is U:
    when T is PixInt:
      if false: #  useRealLayout and T is PixInt:
        ## XXX: WE GENERALLY DON'T WANT THIS.
        ## If we were to use `getPixels` with `useRealLayout` from `likelihood` for the septem veto
        ## this is *NOT* good. We want to perform the clustering logic *on the tight septemboard layout*
        ## and *NOT* on the real layout. Irrespective of whether the user wants the real layout.
        ## If anything this should become a *SEPARATE* option from the `likelihood` septem veto option!
        result = newSeq[PixInt](dat.pixels.len)
        for i in 0 ..< dat.pixels.len:
          ## XXX: use `tightToReal`?
          result[i] = septemPixToRealPix(dat.pixels[i])
      else:
        result = dat.pixels
    else:
      result = dat.pixels
  elif T is PixTpx3:
    doAssert dat.pixels.len == dat.toa.len
    result = newSeq[PixTpx3](dat.pixels.len)
    for i in 0 ..< result.len:
      result[i] = (x: dat.pixels[i].x, y: dat.pixels[i].y, ch: dat.pixels[i].ch,
                   toa: dat.toa[i], toaCombined: dat.toaCombined[i])
  else:
    error("Invalid type : " & $T)

proc recoEvent*[T: SomePix](dat: RecoInputEvent[T],
                            chip, run, searchRadius: int,
                            dbscanEpsilon: float,
                            clusterAlgo: ClusteringAlgorithm,
                            timepixVersion = Timepix1,
                            useRealLayout = false): RecoEvent[T] {.gcsafe, hijackMe.} =
  result.event_number = dat.eventNumber
  result.chip_number = chip

  ## NOTE: The usage of a `PixTpx3` is rather wasteful and is only done, because
  ## it's currently easier to perform the clustering using the default algorithm
  ## with such data. Otherwise we need to keep track which indices end up in what
  ## cluster.
  ## However: it may anyway be smart to avoid the logic of `deleteIntersection` and
  ## instead mark indices in a set (?) and keep track of each index being in what
  ## cluster?
  ## I remember trying *some* kind of set based approach before which turned out
  ## slower, so gotta be careful.

  template recoClusterTmpl(typ, pixels: untyped): untyped {.dirty.} =
    var cluster: seq[Cluster[typ]]
    case clusterAlgo
    of caDefault: cluster = findSimpleCluster(pixels, searchRadius)
    of caDBSCAN:  cluster = findClusterDBSCAN(pixels, dbscanEpsilon)
    result.cluster = newSeq[ClusterObject[T]](cluster.len)
    for i, cl in cluster:
      when typ is PixInt: ## We convert pixels to the real layout *now*
        let mcl = if useRealLayout: septemPixToRealPix(cl) else: cl
        result.cluster[i] = recoCluster(mcl, timepixVersion, T, useRealLayout)
      else:
        result.cluster[i] = recoCluster(cl, timepixVersion, T, useRealLayout)

  if dat[0].len > 0:
    case timepixVersion
    of Timepix1:
      when T is PixInt or T is PixIntTpx3:
        let pixels = getPixels(dat, PixInt, useRealLayout)
        recoClusterTmpl(PixInt, pixels)
      else:
        let pixels = getPixels(dat, Pix, useRealLayout)
        recoClusterTmpl(Pix, pixels)
    of Timepix3:
      when T is PixIntTpx3 or T is PixInt:
        let pixels = getPixels(dat, PixIntTpx3, useRealLayout)
        recoClusterTmpl(PixIntTpx3, pixels)
      else:
        let pixels = getPixels(dat, PixTpx3, useRealLayout)
        recoClusterTmpl(PixTpx3, pixels)
