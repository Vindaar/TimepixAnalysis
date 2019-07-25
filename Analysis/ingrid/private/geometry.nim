import math
import ../ingrid_types
import cdl_cuts

################################################################################
############# Geometry calculation related procs ###############################
################################################################################

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
