import math
import ../ingrid_types

################################################################################
############# Geometry calculation related procs ###############################
################################################################################

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
