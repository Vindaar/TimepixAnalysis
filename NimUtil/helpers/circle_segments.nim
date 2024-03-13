import math

proc circleSegment(r, ϑ: float): float =
  ## Computes the area of a circle segment. That is the area cut off by
  ## a line that cuts the circle at two points.
  ##
  ## `r` is the radius of the circle.
  ##
  ## `ϑ` is the angle between the two points that the line cuts the circle.
  result = r * r / 2.0 * (ϑ - sin(ϑ))

proc areaTriangle(a, b: float): float =
  ## Area of a right triangle for `a, b` the legs.
  result = a * b / 2.0

proc areaCircleTwoLinesCut*(R, r1, r2: float): float =
  ## Computes the area of a circle left after two lines cut the circle, which
  ## are orthogonal at radii `r1, r2` from the center.
  ##
  ## `R` is the radius of the circle.
  template compSegment(r: float): (float, float, float) =
    let rp = R - r
    if rp < 0.0:
      (0.0, 0.0, 0.0)
    else:
      let ϑ = arccos(r / R) * 2.0
      let A = circleSegment(R, ϑ)
      (rp, ϑ, A)

  # compute first and second segment
  let (rp1, ϑ1, A) = compSegment(r1)
  let (rp2, ϑ2, B) = compSegment(r2)
  # compute third circle segment (part that is subtracted twice)
  let α = (ϑ1 + ϑ2 - PI) / 2.0
  let D = circleSegment(R, α)
  # compute triangle that is removed twice
  let C = areaTriangle(rp1, rp2)
  # put it all together
  #dump α
  #dump A
  #dump B
  #dump C
  #dump D
  if α > 0.0:
    result = PI * R * R - A - B + C + D
  else:
    # else the two lines do ``not`` cut part of the same
    result = PI * R * R - A - B

when isMainModule:
  echo areaCircleTwoLinesCut(50.0, 35.0, 35.0)
  echo "Compare : ", PI * 50.0 * 50.0

  echo areaCircleTwoLinesCut(50.0, 49.0, 49.0)
  echo "Compare : ", PI * 50.0 * 50.0
