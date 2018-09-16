import nimpy, math

proc polyaPython(xOb: PyObject, p0, p1, p2: float): seq[float] {.exportpy.} =
  ## Polya function to fit to TOT histogram / charge in electrons of a
  ## run. This is the actual implementation of the polya distribution.
  ## Parameters:
  ## N     = p[0]    scaling factor
  ## G     = p[1]    gas gain
  ## theta = p[2]    parameter, which describes distribution (?! I guess it makes sens
  ##                 since we take its power and it enters gamma)
  var thetaDash = p2 + 1
  result = newSeqOfCap[float](150)
  # TODO: implement this into `nimpy`!
  let coeff1 = (p0 / p1) * pow((thetaDash), thetaDash) / gamma(thetaDash)
  var coeff2: float
  for x in xOb:
    let xNim = x.to(float)
    coeff2 = pow((xNim / p1), p2) * exp(-thetaDash * xNim / p1)
    result.add coeff1 * coeff2
