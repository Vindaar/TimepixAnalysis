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

proc feSpectrumFuncCharge(xOb: PyObject, p0, p1, p2, p3, p4, p5: float): seq[float] {.exportpy.} =
  ## the fitting function used for the Fe spectrum with charge bins, see `XrayCalib.c`
  # fix N_Kbeta/N_Kalpha to theoretical value from XDB
  # TODO: allow conversion of PyObject to seq[float]!!!
  #let x = xOb.to(float)
  result = newSeq[float]()
  for xPy in xOb:
    let x = xPy.to(float)
    let
      p6 = 17.0 / 150.0
      t1 = p0 * exp(-(x - p1) * (x - p1) / (2.0 * p2 * p2))
      t2exp = -(x - (p1 * 3.53 / 2.94)) * (x - (p1 * 3.53/2.94)) / (2.0 * p2 * p2)
      t2 = p0 * p6 * exp(t2exp)
      t3 = p3 * exp(-(x - p4) * (x - p4) / (2.0 * p5 * p5))
      t4exp = -(x - (p4 * 6.49 / 5.90)) * (x - (p4 * 6.49 / 5.90)) / (2.0 * p5 * p5)
      t4 = p3 * p6 * exp(t4exp)
    result.add (t1 + t2 + t3 + t4)

proc linearFunc(xOb: PyObject, a: float): seq[float] {.exportpy.} =
  #let ss = xOb.to(seq[float])
  # TODO: make it possible to convert PyObject to seq
  #result = newSeqOfCap[float](ss.len)
  result = newSeq[float]()
  for xPy in xOb:
    let x = xPy.to(float)
    result.add a * x
