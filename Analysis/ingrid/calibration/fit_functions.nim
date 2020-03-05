import math, seqmath, sequtils
import ../cdlFitting/cdlFitMacro

func sCurveFunc*(p: seq[float], x: float): float =
  ## we fit the complement of a cumulative distribution function
  ## of the normal distribution
  # parameter p[2] == sigma
  # parameter p[1] == x0
  # parameter p[0] == scale factor
  result = normalCdfC(x, p[2], p[1]) * p[0]

func linearFunc*(p: seq[float], x: float): float =
  result = p[0] + x * p[1]

proc thlCalibFunc*(p: seq[float], x: float): float =
  ## we fit a linear function to the charges and mean thl values
  ## of the SCurves
  linearFunc(p, x)

func totCalibFunc*(p: seq[float], x: float): float =
  ## we fit a combination of a linear and a 1 / x function
  ## The function is:
  ## ToT[clock cycles] = a * x + b - (c / (x - t))
  ## where x is the test pulse height in mV and:
  ## p = [a, b, c, t]
  result = p[0] * x + p[1] - p[2] / (x - p[3])

proc expGauss*(p: seq[float], x: float): float =
  # exponential * times (?!) from Christoph's expogaus.c
  if len(p) != 5:
    return Inf
  let p_val = 2.0 * (p[1] * pow(p[4], 2.0) - p[3])
  let q_val = 2.0 * pow(p[4], 2.0) * p[0] + pow(p[3], 2.0) - ln(p[2]) * 2.0 * pow(p[4], 2.0)

  let threshold = - p_val / 2.0 - sqrt( pow(p_val, 2.0) / 4.0 - q_val )

  if x < threshold:
    result = exp(p[0] + p[1] * x)
  else:
    result = p[2] * gauss(x, p[3], p[4])

## Define the `feSpectrumFunc`. Total parameters: 15.
## 14 from 4 ffExpGauss and 1 additional (`p_ar[14]`).
declareFitFunc(feSpectrum):
  ffExpGauss: "Mn-Kalpha-esc"
  ffExpGauss:
    name = "Mn-Kbeta-esc"
    eN = eN("Mn-Kalpha-esc") * p_ar[14] # p_ar[14] is an additional fit parameter
    emu = emu("Mn-Kalpha-esc") * 3.5 / 2.9 # lock to relation to `Mn-Kalpha-esc` arg
    es = es("Mn-Kalpha-esc") # lock to `es` of `Mn-Kalpha`
  ffExpGauss: "Mn-Kalpha"
  ffExpGauss:
    name = "Mn-Kbeta"
    eN = eN("Mn-Kalpha") * p_ar[14]
    emu = emu("Mn-Kalpha") * 6.35 / 5.75 # lock to relation to `Mn-Kalpha` arg
    es = es("Mn-Kalpha") # lock to `es` of `Mn-Kalpha`

func polyaImpl*(p: seq[float], x: float): float =
  ## Polya function to fit to TOT histogram / charge in electrons of a
  ## run. This is the actual implementation of the polya distribution.
  ## Parameters:
  ## N     = p[0]    scaling factor
  ## G     = p[1]    gas gain
  ## theta = p[2]    parameter, which describes distribution (?! I guess it makes sens
  ##                 since we take its power and it enters gamma)
  let
    thetaDash = p[2] + 1
    coeff1 = (p[0] / p[1]) * pow((thetaDash), thetaDash) / gamma(thetaDash)
    coeff2 = pow((x / p[1]), p[2]) * exp(-thetaDash * x / p[1])
  result = coeff1 * coeff2

type
  FitObject* = object
    x*: seq[float]
    y*: seq[float]
    yErr*: seq[float]

template fitForNlopt*(name, funcToCall: untyped): untyped =
  proc `name`*(p: seq[float], fitObj: FitObject): float =
    # TODO: add errors, proper chi^2 intead of unscaled values,
    # which results in huge chi^2 despite good fit
    let x = fitObj.x
    let y = fitObj.y
    let yErr = fitObj.yErr
    var fitY = x.mapIt(`funcToCall`(p, it))
    var diff = newSeq[float](x.len)
    result = 0.0
    for i in 0 .. x.high:
      diff[i] = (y[i] - fitY[i]) / yErr[i]
      result += pow(diff[i], 2.0)
    result = result / (x.len - p.len).float

template fitForNloptLnLikelihood*(name, funcToCall: untyped): untyped =
  proc `name`*(p: seq[float], fitObj: FitObject): float =
    ## Maximum Likelihood Estimator
    ## the Chi Square of a Poisson distributed log Likelihood
    ## Chi^2_\lambda, P = 2 * \sum_i y_i - n_i + n_i * ln(n_i / y_i)
    ## n_i = number of events in bin i
    ## y_i = model prediction for number of events in bin i
    ## derived from likelihood ratio test theorem
    let x = fitObj.x
    let y = fitObj.y
    var fitY = x.mapIt(`funcToCall`(p, it))
    result = 0.0
    for i in 0 ..< x.len:
      if fitY[i].float > 0.0 and y[i].float > 0.0:
        result = result + (fitY[i] - y[i] + y[i] * ln(y[i] / fitY[i]))
      else:
        result = result + (fitY[i] - y[i])
      # ignore empty data and model points
    result = 2 * result

template fitForNloptLnLikelihoodGrad*(name, funcToCall: untyped): untyped =
  proc `name`*(p: seq[float], fitObj: FitObject): (float, seq[float]) =
    ## Maximum Likelihood Estimator
    ## the Chi Square of a Poisson distributed log Likelihood
    ## Chi^2_\lambda, P = 2 * \sum_i y_i - n_i + n_i * ln(n_i / y_i)
    ## n_i = number of events in bin i
    ## y_i = model prediction for number of events in bin i
    ## derived from likelihood ratio test theorem
    # NOTE: do not need last gradients
    let xD = fitObj.x
    let yD = fitObj.y
    var gradRes = newSeq[float](p.len)
    var res = 0.0
    var h: float
    proc fnc(x, y, params: seq[float]): float =
      let fitY = x.mapIt(`funcToCall`(params, it))
      for i in 0 ..< x.len:
        if fitY[i].float > 0.0 and y[i].float > 0.0:
          result = result + (fitY[i] - y[i] + y[i] * ln(y[i] / fitY[i]))
        else:
          result = result + (fitY[i] - y[i])
    res = 2 * fnc(xD, yD, p)
    for i in 0 ..< gradRes.len:
      h = p[i] * sqrt(epsilon(float64))
      var
        modParsUp = p
        modParsDown = p
      modParsUp[i] = p[i] + h
      modParsDown[i] = p[i] - h
      gradRes[i] = (fnc(xD, yD, modParsUp) - fnc(xD, yD, modParsDown)) / (2.0 * h)
      # ignore empty data and model points
    result = (res, gradRes)

#fitForNloptLnLikelihoodGrad(polya, polyaImpl)
#fitForNloptLnLikelihood(polya, polyaImpl)
fitForNlopt(polya, polyaImpl)
