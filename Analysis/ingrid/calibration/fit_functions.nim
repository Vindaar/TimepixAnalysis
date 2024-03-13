import math, seqmath, sequtils, fenv, macros
import ../cdlFitting/cdlFitMacro

type
  NloptFitKind* = enum
    nfChiSq, nfChiSqGrad, nfMle, nfMleGrad

func sCurveFunc*(p: seq[float], x: float): float =
  ## we fit the complement of a cumulative distribution function
  ## of the normal distribution
  # parameter p[2] == sigma
  # parameter p[1] == x0
  # parameter p[0] == scale factor
  result = normalCdfC(x, p[2], p[1]) * p[0]

func linearFunc*(p: seq[float], x: float): float =
  result = p[0] + x * p[1]

func linearFuncNoOffset*(p: seq[float], x: float): float =
  result = x * p[0]

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

## Define the `feSpectrumFuncCharge`. Total parameters: .
## 14 from 4 ffExpGauss and 1 additional (`p_ar[14]`).
declareFitFunc(feSpectrumCharge):
  ffGauss: "Mn-Kalpha-esc"
  ffGauss: "Mn-Kalpha"
  ffGauss:
    name = "Mn-Kbeta-esc"
    gN = gN("Mn-Kalpha-esc") * (17.0 / 150.0)# lock to Kalpha escape peak
    gmu = gmu("Mn-Kalpha-esc") * (3.53 / 2.94)
    gs = gs("Mn-Kalpha-esc") # lock to Kalpha escape peak
  ffGauss:
    name = "Mn-Kbeta"
    gN = gN("Mn-Kalpha") * (17.0 / 150.0) # lock to Kalpha escape peak
    gmu = gmu("Mn-Kalpha") * (6.49 / 5.90)
    gs = gs("Mn-Kalpha") # lock to Kalpha escape peak


func polyaImpl*(p: seq[float], x: float): float =
  ## Polya function to fit to TOT histogram / charge in electrons of a
  ## run. This is the actual implementation of the polya distribution.
  ## Parameters:
  ## N     = p[0]    scaling factor
  ## G     = p[1]    gas gain
  ## theta = p[2]    parameter, which describes distribution (?! I guess it makes sens
  ##                 since we take its power and it enters gamma)
  ##
  ## Description goes back to
  ## `Statisitcs of electron avalanches and ultimate resolution of proportional counters`
  ## by Alkhazov, 1970.
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

import ggplotnim, os
template chiSquareNoGrad(funcToCall: untyped): untyped {.dirty.} =
  ## Chi^2 estimator for no gradient based algorithms
  # TODO: add errors, proper chi^2 intead of unscaled values,
  # which results in huge chi^2 despite good fit
  let x = fitObj.x
  let y = fitObj.y
  let yErr = fitObj.yErr
  var diff = newSeq[float](x.len)
  result = 0.0

  #let df = toDf({ "x" : x,
  #                    "y" : y,
  #                    "yFit" : fitY })
  #ggplot(df, aes("x", "y")) +
  #  geom_histogram(stat = "identity") +
  #  geom_line(aes(y = "yFit")) +
  #  ggsave("/tmp/fit_progress.pdf")

  for i in 0 .. x.high:
    let fitY = funcToCall(p, x[i])
    diff[i] = (y[i] - fitY) / yErr[i]
    result += pow(diff[i], 2.0)
  # return `Chi^2` (and not reduced `Chi^2`!
  result = result

template chiSquareGrad(funcToCall: untyped): untyped {.dirty.} =
  ## Chi^2 estimator for gradient based algorithms
  # TODO: add errors, proper chi^2 intead of unscaled values,
  # which results in huge chi^2 despite good fit
  let xD = fitObj.x
  let yD = fitObj.y
  let yErr = fitObj.yErr
  var gradRes = newSeq[float](p.len)
  var res = 0.0
  var h: float

  proc fnc(x, y, params: seq[float]): float =
    var diff: float
    for i in 0 ..< x.len:
      let fitY = funcToCall(params, x[i])
      diff = (y[i] - fitY) / yErr[i]
      result += diff*diff
  res = fnc(xD, yD, p)
  for i in 0 ..< gradRes.len:
    h = if p[i] != 0.0: p[i] * 1e-12 else: 1e-12 # sqrt(epsilon(float64)) else: sqrt(epsilon(float64))
    var
      modParsUp = p
      modParsDown = p
    modParsUp[i] = p[i] + h
    modParsDown[i] = p[i] - h
    gradRes[i] = (fnc(xD, yD, modParsUp) - fnc(xD, yD, modParsDown)) / (2.0 * h)
    # ignore empty data and model points
  # return `Chi^2` (and not reduced `Chi^2`!
  result = (res, gradRes)

template mleLnLikelihoodNoGrad(funcToCall: untyped): untyped {.dirty.} =
  ## Maximum Likelihood Estimator
  ## the Chi Square of a Poisson distributed log Likelihood
  ## Chi^2_\lambda, P = 2 * \sum_i y_i - n_i + n_i * ln(n_i / y_i)
  ## n_i = number of events in bin i
  ## y_i = model prediction for number of events in bin i
  ## derived from likelihood ratio test theorem
  let x = fitObj.x
  let y = fitObj.y

  #let df = toDf({ "x" : x,
  #                    "y" : y,
  #                    "yFit" : fitY })
  #ggplot(df, aes("x", "y")) +
  #  geom_histogram(stat = "identity") +
  #  geom_line(aes(y = "yFit")) +
  #  ggsave("/tmp/fit_progress_mle_nograd.pdf")
  #sleep(50)

  result = 0.0
  for i in 0 ..< x.len:
    let fitY = funcToCall(p, x[i])
    if fitY.float > 0.0 and y[i].float > 0.0:
      result = result + (fitY - y[i] + y[i] * ln(y[i] / fitY))
    else:
      result = result + (fitY - y[i])
    # ignore empty data and model points
  result = 2 * result

template mleLnLikelihoodGrad(funcToCall: untyped): untyped {.dirty.} =
  ## Maximum Likelihood Estimator for gradient based algorithms
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
    #let df = toDf({ "x" : x,
    #                    "y" : y,
    #                    "yFit" : fitY })
    #ggplot(df, aes("x", "y")) +
    #  geom_histogram(stat = "identity") +
    #  geom_line(aes(y = "yFit")) +
    #  ggsave("/tmp/fit_progress_mle_grad.pdf")
    #sleep(50)
    for i in 0 ..< x.len:
      let fitY = funcToCall(params, x[i])
      if fitY.float > 0.0 and y[i].float > 0.0:
        result = result + (fitY - y[i] + y[i] * ln(y[i] / fitY))
      else:
        result = result + (fitY - y[i])
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

macro fitForNlopt*(name, funcToCall: untyped,
                   nfKind: static NloptFitKind,
                   toExport: static bool = true): untyped =
  let pArg = ident"p"
  let fobj = ident"fitObj"
  var body: NimNode
  var retType: NimNode
  case nfKind
  of nfChiSq:
    body = getAst(chiSquareNoGrad(funcToCall))
    retType = ident"float"
  of nfChiSqGrad:
    body = getAst(chiSquareGrad(funcToCall))
    retType = nnkPar.newTree(ident"float",
                             nnkBracketExpr.newTree(ident"seq",
                                                    ident"float"))
  of nfMle:
    body = getAst(mleLnLikelihoodNoGrad(funcToCall))
    retType = ident"float"
  of nfMleGrad:
    body = getAst(mleLnLikelihoodGrad(funcToCall))
    retType = nnkPar.newTree(ident"float",
                             nnkBracketExpr.newTree(ident"seq",
                                                    ident"float"))
  if toExport:
    result = quote do:
      proc `name`*(`pArg`: seq[float], `fobj`: FitObject): `retType` =
        `body`
  else:
    result = quote do:
      proc `name`(`pArg`: seq[float], `fobj`: FitObject): `retType` =
        `body`
  #echo result.treeRepr
  #echo result.repr

#fitForNloptLnLikelihoodGrad(polya, polyaImpl)
#fitForNloptLnLikelihood(polya, polyaImpl)
fitForNlopt(polya, polyaImpl, nfKind = nfChiSq, toExport = true)
