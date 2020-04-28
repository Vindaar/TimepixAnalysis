#[
Contains all logic that handles the fitting routines.
This module does ``not`` interact with HDF5 files. All routines take
data perform a fit either with Nlopt or mpfit and return data.
The actual functions to be fitted are found in `fit_functions.nim`
]#

import seqmath, sequtils, stats, strformat, fenv
import nlopt, mpfit

import .. / ingrid_types
import fit_functions
import .. / cdlFitting / cdlFitMacro

################################################################################
########## Timepix calibration related fitting routines ########################
################################################################################

const
  NTestPulses = 1000.0
  ScalePulses = 1.01
  CurveHalfWidth = 15

func findDrop(thl: seq[int], count: seq[float]): (float, int, int) =
  ## search for the drop of the SCurve, i.e.
  ##
  ##    ___/\__
  ## __/       \__
  ## 0123456789ABC
  ## i.e. the location of feature A (the center of it)
  ## returns the thl value at the center of the drop
  ## and the ranges to be used for the fit (min and max
  ## as indices for the thl and count seqs)
  # zip thl and count
  let
    mthl = thl.mapIt(it.float)
    mcount = count.mapIt(it.float)
    thlZipped = zip(toSeq(0 .. mthl.high), mthl)
    thlCount = zip(thlZipped, mcount)
    # helper, which is the scaled bound to check whether we need
    # to tighten the view on the valid THLs
    scaleBound = NTestPulses * ScalePulses

  # extract all elements  of thlCount where count larger than 500, return
  # the largest element
  let
    drop = thlCount.filterIt(it[1] > (NTestPulses / 2.0))[^1]
    (pCenterInd, pCenter) = drop[0]

  var minIndex = max(pCenterInd - CurveHalfWidth, 0)
  # test the min index
  if mcount[minIndex.int] > scaleBound:
    # in that case get the first index larger than minIndex, which
    # is smaller than the scaled test pulses
    let
      thlView = thlCount[minIndex.int .. thlCount.high]
      thlSmallerBound = thlView.filterIt(it[1] < scaleBound)
    # from all elements smaller than  the bound, extract the index,
    # [0][0][0]:
    # - [0]: want element at position 0, contains smallest THL values
    # - [0]: above results in ( (`int`, `float`), `float` ), get first tupl
    # - [0]: get `int` from above, which corresponds to indices of THL
    minIndex = thlSmallerBound[0][0][0]
  let maxIndex = min(pCenterInd + CurveHalfWidth, thl.high)

  # index is thl of ind
  result = (pCenter, minIndex, maxIndex)


proc fitToTCalib*(tot: Tot, startFit = 0.0): FitResult =
  var
    # local mutable variables to potentially remove unwanted data for fit
    mPulses = tot.pulses.mapIt(it.float)
    mMean = tot.mean
    mStd = tot.std

  # define the start parameters. Use a good guess...
  let p = [0.4, 64.0, 1000.0, -20.0] # [0.149194, 23.5359, 205.735, -100.0]
  #var pLimitBare: mp_par
  var pLimits = @[(l: -Inf, u: Inf),
                 (l: -Inf, u: Inf),
                 (l: -Inf, u: Inf),
                 (l: -100.0, u: 0.0)]

  if startFit > 0:
    # in this case cut away the undesired parameters
    let ind = mPulses.find(startFit) - 1
    if ind > 0:
      mPulses.delete(0, ind)
      mMean.delete(0, ind)
      mStd.delete(0, ind)

  #let p = [0.549194, 23.5359, 50.735, -1.0]
  let (pRes, res) = fit(totCalibFunc,
                        p,
                        mPulses,
                        mMean,
                        mStd,
                        bounds = pLimits) #@[pLimitBare, pLimitBare, pLimitBare, pLimit])
  echoResult(pRes, res = res)

  # the plot of the fit is performed to the whole pulses range anyways, even if
  result.x = linspace(tot.pulses[0].float, tot.pulses[^1].float, 100)
  result.y = result.x.mapIt(totCalibFunc(pRes, it))
  result.pRes = pRes
  result.pErr = res.error
  result.redChiSq = res.reducedChiSq

proc fitSCurve*(curve: SCurve): FitResult =
  ## performs the fit of the `sCurveFunc` to the given `thl` and `hits`
  ## fields of the `SCurve` object. Returns a `FitScurve` object of the fit result
  const pSigma = 5.0
  let
    (pCenter, minIndex, maxIndex) = findDrop(curve.thl, curve.hits)
    err = curve.thl.mapIt(1.0)
    p = @[NTestPulses, pCenter.float, pSigma]

  let
    thlCut = curve.thl[minIndex .. maxIndex].mapIt(it.float)
    countCut = curve.hits[minIndex .. maxIndex].mapIt(it.float)

    (pRes, res) = fit(sCurveFunc,
                      p,
                      thlCut,
                      countCut,
                      err)
  echoResult(pRes, res = res)
  # create data to plot fit as result
  result.x = linspace(curve.thl[minIndex].float, curve.thl[maxIndex].float, 1000)
  result.y = result.x.mapIt(sCurveFunc(pRes, it))
  result.pRes = pRes
  result.pErr = res.error
  result.redChiSq = res.reducedChiSq

proc fitSCurve*[T](thl, count: seq[T], voltage: int): FitResult =
  when T isnot float:
    let
      mthl = mapIt(thl.float)
      mhits = mapIt(count.float)
    let curve = SCurve(thl: mthl, hits: mhits, voltage: voltage)
  else:
    let curve = SCurve(thl: thl, hits: count, voltage: voltage)
  result = fitSCurve(curve)

proc fitThlCalib*(charge, thl, thlErr: seq[float]): FitResult =

  # determine start parameters
  let p = @[0.0, (thl[1] - thl[0]) / (charge[1] - charge[0])]

  echo "Fitting ", charge, " ", thl, " ", thlErr

  let (pRes, res) = fit(thlCalibFunc, p, charge, thl, thlErr)
  # echo parameters
  echoResult(pRes, res = res)

  result.x = linspace(charge[0], charge[^1], 100)
  result.y = result.x.mapIt(thlCalibFunc(pRes, it))
  result.pRes = pRes
  result.pErr = res.error
  result.redChiSq = res.reducedChiSq

################################################################################
########## Fe calibration spectra related fitting routines #####################
################################################################################

proc getLines(hist, binning: seq[float]): (float, float, float, float, float, float) =
  # define center, std and amplitude of K_alpha line
  # as well as escape peak
  doAssert binning.len > 0
  let peakIdx = argMax(hist)
  let halfIdx = toSeq(0 ..< hist.len).filterIt(hist[it] > (hist[peakIdx] / 2.0))[0]
  let mu_Kalpha = binning[peakIdx]

  let sigma_Kalpha = (muKalpha - binning[halfIdx]) * 0.8
  let n_kalpha = hist[peakIdx]

  let ratio = mu_Kalpha / sigma_Kalpha

  let mu_kalpha_esc = mu_kalpha * 2.9 / 5.75
  let sigma_kalpha_esc = mu_kalpha_esc / ratio
  let n_kalpha_esc = n_kalpha / 10.0
  result = (mu_kalpha, sigma_kalpha, n_kalpha, mu_kalpha_esc, sigma_kalpha_esc, n_kalpha_esc)

func getFeSpectrumParams(hist, binning: seq[float]): seq[float] =
  ## Returns suitable start parameters for the `FeSpectrumFunc` fit function
  ## This includes the 4 peaks we fit with fixed parameters and the additional
  ## 15th parameter for the ratio of the alpha to beta peaks.
  let (mu_kalpha,
       sigma_kalpha,
       n_kalpha,
       mu_kalpha_esc,
       sigma_kalpha_esc,
       n_kalpha_esc) = getLines(hist, binning)
  let params = @[
    FitFuncArgs(
      name: "Mn-Kalpha-esc",
      kind: ffExpGauss,
      #ea: 1e-4,
      ea: 0.0,
      #eb: 1e-5,
      eb: 0.0,
      eN: n_kalpha_esc,
      emu: mu_Kalpha * 2.9 / 5.75,
      es: sigma_kalpha_esc),
    FitFuncArgs(
      name: "Mn-Kbeta-esc",
      kind: ffExpGauss,
      #ea: 1e-4,
      ea: 0.0,
      #eb: 1e-5,
      eb: 0.0,
      eN: fixed,
      emu: fixed,
      es: fixed), # additional parameters fixed, `fixed` is just an overload for `NaN`
    FitFuncArgs(
      name: "Mn-Kalpha",
      kind: ffExpGauss,
      #ea: 1e-4,
      ea: 0.0,
      #eb: 1e-5,
      eb: 0.0,
      eN: n_kalpha,
      emu: mu_Kalpha, # since we count single electrons, index equals number electrons!
      es: sigma_Kalpha), # sigma is approxed to 0.8 times the half width
    FitFuncArgs(
      name: "Mn-Kbeta",
      kind: ffExpGauss,
      #ea: 1e-4,
      ea: 0.0,
      #eb: 1e-5,
      eb: 0.0,
      eN: fixed,
      emu: fixed,
      es: fixed) # additional parameters fixed
  ]
  # serialize them to a `seq[float]` params
  result = params.serialize
  # add the final parameter 15 for the NKalpha/NKbeta ratio
  result.add 17.0 / 150.0

func getFeSpectrumChargeParams(hist, binning: seq[float]): seq[float] =
  ## Returns suitable start parameters for the `FeSpectrumChargeFunc` fit function
  ## This includes the 4 peaks we fit with fixed parameters where 2 peaks are completely
  ## fixed.
  let (mu_kalpha,
       sigma_kalpha,
       n_kalpha,
       mu_kalpha_esc,
       sigma_kalpha_esc,
       n_kalpha_esc) = getLines(hist, binning)
  let params = @[
    FitFuncArgs(
      name: "Mn-Kalpha-esc",
      kind: ffGauss,
      gN: n_kalpha_esc,
      gmu: mu_Kalpha * 2.9 / 5.75,
      gs: sigma_kalpha_esc),
    FitFuncArgs(
      name: "Mn-Kbeta-esc",
      kind: ffGauss,
      gN: fixed,
      gmu: fixed,
      gs: fixed), # additional parameters fixed, `fixed` is just an overload for `NaN`
    FitFuncArgs(
      name: "Mn-Kalpha",
      kind: ffGauss,
      gN: n_kalpha,
      gmu: mu_Kalpha, # since we count single electrons, index equals number electrons!
      gs: sigma_Kalpha), # sigma is approxed to 0.8 times the half width
    FitFuncArgs(
      name: "Mn-Kbeta",
      kind: ffGauss,
      gN: fixed,
      gmu: fixed,
      gs: fixed) # additional parameters fixed
  ]
  # serialize them to a `seq[float]` params
  result = params.serialize

proc getBoundsList(n: int): seq[tuple[l, u: float]] =
  for i in 0 ..< n:
    result.add (l: -Inf, u: Inf)
    # use this (or some other value for lower upper, just not `Inf`
    # to use global fitting routines of NLOPT, which require bounds for
    # all parameters
    #result.add (l: -100.0, u: 100.0)

func getFeSpectrumBounds(hist, binning: seq[float]): seq[tuple[l, u: float]] =
  let (mu_kalpha,
       sigma_kalpha,
       n_kalpha,
       mu_kalpha_esc,
       sigma_kalpha_esc,
       n_kalpha_esc) = getLines(hist, binning)
  result = getBoundsList(15)
  # set bound on paramerters
  # constrain amplitude of K_beta to some positive value
  result[2].l = 0
  result[2].u = 10000
  # constrain amplitude of K_alpha to some positive value
  result[9].l = 0
  result[9].u = 10000

  # location of K_alpha escape peak, little more than half of K_alpha location
  result[3].l = mu_kalpha_esc * 0.8
  result[3].u = mu_kalpha_esc * 1.2

  # result for center of K_alpha peak (should be at around 220 electrons, hits)
  result[10].l = mu_kalpha*0.8
  result[10].u = mu_kalpha*1.2

  # some useful bound for K_alpha escape peak width
  result[4].l = sigma_kalpha_esc*0.5
  result[4].u = sigma_kalpha_esc*1.5

  # some useful bound for K_alpha width
  result[11].l = sigma_kalpha*0.5
  result[11].u = sigma_kalpha*1.5
  # param 14: "N_{K_{#beta}}/N_{K_{#alpha}}"
  # known ratio of two K_alpha and K_beta, should be in some range
  result[14].l = 0.01
  result[14].u = 0.3

func getFeSpectrumChargeBounds(): seq[tuple[l, u: float]] =
  result = getBoundsList(6)
  # NOTE: Bounds for this fit seem unnecessary. Fit converges just fine
  # This way we get the Chi^2/dof
  # result[0].l = 0
  # result[0].u = 10000
  # result[3].l = 0
  # result[3].u = 10000
  #
  # result[1].l = mkalpha_esc*0.8
  # result[1].u = mkalpha_esc*1.2
  # result[4].l = mkalpha*0.8
  # result[4].u = mkalpha*1.2
  #
  # result[2].l = sigma_kalpha_esc*0.5
  # result[2].u = sigma_kalpha_esc*1.5
  # result[5].l = sigma_kalpha*0.5
  # result[5].u = sigma_kalpha*1.5

template fitNlopt*(xData, yData, errData: seq[float],
                   bounds: seq[tuple[l, u: float]],
                   algorithm: nlopt_algorithm,
                   params: seq[float],
                   fn: untyped,
                   body: untyped
              ): untyped =
  doAssert xData.len == yData.len
  echo "START PARAMS : ", params
  var opt = newNloptOpt[FitObject](algorithm, params.len, bounds)
  # get charges in the fit range
  # hand the function to fit as well as the data object we need in it
  # rebin the data
  #var xm = xData.rebin(3)
  #let xm = linspace(xData.min, xData.max, xData.len div 3)
  #var ym = yData.rebin(3).mapIt(it / 3.0)
  #var yErrM = errData.rebin(3)

  let fitObject = FitObject(
    x: xData,
    y: yData,
    yErr: errData
  )
  fitForNlopt(fnNlopt, fn, nfChiSq, toExport = false)
  let varStruct = newVarStruct(fnNlopt, fitObject)
  opt.setFunction(varStruct)
  # set relative precisions of x and y, as well as limit max time the algorithm
  # should take to 5 second (will take much longer, time spent in NLopt lib!)
  # these default values have proven to be working
  opt.xtol_rel = 1e-10
  opt.ftol_rel = 1e-10
  opt.maxtime  = 5.0
  opt.initialStep *= 2.0
  # start actual optimization
  let nloptRes {.inject.} = opt.optimize(params)

  body

  if opt.status < NLOPT_SUCCESS:
    echo opt.status
    echo "nlopt failed!"
  else:
    echo "Nlopt successfully exited with ", opt.status
  # clean up optimizer
  nlopt_destroy(opt.optimizer)

import ggplotnim, os
proc fitFeSpectrumImpl(hist, binning: seq[float]): FeSpecFitData =
  # given our histogram and binning data
  # for fit := (y / x) data
  # fit a double gaussian to the data
  let params = getFeSpectrumParams(hist, binning)
  let bounds = getFeSpectrumBounds(hist, binning)
  # this leaves parameter 7 and 8, as well as 12 and 13 without bounds
  # these describe the exponential factors contributing, since they will
  # be small anyways...
  echo "Len of bounds: ", bounds
  echo "Bounds for pixel fit: ", bounds
  echo "N params for pixel fit: ", len(params)

  # only fit in range up to 350 hits. Can take index 350 on both, since we
  # created the histogram for a binning with width == 1 pixel per hit
  let idx_tofit = toSeq(0 .. binning.high).filterIt(binning[it] < 350)
  let data_tofit = idx_tofit.mapIt(hist[it])
  let bins_tofit = idx_tofit.mapIt(binning[it])
  let err = data_tofit.mapIt(1.0)

  when false:
    let df = seqsToDf({ "x" : bins_to_fit,
                        "y" : data_to_fit,
                        "yFit" : bins_to_fit.mapIt(feSpectrumFunc(params, it))})
    ggplot(df, aes("x", "y")) +
      geom_histogram(stat = "identity") +
      geom_line(aes(y = "yFit")) +
      ggsave("/tmp/start_params.pdf")

    var pResNlopt: seq[float]
    fitNlopt(bins_tofit, data_tofit, err, bounds, LN_COBYLA, params,
             feSpectrumFunc):
      echo nloptRes[0]
      echo nloptRes[1]
      pResNlopt = nloptRes[0]

  #var cfg = MpConfig(xtol: 1e-30,
  #                   ftol: 1e-30,
  #                   gtol: 1e-30,
  #                   maxIter: 10000,
  #                   stepFactor: 100.0)
  #
  let (pRes, res) = fit(feSpectrumfunc,
                        params,
                        bins_tofit,
                        data_tofit,
                        err,
                        bounds = bounds)
                        #config = some(cfg))

  echoResult(pRes, res = res)
  let chiSq = res.chiSq
  let nDof = bins_to_fit.len - params.len
  let yFit = bins_to_fit.mapIt(feSpectrumFunc(pRes, it))
  result = initFeSpecData(data_to_fit,
                          bins_to_fit,
                          idx_kalpha = 10,
                          idx_sigma = 11,
                          pRes = pRes,
                          pErr = res.error,
                          xFit = bins_to_fit,
                          yFit = yFit,
                          chiSq = chiSq,
                          nDof = nDof)

proc fitFeSpectrumChargeImpl(hist, binning: seq[float]): FeSpecFitData =
  # given our histogram and binning data
  # for fit := (y / x) data
  # fit a double gaussian to the data
  let params = getFeSpectrumChargeParams(hist, binning)
  let bounds = getFeSpectrumChargeBounds()
  echo "Len of bounds: ", bounds
  echo "Bounds for pixel fit: ", bounds
  echo "N params for pixel fit: ", len(params)
  # only fit in range up to 350 hits. Can take index 350 on both, since we
  # created the histogram for a binning with width == 1 pixel per hit
  let idx_tofit = toSeq(0 .. binning.high).filterIt(binning[it] >= 200 and binning[it] < 4000)
  let data_tofit = idx_tofit.mapIt(hist[it])
  let bins_tofit = idx_tofit.mapIt(binning[it])
  let err = data_tofit.mapIt(1.0)
  let (pRes, res) = fit(feSpectrumChargeFunc,
                        params,
                        bins_tofit,
                        data_tofit,
                        err)
  echoResult(pRes, res = res)
  let chiSq = res.chiSq
  let nDof = bins_to_fit.len - params.len
  let yFit = bins_to_fit.mapIt(feSpectrumChargeFunc(pRes, it))
  result = initFeSpecData(data_to_fit,
                          bins_to_fit,
                          idx_kalpha = 4,
                          idx_sigma = 5,
                          pRes = pRes,
                          pErr = res.error,
                          xFit = bins_to_fit,
                          yFit = yFit,
                          chiSq = chiSq,
                          nDof = nDof)

proc fitFeSpectrum*[T: SomeInteger](data: seq[T]): FeSpecFitData =
  const binSize = 1.0
  let low = -0.5
  var high = max(data).float + 0.5
  let nbins = (ceil((high - low) / binSize)).int
  # using correct nBins, determine actual high
  high = low + binSize * nbins.float
  let (hist, bin_edges) = data.histogram(bins = nbins, range = (low, high))
  result = fitFeSpectrumImpl(hist.mapIt(it.float), bin_edges[0 .. ^2])

proc fitFeSpectrumCharge*[T](data: seq[T]): FeSpecFitData =
  # divide data (number of electrons) by 1000
  let dataPer1000 = data.mapIt(it.float / 1000.0)
  let (hist, bin_edges) = dataPer1000.histogram(bins = 300)
  result = fitFeSpectrumChargeImpl(hist.mapIt(it.float), bin_edges[0 .. ^2])

proc fitEnergyCalib*(x_ph, x_esc, x_ph_err, x_esc_err: float,
                     energies: array[2, float]): EnergyCalibFitData =
  ## before we plot the Fe spectrum, perform the calculations and fit of
  ## the fit to the spectrum peaks
  ## now we can fit the energy calibration function
  let energies = @energies
  let pixels_peaks = @[x_esc, x_ph]
  var pixels_err = @[x_esc_err, x_ph_err]
  # modify errors due to possibly, increase to 1.0 if they are 0
  pixels_err.applyIt(if it == 0.0: 1.0 else: it)
  let (pRes, res) = fit(linearFuncNoOffset,
                        @[1.0],
                        energies,
                        pixels_peaks,
                        pixels_err)
  echoResult(pRes, res = res)
  let E_calc = linspace(2.0, 7.0, 1000)
  let H_calc = E_calc.mapIt(linearFuncNoOffset(pRes, it))

  let ecData = initEnergyCalibData(energies, pixels_peaks, pixels_err, pRes, res.error,
                                   xFit = E_calc,
                                   yFit = H_calc,
                                   chiSq = res.chiSq,
                                   nDof = res.nfunc - res.nfree)
  let chiSq = res.reducedChiSq
  echo &"a^-1 = {ecData.aInv} +- {ecData.aInvErr}"
  echo "Chi^2 / dof = ", chiSq
  result = ecData

proc fitEnergyCalib*(feSpec: FeSpecFitData,
                     isPixel = true): EnergyCalibFitData =
  ## before we plot the Fe spectrum, perform the calculations and fit of
  ## the fit to the spectrum peaks
  ## now we can fit the energy calibration function
  if isPixel:
    const energies = [2.925, 5.755]
    result = fitEnergyCalib(feSpec.kalpha, feSpec.pRes[3],
                            feSpec.pErr[10], feSpec.pErr[3],
                            energies = energies)
  else:
    const energies = [2.942, 5.899]
    result = fitEnergyCalib(feSpec.kalpha, feSpec.pRes[1],
                            feSpec.pErr[4], feSpec.pErr[1],
                            energies = energies)

################################################################################
######################## Gas gain related fitting routines #####################
################################################################################

template fitPolyaTmpl(charges,
                      counts: seq[float],
                      chipNumber, runNumber: int,
                      actions: untyped): untyped =
  ## template which wraps the boilerplate code around the fitting implementation
  ## chosen
  ## proc to fit a polya distribution to the charge values of the
  ## reconstructed run. Called if `reconstruction` ran with --only_charge.
  ## After charge calc from TOT calib, this proc calculates the gas gain
  ## for this run.
  ## NOTE: charges has 1 element more than counts, due to representing the
  ## bin edges!
  # determine start parameters
  # estimate of 3000 or gas gain
  # TODO: test again with mpfit. Possible to get it working?
  # start parameters
  let
    # typical gain
    gain = 3500.0
    # start parameter for p[0] (scaling) is gas gain * max count value / 2.0
    # as good guess
    scaling = max(counts) * gain / 2.0
    # factor 23.0 found by trial and error! Works
    rms = standardDeviation(counts) / 25.0
    # combine paramters, 3rd arg from `polya.C` ROOT script by Lucian
    p {.inject.} = @[scaling, gain, 1.0] #gain * gain / (rms * rms) - 1.0]

  ## code will be inserted here
  actions

  # set ``x``, ``y`` result and use to create plot
  # TODO: replace by just using charges directly
  result.x = linspace(charges[0], charges[^1], counts.len)
  result.y = result.x.mapIt(polyaImpl(params, it))

  result.pRes = params
  # calc reduced Chi^2 from total Chi^2
  result.redChiSq = minVal / (charges.len - p.len).float

proc filterByCharge(charges, counts: seq[float]): (seq[float], seq[float]) =
  ## filters the given charges and counts seq by the cuts Christoph applied
  ## (see `getGasGain.C` in `resources` directory as reference)
  const
    ChargeLow = 1200.0
    ChargeHigh = 20_000.0
  result[0] = newSeqOfCap[float](charges.len)
  result[1] = newSeqOfCap[float](charges.len)
  for idx, ch in charges:
    # get those indices belonging to charges > 1200 e^- and < 20,000 e^-
    if ch >= ChargeLow and ch <= ChargeHigh:
      result[0].add ch
      result[1].add counts[idx]

proc fitPolyaNim*(charges,
                  counts: seq[float],
                  chipNumber, runNumber: int): FitResult =
  fitPolyaTmpl(charges, counts, chipNumber, runNumber):
    # set ``x``, ``y`` result and use to create plot
    # create NLopt optimizer without parameter bounds
    let bounds = @[(-Inf, Inf), (-Inf, Inf), (0.5, 15.0)]
    var opt = newNloptOpt[FitObject](LN_COBYLA, 3, bounds)
    #var opt = newNloptOpt[FitObject](LD_MMA, 3, bounds)
    # get charges in the fit range
    let (chToFit, countsToFit) = filterByCharge(charges, counts)
    # hand the function to fit as well as the data object we need in it
    let fitObject = FitObject(
      x: chToFit,
      y: countsToFit,
      yErr: countsToFit.mapIt(if it >= 1.0: sqrt(it) else: Inf)
    )
    let varStruct = newVarStruct(polya, fitObject)
    opt.setFunction(varStruct)
    # set relative precisions of x and y, as well as limit max time the algorithm
    # should take to 5 second (will take much longer, time spent in NLopt lib!)
    # these default values have proven to be working
    opt.xtol_rel = 1e-10
    opt.ftol_rel = 1e-10
    opt.maxtime  = 5.0
    # start actual optimization
    let (params, minVal) = opt.optimize(p)
    if opt.status < NLOPT_SUCCESS:
      echo opt.status
      echo "nlopt failed!"
    # clean up optimizer
    nlopt_destroy(opt.optimizer)
    echo &"Result of gas gain opt for chip {chipNumber} and run {runNumber}: "
    echo "\t ", params, " at chisq ", minVal

#proc fitPolyaPython*(charges,
#                     counts: seq[float],
#                     chipNumber, runNumber: int): FitResult =
#  fitPolyaTmpl(charges, counts, chipNumber, runNumber):
#    # set ``x``, ``y`` result and use to create plot
#    #echo "Charges ", charges
#    #let preX = linspace(min(charges), max(charges), 100)
#    #let preY = preX.mapIt(polyaImpl(p, it))
#    #let preTr = getTrace(preX, preY, "gas gain pre fit")
#    #echo preY
#    #plotGasGain(@[trData, preTr], chipNumber, runNumber, false)
#    # create NLopt optimizer without parameter bounds
#    let toFitInds = toSeq(0 ..< counts.len).filterIt(charges[it] > 1200.0)# and charges[it] < 4500.0)
#    let chToFit = toFitInds.mapIt(charges[it])
#    let countsToFit = toFitInds.mapIt(counts[it])
#    # try to fit using scipy.optimize.curve_fit
#    let bPy = @[@[-Inf, -Inf, 0.5], @[Inf, Inf, 15.0]]
#    let scipyOpt = pyImport("scipy.optimize")
#    let pyRes = scipyOpt.curve_fit(polyaPython, chToFit, countsToFit,
#                                   p0=p, bounds = bPy)
#    var params = newSeq[float](p.len)
#    var count = 0
#    for resP in pyRes[0]:
#      params[count] = resP.to(float)
#      inc count
#    let minVal = 0.0
#    echo &"Result of gas gain opt for chip {chipNumber} and run {runNumber}: "
#    echo "\t ", params, " at chisq `not available`"#, minVal

################################################################################
##################### Full dataset fitting routines ############################
################################################################################
# the routines in this block have to make use of full datasets, e.g.
# `CalibrationRuns.h5 + DataRuns.h5` in its entirety

proc fitChargeCalibVsGasGain*(gain, calib, calibErr: seq[float]): FitResult =
  ## fits a linear function to the relation between gas gain and the
  ## calibration factor of the charge Fe spectrum
  # approximate start parameters
  let p = @[30.0,
            (max(calib) - min(calib)) / (max(gain) - min(gain))]
  let (pRes, res) = fit(linearFunc, p, gain, calib, calibErr)
  # echo parameters
  echoResult(pRes, res = res)
  result.x = linspace(min(gain), max(gain), 100)
  result.y = result.x.mapIt(linearFunc(pRes, it))
  result.pRes = pRes
  result.pErr = res.error
  result.redChiSq = res.reducedChiSq
