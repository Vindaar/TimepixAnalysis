import sequtils, strutils, strformat
import ospaths
import future
import seqmath
import nimhdf5
import tables
import mpfit
import zero_functional
import stats
import nlopt
import math
import plotly

import tos_helpers
import helpers/utils
import ingrid_types
import ingridDatabase/databaseRead
import ingridDatabase/databaseDefinitions

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

func polyaImpl(p: seq[float], x: float): float =
  ## Polya function to fit to TOT histogram / charge in electrons of a
  ## run. This is the actual implementation of the polya distribution.
  ## Parameters:
  ## N     = p[0]    scaling factor
  ## G     = p[1]    gas gain
  ## theta = p[2]    parameter, which describes distribution (?! I guess it makes sens
  ##                 since we take its power and it enters gamma)
  # debugecho "Scale ", p[0]
  # debugecho "Gain ", p[1]
  # debugecho "theta ", p[2]
  let
    thetaDash = p[2] + 1
    coeff1 = (p[0] / p[1]) * pow((thetaDash), thetaDash) / gamma(thetaDash)
    coeff2 = pow((x / p[1]), p[2]) * exp(-thetaDash * x / p[1])
  result = coeff1 * coeff2

type
  FitObject = object
    x: seq[float]
    y: seq[float]

func polya(p: seq[float], fitObj: FitObject): float =
  ## Polya function to fit to TOT histogram / charge in electrons of a
  ## run. This is the polya function, which is handed to NLopt to perform
  ## the non-linear optimization of the calc'd chi^2.
  ## Parameters:
  ## N     = p[0]    scaling factor
  ## G     = p[1]    gas gain
  ## theta = p[2]    parameter, which describes distribution (?! I guess it makes sens
  ##                 since we take its power and it enters gamma)
  # debugecho "Scale ", p[0]
  # debugecho "Gain ", p[1]
  # debugecho "theta ", p[2]
  # TODO: add errors, proper chi^2 intead of unscaled values,
  # which results in huge chi^2 despite good fit
  let x = fitObj.x
  let y = fitObj.y
  var fitY = x.mapIt(polyaImpl(p, it))
  var diff = newSeq[float](x.len)
  result = 0.0
  for i in 0 .. x.high:
    diff[i] = y[i] - fitY[i]
    result += pow(diff[i], 2.0)

  result = result / (x.len - p.len).float

func sCurveFunc(p: seq[float], x: float): float =
  ## we fit the complement of a cumulative distribution function
  ## of the normal distribution
  # parameter p[2] == sigma
  # parameter p[1] == x0
  # parameter p[0] == scale factor
  result = normalCdfC(x, p[2], p[1]) * p[0]

func thlCalibFunc(p: seq[float], x: float): float =
  ## we fit a linear function to the charges and mean thl values
  ## of the SCurves
  result = p[0] + x * p[1]

func totCalibFunc(p: seq[float], x: float): float =
  ## we fit a combination of a linear and a 1 / x function
  ## The function is:
  ## ToT[clock cycles] = a * x + b - (c / (x - t))
  ## where x is the test pulse height in mV and:
  ## p = [a, b, c, t]
  result = p[0] * x + p[1] - p[2] / (x - p[3])

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

proc fitToTCalib*(tot: Tot, startFit = 0.0): FitResult =
  var
    # local mutable variables to potentially remove unwanted data for fit
    mPulses = tot.pulses.mapIt(it.float)
    mMean = tot.mean
    mStd = tot.std

  # define the start parameters. Use a good guess...
  let p = [0.149194, 23.5359, 205.735, -100.0]
  var pLimitBare: mp_par
  var pLimit: mp_par
  pLimit.limited = [1.cint, 1]
  pLimit.limits = [-100.0, 0.0]

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
                        bounds = @[pLimitBare, pLimitBare, pLimitBare, pLimit])
  echoResult(pRes, res = res)

  # the plot of the fit is performed to the whole pulses range anyways, even if
  result.x = linspace(tot.pulses[0].float, tot.pulses[^1].float, 100)
  result.y = result.x.mapIt(totCalibFunc(pRes, it))
  result.pRes = pRes
  result.pErr = res.error
  result.redChiSq = res.reducedChiSq

proc plotGasGain*[T](traces: seq[Trace[T]], chipNumber, runNumber: int) =
  ## given a seq of traces (polya distributions for gas gain) plot
  ## the data and the fit, save plots as svg.
  let
    layout = Layout(title: "Polya for gas gain",
                    width: 1200, height: 800,
                    xaxis: Axis(title: "charge / e-"),
                    yaxis: Axis(title: "counts"),
                    autosize: false)
    p = Plot[float](layout: layout, traces: traces)
  # save plots
  let filename = &"out/gas_gain_run_{runNumber}_chip_{chipNumber}.svg"
  p.show(filename)

proc getTrace[T](ch, counts: seq[T], info = "", `type`: PlotType = PlotType.Scatter): Trace[float] =
  result = Trace[float](`type`: `type`)
  # filter out clock cycles larger 300 and assign to `Trace`
  result.xs = ch
  result.ys = counts
  result.name = info

proc fitPolya*(charges,
               counts: seq[float],
               chipNumber, runNumber: int,
               createPlots = true): FitResult =
  ## proc to fit a polya distribution to the charge values of the
  ## reconstructed run. Called if `reconstruction` ran with --only_charge.
  ## After charge calc from TOT calib, this proc calculates the gas gain
  ## for this run.
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
    rms = standardDeviation(counts) / 23.0
    # combine paramters, 3rd arg from `polya.C` ROOT script by Lucian
    p = @[scaling, gain, gain * gain / (rms * rms) - 1.0]

  # data trace
  let trData = getTrace(charges,
                        counts.asType(float64),
                        &"polya data {chipNumber}",
                        PlotType.Bar)
  # create NLopt optimizer without parameter bounds
  var opt = newNloptOpt("LN_COBYLA", 3, @[])
  # hand the function to fit as well as the data object we need in it
  var fitObject: FitObject
  fitObject.x = charges
  fitObject.y = counts
  var varStruct = newVarStruct(polya, fitObject)
  opt.setFunction(varStruct)
  # set relative precisions of x and y, as well as limit max time the algorithm
  # should take to 5 second (will take much longer, time spent in NLopt lib!)
  # these default values have proven to be working
  opt.xtol_rel = 1e-8
  opt.ftol_rel = 1e-8
  opt.maxtime  = 5.0
  # start actual optimization
  let (params, minVal) = opt.optimize(p)
  if opt.status < NLOPT_SUCCESS:
    echo opt.status
    echo "nlopt failed!"
  # clean up optimizer
  nlopt_destroy(opt.optimizer)

  echo "Result of gas gain opt: ", params, " at chisq ", minVal

  # set ``x``, ``y`` result and use to create plot
  result.x = linspace(charges[0], charges[^1], 100)
  result.y = result.x.mapIt(polyaImpl(params, it))

  if createPlots:
    # create plots if desired
    let trFit = getTrace(result.x, result.y, "Gas gain fit")
    plotGasGain(@[trData, trFit], chipNumber, runNumber)

  result.pRes = params
  # calc reduced Chi^2 from total Chi^2
  result.redChiSq = minVal / (charges.len - p.len).float

proc cutFeSpectrum(data: array[4, seq[float64]], eventNum, hits: seq[int64]):
                    seq[(int64, int64, int64)] =
  ## proc which receives the data for the cut, performs the cut and returns tuples of
  ## event numbers, number of hits and the indices of the passing elements
  ## inputs:
  ##    data: array[4, seq[float64]] = array containing 4 sequences
  ##      - pos_x, pos_y, eccentricity, rms_transverse which we need for cuts
  ##    eventNum: seq[int] = sequence containing event numbers of data stored in
  ##        other seqs
  ##    hits: seq[int] = sequence containing the hits of the corresponding event
  ##        which are the data in the final spectrum
  ## outputs:
  ##    seq[int] = hits of the events passing the cuts
  # constants whihc define the cuts
  result = @[]
  const
    cut_x = 7.0
    cut_y = 7.0
    cut_r = 4.5
    cut_ecc_high = 1.3
    cut_rms_trans_high = 1.2

  let
    pos_x = data[0]
    pos_y = data[1]
    ecc = data[2]
    rms_trans = data[3]

  for i in 0 .. pos_x.high:
    let dist = distance( (pos_x[i] - cut_x), (pos_y[i] - cut_y) )
    if dist > cut_r:
      continue
    if ecc[i] > cut_ecc_high:
      continue
    if rms_trans[i] > cut_rms_trans_high:
      continue

    # else we keep these events, hence add event number to output
    result.add (eventNum[i], hits[i], i.int64)

proc createFeSpectrum*(h5f: var H5FileObj, runNumber: int) =
  ## proc which reads necessary reconstructed data from the given H5 file,
  ## performs cuts (as Christoph) and writes resulting spectrum to H5 file
  ## NOTE: currently does not perform a check whether the run is actually a
  ## calibration run or not...
  ## throws:
  ##    HDF5LibraryError = in case a call to the H5 library fails, this may happen
  # spectrum will be calculated for center chip, everything else makes no
  # sense, since we don't have X-rays there
  # obviously means that calibration will only be good for center chip,
  # but that is fine, since we only care about details on that one
  const chip = 3

  # what we need:
  # pos_x, pos_y, eccentricity < 1.3, transverse RMS < 1.2
  var reco_group = recoDataChipBase(runNumber) & $chip
  # get the group from file
  var group = h5f[reco_group.grp_str]
  # get the chip number from the attributes of the group
  let chip_number = group.attrs["chipNumber", int]
  # sanity check:
  assert chip_number == chip
  var
    pos_x_dset = h5f[(group.name / "centerX").dset_str]
    pos_y_dset = h5f[(group.name / "centerY").dset_str]
    ecc_dset   = h5f[(group.name / "eccentricity").dset_str]
    rms_trans_dset = h5f[(group.name / "rmsTransverse").dset_str]
    event_num_dset = h5f[(group.name / "eventNumber").dset_str]
    hits_dset = h5f[(group.name / "hits").dset_str]
  let
    pos_x  = pos_x_dset[float64]
    pos_y  = pos_y_dset[float64]
    ecc    = ecc_dset[float64]
    rms_trans = rms_trans_dset[float64]
    event_num = event_num_dset[int64]
    hits = hits_dset[int64]

  # given this data, filter all events which don't conform
  let hitsSpecTuples = cutFeSpectrum([pos_x, pos_y, ecc, rms_trans], event_num, hits)
  let nEventsPassed = hitsSpecTuples.len
  # with the events to use for the spectrum
  echo "Elements passing cut : ", nEventsPassed
  let
    eventSpectrum = hitsSpecTuples.mapIt(it[0])
    hitsSpectrum = hitsSpecTuples.mapIt(it[1])
    specIndices = hitsSpecTuples.mapIt(it[2])

  # given hits, write spectrum to file
  let
    spectrumDset  = h5f.create_dataset(group.name & "/FeSpectrum",
                                       nEventsPassed,
                                       dtype = int)
    specEventDset = h5f.create_dataset(group.name & "/FeSpectrumEvents",
                                       nEventsPassed,
                                       dtype = int)
    specIndDset = h5f.create_dataset(group.name & "/FeSpectrumIndices",
                                       nEventsPassed,
                                       dtype = int)
  spectrumDset[spectrumDset.all] = hitsSpectrum
  specEventDset[specEventDset.all] = eventSpectrum
  specIndDset[specIndDset.all] = specIndices

# TODO: also does not work I guess, because this won't be declared in case we import
# this module
# Find some way that works!
# when declaredInScope(ingridDatabase):
func calibrateCharge*(totValue: float, a, b, c, t: float): float =
  ## calculates the charge in electrons from the TOT value, based on the TOT calibration
  ## from MarlinTPC:
  ## measured and fitted is ToT[clock cycles] in dependency of TestPuseHeight [mV]
  ## fit function is:
  ##   ToT[clock cycles] = a*TestPuseHeight [mV]  + b - ( c / (TestPuseHeight [mV] -t) )
  ## isolating TestPuseHeight gives:
  ##   TestPuseHeight [mV] = 1/(2*a) * (ToT[clock cycles] - (b - a*t) +
  ##                         SQRT( (ToT[clock cycles] - (b - a*t))^2 +
  ##                               4*(a*b*t + a*c - a*t*ToT[clock cycles]) ) )
  ## conversion from charge to electrons:
  ##   electrons = 50 * testpulse[mV]
  ## so we have:
  ## Charge[electrons] = 50 / (2*a) * (ToT[clock cycles] - (b - a*t) +
  ##                     SQRT( (ToT[clock cycles] - (b - a*t))^2 +
  ##                           4*(a*b*t + a*c - a*t*ToT[clock cycles]) ) )
  # 1.sum term
  let p = totValue - (b - a * t)
  # 2. term of sqrt - neither is exactly the p or q from pq formula
  let q = 4 * (a * b * t  +  a * c  -  a * t * totValue)
  result = (50 / (2 * a))  * (p + sqrt(p * p + q))

proc applyChargeCalibration*(h5f: var H5FileObj, runNumber: int)
  {.raises: [HDF5LibraryError, ref ValueError, Exception]} =
  ## applies the charge calibration to the TOT values of all events of the
  ## given run

  # what we need:
  # TOT values of each run. Run them through calibration function with the TOT calibrated
  # values taken from the `InGridDatabase`
  var chipBase = recoDataChipBase(runNumber)
  # get the group from file
  for grp in keys(h5f.groups):
    if chipBase in grp:
      # now can start reading, get the group containing the data for this chip
      var group = h5f[grp.grp_str]
      # get the chip number from the attributes of the group
      let chipNumber = group.attrs["chipNumber", int]
      let chipName = group.attrs["chipName", string]
      # get dataset of hits
      let
        totDset = h5f[(grp / "ToT").dset_str]
        sumTotDset = h5f[(grp / "sumTot").dset_str]
        vlenInt = special_type(uint16)
        tots = totDset[vlenInt, uint16]
        sumTots = sumTotDset[int64]
      # now calculate charge in electrons for all TOT values
      # need calibration factors from InGrid database for that
      let (a, b, c, t) = getTotCalibParameters(chipName)
      #mapIt(it.mapIt(calibrateCharge(it.float, a, b, c, t)))
      var charge = newSeqWith(tots.len, newSeq[float]())
      var totalCharge = newSeq[float](sumTots.len)
      for i, vec in tots:
        # calculate charge values for individual pixels
        charge[i] = vec.mapIt((calibrateCharge(it.float, a, b, c, t)))
        # and for the sum of all in one cluster
        totalCharge[i] = calibrateCharge(sumTots[i].float, a, b, c, t)
      #let charge = tots --> map(it --> map(it --> calibrateCharge(it.float, a, b, c, t))) --> to(seq[seq[float]])
      # create dataset for charge values
      let vlenFloat = special_type(float64)
      var
        chargeDset = h5f.create_dataset(grp / "charge", charge.len, dtype = vlenFloat)
        totalChargeDset = h5f.create_dataset(grp / "totalCharge", charge.len, dtype = float64)

      template writeDset(dset: H5DataSet, data: untyped) =
        dset[dset.all] = data
        # add attributes for TOT calibration factors used
        dset.attrs["charge_a"] = a
        dset.attrs["charge_b"] = b
        dset.attrs["charge_c"] = c
        dset.attrs["charge_t"] = t
      chargeDset.writeDset(charge)
      totalChargeDset.writeDset(totalCharge)

proc calcGasGain*(h5f: var H5FileObj, runNumber: int) =
  ## fits the polya distribution to the charge values and writes the
  ## fit parameters (including the gas gain) to the H5 file
  var chipBase = recoDataChipBase(runNumber)
  # get the group from file
  echo "Calcing gas gain"
  for grp in keys(h5f.groups):
    if chipBase in grp:
      # now can start reading, get the group containing the data for this chip
      var group = h5f[grp.grp_str]
      # get the chip number from the attributes of the group
      let chipNumber = group.attrs["chipNumber", int]
      let chipName = group.attrs["chipName", string]
      # get dataset of hits
      var chargeDset = h5f[(grp / "charge").dset_str]
      var totDset = h5f[(grp / "ToT").dset_str]
      let vlenFloat = special_type(float64)
      let vlenInt = special_type(uint16)
      # get all charge values as seq[seq[float]] flatten
      let charges = chargeDset[vlenFloat, float64].flatten
      let tots = totDset[vlenInt, uint16].flatten
      # bin the data according to ToT values
      let (a, b, c, t) = getTotCalibParameters(chipName)
      # get bin edges by calculating charge values for all TOT values at TOT's bin edges
      #let bin_edges = mapIt(linspace(-0.5, 249.5, 251), calibrateCharge(it, a, b, c, t))
      let bin_edges = mapIt(linspace(-0.5, 100.5, 101), calibrateCharge(it, a, b, c, t))
      # the histogram counts are the same for ToT values as well as for charge values,
      # so calculate for ToT
      let binned = tots.histogram(bins = 101, range = (0.0, 100.0))
      # given binned histogram, fit polya
      let fitResult = fitPolya(bin_edges, binned.asType(float64), chipNumber, runNumber)
      # now write resulting fit parameters as attributes
      chargeDset.attrs["N"] = fitResult.pRes[0]
      chargeDset.attrs["G"] = fitResult.pRes[1]
      # TODO: Christoph takes the "mean gas gain" by calculating the mean
      # of the `chargePerPixelAssymBin` histogram instead of `G` into
      # account. Why?
      #chargeDset.attrs["G_mean"] = mean(binned.asType(float64))
      chargeDset.attrs["theta"] = fitResult.pRes[2]
      # TODO: get some errors from NLopt?
      #chargeDset.attrs["N_err"] = fitResult.pErr[0]
      #chargeDset.attrs["G_err"] = fitResult.pErr[1]
      #chargeDset.attrs["theta_err"] = fitResutl.pErr[2]
      chargeDset.attrs["redChiSq"] = fitResult.redChiSq

proc applyEnergyCalibration*(h5f: var H5FileObj, runNumber: int, calib_factor: float) =
  ## proc which applies an energy calibration based on the number of hit pixels in an event
  ## using a conversion factor of unit eV / hit pixel to the run given by runNumber contained
  ## in file h5f
  ## throws:
  ##     HDF5LibraryError = in case a call to the H5 library fails, this might be raised

  # what we need:
  # the hits of the clusters is all we need
  var chipBase = recoDataChipBase(runNumber)
  # get the group from file
  for grp in keys(h5f.groups):
    if chipBase in grp:
      # now can start reading, get the group containing the data for this chip
      var group = h5f[grp.grp_str]
      # get the chip number from the attributes of the group
      let chipNumber = group.attrs["chipNumber", int]
      # get dataset of hits
      var hits_dset = h5f[(grp / "hits").dset_str]
      let hits = hits_dset[int64]

      # now calculate energy for all hits
      let energy = mapIt(hits, float(it) * calib_factor)
      # create dataset for energy
      var energy_dset = h5f.create_dataset(grp / "energyFromPixel", energy.len, dtype = float)
      energy_dset[energy_dset.all] = energy
      # attach used conversion factor to dataset
      energy_dset.attrs["conversionFactorUsed"] = calib_factor
