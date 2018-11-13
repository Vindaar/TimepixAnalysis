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
import chroma

import tos_helpers
import helpers/utils
import ingrid_types
import ingridDatabase / [databaseRead, databaseDefinitions]
import procsForPython
from ingridDatabase/databaseWrite import writeCalibVsGasGain

## need nimpy to call python functions
import nimpy

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
  let
    thetaDash = p[2] + 1
    coeff1 = (p[0] / p[1]) * pow((thetaDash), thetaDash) / gamma(thetaDash)
    coeff2 = pow((x / p[1]), p[2]) * exp(-thetaDash * x / p[1])
  result = coeff1 * coeff2

type
  FitObject = object
    x: seq[float]
    y: seq[float]
    yErr: seq[float]

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
  let yErr = fitObj.yErr
  var fitY = x.mapIt(polyaImpl(p, it))
  var diff = newSeq[float](x.len)
  result = 0.0
  for i in 0 .. x.high:
    diff[i] = (y[i] - fitY[i]) / yErr[i]
    result += pow(diff[i], 2.0)
  result = result / (x.len - p.len).float


func sCurveFunc(p: seq[float], x: float): float =
  ## we fit the complement of a cumulative distribution function
  ## of the normal distribution
  # parameter p[2] == sigma
  # parameter p[1] == x0
  # parameter p[0] == scale factor
  result = normalCdfC(x, p[2], p[1]) * p[0]

func linearFunc(p: seq[float], x: float): float =
  result = p[0] + x * p[1]

proc thlCalibFunc(p: seq[float], x: float): float =
  ## we fit a linear function to the charges and mean thl values
  ## of the SCurves
  linearFunc(p, x)

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

proc chargeCalibVsGasGain(gain, calib, calibErr: seq[float]): FitResult =
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

proc plotGasGain*[T](traces: seq[Trace[T]], chipNumber, runNumber: int,
                     toSave = true) =
  ## given a seq of traces (polya distributions for gas gain) plot
  ## the data and the fit, save plots as svg.
  let
    layout = Layout(title: &"Polya for gas gain of chip {chipNumber} " &
                    &"and run {runNumber}",
                    width: 1200, height: 800,
                    xaxis: Axis(title: "charge / e-"),
                    yaxis: Axis(title: "counts"),
                    autosize: false)
    p = Plot[float](layout: layout, traces: traces)
  # save plots
  if toSave:
    let filename = &"out/gas_gain_run_{runNumber}_chip_{chipNumber}.svg"
    p.show(filename)
  else:
    p.show()

proc getTrace[T](ch, counts: seq[T], info = "", `type`: PlotType = PlotType.Scatter): Trace[float] =
  result = Trace[float](`type`: `type`)
  # filter out clock cycles larger 300 and assign to `Trace`
  result.xs = ch
  result.ys = counts
  result.name = info

template fitPolyaTmpl(charges,
                      counts: seq[float],
                      chipNumber, runNumber: int,
                      createPlots = true,
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

  # data trace
  let trData = getTrace(charges,
                        counts.asType(float64),
                        &"polya data {chipNumber}",
                        PlotType.Bar)

  ## code will be inserted here
  actions

  # set ``x``, ``y`` result and use to create plot
  # TODO: replace by just using charges directly
  result.x = linspace(charges[0], charges[^1], counts.len)
  result.y = result.x.mapIt(polyaImpl(params, it))

  if createPlots:
    # create plots if desired
    let trFit = getTrace(result.x, result.y, "Gas gain fit")
    plotGasGain(@[trData, trFit], chipNumber, runNumber)

  result.pRes = params
  # calc reduced Chi^2 from total Chi^2
  result.redChiSq = minVal / (charges.len - p.len).float

proc fitPolyaNim*(charges,
                  counts: seq[float],
                  chipNumber, runNumber: int,
                  createPlots = true): FitResult =
  fitPolyaTmpl(charges, counts, chipNumber, runNumber, createPlots):
    # set ``x``, ``y`` result and use to create plot
    # create NLopt optimizer without parameter bounds
    let bounds = @[(-Inf, Inf), (-Inf, Inf), (0.5, 15.0)]
    var opt = newNloptOpt("LN_COBYLA", 3, bounds)
    # hand the function to fit as well as the data object we need in it
    var fitObject: FitObject
    let toFitInds = toSeq(0 ..< charges.len).filterIt(charges[it] > 1200.0)# and charges[it] < 4500.0)
    let chToFit = toFitInds.mapIt(charges[it])
    let countsToFit = toFitInds.mapIt(counts[it])

    fitObject.x = chToFit#charges
    fitObject.y = countsToFit#counts
    fitObject.yErr = countsToFit.mapIt(if it >= 1.0: sqrt(it) else: Inf)
    var varStruct = newVarStruct(polya, fitObject)
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

proc fitPolyaPython*(charges,
                     counts: seq[float],
                     chipNumber, runNumber: int,
                     createPlots = true): FitResult =
  fitPolyaTmpl(charges, counts, chipNumber, runNumber, createPlots):
    # set ``x``, ``y`` result and use to create plot
    #echo "Charges ", charges
    #let preX = linspace(min(charges), max(charges), 100)
    #let preY = preX.mapIt(polyaImpl(p, it))
    #let preTr = getTrace(preX, preY, "gas gain pre fit")
    #echo preY
    #plotGasGain(@[trData, preTr], chipNumber, runNumber, false)
    # create NLopt optimizer without parameter bounds
    let toFitInds = toSeq(0 ..< counts.len).filterIt(charges[it] > 1200.0)# and charges[it] < 4500.0)
    let chToFit = toFitInds.mapIt(charges[it])
    let countsToFit = toFitInds.mapIt(counts[it])
    # try to fit using scipy.optimize.curve_fit
    let bPy = @[@[-Inf, -Inf, 0.5], @[Inf, Inf, 15.0]]
    let scipyOpt = pyImport("scipy.optimize")
    let pyRes = scipyOpt.curve_fit(polyaPython, chToFit, countsToFit,
                                   p0=p, bounds = bPy)
    var params = newSeq[float](p.len)
    var count = 0
    for resP in pyRes[0]:
      params[count] = resP.to(float)
      inc count
    let minVal = 0.0
    echo &"Result of gas gain opt for chip {chipNumber} and run {runNumber}: "
    echo "\t ", params, " at chisq `not available`"#, minVal

proc cutOnDsets[T](eventNumbers: seq[SomeInteger],
                   region: ChipRegion,
                   posX, posY: seq[T],
                   nPix: seq[SomeInteger],
                   data: varargs[tuple[d: seq[T], lower, upper: T]]):
                     (seq[SomeInteger], seq[SomeInteger], seq[SomeInteger]) =
  ## given tuples of sequences and values, will cut eventNumbers to
  ## those which match all cut criteria
  ## in addition we cut on the chip region and a minimum of 3 pixels
  ## returns a tuple of three seq.
  ## - the event number of the passing events
  ## - the number of pixels in said cluster
  ## - the index of the events, which pass the cut (not necessarily the same
  ##   as the event number!)
  # allow max size of seq to avoid too much reallocation
  for res in fields(result):
    res = newSeqOfCap[SomeInteger](eventNumbers.len)

  for i, ev in eventNumbers:
    # check number of pixels
    if nPix[i] < 3:
      continue
    # now perform region cut
    if not inRegion(posX[i], posY[i], crSilver):
      continue
    var skipThis = false
    for el in data:
      let
        d = el.d[i]
      if d < el.lower or d > el.upper:
        skipThis = true
        break
    if not skipThis:
      # else we add the element
      result[0].add ev
      result[1].add nPix[i]
      result[2].add i.int64

proc cutOnProperties*(h5f: var H5FileObj,
                      group: H5Group,
                      region: ChipRegion,
                      cuts: varargs[tuple[dset: string,
                                          lower, upper: float]]): seq[int] =
  ## applies the cuts from `cuts` and returns a sequence of indices, which pass
  ## the cut.
  ## Any datasets given will be converted to float after reading.
  ## For usage with FADC data, be careful to extract the event numbers using the
  ## indices first, before applying the indices on InGrid data.
  # first get data
  var dsets = newSeqOfCap[seq[float]](cuts.len)
  for c in cuts:
    let dset = h5f[(group.name / c.dset).dset_str]
    # use `convertType` proc from nimhdf5
    let convert = dset.convertType(float)
    dsets.add dset.convert
  # take any sequence for number of events, since all need to be of the same
  # type regardless
  let nEvents = dsets[0].len
  # if chip region not all, get `posX` and `posY`
  var
    posX: seq[float]
    posY: seq[float]
  case region
  of crAll:
    discard
  else:
    posX = h5f[group.name / "centerX", float]
    posY = h5f[group.name / "centerY", float]

  for i in 0 ..< nEvents:
    # cut on region if applicable
    if region != crAll and not inRegion(posX[i], posY[i], region):
      continue
    var skipThis = false
    for j in 0 ..< cuts.len:
      let
        el = cuts[j]
      let
        d = dsets[j][i]
      if d < el.lower or d > el.upper:
        skipThis = true
        break
    if not skipThis:
      # else add this index
      result.add i

proc cutOnProperties*(h5f: var H5FileObj,
                      group: H5Group,
                      cuts: varargs[tuple[dset: string,
                                         lower, upper: float]]): seq[int] {.inline.} =
  ## wrapper around the above for the case of the whole chip as region
  result = h5f.cutOnProperties(group, crAll, cuts)

proc cutFeSpectrum(data: array[4, seq[float64]], eventNum, hits: seq[int64]):
                    (seq[int64], seq[int64], seq[int64]) =
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
  ##    (seq[int], seq[int], seq[int]) = event numbers of passing clusters, number
  ##      of hits of these and indices of these clusters in the input
  # constants which define the cuts
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

  result = cutOnDsets(eventNum, crSilver,
                      pos_x, pos_y, hits,
                      (ecc, -Inf, cut_ecc_high),
                      (rms_trans, -Inf, cut_rms_trans_high))

proc createFeSpectrum*(h5f: var H5FileObj, runNumber, centerChip: int) =
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
  let chip = centerChip

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
  let (eventSpectrum,
       hitsSpectrum,
       specIndices) = cutFeSpectrum([pos_x, pos_y, ecc, rms_trans], event_num, hits)
  let nEventsPassed = eventSpectrum.len
  # with the events to use for the spectrum
  echo "Elements passing cut : ", nEventsPassed

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
  result = (50 / (2 * a)) * (p + sqrt(p * p + q))

proc applyChargeCalibration*(h5f: var H5FileObj, runNumber: int) =
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
      try:
        # `contains` calls `parseChipName`, which might throw a ValueError
        if not inDatabase(chipName):
          raise newException(KeyError, &"No entry for chip {chipName} in InGrid " &
            "database!")
      except ValueError as e:
          raise newException(KeyError, &"No entry for chip {chipName} in InGrid " &
            "database!")

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

proc calcGasGain*(h5f: var H5FileObj, runNumber: int, createPlots = false) =
  ## fits the polya distribution to the charge values and writes the
  ## fit parameters (including the gas gain) to the H5 file

  const
    hitLow = 1.5
    hitHigh = 180.5
    binCount = 181

  var chipBase = recoDataChipBase(runNumber)
  # get the group from file
  echo "Calcing gas gain for run: ", chipBase
  for grp in keys(h5f.groups):
    if chipBase in grp:
      # now can start reading, get the group containing the data for this chip
      var group = h5f[grp.grp_str]
      # get the chip number from the attributes of the group
      let chipNumber = group.attrs["chipNumber", int]
      #if chipNumber != 3:
      #  continue
      let chipName = group.attrs["chipName", string]
      # get dataset of hits
      var chargeDset = h5f[(grp / "charge").dset_str]
      var totDset = h5f[(grp / "ToT").dset_str]
      let vlenFloat = special_type(float64)
      let vlenInt = special_type(uint16)
      # get all charge values as seq[seq[float]] flatten
      let tots = totDset[vlenInt, uint16].flatten
      # bin the data according to ToT values
      let (a, b, c, t) = getTotCalibParameters(chipName)
      # get bin edges by calculating charge values for all TOT values at TOT's bin edges
      #let bin_edges = mapIt(linspace(-0.5, 249.5, 251), calibrateCharge(it, a, b, c, t))
      # skip range from 0 - 1.5 to leave out noisy pixels w/ very low ToT
      let bin_edges = mapIt(linspace(hitLow, hitHigh, binCount + 1), calibrateCharge(it, a, b, c, t))
      # the histogram counts are the same for ToT values as well as for charge values,
      # so calculate for ToT
      let binned = tots.histogram(bins = binCount, range = (hitLow + 0.5, hitHigh + 0.5))
      # given binned histogram, fit polya
      let fitResult = fitPolyaPython(bin_edges,
                                     binned.asType(float64),
                                     chipNumber, runNumber,
                                     createPlots = createPlots)
      # create dataset for polya histogram
      var polyaDset = h5f.create_dataset(group.name / "polya", (binCount + 1, 2), dtype = float64)
      var polyaFitDset = h5f.create_dataset(group.name / "polyaFit", (binCount, 2), dtype = float64)
      let polyaData = block:
        # convert and append one element to bin content array to fill up with 0
        var mbin = binned.asType(float64)
        mbin.add 0.0
        zip(bin_edges, mbin) --> map(@[it[0], it[1]])
      let polyaFitData = zip(fitResult.x, fitResult.y) --> map(@[it[0], it[1]])
      polyaDset[polyaDset.all] = polyaData
      polyaFitDset[polyaFitDset.all] = polyaFitData

      # now write resulting fit parameters as attributes
      template writeAttrs(d: var H5DataSet, fitResult: typed): untyped =
        d.attrs["N"] = fitResult.pRes[0]
        d.attrs["G"] = fitResult.pRes[1]
        # TODO: Christoph takes the "mean gas gain" by calculating the mean
        # of the `chargePerPixelAssymBin` histogram instead of `G` into
        # account. Why?
        #d.attrs["G_mean"] = mean(binned.asType(float64))
        d.attrs["theta"] = fitResult.pRes[2]
        # TODO: get some errors from NLopt?
        #d.attrs["N_err"] = fitResult.pErr[0]
        #d.attrs["G_err"] = fitResult.pErr[1]
        #d.attrs["theta_err"] = fitResutl.pErr[2]
        d.attrs["redChiSq"] = fitResult.redChiSq
      writeAttrs(chargeDset, fitResult)
      writeAttrs(polyaDset, fitResult)
      writeAttrs(polyaFitDset, fitResult)

proc writeFeFitParameters(dset: var H5DataSet,
                          popt, popt_E: seq[float],
                          pcov, pcov_E: seq[seq[float]]) =
  ## writes the fit parameters obtained via a call to `scipy.curve_fit` to
  ## the attributes of the dataset ``dset``, which should be the corresponding
  ## FeSpectrum[Charge] dataset
  proc writeAttrs(dset: var H5DataSet, name: string,
                  p: seq[float], c: seq[seq[float]]) =
    for i in 0 .. p.high:
      dset.attrs[name & $i] = p[i]
      dset.attrs[name & "Err_" & $i] = sqrt(c[i][i])
  dset.writeAttrs("p", popt, pcov)
  dset.writeAttrs("p_E_", popt_E, pcov_E)

proc writeEnergyPerAttrs(dset: var H5DataSet,
                         key: string,
                         scaling: float,
                         popt, pErr: float) =
  ## writes the fit results as `eV per Pixel` / `keV per Electron` to the
  ## given dataset as attributtes
  ## Scaling factor:
  ##   - `eV per Pixel`: 1000.0
  ##   - `keV per electron`: 1e-3
  let aInv = 1 / popt * scaling
  dset.attrs[key] = aInv
  dset.attrs["d_" & key] = aInv * pErr / popt

proc writeTextFields(dset: var H5DataSet,
                     texts: seq[string]) =
  ## writes the text fields, which are printed on the related plot
  ## as attributes
  echo "texts ", texts, " \n\n\n"
  for i in 0 .. texts.high:
    echo "Writing ", texts[i], " to ", dset.name
    dset.attrs[&"text_{i}"] = texts[i]

proc writeFeDset(h5f: var H5FileObj,
                 group, suffix: string,
                 feSpec, ecData: PyObject): (H5DataSet, H5DataSet) =
  ## writes the dataset for the given FeSpectrum (either FeSpectrum
  ## or FeSpectrumCharge), i.e. the actual data points for the plots
  ## created by the Python functions
  let
    feCounts = feSpec.hist.toNimSeq(float)
    feBins = feSpec.binning.toNimSeq(float)
    feFitX = feSpec.x_pl.toNimSeq(float)
    feFitY = feSpec.x_pl.toNimSeq(float)
  template createWriteDset(x, y: seq[float], name: string): untyped =
    var dset = h5f.create_dataset(group / name,
                                  (x.len, 2),
                                  dtype = float)
    let data = zip(x, y) --> map(@[it[0], it[1]]) --> to(seq[seq[float]])
    dset[dset.all] = data
    dset
  result[0] = createWriteDset(feBins, feCounts, "FeSpectrum" & $suffix & "Plot")
  result[1] = createWriteDset(feFitX, feFitY, "FeSpectrum" & $suffix & "FitPlot")

proc fitToFeSpectrum*(h5f: var H5FileObj, runNumber, chipNumber: int,
                      fittingOnly = true, outfiles: seq[string] = @[],
                      writeToFile = true) =
  ## (currently) calls Python functions from `ingrid` Python module to
  ## perform fit to the `FeSpectrum` dataset in the given run number
  ## NOTE: due to calling Python functions, this proc is *extremely*
  ## inefficient!
  ## The optional `writeToFile` flag can be set to `false` if this proc
  ## is to be called if the H5 file is opened read only to reproduce
  ## the results, but not rewrite them to the file (used in `plotData`)

  let pyFitFe = pyImport("ingrid.fit_fe_spectrum")
  # get the fe spectrum for the run
  let groupName = recoDataChipBase(runNumber) & $chipNumber
  var feDset = h5f[(groupName / "FeSpectrum").dsetStr]
  let feData = feDset[int64]
  # call python function with data
  let res = pyFitFe.fitAndPlotFeSpectrum([feData], "", ".", runNumber,
                                         fittingOnly, outfiles)


  # NOTE: this is a workaround for a weird bug we're seeing. If we don't close the
  # library here, we get an error in a call to `deleteAttribute` from within
  # `writeFeFitParameters`, line 738
  # If we just reopen the file we get an error from `existsAttribute` from line
  # 739
  # only if we also revisit the file it works. And this ONLY happens from the
  # `plotData` script, not from the `reconstruction` program.
  # discard h5f.close()
  # h5f = H5file(h5f.name, "rw")
  # h5f.visit_file()

  proc removePref(s: string, pref: string): string =
    result = s
    result.removePrefix(pref)

  proc extractAndWriteAttrs(h5f: var H5FileObj,
                            dset: var H5DataSet,
                            scaling: float,
                            res: PyObject,
                            key: string,
                            suffix = "") =
    let
      popt = res[0].popt.toNimSeq(float)
      pcov = res[0].pcov.toNimSeq(seq[float])
      popt_E = res[1].popt.toNimSeq(float)
      pcov_E = res[1].pcov.toNimSeq(seq[float])
      texts = res[2].toNimSeq(string)
    dset.writeFeFitParameters(popt, popt_E, pcov, pcov_E)

    writeEnergyPerAttrs(dset, key,
                        scaling,
                        popt_E[0],
                        pcov_E[0][0])
    dset.writeTextFields(texts)

    var (grp1, grp2) = h5f.writeFeDset(dset.parent, suffix, res[0], res[1])
    grp1.copy_attributes(dset.attrs)
    grp2.copy_attributes(dset.attrs)

  const eVperPixelScaling = 1e3
  if writeToFile:
    h5f.extractAndWriteAttrs(feDset, eVperPixelScaling, res, "eV_per_pix")

  # run might not have ``totalCharge`` dset, if no ToT calibration is available,
  # but user wishes Fe spectrum fit to # hit pixels
  if h5f.hasTotalChargeDset(runNumber, chipNumber):
    # also fit to the charge spectrum
    let feIdx = h5f[groupName / "FeSpectrumIndices", int64]
    let totChargeData = h5f[groupName / "totalCharge", float64]
    # extract correct clusters from totChargeData using indices
    let totChSpec = feIdx.mapIt(totChargeData[it.int])
    # create and write as a dataset
    var outfilesCharge: seq[string] = @[]
    if outfiles.len > 0:
      let parent = outfiles[0].parentDir
      outfilesCharge = outfiles.mapIt(parent / ("charge_" & it.removePref(parent & "/")))
    var totChDset: H5DataSet
    if writeToFile:
      totChDset = h5f.write_dataset(groupName / "FeSpectrumCharge", totChSpec)
    let resCharge = pyFitFe.fitAndPlotFeSpectrumCharge([totChSpec], "", ".",
                                                       runNumber, fittingOnly,
                                                       outfilesCharge)
    # given resCharge, need to write the result of that fit to H5 file, analogous to
    # `writeFitParametersH5` in Python
    const keVPerElectronScaling = 1e-3
    if writeToFile:
      h5f.extractAndWriteAttrs(totChDset, keVPerElectronScaling,
                               resCharge, "keV_per_electron",
                               "Charge")

  else:
    echo "Warning: `totalCharge` dataset does not exist in file. No fit to " &
      "charge Fe spectrum will be performed!"

proc performChargeCalibGasGainFit*(h5f: var H5FileObj) =
  ## performs the fit of the charge calibration factors vs gas gain fit
  ## Assumes:
  ## - h5f points to a h5 file of `runType == rtCalibration`
  ## - for all runs the Fe spectrum was calculated and fitted
  ## writes the resulting fit data to the ingridDatabase
  # iterate over all runs, extract center chip grou
  var
    calib = newSeq[float64]()
    calibErr = newSeq[float64]()
    gainVals = newSeq[float64]()
    centerChipName: string
    centerChip = 0
    # to differentiate old and new TOS data
    rfKind: RunFolderKind
  var recoGroup = h5f[recoGroupGrpStr]
  if "centerChipName" in recoGroup.attrs:
    centerChipName = recoGroup.attrs["centerChipName", string]
    rfKind = parseEnum[RunFolderKind](recoGroup.attrs["runFolderKind", string])
    centerChip = recoGroup.attrs["centerChip", int]
  else:
    rfKind = rfOldTos
    centerChipName = ""
  for run, grp in runs(h5f):
    var centerChipGrp: H5Group
    # TODO: TAKE OUT once we have ran over old CalibrationRuns again
    # so that we have `centerChipName` attribute there as well!
    # now iterate over chips in this run
    for chpGrp in items(h5f, start_path = grp):
      centerChipGrp = chpGrp
      case rfKind
      of rfOldTos, rfSrsTos:
        if centerChipName.len == 0:
          centerChipName = centerChipGrp.attrs["chipName", string]
      of rfNewTos:
        if ("chip_" & $centerChip) notin centerChipGrp.name:
          # skip this chip, since not the center chip
          echo "Skipping group ", centerChipGrp.name
          continue
        else:
          echo "\t taking group ", centerChipGrp.name
      of rfUnknown:
        echo "Unknown run folder kind. Skipping charge calibration for run " &
          centerChipGrp.name & "!"
        continue
      # read the chip name
      # given correct group, get the `charge` and `FeSpectrumCharge` dsets
      var
        chargeDset = h5f[(centerChipGrp.name / "charge").dset_str]
        feChargeSpec = h5f[(centerChipGrp.name / "FeSpectrumCharge").dset_str]
      let
        keVPerE = feChargeSpec.attrs["keV_per_electron", float64]
        dkeVPerE = feChargeSpec.attrs["d_keV_per_electron", float64]
        gain = chargeDset.attrs["G", float64]
      calib.add keVPerE * 1e6
      calibErr.add dkeVPerE * 1e6
      gainVals.add gain

  # now that we have all, plot them first
  let chGainTrace = Trace[float64](mode: PlotMode.Markers, `type`: PlotType.Scatter)
  chGainTrace.xs = gainVals
  chGainTrace.ys = calib
  chGainTrace.ys_err = newErrorBar(calibErr,
                                   color = Color(r: 0.5, g: 0.5, b: 0.5, a: 1.0))
  chGainTrace.name = "Charge calibration factors vs gas gain"

  # TODO: refactor the following by creating a function which takes care of
  # boilerplate in the whole file here
  let
    fitResult = chargeCalibVsGasGain(gainVals, calib, calibErr)
    fitTrace = Trace[float64](mode: PlotMode.Lines, `type`: PlotType.Scatter,
                              xs: fitResult.x,
                              ys: fitResult.y,
                              name: "ChiSq: " & $fitResult.redChiSq)

  # write results to ingrid database
  writeCalibVsGasGain(gainVals, calib, calibErr, fitResult, centerChipName)

  let
    lo = Layout(title: "Charge calibration factors vs gas gain",
                width: 1200, height: 800,
                xaxis: Axis(title: "Gas gain `G`"),
                yaxis: Axis(title: "Calibration factor `a^{-1}` [1e-6 keV / e]"))
    p = Plot[float64](layout: lo, traces: @[chGainTrace, fitTrace])
  # TODO: change filename to include more info
  p.show("gasgain_vs_calibration_charge.svg")

proc calcEnergyFromPixels*(h5f: var H5FileObj, runNumber: int, calib_factor: float) =
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

proc calcEnergyFromCharge*(h5f: var H5FileObj) =
                             #runType: RunTypeKind) =
  ## proc which applies an energy calibration based on the number of electrons in a cluster
  ## using a conversion factor of unit 1e6 keV / electron to the run given by runNumber contained
  ## in file h5f
  ## The `performChargeCalibGasGainFit` has to be run on the calibration dataset before,
  ## i.e. we need calibration fit results in the `ingridDatabase` for the center chip of
  ## the current detector
  ## throws:
  ##     HDF5LibraryError = in case a call to the H5 library fails, this might be raised
  # get the parameters from
  # TODO: FIX THIS STUPID HACK!!!
  const chip2014 = "D03W63"
  const chip2017 = "H10W69"
  var chipName = ""
  # get the calibration factors for the current center chip
  # TODO: since we always use the same calibration factors for all runs in
  # a single file, make sure we only call this proc once and iterate over
  # all runs in here!
  var
    b: float
    m: float
  # get the group from file

  for num, grp in runs(h5f):
    # now can start reading, get the group containing the data for this chip
    var group = h5f[grp.grp_str]
    echo "Energy from charge calibration for run ", num
    if chipName.len == 0:
      # HACK, only happens in first iteration...
      # only need to set it once
      let nChips = group.attrs["numChips", int]
      if nChips == 1:
        # 2014/15 data
        chipName = chip2014
      elif nChips == 7:
        chipName = chip2017
      else: quit("Invalid number of chips for now! " & $nChips)
      # get parameters during first iter...
      (b, m) = getCalibVsGasGainFactors(chipName)

    # now iterate over chips in this run
    for chpGrp in items(h5f, start_path = group.name):
      if "fadc" in chpGrp.name:
        continue
      # get the chip number from the attributes of the group
      var mchpGrp = chpGrp
      let chipNumber = mchpGrp.attrs["chipNumber", int]
      # get dataset of total charge
      var chargeDset = h5f[(mchpGrp.name / "charge").dset_str]
      var gain: float
      if "G" in chargeDset.attrs:
        gain = chargeDset.attrs["G", float64]
      else:
        # TODO: good idea? or just call calc gas gain?
        discard h5f.close()
        raise newException(Exception, "Gas gain not yet calculated for " &
          "run " & $mchpGrp.name)
      let totCharge = h5f[(mchpGrp.name / "totalCharge"), float64]
      # remove correction factor of 1e6
      let calibFactor = linearFunc(@[b, m], gain) * 1e-6
      # now calculate energy for all hits
      let energy = mapIt(totCharge, it * calibFactor)
      # create dataset for energy
      var energy_dset = h5f.create_dataset(mchpGrp.name / "energyFromCharge", energy.len, dtype = float)
      energy_dset[energy_dset.all] = energy
      # attach used conversion factor to dataset
      energy_dset.attrs["Calibration function"] = "y = m * x + b"
      energy_dset.attrs["calib_b"] = b
      energy_dset.attrs["calib_m"] = m
