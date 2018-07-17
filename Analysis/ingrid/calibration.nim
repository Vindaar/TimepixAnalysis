import sequtils, strutils
import ospaths
import future
import seqmath
import nimhdf5
import tables
import mpfit

import tos_helpers
import ingrid_types

type
  # type to store results of fitting with mpfit
  # defined here instead of ingrid_types, since it's only related
  # to calibration
  FitResult* = object
    x*: seq[float]
    y*: seq[float]
    pRes*: seq[float]
    pErr*: seq[float]
    redChiSq*: float

const
  NTestPulses = 1000.0
  ScalePulses = 1.01
  CurveHalfWidth = 15

func findDrop(thl, count: seq[int]): (float, int, int) =
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

func polya(p: seq[float], x: float): float =
  ## Polya function to fit to TOT histogram / charge in electrons of a
  ## run.
  ## Parameters:
  ## N     = p[0]    scaling factor
  ## G     = p[1]    gas gain
  ## theta = p[2]    parameter, which describes distribution (?! I guess it makes sens
  ##                 since we take its power and it enters gamma)
  let
    thetaDash = p[2] + 1
    coeff1 = (p[0] / p[1]) * pow(thetaDash, thetaDash) / gamma(thetaDash)
    coeff2 = pow(x / p[1], p[2]) * exp(-thetaDash * x / p[1])
  result = coeff1 * coeff2

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


proc cutFeSpectrum(data: array[4, seq[float64]], event_num, hits: seq[int64]): seq[int64] =
  ## proc which receives the data for the cut, performs the cut and returns the
  ## event numbers of the passing elements
  ## inputs:
  ##    data: array[4, seq[float64]] = array containing 4 sequences
  ##      - pos_x, pos_y, eccentricity, rms_transverse which we need for cuts
  ##    event_num: seq[int] = sequence containing event numbers of data stored in
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
    result.add hits[i]

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
  let hits_spectrum = cutFeSpectrum([pos_x, pos_y, ecc, rms_trans], event_num, hits)
  # with the events to use for the spectrum
  echo "Elements passing cut : ", hits_spectrum.len

  # given hits, write spectrum to file
  var spectrum_dset = h5f.create_dataset(group.name & "/FeSpectrum", hits_spectrum.len, dtype = int)
  spectrum_dset[spectrum_dset.all] = hits_spectrum

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
  var chip_base = recoDataChipBase(runNumber)
  # get the group from file
  for grp in keys(h5f.groups):
    if chip_base in grp:
      # now can start reading, get the group containing the data for this chip
      var group = h5f[grp.grp_str]
      # get the chip number from the attributes of the group
      let chip_number = group.attrs["chipNumber", int]
      # get dataset of hits
      var totDset = h5f[(grp / "TOT").dset_str]
      let tots = totDset[int64]
      # now calculate charge in electrons for all TOT values
      let charge = @[0] #mapIt(tots, float(it) * calib_factor)
      # create dataset for charge values 
      var chargeDset = h5f.create_dataset(grp / "charge", charge.len, dtype = float)
      chargeDset[chargeDset.all] = charge

      # add attributes for TOT calibration factors used
      #chargeDset.attrs["conversionFactorUsed"] = calib_factor

proc applyEnergyCalibration*(h5f: var H5FileObj, runNumber: int, calib_factor: float) =
  ## proc which applies an energy calibration based on the number of hit pixels in an event
  ## using a conversion factor of unit eV / hit pixel to the run given by runNumber contained
  ## in file h5f
  ## throws:
  ##     HDF5LibraryError = in case a call to the H5 library fails, this might be raised

  # what we need:
  # the hits of the clusters is all we need
  var chip_base = recoDataChipBase(runNumber)
  # get the group from file
  for grp in keys(h5f.groups):
    if chip_base in grp:
      # now can start reading, get the group containing the data for this chip
      var group = h5f[grp.grp_str]
      # get the chip number from the attributes of the group
      let chip_number = group.attrs["chipNumber", int]
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




  
