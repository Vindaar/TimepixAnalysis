import std / [sequtils, strutils, strformat, tables, math, times]
import std / os except FileInfo
import seqmath, fenv
import nimhdf5
import mpfit
import zero_functional
import chroma

import tos_helpers
import helpers/utils
import ingrid_types
import ingridDatabase / [databaseRead, databaseDefinitions, databaseUtils]
import ingrid / calibration / [fit_functions, calib_fitting, calib_plotting]
from ingridDatabase/databaseWrite import writeCalibVsGasGain

import logging
# set up the logger
var L = newConsoleLogger()
if not dirExists("logs"):
  createDir("logs")
var fL = newFileLogger("logs/calibration.log", fmtStr = verboseFmtStr)
when isMainModule:
  addHandler(L)
  addHandler(fL)

# the filter we use globally in this file
let filter = H5Filter(kind: fkZlib, zlibLevel: 4)

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

proc cutFeSpectrum(pos_x, pos_y, ecc, rms_trans: seq[float], eventNum, hits: seq[int64]):
                    (seq[int64], seq[int64], seq[int64]) =
  ## proc which receives the data for the cut, performs the cut and returns tuples of
  ## event numbers, number of hits and the indices of the passing elements
  ## inputs:
  ##    - pos_x, pos_y, eccentricity, rms_transverse which we need for cuts
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
  result = cutOnDsets(eventNum, crSilver,
                      pos_x, pos_y, hits,
                      (ecc, -Inf, cut_ecc_high),
                      (rms_trans, -Inf, cut_rms_trans_high))

proc cutFeSpectrum*(df: DataFrame): DataFrame =
  ## proc which receives the data for the cut, performs the cut and returns tuples of
  ## event numbers, number of hits and the indices of the passing elements
  ## inputs:
  ##    - df: Needs to contain:
  ##      - centerX, centerY, eccentricity, rmsTransverse which we need for cuts
  ##      - eventNum: seq[int] = sequence containing event numbers of data stored in
  ##        other seqs
  ##      - hits: seq[int] = sequence containing the hits of the corresponding event
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
  doAssert "centerX" in df
  doAssert "centerY" in df
  doAssert "rmsTransverse" in df
  doAssert "eccentricity" in df
  result = df.filter(f{float -> bool:
                        `eccentricity` < cut_ecc_high and
                        `rmsTransverse` < cut_rms_trans_high and
                        inRegion(df["centerX"][idx], df["centerY"][idx], crSilver)})

proc createFeSpectrum*(h5f: H5File, runNumber, centerChip: int) =
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

  info "Plotting facet for variables used for Fe spectrum for run: " & $runNumber
  plotFeSpectrumInfoFacet(pos_x, pos_y, ecc, rms_trans, hits,
                          runNumber = runNumber,
                          chipNumber = chipNumber,
                          pathPrefix = h5f.attrs[PlotDirPrefixAttr, string])

  # given this data, filter all events which don't conform
  let (eventSpectrum,
       hitsSpectrum,
       specIndices) = cutFeSpectrum(pos_x, pos_y, ecc, rms_trans, event_num, hits)
  let nEventsPassed = eventSpectrum.len
  # with the events to use for the spectrum
  info "Elements passing cut for Fe spectrum : ", nEventsPassed, " for run: ", runNumber

  # given hits, write spectrum to file
  let
    spectrumDset  = h5f.create_dataset(group.name & "/FeSpectrum",
                                       nEventsPassed,
                                       dtype = int,
                                       overwrite = true,
                                       filter = filter)
    specEventDset = h5f.create_dataset(group.name & "/FeSpectrumEvents",
                                       nEventsPassed,
                                       dtype = int,
                                       overwrite = true,
                                       filter = filter)
    specIndDset = h5f.create_dataset(group.name & "/FeSpectrumIndices",
                                       nEventsPassed,
                                       dtype = int,
                                       overwrite = true,
                                       filter = filter)
  spectrumDset[spectrumDset.all] = hitsSpectrum
  specEventDset[specEventDset.all] = eventSpectrum
  specIndDset[specIndDset.all] = specIndices

import unchained
func charge*(C: FemtoFarad, U: MilliVolt): UnitLess =
  ## Returns the charge on a capacitance of `C` at a voltage
  ## of `U` in *electrons*
  let e = 1.602176634e-19.C
  result = C * U / e

func capacityToFactor*(C: FemtoFarad): UnitLess =
  ## Returns the number of electrons one `MilliVolt` of voltage corresponds
  ## to on a capacitor of `C`.
  result = charge(C, 1.mV)

func calibrateCharge*(totValue: float, C: FemtoFarad, a, b, c, t: float): float =
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
  let factor = capacityToFactor(C)
  # 1.sum term
  let p = totValue - (b - a * t)
  # 2. term of sqrt - neither is exactly the p or q from pq formula
  let q = 4 * (a * b * t  +  a * c  -  a * t * totValue)
  result = (factor / (2 * a)) * (p + sqrt(p * p + q))

proc invertCharge*(charge: float, C: FemtoFarad, a, b, c, t: float): uint16 =
  ## Inverts the charge calibration for the given `charge` input and returns a
  ## ToT value using the fit parameters.
  ##
  ## See the docstring of `calibrateCharge` above.
  let factor = capacityToFactor(C)
  let eqV = charge / factor # equivalent voltage of the recorded charge
  let tot = (a * eqV + b - (c / (eqV - t) ))
  if tot > 11810: return 11810.uint16
  #doAssert tot < uint16.high.float, " Input was " & $tot & " from charge: " & $charge
  # floor conversion because if not one full clock cycle, it would not record!
  if tot < 0.0: result = 0.uint16
  else: result = tot.floor.uint16

proc applyChargeCalibration*(h5f: H5File, runNumber: int,
                             toDelete: bool = false) =
  ## applies the charge calibration to the TOT values of all events of the
  ## given run
  # what we need:
  # TOT values of each run. Run them through calibration function with the TOT calibrated
  # values taken from the `InGridDatabase`
  var chipBase = recoDataChipBase(runNumber)
  # get the group from file
  info "Applying charge calibration for run: ", runNumber
  let fileInfo = h5f.getFileInfo() # for timepix version
  var runPeriodPrinted: set[uint8] # for each chip number
  for run, chip, grp in chipGroups(h5f):
    doAssert chipBase in grp
    doAssert run == runNumber
    # now can start reading, get the group containing the data for this chip
    var group = h5f[grp.grp_str]
    # get the chip number from the attributes of the group
    let chipNumber = group.attrs["chipNumber", int]
    doAssert chipNumber == chip
    let chipName = group.attrs["chipName", string]
    try:
      # `contains` calls `parseChipName`, which might throw a ValueError
      if not inDatabase(chipName, runNumber):
        raise newException(KeyError, &"No entry for chip {chipName} in InGrid " &
          "database!")
    except ValueError as e:
        raise newException(KeyError, &"No entry for chip {chipName} in InGrid " &
          "database! Internal exception message " & e.msg)
    # get dataset of hits
    let
      totDset = h5f[(grp / "ToT").dset_str]
      sumTotDset = h5f[(grp / "sumTot").dset_str]
      vlenInt = special_type(uint16)
      tots = totDset[vlenInt, uint16]
      sumTots = sumTotDset[int64]
    # now calculate charge in electrons for all TOT values
    # need calibration factors from InGrid database for that
    let ((a, b, c, t), runPeriod) = getTotCalibParameters(chipName, runNumber)
    if chip.uint8 notin runPeriodPrinted:
      runPeriodPrinted.incl chip.uint8
      info "Using calibrations for chip: ", chipName, "(# ", chip, ") from run period: ", runPeriod
    #mapIt(it.mapIt(calibrateCharge(it.float, a, b, c, t)))
    var charge = newSeqWith(tots.len, newSeq[float]())
    var totalCharge = newSeq[float](sumTots.len)
    let capacitance = fileInfo.timepix.getCapacitance()
    for i, vec in tots:
      # calculate charge values for individual pixels
      charge[i] = vec.mapIt((calibrateCharge(it.float, capacitance, a, b, c, t)))
      # and for the sum of all in one cluster
      totalCharge[i] = charge[i].sum
    # create dataset for charge values
    let vlenFloat = special_type(float64)
    if toDelete:
      template ifDelete(name: string): untyped =
        if grp / name in h5f:
          doAssert h5f.delete(grp / name)
      ifDelete("charge")
      ifDelete("totalCharge")
      h5f.flush()

    var
      chargeDset = h5f.create_dataset(grp / "charge", charge.len, dtype = vlenFloat, overwrite = true, filter = filter)
      totalChargeDset = h5f.create_dataset(grp / "totalCharge", charge.len, dtype = float64, overwrite = true, filter = filter)
    template writeDset(dset: H5DataSet, data: untyped) =
      dset[dset.all] = data
      # add attributes for TOT calibration factors used
      dset.attrs["charge_a"] = a
      dset.attrs["charge_b"] = b
      dset.attrs["charge_c"] = c
      dset.attrs["charge_t"] = t
      dset.attrs["runPeriod (charge)"] = runPeriod
    chargeDset.writeDset(charge)
    totalChargeDset.writeDset(totalCharge)

proc calcMeanGainFit*(x, y: seq[float]): float =
  ## Computes the `G_fitmean` value, i.e. the *mean* of the data
  ## described by the fit to the polya distribution.
  result = (zip(x, y) -->
            map(it[0] * it[1]) -->
            fold(0.0, a + it)) / y.sum.float

proc initGasGainIntervalResult(g: GasGainIntervalData,
                               fitResult: FitResult,
                               binned: seq[int], bin_edges: seq[float],
                               sliceIdx: Slice[int],
                               sliceStartEvNum, sliceStopEvNum: int): GasGainIntervalResult =
  result.idx       = g.idx
  result.interval  = g.interval
  result.minInterval  = g.minInterval
  result.tStart    = g.tStart
  result.tStop     = g.tStop
  result.tLength   = (g.tStop - g.tStart).float / 60.0 # to easier keep track of lengths
  result.N         = fitResult.pRes[0]
  result.G_fit     = fitResult.pRes[1]
  result.sliceStart = sliceIdx.a
  result.sliceStop = sliceIdx.b
  result.sliceStartEvNum = sliceStartEvNum
  result.sliceStopEvNum = sliceStopEvNum
  # also calculate the mean of the polya fit. This usually is *not* the same
  # as the `meanGain` (mean of data). `meanGainFit` is what we used errorneously
  # in the past!
  result.G_fitmean = calcMeanGainFit(fitResult.x, fitResult.y)
  # calculate the mean of the data histogram. This is what Krieger calls
  # `G_mean` and uses for the rest of the calculation!
  let meanGain = histMean(binned, bin_edges[0 .. ^2])
  result.G         = meanGain
  result.theta     = fitResult.pRes[2]
  result.redChiSq  = fitResult.redChiSq

proc writePolyaDsets(h5f: H5File, group: H5Group,
                     chargeDset: H5DataSet,
                     binned: seq[int], bin_edges: seq[float],
                     fitResult: FitResult,
                     cutFormula: string,
                     gasGainInterval = none[GasGainIntervalData]()) =
  # all polya datasets are stored in `polya` subgroup
  const
    subgroup = "polyaDsets"
  # create dataset for polya histogram
  let dsetSuffix = if gasGainInterval.isSome:
                     gasGainInterval.get.toDsetSuffix
                   else: ""
  let
    pName = group.name / subgroup / &"polya{dsetSuffix}"
    pFitName = group.name / subgroup / &"polyaFit{dsetSuffix}"
  template ifDelete(name: untyped): untyped =
    if name in h5f:
      doAssert h5f.delete(name)
  ifDelete(group.name / "polya")
  ifDelete(group.name / "polyaFit")
  ifDelete(pName)
  ifDelete(group.name / &"polya{dsetSuffix}")
  ifDelete(pFitName)
  ifDelete(group.name / &"polyaFit{dsetSuffix}")

  var polyaDset = h5f.create_dataset(pName,
                                     (binned.len, 2),
                                     dtype = float64,
                                     overwrite = true,
                                     filter = filter)
  polyaDset[polyaDset.all] = zip(bin_edges, binned) -->
    map(@[it[0], it[1].float])

  # now write resulting fit parameters as attributes
  template writeAttrs(d: H5DataSet, fitResult: FitResult, cutFormula: string): untyped =
    d.attrs["applied Cut for gas gain"] = cutFormula
  writeAttrs(chargeDset, fitResult, cutFormula)
  writeAttrs(polyaDset, fitResult, cutFormula)
  # IMPORTANT: only write these to polya dataset and not to charge dataset. Otherwise the
  # latter will end up with hundreds of attributes which slows everything to a crawl
  polyaDset.attrs["N"] = fitResult.pRes[0]
  polyaDset.attrs["G_fit"] = fitResult.pRes[1]
  polyaDset.attrs["theta"] = fitResult.pRes[2]

proc writeGasGainAttributes(h5f: H5File, group: H5Group, sliceCount: int,
                            interval: int,
                            runPeriod: string) =
  ## Writes the attributes about the number of gas gain slices to the relevant
  ## objects (charge dset & polya group) as well as the "current" gas gain
  ## slice length
  # Note: a `sliceCount` of 0 means that *no* slicing  was used, not that there
  # is no available gas gain! (this also implies interval = 0)
  var gainDsetName = group.name / "gasGainSlices" & $interval
  if gainDsetName notin h5f:
    raise newException(KeyError, "The gas gain dataset : " & $gainDsetName &
      "does not exist in the H5 file: " & $h5f.name & ". Did you call " &
      "`writeGasGainAttributes` before calling `writeGasGainSliceData`?")
  let chargeDset = h5f[(group.name / "charge").dset_str]
  let polyaGrp = h5f[(group.name / "polyaDsets").grp_str]
  let gainDset = h5f[gainDsetName.dset_str]
  template writeAttrs(h5o: untyped): untyped =
    h5o.attrs["numGasGainSlices"] = sliceCount
    h5o.attrs["gasGainInterval"] = interval
    h5o.attrs["runPeriod (gas gain)"] = runPeriod
  writeAttrs(chargeDset)
  writeAttrs(polyaGrp)
  writeAttrs(gainDset)

proc writeGasGainSliceData(h5f: H5File, group: H5Group, slices: seq[GasGainIntervalResult],
                           cutFormula: string) =
  ## writes the information about the gas gain slices, including the fit results
  ## for each polya fit as one composite dataset to the H5 file.
  ##
  ## Dataset name: `gasGainSlices{interval}`
  ## There is a hardlink to `gasGainSlices`, which is overwritten each time a new
  ## gas gain computation is run so the ``latest`` computation is accessible that
  ## way!
  doAssert slices.len > 0, "At least one gas gain time slice is needed!"
  let interval = slices[0].interval
  doAssert slices.allIt(it.interval == interval), "All gas gain time slices must " &
    "be for the same length of " & $interval & " min in this case!"
  let baseName = group.name / "gasGainSlices"
  var dset = h5f.create_dataset(baseName & $(interval.round.int),
                                slices.len,
                                dtype = GasGainIntervalResult,
                                overwrite = true,
                                filter = filter)
  dset[dset.all] = slices
  dset.attrs["applied Cut for gas gain"] = cutFormula
  ## create a
  if baseName in h5f:
    doAssert h5f.delete(baseName), "Could not delete hardlink/dataset " & $baseName
  h5f.create_hardlink(dset.name, baseName)

proc applyGasGainCut(h5f: H5File, group: H5Group): seq[int] =
  ## Performs the cuts, which are used to select the events which we use
  ## to perform the polya fit + gas gain determination. This is based on
  ## C. Krieger's thesis. The cuts are extracted from the `getGasGain.C`
  ## macro found in `resources`.
  const cut_rms_trans_low = 0.1
  const cut_rms_trans_high = 1.5
  const cut_hits_high = 500.0
  result = cutOnProperties(h5f,
                           group,
                           crSilver,
                           ("rmsTransverse", cut_rms_trans_low, cut_rms_trans_high),
                           ("hits", 0.0, cut_hits_high))

proc getGasGainCutFormula(): FormulaNode =
  ## the actual cut formula used for the gas gain input data to clean it
  ## of events, which are deemed too irregular.
  const cut_rms_trans_low = 0.1
  const cut_rms_trans_high = 1.5
  result = f{float -> bool:
    `rmsTransverse` >= cut_rms_trans_low and
    `rmsTransverse` <= cut_rms_trans_high and
    inRegion(df["centerX"][idx], df["centerY"][idx], crSilver) and
    `hits` < 500}

proc applyGasGainCut*(df: DataFrame): DataFrame =
  ## Performs the cuts, which are used to select the events which we use
  ## to perform the polya fit + gas gain determination. This is based on
  ## C. Krieger's thesis. The cuts are extracted from the `getGasGain.C`
  ## macro found in `resources`.
  doAssert "rmsTransverse" in df
  doAssert "centerX" in df
  doAssert "centerY" in df
  result = df
  result["passIdx"] = arange(0, result.len)
  let cutFormula = getGasGainCutFormula()
  result = result.filter(cutFormula)

proc gasGainHistoAndFit(data: seq[float], bin_edges: seq[float],
                        chipNumber, runNumber: int,
                        plotPath: string,
                        gasGainInterval = none[GasGainIntervalData](),
                        useTeX = false):
                          (seq[int], FitResult) =
  ## performs the binning of `data` (should be charges in electron converted from the
  ## ToT values in `calcGasGain`) according to `totBins` and fits the
  # compute histogram according to the unequal width bin edges
  let (binned, _) = data.histogram(bins = bin_edges)
  # ``NOTE: remove last element from bin_edges to have:``
  # ``bin_edges.len == binned.len``
  let bin_edges = bin_edges[0 .. ^2]
  # given binned histogram, fit polya
  info "Fitting polya to run: ", runNumber, " and chip: ", chipNumber
  let fitResult = fitPolyaNim(bin_edges,
                              binned.asType(float64),
                              chipNumber, runNumber)
  # create plots if desired
  info "Plotting polya of run: ", runNumber, " and chip: ", chipNumber, " to ", plotPath
  plotGasGain(bin_edges, binned.asType(float64),
              fitResult.x, fitResult.y,
              fitResult.xMin, fitResult.xMax,
              G_fit = fitResult.pRes[1],
              chiSq = fitResult.redChiSq,
              chipNumber, runNumber,
              pathPrefix = plotPath,
              gasGainInterval = gasGainInterval,
              useTeX = useTeX)
  result = (binned, fitResult)

proc readGasGainDf*(h5f: H5File, grp: string,
                    addDsets: seq[string]): DataFrame =
  var dfChip = newDataFrame()
  var dfAll = newDataFrame()
  h5f.readDsets(dfChip, concat(@["eventNumber"], addDsets), grp)
  h5f.readDsets(dfAll, @["timestamp", "eventNumber"], grp.parentDir)
  result = innerJoin(dfChip, dfAll, "eventNumber")

iterator iterGainSlices*(df: DataFrame,
                         interval, minInterval: float): (GasGainIntervalData, Slice[int]) =
  ## NOTE: the input has to be filtered by the gas gain cuts! (The indices given in
  ## the slice range will then correspond to the indices of that *filtered* DF!)
  let tstamps = df["timestamp"].toTensor(float)
  # determine the start time
  var tStart = tstamps[0]
  let tStop = tstamps[tstamps.size - 1]
  var idxOld = 0
  var idx = 0
  for i in 0 ..< tstamps.size:
    if (tstamps[i] - tStart) >= (interval * 60.0): # convert interval in minutes to seconds
      if (tStop - tstamps[i]) < (minInterval * 60.0):
        # break out of for loop, thereby last slice will be from `tStart` (2nd to last real slice)
        # to end of data
        break

      ## TODO:
      ## We might want to think about making interval adaptive? That is to demand
      ## we have a certain number of entries in `data` to have good enough statistics
      ## (useful for calibration runs for outer chips!)
      # run over bin range
      let g = initInterval(idx, interval, minInterval, tStart.float, tstamps[i])
      yield (g, idxOld ..< min(i, tstamps.size - 1))
      tStart = tstamps[i]
      idxOld = i
      inc idx
  let g = initInterval(idx, interval, minInterval, tStart, tstamps[tstamps.size - 1])
  yield (g, idxOld ..< tstamps.size.int)

iterator iterGainSlicesDF*(df: DataFrame,
                           gainSlices: seq[GasGainIntervalResult]): DataFrame =
  ## given a DF, which must contain the following:
  ## - eventNumber
  ## which is going to be filtered according to each gas gain slice.
  ## Important: Obviously the event numbers are slightly specific to
  ## for each chip, so take care to choose the correct gas gain slices
  ## for each chip
  ##
  ## NOTE: The yielded value's slice range corresponds to the indices of the chip's dataset
  ## after the application of the gas gain cuts!
  var
    lowEv = 0
    highEv: int
  let lastEv = df["eventNumber"][df.high, int]
  for i, gasGainInterval in gainSlices:
    ## NOTE: to calibrate the energy we do ``not`` care about the index of the start
    ## of the slice, but only the event numbers, because these tell us if the
    ## which events start and stop
    lowEv = if i == 0: 0 else: highEv + 1 # start from last high value
    highEv = if i == gainSlices.high: lastEv else: gasGainInterval.sliceStopEvNum
    let sliceDf = df.filter(f{int -> bool: `eventNumber` >= lowEv and
                                           `eventNumber` <= highEv})
    yield sliceDf

proc deleteAllAttrStartingWith(dset: H5DataSet, start: string) =
  ## deletes all attributes starting with string
  dset.attrs.read_all_attributes()
  for key in dset.attrs.attr_tab.keys:
    if key.startsWith(start):
      echo "Deleting ", key
      discard dset.deleteAttribute(key)

proc gasGainDfAndCharges(h5f: H5File, group: string): (DataFrame, seq[seq[float]]) =
  ## Reads the relevant fields to compute the indices participating in the gas gain
  ## calculation, applies the cuts and returns the DF as well as the correct charge
  ## values of all events passing the cuts.
  ## The charge data is *not* flattened, as to allow extracting correct slices for the
  ## gas gain slicing!
  let df = h5f.readGasGainDf(group, @["centerX", "centerY", "rmsTransverse", "hits"])
    .applyGasGainCut()
  let chargeDset = h5f[(group / "charge").dset_str]
  let chFull = chargeDset.readVlen(float)
  var ch = newSeq[seq[float]](df.len) # nested seq of all charge values
  for i in 0 ..< df.len:
    ch[i] = chFull[df["passIdx", i, int]] # extract passing event `i`
  result = (df, ch)

proc handleGasGainFullRun*(h5f: H5File,
                           runNumber, chipNumber: int,
                           bin_edges: seq[float], group, plotPath: string,
                           useTeX = false): seq[GasGainIntervalResult] =
  let cutFormula = $getGasGainCutFormula()
  let (df, chs) = h5f.gasGainDfAndCharges(group)
  let (binned, fitResult) = gasGainHistoAndFit(chs.flatten, bin_edges, #totBins,
                                                                       # chipName,
                                               chipNumber, runNumber,
                                               plotPath,
                                               useTeX = useTeX)
  var grp = h5f[group.grp_str]
  let runGroup = h5f[group.parentDir.grp_str]
  let tstampDset = h5f[(group.parentDir / "timestamp").dset_str]
  let tStart = if "runStart" in grp.attrs: grp.attrs["runStart", string].parseTOSDateString()
               else: tstampDset[0, int].fromUnix()
  let tStop = if "runStop" in grp.attrs: grp.attrs["runStop", string].parseTOSDateString()
               else: tstampDset[tstampDset.shape[0],
                                int].fromUnix()
  let gasGainSingle = GasGainIntervalData(idx: 0, interval: 0.0, tStart: tStart.toUnix().int,
                                          tStop: tStop.toUnix().int)
  let ggRes = initGasGainIntervalResult(gasGainSingle, fitResult,
                                        binned, bin_edges,
                                        0 ..< df["passIdx", int].max,
                                        0, runGroup.attrs["numEventsStored", int])
  let chargeDset = h5f[(group / "charge").dset_str]
  h5f.writePolyaDsets(grp, chargeDset, binned, bin_edges, fitResult, cutFormula)
  result = @[ggRes]

  #h5f.writeGasGainSliceData(group, @[ggRes], cutFormula)
  #
  ## write gas gain slice attribute
  #h5f.writeGasGainAttributes(group, sliceCount = 0, interval = 0)

proc handleGasGainSlice*(h5f: H5File,
                         runNumber, chipNumber: int,
                         interval, minInterval: float,
                         bin_edges: seq[float],
                         group, plotPath: string,
                         useTeX = false): seq[GasGainIntervalResult] =
  # read required data for gas gain cuts & to map clusters to timestamps
  let cutFormula = $getGasGainCutFormula()
  let (df, chs) = h5f.gasGainDfAndCharges(group)
  var sliceCount = 0
  for (gasGainInterval, slice) in iterGainSlices(df, interval, minInterval):
    echo $gasGainInterval
    let chSlice = chs[slice].flatten
    let (binned, fitResult) = gasGainHistoAndFit(
      chSlice, bin_edges, #totBins, # chipName,
      chipNumber, runNumber,
      plotPath,
      gasGainInterval = some(gasGainInterval),
      useTeX = useTeX)

    let grp = h5f[group.grp_str]
    let chargeDset = h5f[(group / "charge").dset_str]
    h5f.writePolyaDsets(grp, chargeDset, binned, bin_edges, fitResult,
                        cutFormula, some(gasGainInterval))

    result.add initGasGainIntervalResult(gasGainInterval, fitResult,
                                         binned, bin_edges,
                                         slice,
                                         df["eventNumber", slice.a, int],
                                         df["eventNumber", slice.b, int])
    inc sliceCount

proc calcGasGain*(h5f: H5File, grp: string,
                  runNumber, chipNumber: int,
                  fileInfo: FileInfo,
                  fullRunGasGain: bool,
                  interval, minInterval: float, plotPath: string,
                  useTeX = false,
                  overwrite = false,
                  printRunPeriod = false) =
  ## Handles gas gain calculation for `runNumber`, given by `grp` (different chips!)
  if isDone(h5f, grp, rfOnlyGasGain, overwrite):
    echo &"INFO Gas gain calculation for run {runNumber} already exists. Skipping. Force via `--overwrite`."
    return

  ## XXX: Probably make the below adjustable!
  const
    hitLow = 2 # start at 2 to avoid noisy pixels w/ 1 hit
    hitHigh = 302
    binWidth = 3
  let totBins = arange(hitLow, hitHigh, binWidth).mapIt(it.float + 0.5)

  # now can start reading, get the group containing the data for this chip
  var group = h5f[grp.grp_str]
  let chipName = group.attrs["chipName", string]
  let ((a, b, c, t), runPeriod) = getTotCalibParameters(chipName, runNumber)
  if printRunPeriod:
    info "Using calibrations for chip: ", chipName, "(# ", chipNumber, ") from run period: ", runPeriod
  # get bin edges by calculating charge values for all TOT values at TOT's bin edges
  let capacitance = fileInfo.timepix.getCapacitance()
  let bin_edges = mapIt(totBins, calibrateCharge(it, capacitance, a, b, c, t))
  block Cleanup: # delete old indices for intervals
    let chargeDset = h5f[(grp / "charge").dset_str]
    chargeDset.deleteAllAttrStartingWith("interval_")
  # get all charge values as seq[seq[float]], ``then`` apply the `passIdx`
  # and only flatten ``after`` that
  var gasGainSliceData: seq[GasGainIntervalResult]
  if not fullRunGasGain:
    gasGainSliceData = h5f.handleGasGainSlice(runNumber, chipNumber,
                                              interval, minInterval,
                                              bin_edges, grp, plotPath,
                                              useTeX = useTeX)
  else:
    gasGainSliceData = h5f.handleGasGainFullRun(runNumber, chipNumber,
                                                bin_edges, grp, plotPath,
                                                useTeX = useTeX)
  # write gas gain slice data and attributes
  let cutFormula = $getGasGainCutFormula()
  h5f.writeGasGainSliceData(group, gasGainSliceData, cutFormula)
  h5f.writeGasGainAttributes(group, gasGainSliceData.len, interval.round.int, runPeriod)

  # write the flag we are done
  writeFlag(h5f, grp, rfOnlyGasGain)

proc calcGasGain*(h5f: H5File, runNumber: int,
                  interval, minInterval: float, fullRunGasGain: bool,
                  plotPath = "",
                  useTeX = false, overwrite = false) =
  ## fits the polya distribution to the charge values and writes the
  ## fit parameters (including the gas gain) to the H5 file
  ## `interval` is the time interval width on which we apply the binning
  ## of each run in minutes
  var chipBase = recoDataChipBase(runNumber)
  let plotPath = if plotPath.len > 0: plotPath else: h5f.attrs[PlotDirPrefixAttr, string]
  # get the group from file
  info "Calulating gas gain for run: ", runNumber
  let fileInfo = h5f.getFileInfo()
  var printRunPeriod = false
  var runPeriodPrinted: set[uint8]
  for run, chip, grp in chipGroups(h5f):
    doAssert chipBase in grp
    if isDone(h5f, grp, rfOnlyGasGain, overwrite):
      echo &"INFO Gas gain calculation for run {run} already exists. Skipping. Force via `--overwrite`."
      continue
    doAssert run == runNumber
    if chip.uint8 notin runPeriodPrinted:
      runPeriodPrinted.incl chip.uint8
      printRunPeriod = true
    else: printRunPeriod = false
    h5f.calcGasGain(grp, run, chip, fileInfo, fullRunGasGain, interval, minInterval, plotPath,
                    useTeX = useTeX, overwrite = overwrite,
                    printRunPeriod = printRunPeriod)

proc writeFeFitParameters(dset: var H5DataSet,
                          popt, popt_E: seq[float],
                          pErr, pErr_E: seq[float]) =
  ## writes the fit parameters obtained via a call to `scipy.curve_fit` to
  ## the attributes of the dataset ``dset``, which should be the corresponding
  ## FeSpectrum[Charge] dataset
  proc writeAttrs(dset: var H5DataSet, name: string,
                  p: seq[float], c: seq[float]) =
    for i in 0 .. p.high:
      dset.attrs[name & $i] = p[i]
      dset.attrs[name & "Err_" & $i] = c[i]
  dset.writeAttrs("p", popt, pErr)
  dset.writeAttrs("p_E_", popt_E, pErr_E)

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

proc writeFeDset(h5f: H5File,
                 group, suffix: string,
                 feSpec: FeSpecFitData): (H5DataSet, H5DataSet) =
  ## writes the dataset for the given FeSpectrum (either FeSpectrum
  ## or FeSpectrumCharge), i.e. the actual data points for the plots
  ## created by the Python functions
  let
    feCounts = feSpec.hist
    feBins = feSpec.binning
    feFitX = feSpec.xFit
    feFitY = feSpec.yFit
  template createWriteDset(x, y: seq[float], name: string): untyped =
    var dset = h5f.create_dataset(group / name,
                                  (x.len, 2),
                                  dtype = float,
                                  overwrite = true,
                                  filter = filter)
    let data = zip(x, y) --> map(@[it[0], it[1]]) --> to(seq[seq[float]])
    dset[dset.all] = data
    dset
  result[0] = createWriteDset(feBins, feCounts, "FeSpectrum" & $suffix & "Plot")
  result[1] = createWriteDset(feFitX, feFitY, "FeSpectrum" & $suffix & "PlotFit")

proc buildTextForFeSpec*(feSpec: FeSpecFitData,
                         ecData: EnergyCalibFitData,
                         isPixel = true,
                         isFadc = false): seq[string] =
  if isPixel:
    result.add &"μ = {feSpec.k_alpha:.1f} pix"
    result.add &"{ecData.aInv:.1f} eV / pix"
  elif isFadc:
    result.add &"μ = {feSpec.k_alpha:.1f} V"
    result.add &"{ecData.aInv:.1f} eV / V"
  else:
    result.add &"μ = {feSpec.k_alpha:.1f}e3 e^-"
    result.add &"{ecData.aInv:.1f} eV / 1000 e^-"
  result.add &"σ = {feSpec.sigma_kalpha / feSpec.k_alpha * 100.0:.2f} %"
  result.add &"χ²/dof = {feSpec.chiSq / feSpec.nDof.float:.2f}"

proc fitToFeSpectrum*(h5f: H5File, runNumber, chipNumber: int,
                      fittingOnly = false, useTeX = false,
                      outfiles: seq[string] = @[],
                      writeToFile = true,
                      plotPath = "",
                      overwrite = false
                     ) =
  ## The optional `writeToFile` flag can be set to `false` if this proc
  ## is to be called if the H5 file is opened read only to reproduce
  ## the results, but not rewrite them to the file (used in `plotData`)

  proc extractAndWriteAttrs(h5f: H5File,
                            dset: var H5DataSet,
                            scaling: float,
                            feSpec: FeSpecFitData,
                            ecData: EnergyCalibFitData,
                            texts: seq[string],
                            key: string,
                            suffix = "") =
    dset.writeFeFitParameters(feSpec.pRes, ecData.pRes, feSpec.pErr, ecData.pErr)
    writeEnergyPerAttrs(dset, key,
                        scaling,
                        ecData.pRes[0],
                        ecData.pErr[0])
    dset.writeTextFields(texts)
    var (grp1, grp2) = h5f.writeFeDset(dset.parent, suffix, feSpec = feSpec)
    grp1.copy_attributes(dset.attrs)
    grp2.copy_attributes(dset.attrs)

  ## ==============================
  ##         Pixel spectrum
  ## ==============================
  let groupName = recoPath(runNumber, chipNumber).string
  if not isDone(h5f, groupName, $rfOnlyFeSpec & "Pixel", overwrite):
    # get the fe spectrum for the run
    var feDset = h5f[(groupName / "FeSpectrum").dsetStr]
    let plotPath = if plotPath.len == 0: h5f.attrs[PlotDirPrefixAttr, string] else: plotPath
    let feData = feDset[int64]
    info "Fit pixel spectrum of run: " & $runNumber & " and chip: " & $chipNumber
    let feSpec = fitFeSpectrum(feData)
    info "Fit energy calibration of run: " & $runNumber & " and chip: " & $chipNumber
    let ecData = fitEnergyCalib(feSpec, isPixel = true)
    let texts = buildTextForFeSpec(feSpec, ecData)
    info "Plot pixel spectrum of run: " & $runNumber & " and chip: " & $chipNumber
    plotFeSpectrum(feSpec, runNumber, chipNumber,
                   texts, isPixel = true,
                   pathPrefix = plotPath,
                   useTeX = useTeX)
    info "Plot energy calibration of run: " & $runNumber & " and chip: " & $chipNumber
    plotFeEnergyCalib(ecData, runNumber, isPixel = true,
                      pathPrefix = plotPath,
                      useTeX = useTeX)

    const eVperPixelScaling = 1e3
    if writeToFile:
      h5f.extractAndWriteAttrs(feDset, eVperPixelScaling, feSpec, ecData, texts, "eV_per_pix")
      writeFlag(h5f, groupName, $rfOnlyFeSpec & "Pixel")

  ## ==============================
  ##         Charge spectrum
  ## ==============================
  # run might not have ``totalCharge`` dset, if no ToT calibration is available,
  # but user wishes Fe spectrum fit to # hit pixels
  if not isDone(h5f, groupName, $rfOnlyFeSpec & "Charge", overwrite) and
     h5f.hasTotalChargeDset(runNumber, chipNumber):
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
      totChDset = h5f.write_dataset(groupName / "FeSpectrumCharge", totChSpec, overwrite = true)
    info "Fit charge spectrum of run: " & $runNumber & " and chip: " & $chipNumber
    let feSpecCharge = fitFeSpectrumCharge(totChSpec)
    info "Fit charge energy calibration of run: " & $runNumber & " and chip: " & $chipNumber
    let ecDataCharge = fitEnergyCalib(feSpecCharge, isPixel = false)
    let textsCharge = buildTextForFeSpec(feSpecCharge, ecDataCharge, isPixel = false)
    info "Plot charge spectrum of run: " & $runNumber & " and chip: " & $chipNumber
    plotFeSpectrum(feSpecCharge, runNumber, chipNumber,
                   textsCharge, isPixel = false,
                   pathPrefix = plotPath,
                   useTeX = useTeX)
    info "Plot charge energy calibration of run: " & $runNumber & " and chip: " & $chipNumber
    plotFeEnergyCalib(ecDataCharge, runNumber, isPixel = false,
                      pathPrefix = plotPath,
                      useTeX = useTeX)

    # given resCharge, need to write the result of that fit to H5 file, analogous to
    # `writeFitParametersH5` in Python
    const keVPerElectronScaling = 1e-3
    if writeToFile:
      h5f.extractAndWriteAttrs(totChDset, keVPerElectronScaling,
                               feSpecCharge,
                               ecDataCharge,
                               textsCharge,
                               "keV_per_electron",
                               "Charge")
      writeFlag(h5f, groupName, $rfOnlyFeSpec & "Charge")
  else:
    echo "Warning: Either already performed 55Fe fit or `totalCharge` dataset does not exist in file. No fit to " &
      "charge Fe spectrum will be performed!"

  ## ==============================
  ##         FADC spectrum
  ## ==============================
  # `--only_fadc` must have been run nowadays. We only compute the `minVals` in the actual FADC reconstruction
  # now!
  let fadcGroup = fadcRecoPath(runNumber)
  if h5f.hasFadc(runNumber) and
     not isDone(h5f, fadcGroup, $rfOnlyFeSpec & "FADC", overwrite) and
     minValBasename(runNumber) in h5f:
    ## Also fit the 55Fe spectrum in the FADC data
    ## XXX: in the future `minvals` may be replaced by a "charge" equivalent based on
    ## an integral of the signal!
    var minValDset = h5f[minValBasename(runNumber).dset_str]
    let minVal = minValDset[float64]
    info "Fit FADC spectrum of run: " & $runNumber
    let feSpecFadc = fitFeSpectrumFadc(minVal)
    info "Fit FADC energy calibration of run: " & $runNumber
    let ecDataFadc = fitEnergyCalib(feSpecFadc, isPixel = false)
    let textsFadc = buildTextForFeSpec(feSpecFadc, ecDataFadc, isPixel = false, isFadc = true)
    info "Plot FADC spectrum of run: " & $runNumber
    plotFeSpectrum(feSpecFadc, runNumber, 3, texts = textsFadc,
                   isPixel = false, isFadc = true, pathPrefix = plotPath,
                   useTeX = useTeX)
    info "Plot FADC energy calibration of run: " & $runNumber
    ## XXX: maybe fit with an offset in this case, due to high activation threshold?
    plotFeEnergyCalib(ecDataFadc, runNumber, isPixel = false, isFadc = true,
                      pathPrefix = plotPath,
                      useTeX = useTeX)
    if writeToFile:
      h5f.extractAndWriteAttrs(minValDset,
                               1e3, ## XXX: fix this number!
                               feSpecFadc,
                               ecDataFadc,
                               textsFadc,
                               "keV_per_V",
                               "Fadc")
      writeFlag(h5f, fadcGroup, $rfOnlyFeSpec & "FADC")

proc fitSpectraBySlices(h5f: H5File,
                        chipGrp: H5Group, chip: int,
                        calib, calibErr, gainVals: var seq[float],
                        gasGainSlices: seq[GasGainIntervalResult]) =
  ## fits the Fe spectra by gas gain time slicesand returns a seq of fit parameters
  ## and gas gain values
  var
    lowEv = 0
    highEv: int
  let df = h5f.readGasGainDf(chipGrp.name,
                             @["totalCharge", "centerX", "centerY", "eccentricity",
                               "rmsTransverse"])
    .cutFeSpectrum()
  let lastEv = df["eventNumber"][df.high, int]
  for i, gasGainInterval in gasGainSlices:
    ## NOTE: to calibrate the energy we do ``not`` care about the index of the start
    ## of the slice, but only the event numbers, because these tell us if the
    ## which events start and stop
    lowEv = if i == 0: 0 else: highEv + 1 # start from last high value
    highEv = if i == gasGainSlices.high: lastEv else: gasGainInterval.sliceStopEvNum
    let sliceDf = df.filter(f{int -> bool: `eventNumber` >= lowEv and
                                           `eventNumber` <= highEv})
    if sliceDf.len < 100:
      echo "INFO: Skipping slice ", gasGainInterval, " due to only ", df.len, " clusters present"
      continue
    let chSlice = sliceDf["totalCharge", float]
    let feSpec = fitFeSpectrumCharge(chSlice.clone.toRawSeq)
    let ecData = fitEnergyCalib(feSpec, isPixel = false)
    let (popt, pErr) = (ecData.pRes[0], ecData.pErr[0])
    const keVPerElectronScaling = 1e-3
    # add fit results
    let aInv = 1.0 / popt * kevPerElectronScaling
    calib.add aInv * 1e6
    calibErr.add (aInv * pErr / popt) * 1e6 ## previoulsy we artificially enlarged this to 1e8
    gainVals.add gasGainInterval.G

func raiseGasGainSliceNumberError() =
  raise newException(ValueError, "Cannot compute the gas gain vs. charge calibration fit for only " &
    "a single gas gain slice interval. You either need more calibration runs for a fit _or_ (if you know " &
    "what you are doing) you can adjust the gas gain slice in the config.toml file via the " &
    "`gasGainInterval` field under `[Calibration]` to a shorter time. However, make sure there are " &
    "at least 10,000 events per slice!\n" &
    "WARNING: After performing a change in that parameter, you MUST rerun `--only_gas_gain` for all " &
    "your data files (with `--overwrite`) before rerunning `--only_gain_fit! Do not forget to adjust " &
    "`minimumGasGainInterval` accordingly.")

proc performChargeCalibGasGainFit*(h5f: H5File,
                                   interval: float,
                                   gcKind: GasGainVsChargeCalibKind = gcMean,
                                   useTeX = false,
                                   plotPath = "",
                                   overwrite = false) =
  ## performs the fit of the charge calibration factors vs gas gain fit
  ## Assumes:
  ## - h5f points to a h5 file of `runType == rtCalibration`
  ## - for all runs the Fe spectrum was calculated and fitted
  ## writes the resulting fit data to the ingridDatabase

  ## NOTE: the `gcKind` only takes effect for input H5 files, which are already
  ## reconstructed using the new sliced gas gain calibration!
  # iterate over all runs, extract center chip grou
  if isDone(h5f, recoGroupGrpStr().string, rfOnlyGainFit, overwrite):
    echo &"INFO Charge calibration vs. gas gain fit already performed. Skipping. Force via `--overwrite`."
    return

  let plotPath = if plotPath.len > 0: plotPath else: h5f.attrs[PlotDirPrefixAttr, string]
  var
    calib = newSeq[float64]()
    calibErr = newSeq[float64]()
    gainVals = newSeq[float64]()
    runPeriods = newSeq[string]()
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
    centerChip = 0
  for run, chip, grp in chipGroups(h5f):
    if chip != centerChip: continue ## only look at the center chip
    let centerChipGrp = h5f[grp.grp_str]
    case rfKind
    of rfOldTos, rfSrsTos:
      if centerChipName.len == 0:
        centerChipName = centerChipGrp.attrs["chipName", string]
    of rfNewTos, rfUnknown: ## XXX: `rfUnknown` corresponds to Tpx3 now!
      if ("chip_" & $centerChip) notin centerChipGrp.name:
        # skip this chip, since not the center chip
        echo "Skipping group ", centerChipGrp.name
        continue
      else:
        echo "\t taking group ", centerChipGrp.name
    let runPeriod = findRunPeriodFor(centerChipName, run)
    # read the chip name
    # given correct group, get the `charge` and `FeSpectrumCharge` dsets
    var
      chargeDset = h5f[(centerChipGrp.name / "charge").dset_str]
      feChargeSpec = h5f[(centerChipGrp.name / "FeSpectrumCharge").dset_str]
    var
      gain: float
    if grp / "gasGainSlices" notin h5f:
      echo "WARNING: input file ", h5f.name, " does not yet contain the `gasGainSlices`" &
        " dataset. This means we use the single `G` value of the `charge` dataset for " &
        " the computation!"
      gain = chargeDset.attrs["G", float64]
    else:
      let gasGainSlices = h5f[grp / "gasGainSlices" & $(interval.round.int),
                              GasGainIntervalResult]
      ## default way for now: calculate the mean
      case gcKind
      of gcMean:
        gain = gasGainSlices.mapIt(it.G).mean
        let
          keVPerE = feChargeSpec.attrs["keV_per_electron", float64]
          dkeVPerE = feChargeSpec.attrs["d_keV_per_electron", float64]
        calib.add keVPerE * 1e6
        calibErr.add dkeVPerE * 1e6 ## previously we artificially enlarged this to 1e8
        gainVals.add gain
      of gcIndividualFits:
        # fit one spectrum and energy calibration per time slice and add data
        h5f.fitSpectraBySlices(centerChipGrp, centerChip,
                               calib, calibErr, gainVals, gasGainSlices)
      of gcNone:
        doAssert false, "You shouldn' use `gcNone` anymore! Use `gcMean`"
    runPeriods.add runPeriod

  if gainVals.len == 1:
    raiseGasGainSliceNumberError()

  # increase smallest errors to lower 1 percentile errors
  let perc1 = calibErr.percentile(1)
  doAssert perc1 < mean(calibErr), "perc1 = " & $perc1 & " vs mean = " & $(mean(calibErr)) &
    "Does your input file only have a single run with one slice?"
  for x in mitems(calibErr):
    if x < perc1:
      x = perc1
  doAssert calibErr.allIt(it >= perc1)

  let runPeriodsUnique = runPeriods.deduplicate
  doAssert runPeriodsUnique.len == 1, "More than one run period found in input file " &
    $h5f.name & " : " & $runPeriodsUnique

  info "Fit charge calibration vs gas gain for file: " & $h5f.name
  let fitResult = fitChargeCalibVsGasGain(gainVals, calib, calibErr)
  writeCalibVsGasGain(gainVals, calib, calibErr, fitResult, centerChipName,
                      runPeriodsUnique[0], suffix = $gcKind)
  # and create the plot
  info "Plot charge calibration vs gas gain for file: " & $h5f.name
  plotGasGainVsChargeCalib(gainVals, calib, calibErr, fitResult,
                           pathPrefix = plotPath, useTeX = useTeX)
  # write the flag it is done
  writeFlag(h5f, recoGroupGrpStr().string, rfOnlyGainFit)

proc calcEnergyFromPixels*(h5f: H5File, runNumber: int, calib_factor: float, overwrite: bool) =
  ## proc which applies an energy calibration based on the number of hit pixels in an event
  ## using a conversion factor of unit eV / hit pixel to the run given by runNumber contained
  ## in file h5f
  ## throws:
  ##     HDF5LibraryError = in case a call to the H5 library fails, this might be raised
  var chipBase = recoDataChipBase(runNumber)
  # get the group from file
  for run, chip, grp in chipGroups(h5f, recoBase()):
    doAssert chipBase in grp
    if isDone(h5f, grp, rfOnlyEnergy, overwrite):
      echo &"INFO Energy calculation for run {run} already exists. Skipping. Force via `--overwrite`."
      continue
    doAssert runNumber == run
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
    var energy_dset = h5f.create_dataset(grp / "energyFromPixel", energy.len, dtype = float,
                                         overwrite = true,
                                         filter = filter)
    energy_dset[energy_dset.all] = energy
    # attach used conversion factor to dataset
    energy_dset.attrs["conversionFactorUsed"] = calib_factor
    # set the flag that this run is finished
    writeFlag(h5f, grp, rfOnlyEnergy)

proc calcEnergyFromCharge*(h5f: H5File, chipGrp: H5Group,
                           interval: float,
                           b, m: float,
                           gcKind: GasGainVsChargeCalibKind,
                           overwrite: bool) =
  ## performs the actual energy calculation, based on the fit parameters of
  ## `b` and `m` on the group `chipGrp`. This includes writing the energy
  ## dataset and attributes!
  # get the chip number from the attributes of the group
  if isDone(h5f, chipGrp.name, rfOnlyEnergyElectrons, overwrite):
    echo &"INFO Energy calculation from charge calibration for run and chip {chipGrp.name}",
          "already exists. Skipping. Force via `--overwrite`."
    return # nothing to do

  var mchpGrp = chipGrp
  let chipNumber = mchpGrp.attrs["chipNumber", int]

  # get dataset of charge (contains all fit parameters for the gas gain fits)
  var chargeDset = h5f[(mchpGrp.name / "charge").dset_str]

  ## TODO
  ## - write / reuse iterator for time slices like for gas gain intervals
  ## - read all total charges (argument for iterator)
  ## - iterator yields slices of total charge and correct gas gain? or just gasGainInterval
  ##   and we use GasGainIntervalData to read correct H5 attributes?
  ##   latter more complex but cleaner
  ## - calc energy on slices, append energy to `energies` sequence
  ## - only write at very end
  let df = h5f.readGasGainDf(chipGrp.name, @["totalCharge"])
  let rawCharge = h5f[mchpGrp.name / "totalCharge", float]
  let chs = df["totalCharge"].toTensor(float).clone.toRawSeq ## NOTE: WARNING raw seq is longer!!!
  var energy = newSeqOfCap[float](chs.len)
  var allIdx = initHashSet[int]()
  let gainSlices = h5f[chipGrp.name / "gasGainSlices" & $(interval.round.int),
                       GasGainIntervalResult]
  var
    lowEv = 0
    highEv: int
  let lastEv = df["eventNumber"][df.high, int]
  for i, gasGainInterval in gainSlices:
    ## NOTE: to calibrate the energy we do ``not`` care about the index of the start
    ## of the slice, but only the event numbers, because these tell us if the
    ## which events start and stop
    lowEv = if i == 0: 0 else: highEv + 1 # start from last high value
    highEv = if i == gainSlices.high: lastEv else: gasGainInterval.sliceStopEvNum
    let sliceDf = df.filter(f{int -> bool: `eventNumber` >= lowEv and
                                           `eventNumber` <= highEv})
    let chSlice = sliceDf["totalCharge", float]
    let gain = gasGainInterval.G
    # remove correction factor of 1e6
    let calibFactor = linearFunc(@[b, m], gain) * 1e-6
    # now calculate energy for all hits
    energy.add mapIt(chSlice, it * calibFactor)
  doAssert energy.len == chs.len, "Invalid calculation, only " & $energy.len &
    " energies for " & $(chs.len) & " cluster!"
  ## TODO: add something similar again?
  #doAssert allIdx.card == energy.len and allIdx.toSeq.max == energy.high

  when false:
    var gain: float
    if "G" in chargeDset.attrs:
      gain = chargeDset.attrs["G", float64]
    else:
      raise newException(Exception, "Gas gain not yet calculated for " &
        "run " & $mchpGrp.name)
    let totCharge = h5f[(mchpGrp.name / "totalCharge"), float64]
    # remove correction factor of 1e6
    let calibFactor = linearFunc(@[b, m], gain) * 1e-6
    # now calculate energy for all hits
    let energy = mapIt(totCharge, it * calibFactor)
    # create dataset for energy
  var energy_dset = h5f.create_dataset(mchpGrp.name / "energyFromCharge",
                                       energy.len, dtype = float,
                                       overwrite = true,
                                       filter = filter)
  energy_dset[energy_dset.all] = energy
  # attach used conversion factor to dataset
  energy_dset.attrs["Calibration function"] = "y = m * x + b"
  energy_dset.attrs["calib_b"] = b
  energy_dset.attrs["calib_m"] = m
  energy_dset.attrs["GasGainVsChargeCalibKind"] = $gcKind

  # store the interval length used to calculate the energies
  energy_dset.attrs["Gas gain interval length"] = interval
  energy_dset.attrs["Number of gas gain intervals"] = gainSlices.len

  writeFlag(h5f, chipGrp.name, rfOnlyEnergy)

proc calcEnergyFromCharge*(h5f: H5File, interval: float,
                           gcKind: GasGainVsChargeCalibKind = gcNone,
                           overwrite = false) =
  ## proc which applies an energy calibration based on the number of electrons in a cluster
  ## using a conversion factor of unit 1e6 keV / electron to the run given by runNumber contained
  ## in file h5f
  ## The `performChargeCalibGasGainFit` has to be run on the calibration dataset before,
  ## i.e. we need calibration fit results in the `ingridDatabase` for the center chip of
  ## the current detector
  ## throws:
  ##     HDF5LibraryError = in case a call to the H5 library fails, this might be raised
  # get the parameters from
  var chipName = ""
  # get the calibration factors for the current center chip
  var
    b: float
    m: float
  # get the group from file
  for run, chip, grp in chipGroups(h5f):
    # get the chip number from the attributes of the group
    if isDone(h5f, grp, rfOnlyEnergyElectrons, overwrite):
      echo &"INFO Energy calculation from charge calibration for run {run} already exists. Skipping. Force via `--overwrite`."
      continue # nothing to do
    # now can start reading, get the group containing the data for this chip
    var group = h5f[grp.parentDir.grp_str]
    echo "Energy from charge calibration for run ", run
    if chipName.len == 0:
      # get center chip name to be able to read fit parameters
      chipName = group.attrs["centerChipName", string]
      # get parameters during first iter...
      (b, m) = getCalibVsGasGainFactors(chipName, run, suffix = $gcKind)
    # now iterate over chips in this run
    let chipGrp = h5f[grp.grp_str]
    h5f.calcEnergyFromCharge(chipGrp, interval, b, m, gcKind, overwrite)
