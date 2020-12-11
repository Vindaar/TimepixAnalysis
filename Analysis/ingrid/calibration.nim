import sequtils, strutils, strformat
import os, ospaths
import future
import seqmath, fenv
import nimhdf5
import tables
import mpfit
import zero_functional
import math
import chroma

import tos_helpers
import helpers/utils
import ingrid_types
import ingridDatabase / [databaseRead, databaseDefinitions]
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

proc cutOnProperties*(h5f: H5FileObj,
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
    try:
      # TODO: this is a workaround for now. `centerX` is the name used for
      # TimepixAnalysis, but the `calibration-cdl.h5` data from Marlin uses
      # PositionX. I don't want to add a `year` field or something to this proc,
      # so for now we just depend on an exception.
      posX = h5f.readAs(group.name / "centerX", float)
      posY = h5f.readAs(group.name / "centerY", float)
    except KeyError:
      posX = h5f.readAs(group.name / "PositionX", float)
      posY = h5f.readAs(group.name / "PositionY", float)

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

proc cutOnProperties*(h5f: H5FileObj,
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

  info "Plotting facet for variables used for Fe spectrum for run: " & $runNumber
  plotFeSpectrumInfoFacet(pos_x, pos_y, ecc, rms_trans, hits,
                          runNumber = runNumber,
                          chipNumber = chipNumber,
                          pathPrefix = h5f.attrs[PlotDirPrefixAttr, string])

  # given this data, filter all events which don't conform
  let (eventSpectrum,
       hitsSpectrum,
       specIndices) = cutFeSpectrum([pos_x, pos_y, ecc, rms_trans], event_num, hits)
  let nEventsPassed = eventSpectrum.len
  # with the events to use for the spectrum
  info "Elements passing cut for Fe spectrum : ", nEventsPassed, " for run: ", runNumber

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

proc applyChargeCalibration*(h5f: var H5FileObj, runNumber: int,
                             toDelete: bool = false) =
  ## applies the charge calibration to the TOT values of all events of the
  ## given run
  # what we need:
  # TOT values of each run. Run them through calibration function with the TOT calibrated
  # values taken from the `InGridDatabase`
  var chipBase = recoDataChipBase(runNumber)
  # get the group from file
  info "Applying charge calibration for run: ", runNumber
  for grp in keys(h5f.groups):
    if chipBase in grp:
      # now can start reading, get the group containing the data for this chip
      var group = h5f[grp.grp_str]
      # get the chip number from the attributes of the group
      let chipNumber = group.attrs["chipNumber", int]
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
      let (a, b, c, t) = getTotCalibParameters(chipName, runNumber)
      #mapIt(it.mapIt(calibrateCharge(it.float, a, b, c, t)))
      var charge = newSeqWith(tots.len, newSeq[float]())
      var totalCharge = newSeq[float](sumTots.len)
      for i, vec in tots:
        # calculate charge values for individual pixels
        charge[i] = vec.mapIt((calibrateCharge(it.float, a, b, c, t)))
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

proc initGasGainIntervalResult(g: GasGainIntervalData,
                               fitResult: FitResult,
                               binned: seq[int], bin_edges: seq[float],
                               sliceIdx: Slice[int],
                               sliceStartEvNum, sliceStopEvNum: int): GasGainIntervalResult =
  result.idx       = g.idx
  result.interval  = g.interval
  result.tStart    = g.tStart
  result.tStop     = g.tStop
  result.N         = fitResult.pRes[0]
  result.G_fit     = fitResult.pRes[1]
  result.sliceStart = sliceIdx.a
  result.sliceStop = sliceIdx.b
  result.sliceStartEvNum = sliceStartEvNum
  result.sliceStopEvNum = sliceStopEvNum
  # also calculate the mean of the polya fit. This usually is *not* the same
  # as the `meanGain` (mean of data). `meanGainFit` is what we used errorneously
  # in the past!
  let meanGainFit = (zip(fitResult.x, fitResult.y) -->
                     map(it[0] * it[1]) -->
                     fold(0.0, a + it)) / fitResult.y.sum.float
  result.G_fitmean = meanGainFit
  # calculate the mean of the data histogram. This is what Krieger calls
  # `G_mean` and uses for the rest of the calculation!
  let meanGain = histMean(binned, bin_edges[0 .. ^2])
  result.G         = meanGain
  result.theta     = fitResult.pRes[2]
  result.redChiSq  = fitResult.redChiSq
  when false:
    var gidx: string
    if gasGainInterval.isSome:
      let g = gasGainInterval.get
      gidx = g.toAttrPrefix()
      d.attrs[g.toSliceStartAttr()] = g.tstart
      d.attrs[g.toSliceStopAttr()] = g.tstop
      d.attrs[gidx & "length"] = g.interval
      d.attrs[gidx[0 ..< ^1]] = $g

    d.attrs[gidx & "N"] = fitResult.pRes[0]
    d.attrs[gidx & "G_fit"] = fitResult.pRes[1]
    d.attrs[gidx & "G"] = meanGain
    d.attrs[gidx & "G_fitmean"] = meanGainFit
    d.attrs[gidx & "theta"] = fitResult.pRes[2]
    # TODO: get some errors from NLopt?
    #d.attrs["N_err"] = fitResult.pErr[0]
    #d.attrs["G_err"] = fitResult.pErr[1]
    #d.attrs["theta_err"] = fitResutl.pErr[2]
    d.attrs[gidx & "redChiSq"] = fitResult.redChiSq

proc writePolyaDsets(h5f: H5FileObj, group: H5Group,
                     chargeDset: H5DataSet,
                     binned: seq[int], bin_edges: seq[float],
                     fitResult: FitResult,
                     cutFormula: string,
                     gasGainInterval = none[GasGainIntervalData]()) =
  # all polya datasets are stored in `polya` subgroup
  const
    subgroup = "polya"
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
  ifDelete(pName)
  ifDelete(group.name / &"polya{dsetSuffix}")
  ifDelete(pFitName)
  ifDelete(group.name / &"polyaFit{dsetSuffix}")

  var polyaDset = h5f.create_dataset(pName,
                                     (binned.len, 2),
                                     dtype = float64,
                                     overwrite = true)
  var polyaFitDset = h5f.create_dataset(pFitName,
                                        (fitResult.x.len, 2),
                                        dtype = float64,
                                        overwrite = true)
  polyaDset[polyaDset.all] = zip(bin_edges, binned) -->
    map(@[it[0], it[1].float])
  polyaFitDset[polyaFitDset.all] = zip(fitResult.x, fitResult.y) -->
    map(@[it[0], it[1]])

  # now write resulting fit parameters as attributes
  template writeAttrs(d: H5DataSet, fitResult: typed): untyped =
    d.attrs["applied Cut for gas gain"] = cutFormula
  writeAttrs(chargeDset, fitResult)
  writeAttrs(polyaDset, fitResult)
  writeAttrs(polyaFitDset, fitResult)

proc writeGasGainSliceData(h5f: H5File, group: H5Group, slices: seq[GasGainIntervalResult]) =
  ## writes the information about the gas gain slices, including the fit results
  ## for each polya fit as one composite dataset to the H5 file.
  ##
  ## Dataset name: `gasGainSlices`
  var dset = h5f.create_dataset(group.name / "gasGainSlices",
                                slices.len,
                                dtype = GasGainIntervalResult,
                                overwrite = true)
  dset[dset.all] = slices

proc applyGasGainCut(h5f: H5FileObj, group: H5Group): seq[int] =
  ## Performs the cuts, which are used to select the events which we use
  ## to perform the polya fit + gas gain determination. This is based on
  ## C. Krieger's thesis. The cuts are extracted from the `getGasGain.C`
  ## macro found in `resources`.
  const cut_rms_trans_low = 0.1
  const cut_rms_trans_high = 1.5
  result = cutOnProperties(h5f,
                           group,
                           crSilver,
                           ("rmsTransverse", cut_rms_trans_low, cut_rms_trans_high))

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
                        gasGainInterval = none[GasGainIntervalData]()):
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
  info "Plotting polya of run: ", runNumber, " and chip: ", chipNumber
  plotGasGain(bin_edges, binned.asType(float64),
              fitResult.x, fitResult.y,
              fitResult.xMin, fitResult.xMax,
              G_fit = fitResult.pRes[1],
              chiSq = fitResult.redChiSq,
              chipNumber, runNumber,
              pathPrefix = plotPath,
              gasGainInterval = gasGainInterval)
  result = (binned, fitResult)

proc readGasGainDf*(h5f: H5FileObj, grp: string, chip: int,
                    addDsets: seq[string]): DataFrame =
  var dfChip = newDataFrame()
  var dfAll = newDataFrame()
  h5f.readDsets(dfChip, concat(@["eventNumber"], addDsets), grp)
  echo dfChip
  #dfChip["clusterIdx"] = arange(0, dfChip.len)
  h5f.readDsets(dfAll, @["timestamp", "eventNumber"], grp.parentDir)
  result = innerJoin(dfChip, dfAll, "eventNumber")

iterator iterGainSlices(df: DataFrame,
                        interval: float): (GasGainIntervalData, Slice[int]) =
  let tstamps = df["timestamp"].toTensor(float)
  # determine the start time
  var tStart = tstamps[0]
  var idxOld = 0
  var idx = 0
  for i in 0 ..< tstamps.size:
    if (tstamps[i] - tStart) >= (interval * 60.0): # convert interval in minutes to seconds
      ## TODO:
      ## We might want to think about making interval adaptive? That is to demand
      ## we have a certain number of entries in `data` to have good enough statistics
      ## (useful for calibration runs for outer chips!)
      # run over bin range
      let g = initInterval(idx, interval, tStart.float, tstamps[i])
      yield (g, idxOld ..< min(i, tstamps.size - 1))
      tStart = tstamps[i]
      idxOld = i
      inc idx
  let g = initInterval(idx, interval, tStart, tstamps[tstamps.size - 1])
  yield (g, idxOld ..< tstamps.size.int)

iterator iterGainSlicesFromAttrs*(h5f: H5File,
                                  group: H5Group): GasGainIntervalResult =
  let gainSlices = h5f[group.name / "gasGainSlices", GasGainIntervalResult]
  for g in gainSlices:
    yield g

proc deleteAllAttrStartingWith(dset: H5DataSet, start: string) =
  ## deletes all attributes starting with string
  dset.attrs.read_all_attributes()
  for key in dset.attrs.attr_tab.keys:
    if key.startsWith(start):
      echo "Deleting ", key
      discard dset.deleteAttribute(key)

proc calcGasGain*(h5f: var H5FileObj, runNumber: int,
                  interval: float) =
  ## fits the polya distribution to the charge values and writes the
  ## fit parameters (including the gas gain) to the H5 file
  ## `interval` is the time interval width on which we apply the binning
  ## of each run in minutes
  const
    hitLow = 2 # start at 2 to avoid noisy pixels w/ 1 hit
    hitHigh = 302
    binWidth = 3
  let totBins = arange(hitLow, hitHigh, binWidth).mapIt(it.float + 0.5)
  var chipBase = recoDataChipBase(runNumber)
  let plotPath = h5f.attrs[PlotDirPrefixAttr, string]

  ## TODO:
  ## - 1. make sure histogram works well with unequal bin widths (write test
  ##      comparing the output using ToT and charge!
  ##     DONE
  ## - 2. compute bins needed for unequal from totBins and applying calibrateCharge
  ##     DONE in test / above proc
  ## - 3. use readDsets to read the totalCharge dataset and hits, divide one another
  ##     cannot really do that, need individual pixel charge for more statistics

  ## - have to filter by gas gain cuts
  ## - need to merge timestamps and charge datasets (at least event number + indices)
  ## -

  ## TODO: make sure to delete all iterGainSlice attributes first? overwriting is super slow!
  # get the group from file
  info "Calulating gas gain for run: ", runNumber
  for grp in keys(h5f.groups):
    if chipBase in grp:
      # now can start reading, get the group containing the data for this chip
      var group = h5f[grp.grp_str]
      # get the chip number from the attributes of the group
      let chipNumber = group.attrs["chipNumber", int]
      let chipName = group.attrs["chipName", string]
      let (a, b, c, t) = getTotCalibParameters(chipName, runNumber)
      # get bin edges by calculating charge values for all TOT values at TOT's bin edges
      let bin_edges = mapIt(totBins, calibrateCharge(it, a, b, c, t))
      # get dataset of hits
      let chargeDset = h5f[(grp / "charge").dset_str]
      chargeDset.deleteAllAttrStartingWith("interval_")
      # get all charge values as seq[seq[float]], ``then`` apply the `passIdx`
      # and only flatten ``after`` that
      let chFull = chargeDset.readVlen(float)
      when not defined(fullRunGasGain):
        # read required data for gas gain cuts & to map clusters to timestamps
        let cutFormula = $getGasGainCutFormula()
        let df = h5f.readGasGainDf(grp, chipNumber, @["centerX", "centerY", "rmsTransverse", "hits"])
          .applyGasGainCut()
        let passIdx = df["passIdx"].toTensor(int).toRawSeq
        let chs = passIdx.mapIt(chFull[it])

        ## TODO
        ## instead of returning sliced data.
        ## - return only valid indices of the full dataset (a HSlice)
        ## - keep passed indices around as a hashset
        ## - walk hslice as a sequence, extract the data from the slice that is in
        ##   passed seq set as the data we will use
        ## CANNOT work, cause we end up with empty last sequence

        ## alternative: we have to just use the timestamps from the attributes instead
        ## and fill up the last timestamp to the end of the run
        ## need iterator which is similar, but returns slices based on attributes
        ## in gas gain calc have to add an attribute for the number of slices in total
        ##

        var sliceCount = 0
        var gasGainSliceData: seq[GasGainIntervalResult]
        for (gasGainInterval, slice) in iterGainSlices(df, interval):
          echo $gasGainInterval
          let chSlice = chs[slice].flatten
          let (binned, fitResult) = gasGainHistoAndFit(
            chSlice, bin_edges, #totBins, # chipName,
            chipNumber, runNumber,
            plotPath,
            gasGainInterval = some(gasGainInterval))

          h5f.writePolyaDsets(group, chargeDset, binned, bin_edges, fitResult,
                              cutFormula,
                              some(gasGainInterval))

          gasGainSliceData.add initGasGainIntervalResult(gasGainInterval, fitResult,
                                                         binned, bin_edges,
                                                         slice,
                                                         df["eventNumber", 0, int],
                                                         df["eventNumber", df.high, int])
          inc sliceCount

        chargeDset.attrs["numGasGainSlices"] = sliceCount
        h5f.writeGasGainSliceData(group, gasGainSliceData)
      else:
        let cutFormula = "No formula available, due to usage of `cutOnProperties`"
        let passIdx = applyGasGainCut(h5f, group)
        let chs = passIdx.mapIt(chFull[it])
        let (binned, fitResult) = gasGainHistoAndFit(chs, bin_edges, #totBins,
                                                                     # chipName,
                                                     chipNumber, runNumber,
                                                     plotPath)

        let tStart = if "runStart" in group.attrs: group.attrs["runStart", string].parseTOSDateString()
                     else: h5f[group.name / "timestamp", 0, int].fromUnix()
        let tStop = if "runStop" in group.attrs: group.attrs["runStop", string].parseTOSDateString()
                     else: h5f[group.name / "timestamp",
                               h5f[(group.name / "timestamp").dset_str].shape[0],
                               int].fromUnix()
        let gasGainSingle = GasGainIntervalData(idx: 0, interval: 0.0, tStart: tStart.toUnix(),
                                                tStop: tStop.toUnix())
        let ggRes = initGasGainIntervalResult(gasGainSingle, fitResult, binned, bin_edges,
                                              0 ..< passIdx.max,
                                              0, group.attrs["numEventsStored"])
        h5f.writePolyaDsets(group, chargeDset, binned, bin_edges, fitResult,
                            cutFormula)
        h5f.writeGasGainSliceData(group, @[ggRes])

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

proc writeFeDset(h5f: var H5FileObj,
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
                                  dtype = float)
    let data = zip(x, y) --> map(@[it[0], it[1]]) --> to(seq[seq[float]])
    dset[dset.all] = data
    dset
  result[0] = createWriteDset(feBins, feCounts, "FeSpectrum" & $suffix & "Plot")
  result[1] = createWriteDset(feFitX, feFitY, "FeSpectrum" & $suffix & "PlotFit")

proc buildTextForFeSpec(feSpec: FeSpecFitData,
                        ecData: EnergyCalibFitData,
                        isPixel = true): seq[string] =
  if isPixel:
    result.add &"mu = {feSpec.k_alpha:.1f} pix"
  else:
    result.add &"mu = {feSpec.k_alpha:.1f}e3 e^-"
  result.add &"{ecData.aInv:.1f} ev / pix"
  result.add &"sigma = {feSpec.sigma_kalpha / feSpec.k_alpha * 100.0:.2f} %"

proc fitToFeSpectrum*(h5f: var H5FileObj, runNumber, chipNumber: int,
                      fittingOnly = false, outfiles: seq[string] = @[],
                      writeToFile = true) =
  ## (currently) calls Python functions from `ingrid` Python module to
  ## perform fit to the `FeSpectrum` dataset in the given run number
  ## NOTE: due to calling Python functions, this proc is *extremely*
  ## inefficient!
  ## The optional `writeToFile` flag can be set to `false` if this proc
  ## is to be called if the H5 file is opened read only to reproduce
  ## the results, but not rewrite them to the file (used in `plotData`)

  # get the fe spectrum for the run
  let groupName = recoDataChipBase(runNumber) & $chipNumber
  var feDset = h5f[(groupName / "FeSpectrum").dsetStr]
  let plotPath = h5f.attrs[PlotDirPrefixAttr, string]
  let feData = feDset[int64]
  info "Fit pixel spectrum of run: " & $runNumber & " and chip: " & $chipNumber
  let feSpec = fitFeSpectrum(feData)
  info "Fit energy calibration of run: " & $runNumber & " and chip: " & $chipNumber
  let ecData = fitEnergyCalib(feSpec, isPixel = true)
  let texts = buildTextForFeSpec(feSpec, ecData)
  info "Plot pixel spectrum of run: " & $runNumber & " and chip: " & $chipNumber
  plotFeSpectrum(feSpec, runNumber, chipNumber,
                 texts, isPixel = true,
                 pathPrefix = plotPath)
  info "Plot energy calibration of run: " & $runNumber & " and chip: " & $chipNumber
  plotFeEnergyCalib(ecData, runNumber, isPixel = true,
                    pathPrefix = plotPath)

  proc extractAndWriteAttrs(h5f: var H5FileObj,
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

  const eVperPixelScaling = 1e3
  if writeToFile:
    h5f.extractAndWriteAttrs(feDset, eVperPixelScaling, feSpec, ecData, texts, "eV_per_pix")

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
    info "Fit charge spectrum of run: " & $runNumber & " and chip: " & $chipNumber
    let feSpecCharge = fitFeSpectrumCharge(totChSpec)
    info "Fit charge energy calibration of run: " & $runNumber & " and chip: " & $chipNumber
    let ecDataCharge = fitEnergyCalib(feSpecCharge, isPixel = false)
    let textsCharge = buildTextForFeSpec(feSpecCharge, ecDataCharge, isPixel = false)
    info "Plot charge spectrum of run: " & $runNumber & " and chip: " & $chipNumber
    plotFeSpectrum(feSpecCharge, runNumber, chipNumber,
                   textsCharge, isPixel = false,
                   pathPrefix = plotPath)
    info "Plot charge energy calibration of run: " & $runNumber & " and chip: " & $chipNumber
    plotFeEnergyCalib(ecDataCharge, runNumber, isPixel = false,
                      pathPrefix = plotPath)

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
  let plotPath = h5f.attrs[PlotDirPrefixAttr, string]
  let runPeriod = h5f.attrs[RunPeriodAttr, string]
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
      calibErr.add dkeVPerE * 1e8 # increase errors by factor 100 for better visibility
      gainVals.add gain

  # increase smallest errors to lower 10 percentile errors
  let perc10 = calibErr.percentile(1)
  doAssert perc10 < mean(calibErr)
  for x in mitems(calibErr):
    if x < perc10:
      x = perc10
  doAssert calibErr.allIt(it >= perc10)

  info "Fit charge calibration vs gas gain for file: " & $h5f.name
  let fitResult = fitChargeCalibVsGasGain(gainVals, calib, calibErr)
  writeCalibVsGasGain(gainVals, calib, calibErr, fitResult, centerChipName, runPeriod)
  # and create the plot
  info "Plot charge calibration vs gas gain for file: " & $h5f.name
  plotGasGainVsChargeCalib(gainVals, calib, calibErr, fitResult,
                           pathPrefix = plotPath)

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

proc calcEnergyFromCharge*(h5f: var H5FileObj, chipGrp: H5Group,
                           interval: float,
                           b, m: float) =
  ## performs the actual energy calculation, based on the fit parameters of
  ## `b` and `m` on the group `chipGrp`. This includes wirting the energy
  ## dataset and attributes!
  # get the chip number from the attributes of the group
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
  let df = h5f.readGasGainDf(chipGrp.name, chipNumber, @["totalCharge"])
  let chs = df["totalCharge"].toTensor(float).clone.toRawSeq ## NOTE: WARNING raw seq is longer!!!
  var energy = newSeqOfCap[float](chs.len)
  var intCount = 0
  var allIdx = initHashSet[int]()
  for gasGainInterval in h5f.iterGainSlicesFromAttrs(chipGrp):
    let slice = gasGainInterval.sliceStart .. gasGainInterval.sliceStop
    block SanityCheck:
      for s in slice:
        doAssert s notin allIdx
        allIdx.incl s
    let chSlice = chs[slice]
    let gain = gasGainInterval.G
    # remove correction factor of 1e6
    let calibFactor = linearFunc(@[b, m], gain) * 1e-6
    # now calculate energy for all hits
    energy.add mapIt(chSlice, it * calibFactor)
    inc intCount
  doAssert energy.len == chs.len, "Invalid calculation, only " & $energy.len &
    " energies for " & $(chs.len) & " cluster!"
  doAssert allIdx.card == energy.len and allIdx.toSeq.max == energy.high


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
                                       overwrite = true)
  energy_dset[energy_dset.all] = energy
  # attach used conversion factor to dataset
  energy_dset.attrs["Calibration function"] = "y = m * x + b"
  energy_dset.attrs["calib_b"] = b
  energy_dset.attrs["calib_m"] = m

  # store the interval length used to calculate the energies
  energy_dset.attrs["Gas gain interval length"] = interval
  energy_dset.attrs["Number of gas gain intervals"] = intCount

proc calcEnergyFromCharge*(h5f: var H5FileObj, interval: float) =
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

  for num, grp in runs(h5f):
    # now can start reading, get the group containing the data for this chip
    var group = h5f[grp.grp_str]
    echo "Energy from charge calibration for run ", num
    if chipName.len == 0:
      # get center chip name to be able to read fit parameters
      chipName = group.attrs["centerChipName", string]
      # get parameters during first iter...
      (b, m) = getCalibVsGasGainFactors(chipName, num.parseInt)

    # now iterate over chips in this run
    for chipGrp in items(h5f, start_path = group.name):
      if "fadc" in chipGrp.name:
        continue
      h5f.calcEnergyFromCharge(chipGrp, interval, b, m)
