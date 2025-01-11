import std / [strutils, os, sequtils, strformat, random]
import pkg / [nimhdf5, seqmath, ggplotnim, arraymancer]

import cdl_cuts, hdf5_utils, geometry, ggplot_utils, cut_utils
import ../ingrid_types
import helpers/utils

import ../projectDefs

import sugar

proc readToAProbabilities*(pathtoToA:string): (DataFrame, seq[int]) =
  let df = readCsv(pathtoToA,sep=' ') 
  let energies = df.getKeys().filterIt(it != "bin").mapIt(it.parseInt)
  result = (df, energies)

proc morph*(df: DataFrame, energy: float, offset = 1): (Tensor[float], Tensor[float]) =
  ## generates a distribution for the appropriate energy `energy` between the
  ## distribution below and above `energy` using linear interpolation
  ## DF needs to have columns:
  ## - "Hist": counts of distributions
  ## - "Bins": bin edges of distributions
  ## - "Dset": name of the target filter combination

  let lineEnergies = getXrayFluorescenceLines()
  let idx = max(lineEnergies.lowerBound(energy) - 1, 0)
  # need idx and idx+offset
  let xrayRef = getXrayRefTable()
  let refLow = xrayRef[idx]
  let refHigh = xrayRef[idx+offset]
  let dfRL = df.filter(f{string -> bool: `Dset` == refLow})
  let dfRH = df.filter(f{string -> bool: `Dset` == refHigh})
  let refLowT = dfRL["Hist", float]
  let refHighT = dfRH["Hist", float]

  doAssert dfRL.len == dfRH.len
  let bins = dfRL["Bins", float]
  var res = zeros[float](refLowT.size.int)
  # walk over each bin and compute linear interpolation between
  let Ediff = abs(lineEnergies[idx] - lineEnergies[idx+offset])
  for i in 0 ..< bins.size:
    res[i] = refLowT[i] * (1 - (abs(lineEnergies[idx] - energy)) / Ediff) +
      refHighT[i] * (1 - (abs(lineEnergies[idx+offset] - energy)) / Ediff)
  result = (bins, res)

proc getInterpolatedDf*(df: DataFrame, num = 1000): DataFrame =
  ## returns a DF with `num` interpolated distributions using next neighbors
  ## for linear interpolation
  let energiesLines = getXrayFluorescenceLines()
  result = df.select(["Bins", "Hist", "Energy", "Dset"])
  let energies = linspace(energiesLines[0], energiesLines[^1], num)
  let xrayRef = getXrayRefTable()
  for idx, E in energies:
    let (bins, res) = morph(df, E, offset = 1)
    var dfMorph = toDf({"Bins" : bins, "Hist" : res})
    dfMorph["Energy"] = E
    dfMorph["Dset"] = "Morph"
    result.add dfMorph

proc getInterpolatedWideDf*(df: DataFrame, num = 1000): DataFrame =
  ## returns a DF with `num` interpolated distributions using next neighbors
  ## for linear interpolation
  let lineEnergies = getXrayFluorescenceLines()
  let energies = linspace(lineEnergies[0], lineEnergies[^1], num)
  let xrayRef = getXrayRefTable()
  result = newDataFrame()
  for tup, subDf in groups(group_by(df, "Variable")):
    var dfLoc = newDataFrame()
    var lastBins = zeros[float](0)
    for idx, E in energies:
      let (bins, res) = morph(subDf, E, offset = 1)
      block Sanity:
        ## really the same bins in all Target/Filter combinations? Should be, but check!
        if lastBins.size > 0:
          doAssert bins == lastBins
        lastBins = bins
      let suffix = "_" & $idx
      dfLoc["Hist" & $suffix] = res
    dfLoc["Bins"] = lastBins
    dfLoc["Variable"] = constantColumn(tup[0][1].toStr, dfLoc.len)
    result.add dfLoc

proc morphToA*(df: DataFrame, lineEnergies: seq[float], energy: float, offset = 1): (Tensor[float], Tensor[float]) =
  ## generates a distribution for the appropriate energy `energy` between the
  ## distribution below and above `energy` using linear interpolation
  ## DF needs to have columns:
  ## - "Hist": counts of distributions
  ## - "Bins": bin edges of distributions
  let idx = max(lineEnergies.lowerBound(energy*1000) - 1, 0)

  # need idx and idx+offset
  let refLow = int(lineEnergies[idx])
  let refHigh = int(lineEnergies[idx+offset])
  let refLowT = df[$refLow, float]
  let refHighT = df[$refHigh, float]
  let bins = df["bin", float]
  var res = zeros[float](refLowT.size.int)
  # walk over each bin and compute linear interpolation between
  let Ediff = abs((lineEnergies[idx]/1000) - (lineEnergies[idx+offset]/1000))
  for i in 0 ..< bins.size:
    res[i] = refLowT[i] * (1 - (abs((lineEnergies[idx]/1000) - energy)) / Ediff) +
      refHighT[i] * (1 - (abs((lineEnergies[idx+offset]/1000) - energy)) / Ediff)
  result = (bins, res)

proc morphTpx3*(df: DataFrame, lineEnergies: seq[float], energy: float, offset = 1): (Tensor[float], Tensor[float]) =
  ## generates a distribution for the appropriate energy `energy` between the
  ## distribution below and above `energy` using linear interpolation
  ## DF needs to have columns:
  ## - "Hist": counts of distributions
  ## - "Bins": bin edges of distributions
  let idx = max(lineEnergies.lowerBound(energy*1000) - 1, 0)
  # need idx and idx+offset
  let refLow = int(lineEnergies[idx])
  let refHigh = int(lineEnergies[idx+offset])
  let dfRL = df.filter(f{string -> bool: `Dset` == $refLow})
  let dfRH = df.filter(f{string -> bool: `Dset` == $refHigh})
  let refLowT = dfRL["Hist", float]
  let refHighT = dfRH["Hist", float]
  doAssert dfRL.len == dfRH.len
  let bins = dfRL["Bins", float]
  var res = zeros[float](refLowT.size.int)
  # walk over each bin and compute linear interpolation between
  let Ediff = abs((lineEnergies[idx]/1000) - (lineEnergies[idx+offset]/1000))
  for i in 0 ..< bins.size:
    res[i] = refLowT[i] * (1 - (abs((lineEnergies[idx]/1000) - energy)) / Ediff) +
      refHighT[i] * (1 - (abs((lineEnergies[idx+offset]/1000) - energy)) / Ediff)
  result = (bins, res)

proc getInterpolatedDfToA*(df: DataFrame, lineEnergies: seq[int],  dftype: string, num = 1000): DataFrame =
  ## returns a DF with `num` interpolated distributions using next neighbors
  ## for linear interpolation
  #echo "enter interpolation"
  let lineEnergies=lineEnergies.mapIt(it.float)
  let energies = linspace((lineEnergies[0]/1000), (lineEnergies[^1]/1000), num)
  var dfLoc = newDataFrame()
  var lastBins = zeros[float](0)
  result = newDataFrame()
  for idx, E in energies:
    let (bins, res) = morphToA(df, lineEnergies, E, offset = 1)
    block Sanity:
      # really the same bins in all Target/Filter combinations? Should be, but check!
      if lastBins.size > 0:
        doAssert bins == lastBins
      lastBins = bins
    let suffix = "_" & $idx
    dfLoc["Hist" & $suffix] = res
  dfLoc["Bins"] = lastBins
  dfLoc["Variable"] = dftype
  #echo dfLoc
  result.add dfLoc

proc getInterpolatedDfToAlong*(df: DataFrame, num = 1000): DataFrame =
  ## returns a DF with `num` interpolated distributions using next neighbors
  ## for linear interpolation
  let ecc_path = TpxDir / "resources/Ecc_P_densitys.csv"
  let (eccdf, eccEnergy_list)= readToAProbabilities(ecc_path)

  let energiesLines = eccEnergy_list.mapIt(it.float)
  result = df.select(["Bins", "Hist", "Energy", "Dset"])
  let energies = linspace((energiesLines[0]/1000), (energiesLines[^1]/1000), num)
  #let xrayRef = getXrayRefTable()
  for idx, E in energies:
    let (bins, res) = morphTpx3(df, energiesLines, E, offset = 1)
    var dfMorph = toDf({"Bins" : bins, "Hist" : res})
    dfMorph["Energy"] = E
    dfMorph["Dset"] = "Morph"
    result.add dfMorph

proc readMorphKind(): MorphingKind
proc calcLikelihoodDataset*(h5f: var H5File,
                            groupName: string,
                            ctx: LikelihoodContext): seq[float]

when false:
  proc calcLikelihoodDatasetIfNeeded*(h5f: var H5File,
                                      grp: string,
                                      ctx: LikelihoodContext) =
    ## Recomputes the likelihood dataset, either if it doesn't exist yet or
    ## if it was computed using a different morphing technique.
    # get the group of this dataset
    let cdlGroup = h5f[grp.grp_str]
    let morphKindUsed = if "MorphingKind" in cdlGroup.attrs:
                          some(parseEnum[MorphingKind](cdlGroup.attrs["MorphingKind", string]))
                        else:
                          none[MorphingKind]()
    let morphKind = readMorphKind() # the morph kind that was used, assign it
    if morphKindUsed.isNone or morphKindUsed.get != morphKind:
      let logLData = h5f.calcLikelihoodDataset(grp, ctx)
      let loglDset = h5f.create_dataset(grp / logLDset,
                                        logLData.len,
                                        float64,
                                        overwrite = true)
      logLDset[logLDset.all] = logLData
      # write used morph kind to group
      cdlGroup.attrs["MorphingKind"] = $morphKind
    # else nothing to do

template withCdlData*(cdlFile, dset: string,
                      year: YearKind, energyDset: InGridDsetKind,
                      igDsets: seq[InGridDsetKind],
                      body: untyped): untyped =
  ## Sorry for the nested templates with a whole bunch of injected variables
  ## future reader... :)
  ##
  ## The `data` variable is the heart of it. An array of the `IngridDsetsKinds` handed
  ## as `igDsets` (plus a few mandatory ones!) of all the data that passes the charge cuts in case of
  ## `fitByRun` and all data in case of `fitByRun = false`.
  ## XXX: This is ugly. Both branches should have the same behavior.
  ##
  ## WARNING: Any dataset that is not part of `igDsets` will not be read, but the `data` array
  ## still has an entry for that enum value! Do not access an enum field without checking the
  ## length of the data or making sure you supply the dataset as `igDsets`!
  ##
  ## `withCdlData` should ``not`` be used on its own ``unless`` your are fine with
  ## the fact that it pre-applies the charge cuts to the CDL data around the main
  ## fitted peak of each target/filter dataset (by run)!
  var grp_name {.inject.} = cdlPrefix($year) & dset
  # create global vars for xray and normal cuts table to avoid having
  # to recreate them each time
  let xrayCutsTab {.global.} = getXrayCleaningCuts()
  var cutsTab {.global.}: OrderedTable[string, Cuts]
  case year
  of yr2014:
    cutsTab = getEnergyBinMinMaxVals2014()
  of yr2018:
    cutsTab = getEnergyBinMinMaxVals2018()

  var frameworkKind {.inject.} = fkMarlin
  var fitByRun {.inject.} = false
  echo "Opening file to build LogL from ", cdlFile, " for dset: ", dset
  withH5(cdlFile, "r"):
    if "FrameworkKind" in h5f.attrs:
      frameworkKind = parseEnum[FrameworkKind](h5f.attrs["FrameworkKind", string])
    if "fitByRun" in h5f.attrs:
      fitByRun = parseBool(h5f.attrs["fitByRun", string])
    # open h5 file using template
    ## XXX: if `fitByRun` instead read by each run!
    let chargeStr = igTotalCharge.toDset(frameworkKind)
    var cuts {.inject.} = cutsTab[dset]
    let xrayCuts {.inject.} = xrayCutsTab[dset]
    var data {.inject.}: array[InGridDsetKind, seq[float]]

    ## XXX: make gain an InGridDsetKind and insert it here based on the gas gain slices?

    for s in mitems(data):
      s = newSeqOfCap[float](50_000)
    if fitByRun:
      ## In cas of fitByRun we need to read data from all runs and stack the data into the
      ## variables
      ## for run
      for runGrp in items(h5f, start_path = grp_name, depth = 1):
        let chargeFull = h5f.readAs(runGrp.name / chargeStr, float64)
        # read charge cut values, then read rest and filter to charge
        if "ChargeLow" notin runGrp.attrs or "ChargeHigh" notin runGrp.attrs:
          raise newException(KeyError, "The group " & $runGrp.name & " does not contain the charge " &
            "cut bounds. This is likely because you use an old `calibration-cdl*.h5` file.")
        let dataFull = h5f.readInGridDsetKind(runGrp.name, igDsets)
        cuts.minCharge  = runGrp.attrs["ChargeLow", float]
        cuts.maxCharge = runGrp.attrs["ChargeHigh", float]

        ## XXX: just get passing indices, then do rest in one?
        var passIdx = newSeqOfCap[int](50_000)
        for i, c in chargeFull:
          # now add all data that passes the charge cut. All other cuts will be used later
          ## IMPORTANT: This here implies that we _cannot_ get the CDL data _without_ cuts to
          ## the charge. However, we don't *do* that anywhere anyway, unless the user calls
          ## `withCdlData` manually.
          if c >= cuts.minCharge and c <= cuts.maxCharge:
            passIdx.add i
        for d in igDsets:
          data[d].add dataFull[d][passIdx]
          #for i, c in chargeFull:
          #  # now add all data that passes the charge cut. All other cuts will be used later
          #  ## IMPORTANT: This here implies that we _cannot_ get the CDL data _without_ cuts to
          #  ## the charge. However, we don't *do* that anywhere anyway, unless the user calls
          #  ## `withCdlData` manually.
          #  if c >= cuts.minCharge and c <= cuts.maxCharge:
          #    data[d].add dataFull[d][i]
    else:
      let grp = h5f[grp_name.grp_str]
      if "ChargeLow" notin grp.attrs or "ChargeHigh" notin grp.attrs:
        raise newException(KeyError, "The group " & $grp.name & " does not contain the charge " &
          "cut bounds. This is likely because you use an old `calibration-cdl*.h5` file.")
      cuts.minCharge  = grp.attrs["ChargeLow", float]
      cuts.maxCharge = grp.attrs["ChargeHigh", float]
      # the charge cuts will be applied in the regular `minCharge`/`maxCharge` filter
      # now read all data
      data = h5f.readInGridDsetKind(grp_name, igDsets)
    body

template withLogLFilterCuts*(cdlFile, dset: string,
                             year: YearKind, energyDset: InGridDsetKind,
                             igDsets: seq[InGridDsetKind],
                             body: untyped): untyped =
  ## Note: see the injected variables in the `withCdlData` template!
  for d in LogLCutDsets:
    if d notin igDsets:
      raise newException(ValueError, "The dataset: " & $d & " is not in `igDsets`, but is " &
        "required to perform the LogL filter cuts!")

  withCdlData(cdlFile, dset, year, energyDset, igDsets):
    for i {.inject.} in 0 ..< data[energyDset].len:
      let
        # first apply Xray cuts (see C. Krieger PhD Appendix B & C)
        regionCut  = inRegion(data[igCenterX][i], data[igCenterY][i], crSilver)
        xRmsCut    = data[igRmsTransverse][i].float >= xrayCuts.minRms and data[igRmsTransverse][i] <= xrayCuts.maxRms
        xLengthCut = data[igLength][i].float        <= xrayCuts.maxLength
        xEccCut    = data[igEccentricity][i].float  <= xrayCuts.maxEccentricity
        # then apply reference cuts (charge cut already applied if `fitByRun` in use!)
        chargeCut  = if fitByRun: true
                     else: data[igTotalCharge][i].float > cuts.minCharge and data[igTotalCharge][i] < cuts.maxCharge
        rmsCut     = data[igRmsTransverse][i].float > cuts.minRms    and data[igRmsTransverse][i] < cuts.maxRms
        lengthCut  = data[igLength][i].float        < cuts.maxLength
        pixelCut   = data[igHits][i].float          > cuts.minPix
      # add event to likelihood if all cuts passed
      if allIt([regionCut, xRmsCut, xLengthCut, xEccCut, chargeCut, rmsCut, lengthCut, pixelCut], it):
        body

template withXrayCleaningCuts*(cdlFile, dset: string,
                               year: YearKind, energyDset: InGridDsetKind,
                               igDsets: seq[InGridDsetKind],
                               body: untyped): untyped =
  ## Note: see the injected variables in the `withCdlData` template!
  for d in [igCenterX, igCenterY, igTotalCharge, igRmsTransverse, igLength, igHits]:
    if d notin igDsets:
      raise newException(ValueError, "The dataset: " & $d & " is not in `igDsets`, but is " &
        "required to perform the X-ray reference cuts!")

  withCdlData(cdlFile, dset, year, energyDset, igDsets):
    for i {.inject.} in 0 ..< data[energyDset].len:
      let
        # first apply Xray cuts (see C. Krieger PhD Appendix B & C)
        regionCut  = inRegion(data[igCenterX][i], data[igCenterY][i], crSilver)
        xRmsCut    = data[igRmsTransverse][i].float >= xrayCuts.minRms and data[igRmsTransverse][i] <= xrayCuts.maxRms
        xLengthCut = data[igLength][i].float        <= xrayCuts.maxLength
        xEccCut    = data[igEccentricity][i].float  <= xrayCuts.maxEccentricity
      # add event to likelihood if all cuts passed
      if allIt([regionCut, xRmsCut, xLengthCut, xEccCut], it):
        body

template withXrayRefCuts*(cdlFile, dset: string,
                          year: YearKind, energyDset: InGridDsetKind,
                          igDsets: seq[InGridDsetKind],
                          body: untyped): untyped =
  ## Note: see the injected variables in the `withCdlData` template!
  for d in [igCenterX, igCenterY, igTotalCharge, igRmsTransverse, igLength, igHits]:
    if d notin igDsets:
      raise newException(ValueError, "The dataset: " & $d & " is not in `igDsets`, but is " &
        "required to perform the X-ray reference cuts!")

  withCdlData(cdlFile, dset, year, energyDset, igDsets):
    for i {.inject.} in 0 ..< data[energyDset].len:
      let
        # first apply Xray cuts (see C. Krieger PhD Appendix B & C)
        regionCut  = inRegion(data[igCenterX][i], data[igCenterY][i], crSilver)
        # then apply reference cuts (charge cut already applied if `fitByRun` in use!)
        chargeCut  = if fitByRun: true
                     else: data[igTotalCharge][i].float > cuts.minCharge and data[igTotalCharge][i] < cuts.maxCharge
        rmsCut     = data[igRmsTransverse][i].float > cuts.minRms    and data[igRmsTransverse][i] < cuts.maxRms
        lengthCut  = data[igLength][i].float        < cuts.maxLength
        pixelCut   = data[igHits][i].float          > cuts.minPix
      # add event to likelihood if all cuts passed
      if allIt([regionCut, chargeCut, rmsCut, lengthCut, pixelCut], it):
        body

proc readRawRefData*(
  cdlFile, dset: string,
  year: YearKind,
  energyDset: InGridDsetKind
                   ): tuple[eccs, ldiv, frac: seq[float]] =
  var eccs = newSeq[float]()
  var ldiv = newSeq[float]()
  var frac = newSeq[float]()
  #withXrayRefCuts(cdlFile, dset, year, energyDset):
  withLogLFilterCuts(cdlFile, dset, year, energyDset, LogLCutDsets):
    # read filtered data
    eccs.add data[igEccentricity][i]
    ldiv.add data[igLengthDivRmsTrans][i]
    frac.add data[igFractionInTransverseRms][i]
  result = (eccs: eccs, ldiv: ldiv, frac: frac)

proc readRawLogLFilteredData*(
  cdlFile, dset: string,
  year: YearKind,
  energyDset: InGridDsetKind
                   ): tuple[eccs, ldiv, frac, energy: seq[float]] =
  var eccs = newSeq[float]()
  var ldiv = newSeq[float]()
  var frac = newSeq[float]()
  var E = newSeq[float]()
  withLogLFilterCuts(cdlFile, dset, year, energyDset, LogLCutDsets):
    # read filtered data
    eccs.add data[igEccentricity][i]
    ldiv.add data[igLengthDivRmsTrans][i]
    frac.add data[igFractionInTransverseRms][i]
    E.add data[energyDset][i]
  result = (eccs: eccs, ldiv: ldiv, frac: frac, energy: E)

proc buildRefHistos(year: YearKind, eccs, ldiv, frac: seq[float]): tuple[ecc, ldivRms, fracRms: HistTuple] =
  proc assignRes(x: seq[float], name: string): HistTuple =
    var
      numBins: int
      minVal: float
      maxVal: float
    case year
    of yr2014:
      (numBins, minVal, maxVal) = cdlToXrayBinning2014(name)
    of yr2018:
      (numBins, minVal, maxVal) = cdlToXrayBinning2018(name)
    let (h, b) = histogram(x,
                           numBins,
                           range = (minVal, maxVal),
                           upperRangeBinRight = false)
    # combine the hist bins data to a seq2D
    result = (bins: b[0 .. ^2], hist: h.mapIt(it.float))
  result[0] = assignRes(eccs, igEccentricity.toDset(fkTpa))
  result[1] = assignRes(ldiv, igLengthDivRmsTrans.toDset(fkTpa))
  result[2] = assignRes(frac, igFractionInTransverseRms.toDset(fkTpa))

proc applyMutation(
  fns: seq[FormulaNode], eccs, ldiv, frac: seq[float]
     ): tuple[eccs, ldiv, frac: seq[float]] =
  var df = toDf({"eccentricity" : eccs, "lengthDivRmsTrans" : ldiv, "fractionInTransverseRms" : frac})
  when false:
    ggplot(df.filter(f{`eccentricity` < 2.0}), aes("eccentricity")) +
      geom_histogram(bins = 100) +
      ggsave("/t/histo_ecc_original_" & $dset & ".pdf")

  df = df.mutate(fns)
  let eccs = df["eccentricity", float].toSeq1D
  let ldiv = df["lengthDivRmsTrans", float].toSeq1D
  let frac = df["fractionInTransverseRms", float].toSeq1D

  when false:
    ggplot(df.filter(f{`eccentricity` < 2.0}), aes("eccentricity")) +
      geom_histogram(bins = 100) +
      ggsave("/t/histo_ecc_" & $dset & ".pdf")
  result = (eccs, ldiv, frac)

proc genRefData*(cdlFile, dset: string,
                 year: YearKind, energyDset: InGridDsetKind,
                 fns: seq[FormulaNode] = @[]):
                   tuple[ecc, ldivRms, fracRms: HistTuple] =
  # now compute histograms
  var (eccs, ldiv, frac) = readRawRefData(cdlFile, dset, year, energyDset)
  if fns.len > 0:
    (eccs, ldiv, frac) = fns.applyMutation(eccs, ldiv, frac)
  result = year.buildRefHistos(eccs, ldiv, frac)

proc calcFns(ctx: LikelihoodContext, energy: float, eccs, ldiv, frac: seq[float]): seq[FormulaNode] =
  ## Computes the formulas used to modifiy the CDL data based on the possible
  ## `CDLStretcher`
  if ctx.stretch.isSome:
    let st = ctx.stretch.get
    proc toFn(dset: string, cdl: seq[float]): FormulaNode =
      # compute interpolated data min / max
      let
        lp = st.lines[dset] # get line param tuple
        cdlMin = cdl.percentile(1) #min
        cdlMax = cdl.percentile(99) # max
      var
        dataMin = (lp.mins.m * energy + lp.mins.b) * cdlMin
        dataMax = (lp.maxs.m * energy + lp.maxs.b) * cdlMax
      let dS = dset # local copy due to `f{}` capturing it

      when false:
        var f = open("/t/data.csv", fmAppend)
        f.write(&"{cdlMin},{cdlMax},{dataMin},{dataMax},{dset},{energy}\n")
        f.close()
        echo "Min ", cdlMin, " to ", cdlMax, " for ", dataMin, " to ", dataMax, " for dset ", dset, " and energy ", energy

      if dataMin > cdlMin:
        dataMin = cdlMin
      if dataMax < cdlMax:
        dataMax = cdlMax

      result = f{float: dS ~ (idx(dS) - cdlMin) / (cdlMax - cdlMin) * (dataMax - dataMin) + dataMin}
    result.add toFn("eccentricity", eccs)
    result.add toFn("lengthDivRmsTrans", ldiv)
    result.add toFn("fractionInTransverseRms", frac)

proc genRefData*(dset: string, ctx: LikelihoodContext): tuple[ecc, ldivRms, fracRms: HistTuple] =
  const xrayEnergies = getXrayFluorescenceLines()
  const invXrayTab = getInverseXrayRefTable()
  # now compute histograms
  var (eccsR, ldivR, fracR) = readRawRefData(ctx.cdlFile, dset, ctx.year, ctx.energyDset)
  # compute the correct functions for the reference data if needed
  let dsetEnergy = xrayEnergies[invXrayTab[dset]]
  let fns = calcFns(ctx, dsetEnergy, eccsR, ldivR, fracR)
  (eccsR, ldivR, fracR) = fns.applyMutation(eccsR, ldivR, fracR)
  result = ctx.year.buildRefHistos(eccsR, ldivR, fracR)

proc readRefDsetsDF*(ctx: LikelihoodContext): DataFrame =
  ## reads the reference datasets from the `refFile` and returns them.

  ## TODO: implement proper handling of 55Fe distribution matching by replacing
  ## `genRefData` with the version taking `ctx`
  const xray_ref = getXrayRefTable()
  var
    ecc_ref = initTable[string, HistTuple]()
    lengthDivRmsTrans_ref = initTable[string, HistTuple]()
    fracRmsTrans_ref = initTable[string, HistTuple]()
  result = newDataFrame()
  let energies = getXrayFluorescenceLines()
  for idx in 0 ..< xray_ref.len:
    ## read the data from the CDL file and generate the reference data using cuts
    let dsetName = xray_ref[idx]
    let E = energies[idx]
    let (ecc, ldiv, frac) = dsetName.genRefData(ctx)
    template addDf(h: HistTuple, key: string, E: float, dkKind: InGridDsetKind) =
      var dfDset = toDf({ "Bins" : h.bins, "Hist" : h.hist,
                          "Dset" : key, "Energy" : E,
                          "Variable" : dkKind.toDset(fkTpa) })
      dfDset["Dset"] = dset_name
      result.add dfDset
    addDf(ecc, dsetName, E, igEccentricity)
    addDf(ldiv, dsetName, E, igLengthDivRmsTrans)
    addDf(frac, dsetName, E, igFractionInTransverseRms)

proc calcLogL*(e, l, f: float, eccs, ldiv, frac: HistTuple): float =
  ## Computes the log likelihood value for the given eccentricity, ldivT and fT
  ## based on the given `HistTuple`
  ##
  ## XXX: Can we do better than the current approach based on the binned histograms?
  ## Maybe better use a spline interpolation of a normalized KDE?
  result += logLikelihood(eccs[1], e, eccs[0])
  result += logLikelihood(ldiv[1], l, ldiv[0])
  result += logLikelihood(frac[1], f, frac[0])
  result *= -1.0

proc calcLogLwithToA*(e, l, f, t: float, eccs, ldiv, frac, toa: HistTuple): float =
  ## Computes the log likelihood value for the given eccentricity, ldivT and fT
  ## based on the given `HistTuple`
  ##
  ## XXX: Can we do better than the current approach based on the binned histograms?
  ## Maybe better use a spline interpolation of a normalized KDE?
  result += logLikelihood(eccs[1], e, eccs[0])
  result += logLikelihood(ldiv[1], l, ldiv[0])
  result += logLikelihood(frac[1], f, frac[0])
  result += logLikelihood(toa[1], t, toa[0])
  result *= -1.0

proc buildLogLHist*(dset: string, ctx: LikelihoodContext): tuple[logL, energy: seq[float]] =
  ## given a file `h5file` containing a CDL calibration dataset
  ## `dset` apply the cuts on all events and build the logL distribution
  ## for the energy range.
  ## `dset` needs to be of the elements contained in the returned map
  ## of tos_helpers.`getXrayRefTable`
  ## Default `region` is the gold region
  ## Returns a tuple of the actual `logL` values and the corresponding `energy`.
  # decent size that is O(correct)
  result[0] = newSeqOfCap[float](100_000)
  result[1] = newSeqOfCap[float](100_000)
  when false:
    withLogLFilterCuts(ctx.cdlFile, dset, ctx.year, ctx.energyDset):
      result[0].add data[igLikelihood][i]
      result[1].add data[ctx.energyDset][i]
    # now create plots of all ref likelihood distributions
    echo "max is inf ? ", min(result[0]), " ", max(result[0])
  else:
    let (eccs, ldiv, frac) = genRefData(dset, ctx) # , fns)
    let (eccsR, ldivR, fracR, energy) = readRawLogLFilteredData(ctx.cdlFile, dset, ctx.year, ctx.energyDset)
    for i in 0 ..< eccsR.len:
      result[0].add calcLogL(eccsR[i], ldivR[i], fracR[i], eccs, ldiv, frac)
      result[1].add energy[i]

proc computeLogLDistributions*(ctx: LikelihoodContext): DataFrame =
  ## Computes the LogL distributions from thc CDL data file (`cdlFile`) by applying
  ## both sets of cuts (`getXraySpetrcumCuts` and `getEnergyBinMinMaxVals201*`) to the
  ## data in `buildLogLHist` and binning it according to the number and bin width
  ## that we use for the logL distributions.
  const
    xrayRef = getXrayRefTable()
    # logL binning range
    nbins = 300 # TODO: Increase number of bins for logL cut value closer to target?
    # range of histogram in logL
    logLrange = (0.0, 30.0)

  # get the correct binning for the histograms
  let bins = linspace(logLrange[0], logLrange[1], nbins + 1, endpoint = true)
  let energies = getXrayFluorescenceLines()
  for idx, dset in xrayRef:
    # compute the histogram of the CDL data
    let hist = buildLogLHist(dset, ctx)[0]
    var df = toDf( {"Bins" : bins[0 .. ^2], "Hist" : histogram(hist, nbins, logLrange)[0] })
    df["Dset"] = dset
    df["Energy"] = energies[idx]
    result.add df

proc readRawSimData*(energy: string): tuple[eccs, ldiv, frac, toal, energy: seq[float]] =
  ## maybe adding cuts could improve this should try?
  const path = "/reconstruction/run_0/chip_0"
  const ecc = "eccentricity"
  const ldivs = "lengthDivRmsTrans"
  const ftrans = "fractionInTransverseRms"
  const toa = "toaLength"
  const en = "energyFromCharge"
  var fname= TpxDir / "resources/sim_cdl_refs/" & energy & "_3cm_Ar_Isobutane_977_23_787.h5"
  var h5f = H5open(fname, "r")

  let eccs = h5f[(path / ecc), float]
  let ldiv = h5f[(path / ldivs), float]
  let frac = h5f[(path / ftrans), float]
  let toal = h5f[(path / toa), float]
  let E = h5f[(path / en), float]

  result = (eccs: eccs, ldiv: ldiv, frac: frac, toal: toal, energy: E)

proc buildLogLHistusingsim*(dset: string, ctx: LikelihoodContext): tuple[logL, energy: seq[float]] =
  ## given a file `h5file` containing a CDL calibration dataset
  ## `dset` apply the cuts on all events and build the logL distribution
  ## for the energy range.
  ## `dset` needs to be of the elements contained in the returned map
  ## of tos_helpers.`getXrayRefTable`
  ## Default `region` is the gold region
  ## Returns a tuple of the actual `logL` values and the corresponding `energy`.

  #get dataframes
  let ecc_path = TpxDir / "resources/Ecc_P_densitys.csv"
  let ldiv_path = TpxDir / "resources/ldiv_P_densitys.csv"
  let ftrans_path = TpxDir / "resources/ftrans_P_densitys.csv"
  let toa_path = TpxDir / "resources/ToA_P_densitys.csv"
  let (eccdf, eccEnergy_list)= readToAProbabilities(ecc_path)
  let (ldivdf, ldivEnergy_list)= readToAProbabilities(ldiv_path)
  let (fracdf, ftransEnergy_list)= readToAProbabilities(ftrans_path)
  let (toadf, toaEnergy_list)= readToAProbabilities(toa_path)
  let eccs:HistTuple = (toseq(eccdf["bin", float64]),toseq(eccdf[dset, float64]))
  let ldiv:HistTuple = (toseq(ldivdf["bin", float64]),toseq(ldivdf[dset, float64]))
  let frac:HistTuple = (toseq(fracdf["bin", float64]),toseq(fracdf[dset, float64]))
  let toa:HistTuple = (toseq(toadf["bin", float64]),toseq(toadf[dset, float64]))

  # decent size that is O(correct)
  result[0] = newSeqOfCap[float](100_000)
  result[1] = newSeqOfCap[float](100_000)

  let (eccsR, ldivR, fracR, toaR, energy) = readRawSimData(dset)
  if ctx.vetoCfg.useToAlnLCut:
    for i in 0 ..< eccsR.len:
      ##Add 2 to the simulated ref data since its missing Timewalk, 
      ##if at some point a Timewalk calibration is used this might be not necessary any more
      result[0].add calcLogLwithToA(eccsR[i], ldivR[i], fracR[i],(toaR[i]+2), eccs, ldiv, frac, toa)
      result[1].add energy[i]
  else:
    for i in 0 ..< eccsR.len:
      result[0].add calcLogL(eccsR[i], ldivR[i], fracR[i], eccs, ldiv, frac)
      result[1].add energy[i]
  
proc computeLogLDistributionsusingsim*(ctx: LikelihoodContext): DataFrame =
  ## Computes the LogL distributions from thc CDL data file (`cdlFile`) by applying
  ## both sets of cuts (`getXraySpetrcumCuts` and `getEnergyBinMinMaxVals201*`) to the
  ## data in `buildLogLHist` and binning it according to the number and bin width
  ## that we use for the logL distributions.
  const
    # logL binning range
    nbins = 300 # TODO: Increase number of bins for logL cut value closer to target?
    # range of histogram in logL
    logLrange = (0.0, 30.0)

  # get the correct binning for the histograms
  let bins = linspace(logLrange[0], logLrange[1], nbins + 1, endpoint = true)

  #get data
  let ecc_path = TpxDir / "resources/Ecc_P_densitys.csv"
  let (eccdf, eccEnergy_list)= readToAProbabilities(ecc_path)

  let energies = eccEnergy_list

  # compute the histogram of the CDL data
  for energy in energies:
    let hist = buildLogLHistusingsim($energy, ctx)[0]
    var df = toDf( {"Bins" : bins[0 .. ^2], "Hist" : histogram(hist, nbins, logLrange)[0] })
    df["Dset"] = energy
    df["Energy"] = (energy/1000)
    result.add df

proc getLogLData*(ctx: LikelihoodContext): DataFrame =
  ## Returns the raw logL data for all target/filter datasets
  const xrayRef = getXrayRefTable()
  let energies = getXrayFluorescenceLines()
  for idx, dset in xrayRef:
    let (logL, energy) = buildLogLHist(dset, ctx)
    let df = toDf({"logL" : logL, "Energies" : energy, "Energy" : energies[idx], "Dset" : dset})
    result.add df

proc readLogLVariableData*(h5f: var H5File,
                           groupName: string,
                           energyDset: InGridDsetKind,
                           ctx: LikelihoodContext
                          ):
                             (seq[float], seq[float], seq[float], seq[float], seq[float]) =
  # get the datasets needed for LogL
  doAssert energyDset in {igEnergyFromPixel, igEnergyFromCharge}
  let
    ecc = h5f[(groupName / "eccentricity"), float64]
    lengthDivRmsTrans = h5f[(groupName / "lengthDivRmsTrans"), float64]
    fracRmsTrans = h5f[(groupName / "fractionInTransverseRms"), float64]
    # energy to choose correct bin
    energies = h5f[(groupName / energyDset.toDset), float64]
  var toa: seq[float64]
  if ctx.vetoCfg.useToAlnLCut:
    toa = h5f[(groupName / "toaLength"), float64]
  else:
    toa = toSeq(0 ..< ecc.len).mapIt(0.0)
  result = (ecc, lengthDivRmsTrans, fracRmsTrans, toa, energies)

proc readRefDsets*(ctx: LikelihoodContext): tuple[ecc, ldivRms, fracRms: Table[string, HistTuple]] =
  ## Reads the data from the CDL file and generates the reference distributions based
  ## on the cuts defined in `cdl_cuts.nim`
  # create a table, which stores the reference datasets from the ref file
  const xrayRef = getXrayRefTable()
  var
    ecc_ref = initTable[string, HistTuple]()
    lengthDivRmsTrans_ref = initTable[string, HistTuple]()
    fracRmsTrans_ref = initTable[string, HistTuple]()
    df = newDataFrame()

  ## XXX: is it really sensible to compute logL values first? I mean having the option to *only*
  ## compute them seems fine I guess, but in the det effective ε tool we have seen how easy it
  ## is to compute & cut at the same time.

  for idx, dsetName in xrayRef:
    ## read the data from the CDL file and generate the reference data using cuts
    let (eccs, ldiv, frac) = dsetName.genRefData(ctx)
    var dfDset = toDf({ "Eccentricity" : eccs[0],
                        "Ecc #" : eccs[1],
                        "L / RMS_trans" : ldiv[0],
                        "L / RMS_trans #" : ldiv[1],
                        "fracRmsTrans" : frac[0],
                        "fracRmsTrans #" : frac[1] })
    dfDset["Dset"] = dset_name
    df.add dfDset

    ecc_ref[dsetName] = eccs
    lengthDivRmsTrans_ref[dsetName] = ldiv
    fracRmsTrans_ref[dsetName] = frac

  block RefPlots:
    var labelOrder = initTable[Value, int]()
    for idx, el in xrayRef:
      labelOrder[%~ el] = idx
    ggplot(df.filter(f{`Eccentricity` < 5.0}), aes("Eccentricity", "Ecc #", fill = "Dset")) +
      ggridges("Dset", overlap = 1.75, labelOrder = labelOrder) +
      geom_histogram(
        stat = "identity", position = "identity", hdKind = hdOutline,
        color = "black", lineWidth = 1.0
        #alpha = 0.5,
      ) +
      ggtitle(&"Eccentricity of reference file, year: {ctx.year}") +
      ggsave(&"out/eccentricity_ridgeline_{ctx.cdlFile.extractFilename}_{ctx.year}.pdf",
              width = 800, height = 480)
    ggplot(df, aes("L / RMS_trans", "L / RMS_trans #", fill = "Dset")) +
      ggridges("Dset", overlap = 1.75, labelOrder = labelOrder) +
      geom_histogram(
        stat = "identity", position = "identity", hdKind = hdOutline,
        color = "black", lineWidth = 1.0
        #alpha = 0.5,
      ) +
      ggtitle(&"L / RMS_trans of reference file, year: {ctx.year}") +
      ggsave(&"out/lengthDivRmsTrans_ridgeline_{ctx.cdlFile.extractFilename}_{ctx.year}.pdf",
              width = 800, height = 480)
    ggplot(data = df.filter(f{Value: isNull(df["fracRmsTrans"][idx]) == (%~ false)}),
           aes("fracRmsTrans", "fracRmsTrans #", fill = "Dset")) +
      ggridges("Dset", overlap = 1.75, labelOrder = labelOrder) +
      geom_histogram(
        stat = "identity", position = "identity", hdKind = hdOutline,
        color = "black", lineWidth = 1.0
        #alpha = 0.5,
      ) +
      ggtitle(&"fracRmsTrans of reference file, year: {ctx.year}") +
      ggsave(&"out/fracRmsTrans_ridgeline_{ctx.cdlFile.extractFilename}_{ctx.year}.pdf",
              width = 800, height = 480)
  result = (ecc: ecc_ref, ldivRms: lengthDivRmsTrans_ref, fracRms: fracRmsTrans_ref)

func calcLikelihoodForEvent*(energy, eccentricity, lengthDivRmsTrans, fracRmsTrans: float,
                             refSetTuple: tuple[ecc,
                                                ldivRms,
                                                fracRms: Table[string, HistTuple]]): float =
  ## XXX: Manual tuple unpacking due to: https://github.com/nim-lang/Nim/issues/19364
  let ecc_ref = refSetTuple[0]
  let lengthDivRmsTrans_ref = refSetTuple[1]
  let fracRmsTrans_ref = refSetTuple[2]
  # try simple logL calc
  let refset = toRefDset(energy) # / 1000.0) division needed for E from Pix,
                                 # since that is in eV inst of keV
  result = calcLogL(
    eccentricity, lengthDivRmsTrans, fracRmsTrans,
    ecc_ref[refset], lengthDivRmsTrans_ref[refset], fracRmsTrans_ref[refset]
  )

proc calcMorphedLikelihoodForEvent*(eccentricity, lengthDivRmsTrans, fracRmsTrans: float,
                                    refDf: DataFrame, idx: int, ctx: LikelihoodContext): float =
  # try simple logL calc
  ## XXX: replace usages of `{.global.}` here!
  var
    eccDf {.global.}, ldivDf {.global.}, fracDf {.global.}: DataFrame
  once:
    eccDf = refDf.filter(f{string -> bool: `Variable` == "eccentricity"})
    ldivDf = refDf.filter(f{string -> bool: `Variable` == "lengthDivRmsTrans"})
    fracDf = refDf.filter(f{string -> bool: `Variable` == "fractionInTransverseRms"})
  proc addLog(arg: InGridDsetKind, val: float, res: var float, df: DataFrame) =
    let bins = df["Bins", float].toRawSeq
    let hist = df["Hist_" & $idx, float].toRawSeq
    res += logLikelihood(hist,
                         val,
                         bins)

  addLog(igEccentricity, eccentricity, result, eccDf)
  addLog(igLengthDivRmsTrans, lengthDivRmsTrans, result, ldivDf)
  addLog(igFractionInTransverseRms, fracRmsTrans, result, fracDf)
  result *= -1.0

proc calcMorphedLikelihoodForEventToA*(eccentricity, lengthDivRmsTrans, fracRmsTrans, toa: float,
                                    refDf: DataFrame, idx: int, ctx: LikelihoodContext): float =
  # try simple logL calc
  ## XXX: replace usages of `{.global.}` here!

  var
    eccDf {.global.}, ldivDf {.global.}, fracDf {.global.}, ToADf {.global.}: DataFrame
  once:
    eccDf = refDf.filter(f{string -> bool: `Variable` == "eccentricity"})
    ldivDf = refDf.filter(f{string -> bool: `Variable` == "lengthDivRmsTrans"})
    fracDf = refDf.filter(f{string -> bool: `Variable` == "fractionInTransverseRms"})
    ToADf = refDf.filter(f{string -> bool: `Variable` == "ToAlength"})
  proc addLog(arg: InGridDsetKind, val: float, res: var float, df: DataFrame) =
    let bins = df["Bins", float].toRawSeq
    let hist = df["Hist_" & $idx, float].toRawSeq
    res += logLikelihood(hist,
                         val,
                         bins)
  addLog(igEccentricity, eccentricity, result, eccDf)
  addLog(igLengthDivRmsTrans, lengthDivRmsTrans, result, ldivDf)
  addLog(igFractionInTransverseRms, fracRmsTrans, result, fracDf)
  addLog(igToaLength, toa, result, ToADf)
  result *= -1.0

proc calcLikelihoodForEvent*(ctx: LikelihoodContext,
                             energy, ecc, lengthDivRmsTrans, fracRmsTrans: float): float =
  case ctx.morph
  of mkNone:
    result = calcLikelihoodForEvent(energy, ecc, lengthDivRmsTrans, fracRmsTrans,
                                    ctx.refSetTuple)
  of mkLinear:
    let idx = min(ctx.refDfEnergy.lowerBound(energy), ctx.numMorphedEnergies - 1)
    result = calcMorphedLikelihoodForEvent(ecc, lengthDivRmsTrans, fracRmsTrans,
                                           ctx.refDf, idx, ctx)
  #echo "result ? ", result, " based on ", ecc, " ", lengthDivRmsTrans, " ", fracRmsTrans

proc calcLikelihoodForEventToA*(ctx: LikelihoodContext,
                             energy, ecc, lengthDivRmsTrans, fracRmsTrans, toa: float): float =
  case ctx.morph
  of mkNone:
    discard
  of mkLinear:
    let idx = min(ctx.refDfEnergy.lowerBound(energy), ctx.numMorphedEnergies - 1)
    result = calcMorphedLikelihoodForEventToA(ecc, lengthDivRmsTrans, fracRmsTrans, toa,
                                           ctx.refDf, idx, ctx)

proc calcLikelihood*(ctx: LikelihoodContext, df: DataFrame): seq[float] =
  ## Convenience helper that calculates the likelihood values for all events stored in `df`.
  ## `df` must contain the energyFromCharge, eccentricity, lengthDivRmsTrans and fractionInTransverseRms
  ## columns (under that name!).
  let E = igEnergyFromCharge.toDset()
  let ε = igEccentricity.toDset()
  let l = igLengthDivRmsTrans.toDset()
  let f = igFractionInTransverseRms.toDset()
  doAssert E in df
  doAssert ε in df
  doAssert l in df
  doAssert f in df

  result = df.mutate(f{float: "logL" ~ ctx.calcLikelihoodForEvent(idx(E), idx(ε), idx(l), idx(f))})["logL", float].toSeq1D

proc calcLikelihoodDataset*(h5f: var H5File, groupName: string, ctx: LikelihoodContext): seq[float] =
  let (ecc,
       lengthDivRmsTrans,
       fracRmsTrans, toa,
       energies) = h5f.readLogLVariableData(groupName, ctx.energyDset, ctx)


  # create seq to store data logL data for this chip
  ## create a crazy man's plot
  proc crazyPlot(name: string) =
    let dfFiltered = ctx.refDf.filter(f{string: `Variable` == name})
    let histKeysName = dfFiltered.getKeys.filterIt("Hist" in it)
    var df = dfFiltered
      .gather(histKeysName, key = "Idx", value = "Hist")
    var labelOrder = initTable[Value, int]()
    for i in 0 ..< ctx.numMorphedEnergies:
      labelOrder[%~ ("Hist_" & $i)] = i
    #echo df
    ggplot(df, aes("Bins", "Hist", fill = "Idx")) +
      facet_wrap("Idx") +
      geom_histogram(stat = "identity", position = "identity", hdKind = hdOutline) +
      ggsave("/tmp/crazy_" & $name & "_facet_100.pdf", width = 6000, height = 4000)
    #ggplot(df, aes("Bins", "Hist", fill = "Idx")) +
    #  ggridges("Idx", overlap = 1.5, labelOrder = labelOrder) +
    #  geom_histogram(stat = "identity", position = "identity",
    #                 color = some(color(0.0, 0.0, 0.0)),
    #                 hdKind = hdOutline) +
    #  ggsave("/tmp/crazy_" & $name & "_ridge_100.pdf", width = 1900, height = 8000)
  #crazyPlot("igEccentricity")
  #crazyPlot("igLengthDivRmsTrans")
  #crazyPlot("igFractionInTransverseRms")
  result = newSeq[float64](ecc.len)
  echo "[INFO]: Performing likelihood compute using morph kind: ", ctx.morph, " for chip: ", groupName
  for i in 0 .. ecc.high:
    if ctx.vetoCfg.useToAlnLCut:
      result[i] = ctx.calcLikelihoodForEventToA(energies[i], ecc[i], lengthDivRmsTrans[i], fracRmsTrans[i], toa[i])
    else:
      result[i] = ctx.calcLikelihoodForEvent(energies[i], ecc[i], lengthDivRmsTrans[i], fracRmsTrans[i])

import parsetoml
template withConfig(body: untyped): untyped =
  const sourceDir = currentSourcePath().parentDir
  let config {.inject.} = parseToml.parseFile(sourceDir / "../config.toml")
  body

proc readMorphKind(): MorphingKind =
  ## reads the `morphingKind` field from the TOML file
  withConfig:
    result = parseEnum[MorphingKind](config["Likelihood"]["morphingKind"].getStr)

proc readNeuralNetCutKind(): NeuralNetCutKind =
  ## reads the `neuralNetCutMethod` field from the TOML file
  withConfig:
    result = parseEnum[NeuralNetCutKind](config["Likelihood"]["neuralNetCutKind"].getStr)

proc readSignalEff(): float =
  ## reads the `signalEfficiency` field from the TOML file
  withConfig:
    result = config["Likelihood"]["signalEfficiency"].getFloat

proc readConfigField(field: string): string =
  withConfig:
    result = config["Likelihood"][field].getStr

proc determineCutValue*[T: seq | Tensor](hist: T, eff: float): int =
  ## given a histogram `hist`, determine the correct bin to cut at to achieve
  ## a software efficiency of `eff`
  var
    cur_eff = 0.0
    last_eff = 0.0
  let hist_sum = hist.sum.float
  var curSum = 0.0
  while cur_eff < eff:
    last_eff = cur_eff
    curSum += hist[result].float
    cur_eff = curSum / hist_sum
    inc result
  echo "Efficiency is at ", cur_eff, " and last Eff was ", last_eff
  echo "Total elements are: ", histSum, " relative are ", hist[0 .. result].sum

proc determineCutValueData*[T: seq | Tensor](data: T, eff: float): int =
  ## Determine the cut value simply based on the quantile for the given `eff`
  result = (data.len.float * eff).round.int

proc calcCutValueTab*(ctx: LikelihoodContext): CutValueInterpolator =
  ## returns a table mapping the different CDL datasets to the correct cut values
  ## based on the chip center region
  # read signal efficiency (default 80%) from TOML file

  ## TODO: finish implementation of 55Fe for morphed data!
  let morphKind = ctx.morph
  let xray_ref = getXrayRefTable()
  case morphKind
  of mkNone:
    let logLData = getLogLData(ctx)
    # get the cut value for a software efficiency of 80%
    result = initCutValueInterpolator(morphKind)
    for tup, subDf in groups(logLData.group_by("Dset")):
      echo "Starting FOR LOOP--------------------------------------------------"
      let dset = tup[0][1].toStr
      let logL = subDf["logL", float]
      let cutValIdx = determineCutValueData(logL, ctx.vetoCfg.signalEfficiency)
      # incl the correct values for the logL cut values
      result[dset] = logL.sorted(SortOrder.Ascending)[cutValIdx]
      echo "Cut value of ", dset, " is ", result[dset]
  of mkLinear:
    ## given logLHist compute interpolated DF from it
    ## first need a logLHist in DF format
    #
    ## XXX: Can we update this to also do better than the current binned approach?
    result = initCutValueInterpolator(morphKind)
    var logHists: DataFrame
    var dfInterp: DataFrame
    if ctx.vetoCfg.usesimref:
      logHists = computeLogLDistributionsusingsim(ctx)
      dfInterp = logHists.getInterpolatedDfToAlong()
        .filter(f{string: `Dset` == "Morph"})
    else:
      logHists = computeLogLDistributions(ctx)
      dfInterp = logHists.getInterpolatedDf()
        .filter(f{string: `Dset` == "Morph"})
    var
      energies = newSeq[float]()
      cutVals = newSeq[float]()
    for tup, subDf in groups(dfInterp.group_by("Energy")):
      let energy = tup[0][1].toFloat
      when false:
        var efficiency = 0.8
        echo energy
        if energy < 1.0:
          efficiency = 0.6
      let cutVal = determineCutValue(subDf["Hist", float], ctx.vetoCfg.signalEfficiency)
      # incl the correct values for the logL cut values
      energies.add energy
      cutVals.add subDf["Bins", float][cutVal]
    let sorted = zip(energies, cutVals).sortedByIt(it[0])
    result.lnLCutEnergies = sorted.mapIt(it[0])
    result.lnLCutValues = sorted.mapIt(it[1]).toTensor
  else:
    result = getChristophCutVals()
  when false: # not defined(release) or defined(DEBUG):
    #echo logHists
    echo "logHists ", logHists
    #echo "Bins are ", bins
    echo "Cut values are ", result
    #echo mapIt(logHists, it.sum)
    #echo "Corresponding to logL values of ", result

proc calcLowHigh(df: DataFrame, dset: string, percentile, scaleCutoff: float): (float, float) =
  ## Compute the low and high cut value based on the pre defined percentile
  ## and scale factor.
  let data = df[dset, float]
  let samples = linspace(data.min, data.max, 1000).toTensor
  # get the maximum by computing the KDE & getting the sample of the maximum index
  let kdeData = kde(data, samples = samples, normalize = true)
  let kdeArgMax = kdeData.argmax(0)
  let maxVal = samples[kdeArgMax[0]]
  let cutoff = maxVal * scaleCutoff
  # now compute the desired percentile of the rise time data after cutting
  # to `cutoff`
  let dfF = df.filter(f{float -> bool: idx(dset) < cutoff})
  let dataCut = dfF[dset, float]
  let perc = (percentile * 100.0).round.int # convert to integer (0, 100)
  result = (percentile(dataCut, 100 - perc), percentile(dataCut, perc))
  when true:
    echo "For dset: ", dset, " cutting at ", result[0], " to ", result[1]
    # plot the dataset including the percentile cutoffs & hard cutoff
    let dfProps = toDf({ "y" : @[0.0, kdeData.max],
                         "Low" : @[result[0], result[0]],
                         "High" : @[result[1], result[1]],
                         "HardCut" : @[cutoff, cutoff] })
    let setting = df["Settings", string][0]
    ggplot(df.filter(f{idx(dset) < 2.0 * cutoff}), aes(dset)) +
      geom_density(normalize = true) +
      geom_line(data = dfProps, aes = aes("Low", "y"),
                       color = "blue") +
      geom_line(data = dfProps, aes = aes("High", "y"),
                       color = "red") +
      geom_line(data = dfProps, aes = aes("HardCut", "y"),
                       color = "black") +
      ggtitle(&"fadc {dset} {setting} percentile {percentile} scaleCut {scaleCutoff}") +
      ggsave(&"/tmp/fadc_{dset}_{setting}_perc_{percentile}_scaleCut_{scaleCutoff}.pdf")

proc determineFadcVetoCutoff(fname: string, vetoPercentile, fadcScaleCutoff: float): array[FadcSetting, FadcCuts] =
  # 1. get rise time & fall time data, cleaned with some cuts
  # 2. determine correct subsets correlating to FADC settings
  # 3. determine peak of distribution. find hard cut for data
  #   based on peak * scale
  # 4. determine cut based on percentile of data
  # 5. store result in a mapping of valid runs -> `FadcCuts`
  if fname.len == 0: return # nothing to do
  let h5f = H5open(fname, "r")
  var df = readFilteredFadc(h5f)
  # now for each FADC setting:
  for (tup, subDf) in groups(df.group_by("Settings")):
    let (riseLow, riseHigh) = calcLowHigh(subDf, "riseTime", vetoPercentile, fadcScaleCutoff)
    let (fallLow, fallHigh) = calcLowHigh(subDf, "fallTime", vetoPercentile, fadcScaleCutoff)
    # parse string of group label to setting
    let fSetting = parseEnum[FadcSetting](tup[0][1].toStr)
    result[fSetting] = (active: true, # any setting we don't see from this file will be inactive.
                        riseLow: riseLow, riseHigh: riseHigh,
                        fallLow: fallLow, fallHigh: fallHigh,
                        skewness: -0.4)
  discard h5f.close()

proc initLikelihoodContext*(
  cdlFile: string, year: YearKind, region: ChipRegion,
  energyDset: InGridDsetKind,
  timepix: TimepixVersion,
  morphKind: MorphingKind,
  stretch = none[CdlStretch](),
  numMorphedEnergies = 1000,
  centerChip = 3,
  numChips = 7,
  useTeX: bool = false,
  # NN cut
  useNeuralNetworkCut: bool = false,
  nnModelPath = "",
  nnSignalEff: float = 0.0,
  nnCutKind: NeuralNetCutKind = nkNone,
  # lnL cut & settings
  useLnLCut: bool = false,
  signalEfficiency: float = 0.0,
  #usesim
  usesimref: bool = false,
  #ToACut
  useToACut: bool = false,
  ToAcutValue: int = 12,
  #ToAlnLCut
  useToAlnLCut: bool = false,
  # Septem veto related
  septemVeto: bool = false,
  clusterAlgo: ClusteringAlgorithm = caDBSCAN,
  searchRadius: int = 50,
  dbscanEpsilon: float = 65.0,
  useRealLayout: bool = true,
  # line veto related
  lineVetoKind: LineVetoKind = lvNone,
  eccLineVetoCut: float = 0.0,
  # FADC veto related
  fadcVeto: bool = false,
  calibFile: string = "",
  vetoPercentile: float = 0.0,
  fadcScaleCutoff: float = 0.0,
  septemLineVetoEfficiencyFile = "/tmp/septem_veto_before_after.txt",
  rngSeed: int = 299_792_458,
  flags: set[LogLFlagKind] = {},
  readLogLData = false,
  plotPath = "",
                           ): LikelihoodContext =
  ## The configuration elements are generally picked according to the following priority:
  ## 1. command line argument
  ## 2. environment variable (if supported for the option)
  ## 3. config file
  ## 4. default
  #
  ## line veto related
  ## The environment variable `PLOT_SEPTEM_E_CUTOFF` can be used to adjust the energy
  ## cutoff for which events to plot when running with `--plotSeptem`. In addition the
  ## variable `USE_TEX` can be adjusted to generate TikZ TeX plots.


  let useTeX = if useTeX: useTeX
               else: getEnv("USE_TEX", "false").parseBool
  ## XXX: add these to config.toml and as a cmdline argument in addition!
  let PlotCutEnergy = getEnv("PLOT_SEPTEM_E_CUTOFF", "5.0").parseFloat
  let useRealLayout = if useRealLayout: useRealLayout
                      else: parseBool(getEnv("USE_REAL_LAYOUT", "true"))
  ## XXX: Add config.toml field for these!
  let lvKindEnv = parseEnum[LineVetoKind](getEnv("LINE_VETO_KIND", "lvNone"))
  let lvKind = if lineVetoKind != lvNone: lineVetoKind
               elif lvKindEnv != lvNone: lvKindEnv
               elif septemVeto: lvRegularNoHLC # in this case don't need HLC
               else: lvRegular # in case line veto only (or ignored anyway)
  let eccLvCutEnv = parseFloat(getEnv("ECC_LINE_VETO_CUT", "0.0"))
  let eccLineVetoCut = if eccLineVetoCut > 0.0: eccLineVetoCut
                       elif eccLvCutEnv > 0.0: eccLvCutEnv
                       else: 1.0 # default

  let nnCutKind = if nnCutKind != nkNone: nnCutKind
                  else: readNeuralNetCutKind()

  let signalEff = if signalEfficiency > 0.0: signalEfficiency
                  else: readSignalEff() # else read from `config.toml` file

  if fadcVeto and calibFile.len == 0:
    raise newException(ValueError, "When using the FADC veto an H5 file containing the ⁵⁵Fe calibration data " &
      "of the same run period is required.")
  let fadcVetoes = determineFadcVetoCutoff(calibFile, vetoPercentile, fadcScaleCutoff)
  let vetoCfg = VetoSettings(
    # NN cut
    useNeuralNetworkCut: useNeuralNetworkCut,
    nnModelPath: nnModelPath,
    nnSignalEff: nnSignalEff,
    nnCutKind: nnCutKind,
    # lnL cut
    useLnLCut: useLnLCut,
    signalEfficiency: signalEff,
    # septem veto related
    centerChip: centerChip,
    numChips: numChips,
    clusterAlgo: clusterAlgo,
    searchRadius: searchRadius,
    dbscanEpsilon: dbscanEpsilon,
    useRealLayout: useRealLayout,
    #use sim
    usesimref: usesimref,
    # ToACut
    useToACut: useToACut,
    ToAcutValue: ToAcutValue,
    #ToAlnlcut
    useToAlnLCut:useToAlnLCut,
    # line veto
    lineVetoKind: lvKind,
    eccLineVetoCut: eccLineVetoCut,
    # FADC veto
    calibFile: calibFile,
    vetoPercentile: vetoPercentile,
    fadcScaleCutoff: fadcScaleCutoff,
    fadcVetoes: fadcVetoes
  )
  result = LikelihoodContext(
    cdlFile: cdlFile,
    year: year,
    region: region,
    morph: morphKind,
    energyDset: energyDset,
    timepix: timepix,
    stretch: stretch,
    numMorphedEnergies: numMorphedEnergies,
    # misc
    useTeX: useTeX,
    flags: flags,
    vetoCfg: vetoCfg,
    septemLineVetoEfficiencyFile: septemLineVetoEfficiencyFile,
    rngSeed: rngSeed,
    rnd: initRand(rngSeed),
    plotPath: plotPath,
  )

  if useLnLCut or readLogLData: # if user wants the data, be explicit about `readLogLData`
    case result.morph
    of mkNone:
      result.refSetTuple = result.readRefDsets()
    of mkLinear:
      if vetoCfg.usesimref:
        echo "use sim"
        var redf: DataFrame
        let ecc_path = TpxDir / "resources/Ecc_P_densitys.csv"
        let ldiv_path = TpxDir / "resources/ldiv_P_densitys.csv"
        let ftrans_path = TpxDir / "resources/ftrans_P_densitys.csv"
        let toa_path = TpxDir / "resources/ToA_P_densitys.csv"
        let (eccdf, eccEnergy_list)= readToAProbabilities(ecc_path)
        let eccref = getInterpolatedDfToA(eccdf,eccEnergy_list, "eccentricity",num = result.numMorphedEnergies)
        let (ldivdf, ldivEnergy_list)= readToAProbabilities(ldiv_path)
        let ldivref = getInterpolatedDfToA(ldivdf,ldivEnergy_list, "lengthDivRmsTrans",num = result.numMorphedEnergies)
        let (ftransdf, ftransEnergy_list)= readToAProbabilities(ftrans_path)
        let ftransref = getInterpolatedDfToA(ftransdf,ftransEnergy_list, "fractionInTransverseRms",num = result.numMorphedEnergies)
        let (toadf, toaEnergy_list)= readToAProbabilities(toa_path)
        let toaref = getInterpolatedDfToA(toadf,toaEnergy_list, "ToAlength",num = result.numMorphedEnergies)
        redf.add eccref
        redf.add ldivref
        redf.add ftransref
        redf.add toaref
        result.refDf = redf
        let lineEnergies = eccEnergy_list.mapIt(it.float)
        result.refDfEnergy = linspace((lineEnergies[0]/1000), (lineEnergies[^1]/1000), result.numMorphedEnergies)
      else:
        result.refDf = result.readRefDsetsDF()
          .getInterpolatedWideDf(num = result.numMorphedEnergies)
        let lineEnergies = getXrayFluorescenceLines()
        result.refDfEnergy = linspace(lineEnergies[0], lineEnergies[^1], result.numMorphedEnergies)
