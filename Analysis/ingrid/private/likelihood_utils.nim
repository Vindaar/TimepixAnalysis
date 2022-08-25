import strutils, os, sequtils, strformat
import nimhdf5, seqmath, ggplotnim

import cdl_cuts, hdf5_utils, geometry
import ../ingrid_types
import helpers/utils

proc `[]`*(cv: CutValueInterpolator, e: float): float =
  case cv.kind
  of mkNone: result = cv.cutTab[e.toRefDset()]
  of mkLinear:
    let idx = min(cv.cutEnergies.lowerBound(e), cv.cutEnergies.high)
    result = cv.cutValues[idx]

import sugar
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
  let refLowT = df.filter(f{string -> bool: `Dset` == refLow})["Hist", float]
  let refHighT = df.filter(f{string -> bool: `Dset` == refHigh})["Hist", float]
  let bins = df.filter(f{string -> bool: `Dset` == refLow})["Bins", float]
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
    dfMorph["Energy"] = constantColumn(E, dfMorph.len)
    dfMorph["Dset"] = constantColumn("Morph", dfMorph.len)
    result.add dfMorph

proc getInterpolatedWideDf*(df: DataFrame, num = 1000): DataFrame =
  ## returns a DF with `num` interpolated distributions using next neighbors
  ## for linear interpolation
  let lineEnergies = getXrayFluorescenceLines()
  let energies = linspace(lineEnergies[0], lineEnergies[^1], num)
  let xrayRef = getXrayRefTable()
  result = newDataFrame()
  for tup, subDf in groups(group_by(df, "Variable")):
    # echo subDf
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

proc readMorphKind(): MorphingKind
proc calcLikelihoodDataset*(h5f: var H5File,
                            cdlFile: string,
                            groupName: string, year: YearKind,
                            morphKind: MorphingKind,
                            energyDset: InGridDsetKind): seq[float]

proc calcLikelihoodDatasetIfNeeded*(h5f: var H5File,
                                    cdlFile: string,
                                    grp: string,
                                    logLDset: string,
                                    year: YearKind,
                                    energyDset: InGridDsetKind) =
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
    let logLData = h5f.calcLikelihoodDataset(cdlFile, grp, year, morphKind, energyDset)
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
                      body: untyped): untyped =
  ## Sorry for the nested templates with a whole bunch of injected variables
  ## future reader... :)
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
  echo "Opening file to build LogL from ", cdlFile, " for dset: ", dset
  withH5(cdlFile, "rw"):
    if "FrameworkKind" in h5f.attrs:
      frameworkKind = parseEnum[FrameworkKind](h5f.attrs["FrameworkKind", string])
    # open h5 file using template
    let
      energyStr = igEnergyFromCharge.toDset(frameworkKind)
      centerXStr = igCenterX.toDset(frameworkKind)
      centerYStr = igCenterY.toDset(frameworkKind)
      eccStr = igEccentricity.toDset(frameworkKind)
      lengthStr = igLength.toDset(frameworkKind)
      chargeStr = igTotalCharge.toDset(frameworkKind)
      rmsTransStr = igRmsTransverse.toDset(frameworkKind)
      ldivStr = igLengthDivRmsTrans.toDset(frameworkKind)
      fracRmsStr = igFractionInTransverseRms.toDset(frameworkKind)
      npixStr = igHits.toDset(frameworkKind)
      # read the datasets
      energy {.inject.} = h5f.readAs(grp_name / energyStr, float64)
      centerX {.inject.} = h5f.readAs(grp_name / centerXStr, float64)
      centerY {.inject.} = h5f.readAs(grp_name / centerYStr, float64)
      ecc {.inject.} = h5f.readAs(grp_name / eccStr, float64)
      length {.inject.} = h5f.readAs(grp_name / lengthStr, float64)
      charge {.inject.} = h5f.readAs(grp_name / chargeStr, float64)
      rmsTrans {.inject.} = h5f.readAs(grp_name / rmsTransStr, float64)
      ldivRms {.inject.} = h5f.readAs(grp_name / ldivStr, float64)
      fracRms {.inject.} = h5f.readAs(grp_name / fracRmsStr, float64)
      npix {.inject.} = h5f.readAs(grp_name / npixStr, float64)
      # get the cut values for this dataset
      cuts {.inject.} = cutsTab[dset]
      xrayCuts {.inject.} = xrayCutsTab[dset]
    body

template withLogLFilterCuts*(cdlFile, dset: string,
                              year: YearKind, energyDset: InGridDsetKind,
                              body: untyped): untyped =
  ## Note: see the injected variables in the `withCdlData` template!
  withCdlData(cdlFile, dset, year, energyDset):
    # for likelihood dataset: aside from `resources/calibration-cdl.h5`, every other file
    # may not yet have access to the likelihood dataset. So we have to check for that and
    # if it does not exist yet, it has to be calculated.
    let logLStr = igLikelihood.toDset(frameworkKind)
    if grp_name / logLStr notin h5f:
      raise newException(ValueError, "Likelihood dataset does not yet exist in cdlFile " & $cdlFile & ".")

    let logL {.inject.} = h5f.readAs(grp_name / logLStr, float64)
    for i {.inject.} in 0 ..< energy.len:
      let
        # first apply Xray cuts (see C. Krieger PhD Appendix B & C)
        regionCut  = inRegion(centerX[i], centerY[i], crSilver)
        xRmsCut    = rmsTrans[i].float >= xrayCuts.minRms and rmsTrans[i] <= xrayCuts.maxRms
        xLengthCut = length[i].float   <= xrayCuts.maxLength
        xEccCut    = ecc[i].float      <= xrayCuts.maxEccentricity
        # then apply reference cuts
        chargeCut  = charge[i].float   > cuts.minCharge and charge[i]   < cuts.maxCharge
        rmsCut     = rmsTrans[i].float > cuts.minRms    and rmsTrans[i] < cuts.maxRms
        lengthCut  = length[i].float   < cuts.maxLength
        pixelCut   = npix[i].float     > cuts.minPix
      # add event to likelihood if all cuts passed
      if allIt([regionCut, xRmsCut, xLengthCut, xEccCut, chargeCut, rmsCut, lengthCut, pixelCut], it):
        if logL[i] != Inf:
          body

template withXrayRefCuts*(cdlFile, dset: string,
                          year: YearKind, energyDset: InGridDsetKind,
                          body: untyped): untyped =
  ## Note: see the injected variables in the `withCdlData` template!
  withCdlData(cdlFile, dset, year, energyDset):
    for i {.inject.} in 0 ..< energy.len:
      let
        # first apply Xray cuts (see C. Krieger PhD Appendix B & C)
        regionCut  = inRegion(centerX[i], centerY[i], crSilver)
        # then apply reference cuts
        chargeCut  = charge[i].float   > cuts.minCharge and charge[i]   < cuts.maxCharge
        rmsCut     = rmsTrans[i].float > cuts.minRms    and rmsTrans[i] < cuts.maxRms
        lengthCut  = length[i].float   < cuts.maxLength
        pixelCut   = npix[i].float     > cuts.minPix
      # add event to likelihood if all cuts passed
      if allIt([regionCut, chargeCut, rmsCut, lengthCut, pixelCut], it):
        body

proc genRefData*(cdlFile, dset: string,
                 year: YearKind, energyDset: InGridDsetKind):
                   (seq[seq[float]], seq[seq[float]], seq[seq[float]]) =
  var eccs = newSeq[float]()
  var ldiv = newSeq[float]()
  var frac = newSeq[float]()
  withXrayRefCuts(cdlFile, dset, year, energyDset):
    # read filtered data
    eccs.add ecc[i]
    ldiv.add ldivRms[i]
    frac.add fracRms[i]
  # now compute histograms
  proc assignRes(x: seq[float], name: string): seq[seq[float]] =
    var
      numBins: int
      minVal: float
      maxVal: float
    case year
    of yr2014:
      (numBins, minVal, maxVal) = cdlToXrayBinning2014(name)
    of yr2018:
      (numBins, minVal, maxVal) = cdlToXrayBinning2018(name)
    let (hist, bins) = histogram(x,
                                 numBins,
                                 range = (minVal, maxVal),
                                 upperRangeBinRight = false)
    # combine the hist bins data to a seq2D
    result = @[bins[0 .. ^2], hist.mapIt(it.float)]
  result[0] = assignRes(eccs, igEccentricity.toDset(fkTpa))
  result[1] = assignRes(ldiv, igLengthDivRmsTrans.toDset(fkTpa))
  result[2] = assignRes(frac, igFractionInTransverseRms.toDset(fkTpa))

proc readRefDsetsDF(cdlFile: string,
                    yearKind: YearKind): DataFrame =
  ## reads the reference datasets from the `refFile` and returns them.

  const xray_ref = getXrayRefTable()
  var
    ecc_ref = initTable[string, histTuple]()
    lengthDivRmsTrans_ref = initTable[string, histTuple]()
    fracRmsTrans_ref = initTable[string, histTuple]()
    result = newDataFrame()
  let energies = getXrayFluorescenceLines()
  for idx in 0 ..< xray_ref.len:
    ## read the data from the CDL file and generate the reference data using cuts
    let dsetName = xray_ref[idx]
    let E = energies[idx]
    let (ecc, ldiv, frac) = genRefData(cdlFile, dsetName, yearKind, igEnergyFromCharge)
    template addDf(h: histTuple, key: string, E: float, dkKind: InGridDsetKind) =
      var dfDset = toDf({ "Bins" : h.bins, "Hist" : h.hist,
                          "Dset" : key, "Energy" : E,
                          "Variable" : dkKind.toDset(fkTpa) })
      dfDset["Dset"] = dset_name
      result.add dfDset
    addDf((bins: ecc[0], hist: ecc[1]), dsetName, E, igEccentricity)
    addDf((bins: ldiv[0], hist: ldiv[1]), dsetName, E, igLengthDivRmsTrans)
    addDf((bins: frac[0], hist: frac[1]), dsetName, E, igFractionInTransverseRms)

proc buildLogLHist*(cdlFile, dset: string,
                    year: YearKind,
                    energyDset: InGridDsetKind,
                    region: ChipRegion = crGold): tuple[logL, energy: seq[float]] =
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
  withLogLFilterCuts(cdlFile, dset, year, energyDset):
    result[0].add logL[i]
    result[1].add energy[i]
  # now create plots of all ref likelihood distributions
  echo "max is inf ? ", min(result[0]), " ", max(result[0])

proc computeLogLDistributions*(cdlFile: string, yearKind: YearKind,
                               energyDset: InGridDsetKind,
                               region: ChipRegion = crGold): DataFrame =
  ## Computes the LogL distributions from thc CDL data file (`cdlFile`) by applying
  ## both sets of cuts (`getXraySpetrcumCuts` and `getEnergyBinMinMaxVals201*`) to the
  ## data in `buildLogLHist` and binning it according to the number and bin width
  ## that we use for the logL distributions.
  const
    xrayRef = getXrayRefTable()
    # logL binning range
    nbins = 200 # TODO: Increase number of bins for logL cut value closer to target?
    # range of histogram in logL
    logLrange = (0.0, 30.0)

  # get the correct binning for the histograms
  let bins = linspace(logLrange[0], logLrange[1], nbins + 1, endpoint = true)
  let energies = getXrayFluorescenceLines()
  for idx, dset in xrayRef:
    # compute the histogram of the CDL data
    let hist = buildLogLHist(cdlFile, dset, yearKind, energyDset, region)[0]
    var df = toDf( {"Bins" : bins[0 .. ^2], "Hist" : histogram(hist, nbins, logLrange)[0] })
    df["Dset"] = constantColumn(dset, df.len)
    df["Energy"] = constantColumn(energies[idx], df.len)
    result.add df

proc readLogLVariableData*(h5f: var H5File,
                           groupName: string,
                           energyDset: InGridDsetKind
                          ):
                             (seq[float], seq[float], seq[float], seq[float]) =
  # get the datasets needed for LogL
  doAssert energyDset in {igEnergyFromPixel, igEnergyFromCharge}
  let
    ecc = h5f[(groupName / "eccentricity"), float64]
    lengthDivRmsTrans = h5f[(groupName / "lengthDivRmsTrans"), float64]
    fracRmsTrans = h5f[(groupName / "fractionInTransverseRms"), float64]
    # energy to choose correct bin
    energies = h5f[(groupName / energyDset.toDset), float64]
  result = (ecc, lengthDivRmsTrans, fracRmsTrans, energies)

proc readRefDsets*(cdlFile: string, yearKind: YearKind): tuple[ecc, ldivRms, fracRms: Table[string, histTuple]] =
  ## Reads the data from the CDL file and generates the reference distributions based
  ## on the cuts defined in `cdl_cuts.nim`
  # create a table, which stores the reference datasets from the ref file
  const xray_ref = getXrayRefTable()
  var
    ecc_ref = initTable[string, histTuple]()
    lengthDivRmsTrans_ref = initTable[string, histTuple]()
    fracRmsTrans_ref = initTable[string, histTuple]()
    df = newDataFrame()
  for dset_name in values(xray_ref):
    ## read the data from the CDL file and generate the reference data using cuts
    let (ecc, ldiv, frac) = genRefData(cdlFile, dsetName, yearKind, igEnergyFromCharge)
    var dfDset = toDf({ "Eccentricity" : ecc[0],
                        "Ecc #" : ecc[1],
                        "L / RMS_trans" : ldiv[0],
                        "L / RMS_trans #" : ldiv[1],
                        "fracRmsTrans" : frac[0],
                        "fracRmsTrans #" : frac[1] })
    dfDset["Dset"] = dset_name
    df.add dfDset

  block RefPlots:
    var labelOrder = initTable[Value, int]()
    for idx, el in xrayRef:
      labelOrder[%~ el] = idx
    ggplot(df, aes("Eccentricity", "Ecc #", fill = "Dset")) +
      ggridges("Dset", overlap = 1.75, labelOrder = labelOrder) +
      geom_histogram(stat = "identity", position = "identity", alpha = 0.5, hdKind = hdOutline) +
      ggtitle(&"Eccentricity of reference file, year: {yearKind}") +
      ggsave(&"out/eccentricity_ridgeline_{cdlFile.extractFilename}_{yearKind}.pdf",
              width = 800, height = 480)
    ggplot(df, aes("L / RMS_trans", "L / RMS_trans #", fill = "Dset")) +
      ggridges("Dset", overlap = 1.75, labelOrder = labelOrder) +
      geom_histogram(stat = "identity", position = "identity", alpha = 0.5, hdKind = hdOutline) +
      ggtitle(&"L / RMS_trans of reference file, year: {yearKind}") +
      ggsave(&"out/lengthDivRmsTrans_ridgeline_{cdlFile.extractFilename}_{yearKind}.pdf",
              width = 800, height = 480)
    ggplot(data = df.filter(f{Value: isNull(df["fracRmsTrans"][idx]) == (%~ false)}),
           aes("fracRmsTrans", "fracRmsTrans #", fill = "Dset")) +
      ggridges("Dset", overlap = 1.75, labelOrder = labelOrder) +
      geom_histogram(stat = "identity", position = "identity", alpha = 0.5, hdKind = hdOutline) +
      ggtitle(&"fracRmsTrans of reference file, year: {yearKind}") +
      ggsave(&"out/fracRmsTrans_ridgeline_{cdlFile.extractFilename}_{yearKind}.pdf",
              width = 800, height = 480)
  result = (ecc: ecc_ref, ldivRms: lengthDivRmsTrans_ref, fracRms: fracRmsTrans_ref)

func calcLikelihoodForEvent*(energy, eccentricity, lengthDivRmsTrans, fracRmsTrans: float,
                             refSetTuple: tuple[ecc,
                                                ldivRms,
                                                fracRms: Table[string, histTuple]]): float =
  ## XXX: Manual tuple unpacking due to: https://github.com/nim-lang/Nim/issues/19364
  let ecc_ref = refSetTuple[0]
  let lengthDivRmsTrans_ref = refSetTuple[1]
  let fracRmsTrans_ref = refSetTuple[2]
  # try simple logL calc
  let refset = toRefDset(energy) # / 1000.0) division needed for E from Pix,
                                 # since that is in eV inst of keV
  result += logLikelihood(ecc_ref[refset][1],
                          eccentricity,
                          ecc_ref[refset][0])
  result += logLikelihood(lengthDivRmsTrans_ref[refset][1],
                          lengthDivRmsTrans,
                          lengthDivRmsTrans_ref[refset][0])
  result += logLikelihood(fracRmsTrans_ref[refset][1],
                          fracRmsTrans,
                          fracRmsTrans_ref[refset][0])
  result *= -1.0

proc calcMorphedLikelihoodForEvent*(eccentricity, lengthDivRmsTrans, fracRmsTrans: float,
                                    refDf: DataFrame, idx: int): float =
  # try simple logL calc
  var
    eccDf {.global.}, ldivDf {.global.}, fracDf {.global.}: DataFrame
  once:
    eccDf = refDf.filter(f{string -> bool: `Variable` == "igEccentricity"})
    ldivDf = refDf.filter(f{string -> bool: `Variable` == "igLengthDivRmsTrans"})
    fracDf = refDf.filter(f{string -> bool: `Variable` == "igFractionInTransverseRms"})
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

proc calcLikelihoodDataset*(h5f: var H5File,
                            cdlFile: string,
                            groupName: string,
                            year: YearKind,
                            morphKind: MorphingKind,
                            energyDset: InGridDsetKind # the energy dataset to use (pixel / charge)
                           ): seq[float] =
  const num = 1000
  let (ecc,
       lengthDivRmsTrans,
       fracRmsTrans,
       energies) = h5f.readLogLVariableData(groupName, energyDset)


  var refSetTuple {.global.}: tuple[ecc, ldivRms, fracRms: Table[string, histTuple]]
  var refDf {.global.}: DataFrame
  var refDfEnergy {.global.}: seq[float]
  case morphKind
  of mkNone:
    once:
      refSetTuple = cdlFile.readRefDsets(year)
  of mkLinear:
    ## XXX: `refDf` will `only` be defined the first time. However, this is
    ## fine, because it's only used the first time in `calcMorphedLikelihoodForEvent`
    ## as well! However, this *must* be fixed.
    once:
      refDf = cdlFile.readRefDsetsDF(year)
        .getInterpolatedWideDf(num = num)
      let lineEnergies = getXrayFluorescenceLines()
      refDfEnergy = linspace(lineEnergies[0], lineEnergies[^1], num)

  # create seq to store data logL data for this chip
  ## create a crazy man's plot
  proc crazyPlot(name: string) =
    let dfFiltered = refDf.filter(f{string: `Variable` == name})
    let histKeysName = dfFiltered.getKeys.filterIt("Hist" in it)
    var df = dfFiltered
      .gather(histKeysName, key = "Idx", value = "Hist")
    var labelOrder = initTable[Value, int]()
    for i in 0 ..< num:
      labelOrder[%~ ("Hist_" & $i)] = i
    echo df
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
  echo "[INFO]: Performing likelihood compute using morph kind: ", morphKind, " for chip: ", groupName
  for i in 0 .. ecc.high:
    case morphKind
    of mkNone:
      let logL = calcLikelihoodForEvent(energies[i], ecc[i], lengthDivRmsTrans[i], fracRmsTrans[i],
                                        refSetTuple)
      # add logL to the sequence. May be Inf though
      result[i] = logL
    of mkLinear:
      let idx = min(refDfEnergy.lowerBound(energies[i]), num - 1)
      let logL = calcMorphedLikelihoodForEvent(ecc[i], lengthDivRmsTrans[i], fracRmsTrans[i],
                                               refDf, idx)
      # add logL to the sequence. May be Inf though
      result[i] = logL

import parsetoml
template withConfig(body: untyped): untyped =
  const sourceDir = currentSourcePath().parentDir
  let config {.inject.} = parseToml.parseFile(sourceDir / "../config.toml")
  body

proc readMorphKind(): MorphingKind =
  ## reads the `morphingKind` field from the TOML file
  withConfig:
    result = parseEnum[MorphingKind](config["Likelihood"]["morphingKind"].getStr)

proc readSignalEff(): float =
  ## reads the `signalEfficiency` field from the TOML file
  withConfig:
    result = config["Likelihood"]["signalEfficiency"].getFloat

proc determineCutValue*[T: seq | Tensor](hist: T, eff: float): int =
  ## given a histogram `hist`, determine the correct bin to cut at to achieve
  ## a software efficiency of `eff`
  var
    cur_eff = 0.0
    last_eff = 0.0
  let hist_sum = hist.sum.float
  while cur_eff < eff:
    inc result
    last_eff = cur_eff
    cur_eff = hist[0 .. result].sum.float / hist_sum
  echo "Efficiency is at ", cur_eff, " and last Eff was ", last_eff

proc calcCutValueTab*(cdlFile: string, yearKind: YearKind,
                      energyDset: InGridDsetKind): CutValueInterpolator =
  ## returns a table mapping the different CDL datasets to the correct cut values
  ## based on the chip center region
  # read signal efficiency (default 80%) from TOML file
  let efficiency = readSignalEff()
  let morphKind = readMorphKind()
  let xray_ref = getXrayRefTable()
  let logHists = computeLogLDistributions(cdlFile, yearKind, energyDset, crGold)
  case morphKind
  of mkNone:
    # get the cut value for a software efficiency of 80%
    result = initCutValueInterpolator(morphKind)
    for tup, subDf in groups(logHists.group_by("Dset")):
      echo "Starting FOR LOOP--------------------------------------------------"
      let dset = tup[0][1].toStr
      let cutVal = determineCutValue(subDf["Hist", float], efficiency)
      # incl the correct values for the logL cut values
      result[dset] = subDf["Bins", float][cutVal]
      echo "Cut value of ", dset, " is ", result[dset]
    #if true: quit()
  of mkLinear:
    ## given logLHist compute interpolated DF from it
    ## first need a logLHist in DF format
    result = initCutValueInterpolator(morphKind)
    let dfInterp = logHists.getInterpolatedDf()
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
      let cutVal = determineCutValue(subDf["Hist", float], efficiency)
      # incl the correct values for the logL cut values
      energies.add energy
      cutVals.add subDf["Bins", float][cutVal]
    let sorted = zip(energies, cutVals).sortedByIt(it[0])
    result.cutEnergies = sorted.mapIt(it[0])
    result.cutValues = sorted.mapIt(it[1]).toTensor
  else:
    result = getChristophCutVals()
  when not defined(release) or defined(DEBUG):
    #echo logHists
    echo "logHists ", logHists
    #echo "Bins are ", bins
    echo "Cut values are ", result
    #echo mapIt(logHists, it.sum)
    #echo "Corresponding to logL values of ", result
