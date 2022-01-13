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
    var dfMorph = seqsToDf({"Bins" : bins, "Hist" : res})
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
    echo subDf
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

proc readRefDsetsDF(refFile: string,
                    yearKind: YearKind): DataFrame =
  ## reads the reference datasets from the `refFile` and returns them.
  var h5ref = H5open(refFile, "r")

  const xray_ref = getXrayRefTable()
  for dkKind in [igEccentricity, igLengthDivRmsTrans, igFractionInTransverseRms]:
    # create a table, which stores the reference datasets from the ref file
    let dset = dkKind.toDset(fkTpa)
    var tab = initTable[string, histTuple]()
    let xrayRef = getXrayRefTable()
    let energies = getXrayFluorescenceLines()
    let ranges = getEnergyBinning()
    for idx in 0 ..< xray_ref.len:
      let dset_name = xray_ref[idx]
      var data = h5ref[(dset_name / dset).dset_str]
      tab[dset_name] = data.readAs(float64).reshape2D(data.shape).splitSeq(float64)
      var dfDset = seqsToDf({ "Bins" : tab[dset_name].bins, "Hist" : tab[dset_name].hist })
      dfDset["Dset"] = constantColumn(xray_ref[idx], dfDset.len)
      dfDset["Energy"] = constantColumn(energies[idx], dfDset.len)
      dfDset["Variable"] = constantColumn($dkKind, dfDset.len)
      result.add dfDset

proc calcLikelihoodDataset*(h5f: var H5File, refFile: string,
                           groupName: string, year: YearKind): seq[float]

proc buildLogLHist*(cdlFile, refFile, dset: string,
                    year: YearKind,
                    region: ChipRegion = crGold): tuple[logL, energy: seq[float]] =
  ## given a file `h5file` containing a CDL calibration dataset
  ## `dset` apply the cuts on all events and build the logL distribution
  ## for the energy range.
  ## `dset` needs to be of the elements contained in the returned map
  ## of tos_helpers.`getXrayRefTable`
  ## Default `region` is the gold region
  ## Returns a tuple of the actual `logL` values and the corresponding `energy`.
  var grp_name = cdlPrefix($year) & dset
  # create global vars for xray and normal cuts table to avoid having
  # to recreate them each time
  let xrayCutsTab {.global.} = getXrayCleaningCuts()
  var cutsTab {.global.}: Table[string, Cuts]
  case year
  of yr2014:
    cutsTab = getEnergyBinMinMaxVals2014()
  of yr2018:
    cutsTab = getEnergyBinMinMaxVals2018()

  var frameworkKind = fkMarlin

  echo "Opening file to build LogL from ", cdlFile
  withH5(cdlFile, "rw"):
    if "FrameworkKind" in h5f.attrs:
      frameworkKind = parseEnum[FrameworkKind](h5f.attrs["FrameworkKind", string])
    # open h5 file using template
    let
      energyStr = igEnergyFromCharge.toDset(frameworkKind)
      logLStr = igLikelihood.toDset(frameworkKind)
      centerXStr = igCenterX.toDset(frameworkKind)
      centerYStr = igCenterY.toDset(frameworkKind)
      eccStr = igEccentricity.toDset(frameworkKind)
      lengthStr = igLength.toDset(frameworkKind)
      chargeStr = igTotalCharge.toDset(frameworkKind)
      rmsTransStr = igRmsTransverse.toDset(frameworkKind)
      npixStr = igHits.toDset(frameworkKind)

    # for likelihood dataset: aside from `resources/calibration-cdl.h5`, every other file
    # may not yet have access to the likelihood dataset. So we have to check for that and
    # if it does not exist yet, it has to be calculated.
    ## TODO: make running this dependent on whether the given `MorphingKind` is the same
    ## as the one used in the file currently! `if name in h5f or h5f.attrs["MorphingKind", string] == morphKind` ish
    if true: #grp_name / logLStr notin h5f:
      let logLData = h5f.calcLikelihoodDataset(refFile, grp_name, year)
      let loglDset = h5f.create_dataset(grp_name / logLStr,
                                        (logLData.len, 1),
                                        float64,
                                        overwrite = true)
      logLDset[logLDset.all] = logLData

    let
      energy = h5f.readAs(grp_name / energyStr, float64)
      logL = h5f.readAs(grp_name / logLStr, float64)
      centerX = h5f.readAs(grp_name / centerXStr, float64)
      centerY = h5f.readAs(grp_name / centerYStr, float64)
      ecc = h5f.readAs(grp_name / eccStr, float64)
      length = h5f.readAs(grp_name / lengthStr, float64)
      charge = h5f.readAs(grp_name / chargeStr, float64)
      rmsTrans = h5f.readAs(grp_name / rmsTransStr, float64)
      npix = h5f.readAs(grp_name / npixStr, float64)
      # get the cut values for this dataset
      cuts = cutsTab[dset]
      xrayCuts = xrayCutsTab[dset]
    result[0] = newSeqOfCap[float](energy.len)
    result[1] = newSeqOfCap[float](energy.len)
    for i in 0 .. energy.high:
      let
        # first apply Xray cuts (see C. Krieger PhD Appendix B & C)
        regionCut  = inRegion(centerX[i], centerY[i], crSilver)
        xRmsCut    = rmsTrans[i] >= xrayCuts.minRms and rmsTrans[i] <= xrayCuts.maxRms
        xLengthCut = length[i]   <= xrayCuts.maxLength
        xEccCut    = ecc[i]      <= xrayCuts.maxEccentricity
        # then apply reference cuts
        chargeCut  = charge[i]   > cuts.minCharge and charge[i]   < cuts.maxCharge
        rmsCut     = rmsTrans[i] > cuts.minRms    and rmsTrans[i] < cuts.maxRms
        lengthCut  = length[i]   < cuts.maxLength
        pixelCut   = npix[i]     > cuts.minPix
      # add event to likelihood if all cuts passed
      if allIt([regionCut, xRmsCut, xLengthCut, xEccCut, chargeCut, rmsCut, lengthCut, pixelCut], it):
        if logL[i] != Inf:
          result[0].add logL[i]
          result[1].add energy[i]

  # now create plots of all ref likelihood distributions
  echo "max is inf ? ", min(result[0]), " ", max(result[0])

proc computeLogLDistributions*(cdlFile, refFile: string, yearKind: YearKind,
                               region: ChipRegion = crGold): DataFrame =
  ## Computes the LogL distributions from thc CDL data file (`cdlFile`) by applying
  ## both sets of cuts (`getXraySpetrcumCuts` and `getEnergyBinMinMaxVals201*`) to the
  ## data in `buildLogLHist` and binning it according to the number and bin width
  ## that we use for the logL distributions.
  const
    xray_ref = getXrayRefTable()
    # logL binning range
    nbins = 200 # NO CHANGE IF SET TO 200
    # range of histogram in logL
    logLrange = (0.0, 30.0)

  let
    # get the correct binning for the histograms
    bins = linspace(logLrange[0], logLrange[1], nbins + 1, endpoint = true)
    # get the raw log likelihood histograms (seq of individual values). Takes into account
    # cuts on properties and chip region
    rawLogHists = mapIt(toSeq(values(xray_ref)),
                        # build hist and get `[0]` to only get `logL` values
                        buildLogLHist(cdlFile, refFile, it, yearKind, region)[0]
    )
  # given raw log histograms, create the correctly binned histograms from it
  let energies = getXrayFluorescenceLines()
  for idx, val in xrayRef:
    var df = seqsToDf( {"Bins" : bins[0 .. ^2], "Hist" : histogram(rawLogHists[idx], nbins, logLrange)[0] })
    df["Dset"] = constantColumn(val, df.len)
    df["Energy"] = constantColumn(energies[idx], df.len)
    result.add df

proc readLogLVariableData*(h5f: var H5File,
                           groupName: string):
                             (seq[float], seq[float], seq[float], seq[float]) =
  # get the datasets needed for LogL
  let
    ecc = h5f[(groupName / "eccentricity"), float64]
    lengthDivRmsTrans = h5f[(groupName / "lengthDivRmsTrans"), float64]
    fracRmsTrans = h5f[(groupName / "fractionInTransverseRms"), float64]
    # energy to choose correct bin
    energies = h5f[(groupName / "energyFromCharge"), float64]
  result = (ecc, lengthDivRmsTrans, fracRmsTrans, energies)

proc readRefDsets*(refFile: string, yearKind: YearKind): tuple[ecc, ldivRms, fracRms: Table[string, histTuple]] =
  ## reads the reference datasets from the `refFile` and returns them.
  var h5ref = H5open(refFile, "r")
  # create a table, which stores the reference datasets from the ref file
  const xray_ref = getXrayRefTable()

  # check if `FrameworkKind` defined in root of `h5ref` and if it is use the naming
  # scheme from there, otherwise use fkMarlin
  var frameworkKind: FrameworkKind = fkMarlin
  if "FrameworkKind" in h5ref.attrs:
    frameworkKind = parseEnum[FrameworkKind](h5ref.attrs["FrameworkKind", string])

  var
    ecc_ref = initTable[string, histTuple]()
    lengthDivRmsTrans_ref = initTable[string, histTuple]()
    fracRmsTrans_ref = initTable[string, histTuple]()

  var
    ## TODO: the following uses the `toDset` proc, which does not make sense for the original
    ## `XrayReferenceFile.h5` file, since that uses notation different from both normal Marlin
    ## and TPA. Also the `toDset` proc correctly fails for `igLengthDivRmsTrans` field, because
    ## it does not exist in Marlin. That of course is not the case for the reference file, where
    ## it is stored as ``lengthdivbyrmsy``.
    eccStr: string
    ldivrmsStr: string
    frmstStr: string
  case frameworkKind
  of fkMarlin:
    eccStr = "excentricity"
    ldivrmsStr = "lengthdivbyrmsy"
    frmstStr = "fractionwithinrmsy"
  of fkTpa:
    eccStr = igEccentricity.toDset(frameworkKind)
    ldivrmsStr = igLengthDivRmsTrans.toDset(frameworkKind)
    frmstStr = igFractionInTransverseRms.toDset(frameworkKind)

  var df: DataFrame
  for dset_name in values(xray_ref):
    # naming scheme does not depend on the actual data being processed, but only on what was used to
    # generate the `XrayReferenceFile.h5`
    var
      ecc = h5ref[(dset_name / eccStr).dset_str]
      ldivrms = h5ref[(dset_name / ldivrmsStr).dset_str]
      frmst = h5ref[(dset_name / frmstStr).dset_str]

    # to get the reference datasets, we read from the H5DataSet, reshape it to
    # a (N, 2) nested seq and finally split the two columns into two individual
    # sequences converted to float64
    ecc_ref[dset_name] = ecc.readAs(float64).reshape2D(ecc.shape).splitSeq(float64)
    lengthDivRmsTrans_ref[dset_name] = ldivrms.readAs(float64).reshape2D(ldivrms.shape).splitSeq(float64)
    fracRmsTrans_ref[dset_name] = frmst.readAs(float64).reshape2D(frmst.shape).splitSeq(float64)

    var dfDset = seqsToDf({ "Eccentricity" : ecc_ref[dset_name].bins, "Ecc #" : ecc_ref[dset_name].hist,
                            "L / RMS_trans" : lengthDivRmsTrans_ref[dset_name].bins,
                            "L / RMS_trans #" : lengthDivRmsTrans_ref[dset_name].hist,
                            "fracRmsTrans" : fracRmsTrans_ref[dset_name].bins,
                            "fracRmsTrans #" : fracRmsTrans_ref[dset_name].hist })
    dfDset["Dset"] = constantColumn(dset_name, dfDset.len)
    df.add dfDset

  block RefPlots:
    ggplot(df, aes("Eccentricity", "Ecc #", fill = "Dset")) +
      geom_histogram(stat = "identity", position = "identity", alpha = some(0.5)) +
      ggtitle(&"Eccentricity of reference file, year: {yearKind}") +
      ggsave(&"out/eccentricity_{refFile.extractFilename}_{yearKind}.pdf",
              width = 800, height = 480)
    ggplot(df, aes("L / RMS_trans", "L / RMS_trans #", fill = "Dset")) +
      geom_histogram(stat = "identity", position = "identity", alpha = some(0.5)) +
      ggtitle(&"L / RMS_trans of reference file, year: {yearKind}") +
      ggsave(&"out/lengthDivRmsTrans_{refFile.extractFilename}_{yearKind}.pdf",
              width = 800, height = 480)
    ggplot(data = df.filter(f{Value: isNull(df["fracRmsTrans"][idx]) == (%~ false)}),
           aes("fracRmsTrans", "fracRmsTrans #", fill = "Dset")) +
      geom_histogram(stat = "identity", position = "identity", alpha = some(0.5)) +
      ggtitle(&"fracRmsTrans of reference file, year: {yearKind}") +
      ggsave(&"out/fracRmsTrans_{refFile.extractFilename}_{yearKind}.pdf",
              width = 800, height = 480)

    let xrayRef = getXrayRefTable()
    var labelOrder = initTable[Value, int]()
    for idx, el in xrayRef:
      labelOrder[%~ el] = idx

    ggplot(df, aes("Eccentricity", "Ecc #", fill = "Dset")) +
      ggridges("Dset", overlap = 1.75, labelOrder = labelOrder) +
      geom_histogram(stat = "identity", position = "identity") +
      ggtitle(&"Eccentricity of reference file, year: {yearKind}") +
      ggsave(&"out/eccentricity_ridgeline_{refFile.extractFilename}_{yearKind}.pdf",
              width = 800, height = 480)
    ggplot(df, aes("L / RMS_trans", "L / RMS_trans #", fill = "Dset")) +
      ggridges("Dset", overlap = 1.75, labelOrder = labelOrder) +
      geom_histogram(stat = "identity", position = "identity") +
      ggtitle(&"L / RMS_trans of reference file, year: {yearKind}") +
      ggsave(&"out/lengthDivRmsTrans_ridgeline_{refFile.extractFilename}_{yearKind}.pdf",
              width = 800, height = 480)
    ggplot(data = df.filter(f{Value: isNull(df["fracRmsTrans"][idx]) == (%~ false)}),
           aes("fracRmsTrans", "fracRmsTrans #", fill = "Dset")) +
      ggridges("Dset", overlap = 1.75, labelOrder = labelOrder) +
      geom_histogram(stat = "identity", position = "identity") +
      ggtitle(&"fracRmsTrans of reference file, year: {yearKind}") +
      ggsave(&"out/fracRmsTrans_ridgeline_{refFile.extractFilename}_{yearKind}.pdf",
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

proc readMorphKind(): MorphingKind
proc calcLikelihoodDataset*(h5f: var H5File,
                            refFile: string,
                            groupName: string,
                            year: YearKind): seq[float] =
  const num = 1000
  let (ecc,
       lengthDivRmsTrans,
       fracRmsTrans,
       energies) = h5f.readLogLVariableData(groupName)

  var refSetTuple: tuple[ecc, ldivRms, fracRms: Table[string, histTuple]]
  var refDf: DataFrame
  var refDfEnergy: seq[float]
  var morphKind: MorphingKind
  once:
    morphKind = readMorphKind()
  case morphKind
  of mkNone: refSetTuple = readRefDsets(refFile, year)
  of mkLinear:
    once:
      refDf = readRefDsetsDF(refFile, year)
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
  #if true: quit()
  result = newSeq[float64](ecc.len)
  for i in 0 .. ecc.high:
    case morphKind
    of mkNone:
      let logL = calcLikelihoodForEvent(energies[i], ecc[i], lengthDivRmsTrans[i], fracRmsTrans[i],
                                        refSetTuple)
      # add logL to the sequence. May be Inf though
      result[i] = logL
      # If logL != Inf:
      #   discard
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
