# std
import std / [options, strformat, strutils, sequtils]

# nimble
import pkg / [nimhdf5, ggplotnim, arraymancer]

# local
import hdf5_utils, ggplot_utils, cdl_cuts, geometry, likelihood_utils
import ../ingrid_types


#[
This file contains the relevant code dealing with the stretching of CDL data to morph it
into a form more similar to a reference 55Fe dataset given to the `likelihood` call.

To better understand the idea, look at `tools/determineEffectiveEfficiency` and corresponding
notes in StatusAndProgress / journal / (resources/tpx3_effective_efficiency.org/pdf).
]#

## XXX: we have hardcoded the energies here. This is far from ideal!
const EscapepeakEnergy = 2.9
const PhotopeakEnergy = 5.9
const RmsCleaningCut = 1.5

const EscapepeakStr = "Escapepeak"
const PhotopeakStr = "Photopeak"

proc readRunData(h5f: H5File, chip: int): DataFrame =
  ## Simply reads the required data from the 55Fe file for the likelihood
  ## relevant parts.
  result = h5f.readDsets(chipDsets =
    some((chip: chip,
          dsets: @[igEnergyFromCharge.toDset(),
                   igRmsTransverse.toDset(),
                   igLengthDivRmsTrans.toDset(),
                   igFractionInTransverseRms.toDset(),
                   igEccentricity.toDset(),
                   igCenterX.toDset(),
                   igCenterY.toDset(),
                   igLength.toDset(),
                   igLikelihood.toDset()]))
  )

proc filterEvents(df: DataFrame, eccFilters: tuple[e, p: float]): DataFrame =
  ## This is equivalent to
  let xrayCutsTab {.global.} = getXrayCleaningCuts()
  template applyFilters(dfI: typed, eccCut: float): untyped {.dirty.} =
    let minRms = xrayCuts.minRms
    let maxRms = xrayCuts.maxRms
    let maxLen = xrayCuts.maxLength
    let maxEcc = eccCut
    dfI.filter(f{float -> bool: idx(igRmsTransverse.toDset()) < RmsCleaningCut and
      inRegion(idx("centerX"), idx("centerY"), crSilver) and
      idx("rmsTransverse") >= minRms and
      idx("rmsTransverse") <= maxRms and
      idx("length") <= maxLen and
      idx("eccentricity") <= maxEcc
    })
  if "Peak" notin df:
    raise newException(KeyError, "The given input DF does not contain a `Peak` column. Call the " &
      "`splitPeaks` procedure before calling this procedure!")
  result = newDataFrame()
  for (tup, subDf) in groups(df.group_by("Peak")):
    case tup[0][1].toStr
    of EscapepeakStr:
      let dset = EscapepeakEnergy.toRefDset()
      let xrayCuts = xrayCutsTab[dset]
      let dfF = applyFilters(subDf, eccFilters.e)
      result.add dfF
    of PhotopeakStr:
      let dset = PhotopeakEnergy.toRefDset()
      let xrayCuts = xrayCutsTab[dset]
      let dfF = applyFilters(subDf, eccFilters.p)
      result.add dfF
    else: doAssert false, "Invalid name"

proc splitPeaks(df: DataFrame): DataFrame =
  ## Split the input DF based on the energy column such that it contains a classification
  ## of clusters corresponding to the escape and photo peak each. Filter out every other cluster.
  let eD = igEnergyFromCharge.toDset()
  echo df
  result = df.mutate(f{float -> string: "Peak" ~ (
    if idx(eD) < 3.5 and idx(eD) > 2.5:
      EscapepeakStr
    elif idx(eD) > 4.5 and idx(eD) < 7.5:
      PhotopeakStr
    else:
      "None")})
    .filter(f{`Peak` != "None"})

proc readFe55File(fname: string): DataFrame =
  ## Given a single input file, performs application of the likelihood cut for all
  ## runs in it, split by photo & escape peak. Returns a DF with column indicating
  ## the peak, energy of each event & a column whether it passed the likelihood cut.
  ## Only events that are pass the input cuts are stored.
  let h5f = H5open(fname, "r")
  let fileInfo = h5f.getFileInfo()
  result = h5f.readRunData(fileInfo.centerChip)
    .splitPeaks()
    .filterEvents(eccFilters = (e: 5.0, p: 5.0)) # an eccentricity of 5 ought to be
                                                 # enough for everyone...
  discard h5f.close()

proc readCdlData(dfFe: DataFrame, energy: float, cdlFile: string): DataFrame =
  # map input fake energy to reference dataset
  let grp = energy.toRefDset()
  let passedInds = block:
    var res = newSeqOfCap[int](100_000)
    withLogLFilterCuts(cdlFile, grp, yr2018, igEnergyFromCharge, LogLCutDsets):
      res.add i
    res
  const xray_ref = getXrayRefTable()
  result = newDataFrame()
  let h5f = H5open(cdlFile, "r")
  for dset in IngridDsetKind:
    try:
      let d = dset.toDset()
      if d notin dfFe: continue # skip things not in input
      var cdlFiltered = newSeq[float](passedInds.len)
      let cdlRaw = h5f[cdlGroupName(grp, "2019", d), float]
      for i, idx in passedInds:
        cdlFiltered[i] = cdlRaw[idx]
      result[d] = cdlFiltered
    except AssertionError:
      continue
  discard h5f.close()

proc toEdf(x: seq[float], bins: seq[float]): seq[float] =
  ## Computes the EDF of the input data `x` given some `bins`.
  ##
  ## The bins are used as boundaries to count elements below each bin edge. The
  ## result contains `bins.len` elements, all in [0, 1]
  let xS = x.sorted
  var i = 0
  result = newSeq[float](bins.len)
  for j, b in bins:
    while i < xS.len and xS[i] <= b:
      inc i
    result[j] = i.float / x.len.float
  doAssert result.allIt(it <= 1.0)

proc kolmogorovSmirnov(x1, x2: seq[float]): float =
  ## Compute the Kolmogorov-Smirnov test for two datasets.
  ##
  ## The data is binned first to min & max of the combined data range and based on the
  ## associated EDF the KS test is performed.
  ##
  ## ``KS(x) = sup_x | EDF₁(x) - EDF₂(x) |``
  ##
  ## or in ~code
  ##
  ## ``KS(x) = min(abs(EDF₁ -. EDF₂(x)))``
  let range = (min: min(x1.min, x2.min), max: max(x1.max, x2.max))
  ## XXX: choose sensible binning & likely not hardcode it here!
  let bins = linspace(range[0], range[1], 200)
  let x1Edf = x1.toEdf(bins)
  let x2Edf = x2.toEdf(bins)
  result = 0.0
  for i, b in bins:
    result = max( result, abs(x1Edf[i] - x2Edf[i]) )

proc plotEccs(ecVal, ksVal: seq[float], peak: string) =
  ## Plots the walk of the iterative algorithm through the Kolmogorov-Smirnov / eccentricity space
  let ecDf = toDf({"ecc" : ecVal, "b" : toSeq(0 ..< ecVal.len), "ks" : ksVal})
  ggplot(ecDf, aes("ecc", "ks", color = "b")) +
    geom_point() + geom_line() +
    ggsave("/t/eccs_tested_" & $peak & ".pdf")

proc plotFeCdlDistributions(dfFe, dfCdl: DataFrame, dset, typ: string,
                            energy: float, suffix = "") =
  ## Creates a plot showing the current state of the 55Fe data and CDL data
  ## including the mean and a KDE of the distributions.
  let df = bind_rows([("55Fe", dfFe), ("CDL", dfCdl)], "Type")
  let dM = dfFe[dset, float].mean()
  let cM = dfCdl[dset, float].mean()
  let yM = 1.5
  let dfLine = toDf({ "x" : [dM, dM, cM, cM], "y" : [0.0, yM, 0.0, yM],
                      "Type" : ["55Fe", "55Fe", "CDL", "CDL"]})
  let fname = "/t/" & dset & "_" & $typ & $suffix & "_compare.pdf"
  echo "Saving plot: ", fname
  var dfDens = newDataFrame()
  for t, s in groups(df.group_by("Type")):
    let x = s[dset, float]
    let xs = linspace(x.min, x.max, 1000)
    let ys = kde(x, normalize = true, adjust = 2.5)
    dfDens.add toDf({"x" : xs, "y" : ys, "Type" : t[0][1]})
  ggplot(df, aes(dset, fill = "Type")) +
    geom_histogram(bins = 40, hdKind = hdOutline, alpha = 0.5, position = "identity", density = true) +
    geom_line(data = dfDens, aes = aes(x = "x", y = "y", color = "Type"), fillColor = color(0,0,0,0)) +
    geom_line(data = dfLine, aes = aes(x = "x", y = "y")) +
    ggtitle("Property: " & dset & ", Energy: " & $typ) +
    ggsave(fname)

proc applyMutation(fnTab: var Table[string, FormulaNode],
                   cdlDf, peakDf: DataFrame,
                   peak: string, energy: float
                  ): DataFrame =
  result = cdlDf.clone()

  for dset in keys(fnTab):
    let data = peakDf[dset, float]
    let cdl = cdlDf[dset, float]
    let
      dataMin = data.percentile(1) #.min
      dataMax = data.percentile(99) # max
      cdlMin = cdl.percentile(1) #min
      cdlMax = cdl.percentile(99) # max
    let dS = dset # local copy due to `f{}` capturing it
    let fn = f{float: dS ~ (idx(dS) - cdlMin) / (cdlMax - cdlMin) * (dataMax - dataMin) + dataMin}
    fnTab[dset] = fn
    result = result.mutate(fn)

## FIX UP arguments and clean up.
proc determineEccentricityCutoff(
  data, dataCdl: DataFrame,
  showPlots = false
     ): tuple[ecc_esc, ecc_pho: float] =
  ## Performs the iterative approach to the best possible eccentricity cutoff & associated
  ## CDL data stretching range based on the given 55Fe data (`data`) and the CDL data (`dataCdl`).

  ## XXX: make the following parameters?
  const Δε = 0.1
  const α = 0.1
  const absKS = 0.005
  const ΔKS = 1e-3
  for (tup, peakDf) in groups(data.group_by("Peak")):
    let peak = tup[0][1].toStr
    let energy = if peak == EscapepeakStr: EscapepeakEnergy else: PhotopeakEnergy
    var ecc = 1.9 ## Turn this into a parameter
    var lastKs = Inf
    var ks = Inf
    var fnTab = { "eccentricity" : f{""}, "lengthDivRmsTrans" : f{""},
                  "fractionInTransverseRms" : f{""} }.toTable

    # to plot the progress
    var ecVal = newSeq[float]()
    var ksVal = newSeq[float]()

    var n = 0
    var sign = 1.0
    while ks > absKS:
      let peakDf = peakDf.filter(f{`eccentricity` < ecc})
      var dfCdl = dataCdl.filter(f{`Peak` == peak})
      dfCdl = fnTab.applyMutation(dfCdl, peakDf, peak, energy) #, genPlots = false)
      ks = kolmogorovSmirnov(peakDf["eccentricity", float].toSeq1D,
                             dfCdl["eccentricity", float].toSeq1D)
      if showPlots:
        plotFeCdlDistributions(peakDf, dfCdl, "eccentricity", peak, energy, "_stretch")

      ecVal.add ecc
      ksVal.add ks

      # stop early and not adjust `ecc` if final conditions met
      if ks < absKS or abs(ks - lastKs) < ΔKS: break

      if ks > lastKs: # only possibly change sign if we're worse than before!
        sign = -sign
      let adj = sign * Δε * exp(- n.float * α)
      ecc = ecc - adj

      lastKs = ks
      inc n
      if showPlots:
        plotEccs(ecVal, ksVal, peak)

    if peak == EscapepeakStr:
      result.ecc_esc = ecc
    else:
      result.ecc_pho = ecc

    if showPlots:
      plotEccs(ecVal, ksVal, peak)

proc lineParams(x: seq[float], ys: seq[(float, float)]): LineParams =
  let y1 = ys[1]
  let y0 = ys[0]
  let m = (y1[0] / y1[1] - y0[0] / y0[1]) / (x[1] - x[0])
  let b = y0[0] / y0[1] - m * x[0]
  result = (m: m, b: b)

proc initCdlStretch*(fe55, cdlFile: string): Option[CdlStretch] =
  if fe55.len == 0: return none[CdlStretch]()

  let data = readFe55File(fe55)
  var dataCdl = newDataFrame()
  dataCdl.add(data.readCdlData(EscapepeakEnergy, cdlFile).mutate(f{"Peak" <- "Escapepeak"}))
  dataCdl.add(data.readCdlData(PhotopeakEnergy, cdlFile).mutate(f{"Peak" <- "Photopeak"}))

  let (eccEsc, eccPho) = determineEccentricityCutoff(data, dataCdl) #, showPlots = true)

  var lineTab = initTable[string, tuple[mins, maxs: LineParams]]()
  for dset in ["eccentricity", "lengthDivRmsTrans", "fractionInTransverseRms"]:
    var mins = newSeq[(float, float)]()
    var maxs = newSeq[(float, float)]()
    for (tup, subDf) in groups(data.group_by("Peak")):
      let peak = tup[0][1].toStr
      let ecc = if peak == EscapepeakStr: eccEsc else: eccPho
      let pdf = subDf.filter(f{`eccentricity` < ecc})
      let cdl = dataCdl.filter(f{`Peak` == peak})
      mins.add((pdf[dset, float].percentile(1), cdl[dset, float].percentile(1)))
      maxs.add((pdf[dset, float].percentile(99), cdl[dset, float].percentile(99)))

    let energies = @[EscapepeakEnergy, PhotopeakEnergy]
    lineTab[dset] = (mins: lineParams(energies, mins), maxs: lineParams(energies, maxs))


    #let df = toDf({"energy" : [2.9, 5.9], "mins" : mins, "maxs" : maxs})
    #ggplot(df, aes("energy", "mins")) +
    #  geom_point() +
    #  ggsave("/t/mins_" & $dset & ".pdf")
    #ggplot(df, aes("energy", "maxs")) +
    #  geom_point() +
    #  ggsave("/t/maxs_" & $dset & ".pdf")

  ## also do linear interpolation for 1 and 99 percentile of fe55 data? so that we have the
  ## correct

  ## what we finally need is not the effective efficiency cutoffs, but rather the scaling
  ## function that depends on the percentiles 1 and 99 of the data. For the CDL data we compute
  ## it, fine. But for the 55Fe data by definition we have to compute it.

  when false:
    let xrayCuts = getXrayCleaningCuts()
    let energies = getXrayFluorescenceLines()
    let eccCuts = toSeq(values(xrayCuts)).mapIt(it.maxEccentricity)
    echo eccCuts
    var df = toDf({"energies" : energies, "cuts" : eccCuts, "Type" : "CDL"})
    df.add toDf({"energies" : [2.9, 5.9], "cuts" : [ecc_esc, ecc_pho], "Type" : "55Fe"})
    let m = (ecc_pho / 1.3 - ecc_esc / 1.4) / (5.9 - 2.9)
    let b = ecc_pho / 1.3  - m * 5.9
    var extraCuts = newSeq[float]()
    for i, cut in eccCuts:
      extraCuts.add (energies[i] * m + b) * cut

    ## XXX generalize the above for each property in such a way as just using the actual limits
    ## of the data as the foundation?

  result = some(CdlStretch(fe55: fe55, eccEsc: eccEsc, eccPho: eccPho,
                           lines: lineTab))
