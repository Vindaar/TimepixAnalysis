import ggplotnim, seqmath, sequtils, os, sugar, strscans, strformat, strutils, sugar, sets
import ingrid / [tos_helpers, ingrid_types]
from arraymancer import tensor

import nimhdf5, numericalnim, unchained

import ggplotnim / ggplot_vegatex

import cligen

# get serialization support for DataFrame and Tensor
import datamancer / serialize

const Data2014 = "../../resources/background-rate-gold.2014+2015.dat"
const Ecol = "Energy"
const Ccol = "Counts"
const Rcol = "Rate"

type
  LogLFile = object
    name: string
    year: int
    totalTime: float
    df: DataFrame

var DsetNames = newSeq[string]()
for dkKind in InGridDsetKind:
  if dkKind notin {igNumClusters, igFractionInHalfRadius, igRadiusDivRmsTrans,
                    igRadius, igBalance, igLengthDivRadius, igInvalid, igHits, igTotalCharge, igEventNumber}:
    DsetNames.add dkKind.toDset(fkTpa)

var noisyPixels = newSeq[(int, int)]()
import unchained
defUnit(cm²)
proc scaleDset(data: Column, totalTime, factor: float): Column =
  ## scales the data in `data` according to the area of the gold region,
  ## total time and bin width. The data has to be pre binned of course!
  # XXX: make sure we never use wrong area if input data makes use of `--chipRegion` feature!
  var area = pow(0.95 - 0.45, 2).cm² # area of gold region!

  let pixSize = 55.MicroMeter * 55.MicroMeter
  let removedSize = noisyPixels.len * pixSize
  area = area - removedSize.to(cm²)

  const bin_width = 0.2 # 0.392
  const shutter_open = 1.0 # Does not play any role, because totalDuration takes
                           # this into account already!
  let scale = factor / (totalTime * shutter_open * area.float * bin_width) #* 1e5
  result = toColumn data.toTensor(float).map_inline(x * scale)

proc readTime(h5f: H5File): float =
  let lhGrp = h5f["/likelihood".grp_str]
  result = lhGrp.attrs["totalDuration", float]

proc readVetoCfg(h5f: H5File): (set[LogLFlagKind], VetoSettings) =
  let lhGrp = h5f[likelihoodGroupGrpStr]
  let ctx = deserializeH5[LikelihoodContext](h5f, "logLCtx", lhGrp.name,
                                             exclude = @["rnd", "refSetTuple", "refDf", "refDfEnergy"])
  result = (ctx.flags, ctx.vetoCfg)

proc calcVetoEfficiency(flags: set[LogLFlagKind], vetoCfg: VetoSettings): float =
  ## Calculates the combined efficiency of the detector given the used vetoes, the
  ## random coincidence rates of each and the FADC veto percentile used to compute
  ## the veto cut.
  const
    septemVetoRandomCoinc = 0.7841029411764704    # only septem veto random coinc based on bootstrapp
    lineVetoRandomCoinc = 0.8601764705882353      # lvRegular based on bootstrapped fake data
    septemLineVetoRandomCoinc = 0.732514705882353 # lvRegularNoHLC based on bootstrapped fake data
  result = 1.0 # starting efficiency
  if fkFadc in flags:
    doAssert vetoCfg.vetoPercentile > 0.5, "Veto percentile was below 0.5: " &
      $vetoCfg.vetoPercentile # 0.5 would imply nothing left & have upper percentile
    # input is e.g. `0.99` so: invert, multiply (both sided percentile cut), invert again
    result = 1.0 - (1.0 - vetoCfg.vetoPercentile) * 2.0
  if {fkSeptem, fkLineVeto} - flags == {}: # both septem & line veto
    result *= septemLineVetoRandomCoinc
  elif fkSeptem in flags:
    result *= septemVetoRandomCoinc
  elif fkLineVeto in flags:
    result *= lineVetoRandomCoinc

  # now multiply by lnL or MLP cut efficiency. LnL efficiency is likely always defined, but NN
  # only if NN veto was used
  if vetoCfg.nnEffectiveEff > 0.0:
    result *= vetoCfg.nnEffectiveEff
  elif vetoCfg.signalEfficiency > 0.0:
    result *= vetoCfg.signalEfficiency
  else:
    raise newException(ValueError, "Input has neither a well defined NN signal efficiency nor a regular likelihood " &
      "signal efficiency. Stopping.")

proc scaleDset(h5f: H5File, data: Column, factor: float): Column =
  result = scaleDset(data, h5f.readTime(), factor)

import macros
macro log(verbose: bool, body: untyped): untyped =
  result = newStmtList()
  for arg in body:
    let b = quote do:
      if `verbose`:
        let argStr = $`arg`
        echo "[INFO]:", argStr
    result.add b

proc extractYear(f: string, verbose: bool): int =
  ## We assume the filenames
  ## - always contain a year
  ## - always separated by underscores
  let fname = f.extractFilename
  var dummy: string
  if fname.scanf("$*_$i_$*", dummy, result, dummy):
    log(verbose):
      &"found a year {result}"

proc toIdx(x: float): int = (x / 14.0 * 256.0).round.int.clamp(0, 255)

proc readFiles(files: seq[string], names: seq[string], region: ChipRegion,
               energyDset: string,
               centerChip: int,
               energyMin, energyMax: float,
               readToA: bool,
               toaCutoff: seq[float],
               applyEfficiencyNormalization: bool,
               verbose: bool): seq[LogLFile] =
  ## reads all H5 files given and stacks them to a single
  ## DF. An additional column is added, which corresponds to
  ## the filename. That way one can later differentiate which
  ## events belong to what and decide if data is supposed to
  ## be accumulated or not
  ##
  ## After reading we will cut to the given `region`.
  doAssert names.len == 0 or names.len == files.len, "Need one name for each input file!"
  let noiseSet = noisyPixels.toHashSet
  for idx, file in files:
    log(verbose):
      &"Opening file: {file} of files : {files}"
    let h5f = H5open(file, "r")
    var df = h5f.readDsets(likelihoodBase(), some((centerChip, DsetNames)), verbose = verbose)
    df = df
      .rename(f{Ecol <- energyDset})
      .filter(f{float -> bool: (`centerX`.toIdx, `centerY`.toIdx) notin noiseSet })
      .filter(f{float: idx(Ecol) >= energyMin and idx(Ecol) <= energyMax})
    if readToA:
      let toaC = toaCutoff[idx]
      log(verbose):
        &"Cutting toa length to {toaC}, input df length: {df.len}"
      df = df.filter(f{`toaLength` < toaC})
      log(verbose):
        &"df len now {df.len}"

    let (flags, vetoCfg) = readVetoCfg(h5f)
    var totalTime = readTime(h5f)
    # apply the total efficiency based on all vetoes to the data
    let totalEff = calcVetoEfficiency(flags, vetoCfg)
    if applyEfficiencyNormalization:
      # just multiply total time by efficiency to get 'effective background time'
      totalTime *= totalEff

    doAssert not df.isNil, "Read DF is nil. Likely you gave a non existant chip number. Is " &
      $centerChip & " really the center chip in your input file?"
    let fname = if names.len > 0: names[idx]
                else: file.extractFilename
    df["File"] = constantColumn(fname, df.len)
    df["ε_total"] = totalEff
    df["ε_eff"] = if fkLogL in flags: vetoCfg.signalEfficiency else: vetoCfg.nnEffectiveEff
    df["Classifier"] = if fkLogL in flags: "LnL" else: "MLP"
    df["Scinti"] = fkScinti in flags
    df["FADC"] = fkFadc in flags
    df["Septem"] = fkSeptem in flags
    df["Line"] = fkLineVeto in flags
    log(verbose):
      &"File: {fname}"
      &"Elements before cut: {df.len}"
    df = df.filter(f{float: inRegion(`centerX`, `centerY`, region) == true})
    log(verbose):
      &"Elements after cut: {df.len}"
    result.add LogLFile(name: fname,
                        totalTime: totalTime,
                        df: df,
                        year: extractYear(fname, verbose))
    discard h5f.close()

proc read2014Data(log: bool): DataFrame =
  ## Reads the 2014 data from the CSV file and prepares it to be added to the
  ## DFs from `readFiles`.
  result = readCsv(Data2014, sep = ' ', header = "#")
    .rename(f{Rcol <- "Rate[/keV/cm²/s]"}, f{"yMax" <- "dRateUp"},
            f{"yMin" <- "dRateDown"}, f{Ecol <- "E[keV]"})
  let lastLine = toDf({ Ecol : @[10.1], Rcol : @[0.0], "yMin" : @[0.0], "yMax" : @[0.0] })
  result.add lastLine
  if not log:
    result = result.mutate(f{Rcol ~ 1e5 * `Rate`}, f{"yMin" ~ 1e5 * `yMin`},
                           f{"yMax" ~ 1e5 * `yMax`})
  result = result.mutate(f{"yMin" ~ `Rate` - `yMin`}, f{"yMax" ~ `Rate` + `yMax`},
                         f{Ecol ~ `Energy` - 0.1})
  result["Dataset"] = constantColumn("2014/15", result.len)
  result["Classifier"] = "LnL"
  result["ε_eff"] = 0.8
  result["Scinti"] = false
  result["FADC"] = false
  result["Septem"] = false
  result["Line"] = false
  result["ε_total"] = 0.8

proc histogram(df: DataFrame): DataFrame =
  ## Calculates the histogam of the energy data in the `df` and returns
  ## a histogram of the binned data
  ## TODO: allow to do this by combining different `File` values
  let (hist, bins) = histogram(df[Ecol].toTensor(float).toSeq1D,
                               range = (0.0, 20.0), bins = 100)
  result = toDf({ Ecol : bins, Ccol : concat(hist, @[0]) })

template sumIt(s: seq[typed], body: untyped): untyped =
  ## why the heck did I write this template?
  ## The code may be correct, but the usage was utterly wrong (unless I'm *now* extremely
  ## wrong, but I'm pretty sure I am not Jan 2022)
  var res: float
  for it {.inject.} in s:
    res += body
  res

proc copyScalars(to: var DataFrame, frm: DataFrame) =
  for k in getKeys(frm):
    let col = frm[k].unique
    if col.len == 1:
      withNativeTensor(col, x):
        to[k] = x[0]
    else:
      if k in ["ε_total", "ε_eff"]: # one efficiency value for each run period
        # compute the mean instead! (technically weighted by time would be better)
        to[k] = col.toTensor(float).mean

proc flatScale(files: seq[LogLFile], factor: float, verbose: bool, dropCounts = true): DataFrame =
  #var count {.global.} = 0
  var df: DataFrame
  for f in files:
    if not f.df.isNil:
      df.add f.df
  result = newDataFrame()

  proc getTotalTime(f: string): float =
    ## retrieves the correct `totalTime` for the given file `f` from the input `files`
    for file in files:
      if file.name == f:
        ## Accumulate all times of files that match `f` (either same input file or
        ## same `--names` argument
        result += file.totalTime

  when false: ## XXX: hack to support Tpx3. Implement properly!
    let fname = "tpx3" #tup[0][1].toStr
    let totalTime = fname.getTotalTime() # before: files.sumIt(it.totalTime)
  for tup, subDf in groups(df.group_by("File")):
    let fname = tup[0][1].toStr
    let totalTime = fname.getTotalTime() # before: files.sumIt(it.totalTime)
    #subDf.showBrowser()
  #for tup, subDf in groups(df.group_by("L<L_median")):
    var dfLoc = newDataFrame()
    #dfLoc = df.histogram() ## XXX: Tpx3
    dfLoc = subDf.histogram()
    dfLoc = dfLoc.mutate(f{float: "CountErr" ~ sqrt(`Counts`)})
    dfLoc[Rcol] = dfLoc[Ccol].scaleDset(totalTime, factor)
    dfLoc["totalTime"] = totalTime.Second.to(Hour).float
    dfLoc["RateErr"] = dfLoc["CountErr"].scaleDset(totalTime, factor)
    dfLoc["Dataset"] = constantColumn(fname, dfLoc.len) #"2017/18_" & $count, dfLoc.len)
    dfLoc = dfLoc.mutate(f{"yMin" ~ `Rate` - `RateErr`}, f{"yMax" ~ `Rate` + `RateErr`})
    dfLoc["yMin"] = dfLoc["yMin", float].map_inline:
      if x >= 0.0: x
      else: 0.0
    if dropCounts:
      dfLoc = dfLoc.drop(["Counts", "CountErr"])
    # copy over all the scalar fields in the data (e.g. `Classifier` and `ε`)
    dfLoc.copyScalars(subDf)
    result.add dfLoc
  log(verbose):
    result.pretty(-1)

proc flatScaleMedian(df: DataFrame, factor: float, totalTime: float, verbose: bool): DataFrame =
  result = newDataFrame()
  for tup, subDf in groups(df.group_by("Variable")):
    var dfVal = newDataFrame()
    for tupVal, subDfVal in groups(subDf.group_by("Value")):
      var dfLoc = newDataFrame()
      dfLoc = subDfVal.histogram()
      dfLoc = dfLoc.mutate(f{float: "CountErr" ~ sqrt(`Counts`)})
      dfLoc[Rcol] = dfLoc[Ccol].scaleDset(totalTime, factor)
      dfLoc["RateErr"] = dfLoc["CountErr"].scaleDset(totalTime, factor)
      dfLoc["Value"] = constantColumn(tupVal[0][1].toBool, dfLoc.len) #"2017/18_" & $count, dfLoc.len)
      #inc count
      dfLoc = dfLoc.mutate(f{"yMin" ~ `Rate` - `RateErr`}, f{"yMax" ~ `Rate` + `RateErr`})
      dfLoc["yMin"] = dfLoc["yMin"].toTensor(float).map_inline:
        if x >= 0.0: x
        else: 0.0
      dfLoc.drop("Counts")
      dfVal.add dfLoc
    dfVal["Variable"] = constantColumn(tup[0][1].toStr, dfVal.len)
    result.add dfVal
  log(verbose):
    result.pretty(-1)

proc calcIntegratedBackgroundRate(df: DataFrame, factor: float,
                                  energyRange: Slice[float] = 0.0 .. 10.0): float =
  ## returns the integrated background rate given by the Rate
  ## stored in the given DataFrame, integrated over the energy
  ## range stored in the DF (typically that should be 0 - 10 keV)
  ##
  ## It is assumed that the energy is given in keV and the rate in
  ## keV⁻¹ cm⁻² s⁻¹.
  let df = df.filter(f{float: `Energy` >= energyRange.a and `Energy` <= energyRange.b})
  let energies = df[Ecol].toTensor(float)
  let rate = df[Rcol].toTensor(float)
  result = trapz(rate.toSeq1D, energies.toSeq1D) / factor

proc computeMedianBools(df: DataFrame): DataFrame =
  var df = df
  var medianNames = newSeq[string]()
  proc computeMedianBool(dset: string) =
    let medianVal = df[dset, float].toSeq1D.median()
    let nameStr = $dset & "<" & $dset & "_median"
    medianNames.add nameStr
    df = df.mutate(f{float -> bool: nameStr ~ df[dset][idx] < medianVal})
  for dset in DsetNames:
    if dset != "energyFromCharge":
      computeMedianBool(dset)

  result = df.gather(medianNames, key = "Variable", value = "Value")

proc addNoisyPixels() =
  for x in 150 .. 175:
    for y in 130 .. 162:
      noisyPixels.add (x, y)

  for x in 125 .. 135:
    for y in 110 .. 120:
      noisyPixels.add (x, y)

proc plotMedianBools(df: DataFrame, fnameSuffix, title: string,
                     suffix, outpath: string, verbose: bool) =
  let transparent = color(0, 0, 0, 0)
  let fname = &"{outpath}/background_rate_median_bools_{fnameSuffix}.pdf"
  log(true): # always
    &"Storing plot in {fname}"
  #let df = df.select([Ecol, Rcol, "Variable", "Value"])
  log(verbose):
    df.pretty(-1)

  ggplot(df, aes(Ecol, Rcol, fill = "Value")) +
    facet_wrap("Variable", scales = "free") +
    geom_histogram(stat = "identity", position = "identity", alpha = 0.5,
                   color = transparent,
                   hdKind = hdOutline) +
    # too busy in small plot
    #geom_point(binPosition = "center", position = "identity") +
    #geom_errorbar(binPosition = "center",
    #              aes = aes(yMin = "yMin", yMax = "yMax"),
    #              errorBarKind = ebLines) +
    scale_y_continuous() +
    xlab("Energy [keV]") +
    ylab("Rate [10⁻⁵ keV⁻¹ cm⁻² s⁻¹]") +
    #xlim(0, 12.0) +
    ggtitle(&"{title}{suffix}") +
    ggsave(fname, width = 1920, height = 1080)

proc plotBackgroundRate(df: DataFrame, fnameSuffix, title: string,
                        outpath, outfile: string,
                        show2014: bool, suffix: string,
                        hidePoints, hideErrors, fill: bool,
                        useTeX, showPreliminary, genTikZ: bool,
                        showNumClusters, showTotalTime: bool,
                        topMargin,
                        yMax: float,
                        energyMin, energyMax: float,
                        logPlot: bool,
                        applyEfficiencyNormalization: bool
                       ) =
  var df = df # mutable copy
  if logPlot:
    # make sure to remove 0 entries if we do a log plot
    df = df.filter(f{idx(Rcol) > 0.0})

  var titleSuff = if title.len > 0: title
                  elif show2014:
                    "Background rate of Run 2 & 3 (2017/18) compared to 2014/15"
                  else:
                    "Background rate of Run 2 & 3 (2017/18)"
  let transparent = color(0, 0, 0, 0)
  let fname = if outfile.len > 0: outpath / outfile
              else: &"{outpath}/background_rate_{fnameSuffix}.pdf"
  log(true): # always write this!
    df.pretty(-1, precision = 8)
    &"INFO: storing plot in {fname}"
  var plt: GgPlot
  let numDsets = df.unique("Dataset").len

  ## XXX: fix me, it's wrong because here we already have the binned data!
  ##if showNumClusters:
  ##  titleSuff.add &" #clusters={}"


  if showTotalTime:
    if numDsets > 1:
      echo "[WARNING]: Printing total background time currently only supported " &
        "for single datasets."
    else:
      let time = df["totalTime", float][0]
      titleSuff.add &" background time={time.Hour}"

  if applyEfficiencyNormalization:
    titleSuff.add &", normalized by efficiency"

  if numDsets > 1 and fill:
    plt = ggplot(df, aes(Ecol, Rcol, fill = "Dataset")) +
      geom_histogram(stat = "identity", position = "identity", alpha = 0.5,
                     color = transparent,
                     hdKind = hdOutline)
  elif numDsets > 1:
    plt = ggplot(df, aes(Ecol, Rcol, color = "Dataset")) +
      geom_histogram(
        stat = "identity", position = "identity", alpha = 0.0,
        lineWidth = 2.0,
        hdKind = hdOutline
      )
  else:
    plt = ggplot(df, aes(Ecol, Rcol)) +
      geom_histogram(stat = "identity", position = "identity",
                     alpha = 0.5, color = transparent, hdKind = hdOutline) +
      minorGridLines()
      #scale_x_continuous(breaks = 20)
  if not hidePoints:
    plt = plt + geom_point(binPosition = "center", position = "identity")
  if not hideErrors:
    plt = plt + geom_errorbar(binPosition = "center",
                  aes = aes(yMin = "yMin", yMax = "yMax"),
                  errorBarKind = ebLines)
  # force y continuous scales & set limit
  if not logPlot:
    plt = plt + scale_y_continuous() +
      xlim(0, energyMax)
  else:
    plt = plt + scale_y_log10()

  if yMax > 0.0:
    plt = plt + ylim(0.0, yMax)
  if useTeX and showPreliminary:
    plt = plt + annotate("GridPix preliminary",
                         x = 4.3, y = 2.0,
                         rotate = -30.0,
                         font = font(16.0, color = color(0.92, 0.92, 0.92)),
                         backgroundColor = transparent)
  if useTeX:
    plt = plt +
    xlab(r"Energy [\si{keV}]") +
    ylab(r"Rate [\SI{1e-5}{keV⁻¹ cm⁻² s⁻¹}]", margin = 1.6) +
    #minorGridLines() +
    ggtitle(titleSuff) +
    themeLatex(fWidth = 0.9, width = 600, baseTheme = singlePlot)
    if genTikZ:
      plt + ggsave(fname.replace(".pdf", ".tex"), width = 600, height = 360, useTeX = true, onlyTikZ = true)
    else:
      plt + ggsave(fname, width = 600, height = 360, useTex = true, standalone = true)
  else:
    plt + xlab("Energy [keV]") +
    ylab("Rate [10⁻⁵ keV⁻¹ cm⁻² s⁻¹]") +
    ggtitle(titleSuff) +
    margin(top = topMargin) +
    ggsave(fname, width = 800, height = 480)

  df = df.drop(["Dataset", "yMin", "yMax"])
  df.writeCsv("/tmp/background_rate_data.csv")

proc plotEfficiencyComparison(files: seq[LogLFile], outpath: string, verbose: bool) =
  ## creates a plot comparing different logL cut efficiencies. Each filename should contain the
  ## efficiency in the format
  ## `_eff_<eff>.h5`
  ## - for each file:
  ##   - for each bin:
  ##     - eff / sqrt(counts)
  ## plot geom_line + geom_point of these values
  ##
  ## XXX: clean this up?
  var df = newDataFrame()
  for f in files:
    doAssert "_eff_" in f.name
    # ugly for now ## <-- this can come from veto efficiency now!
    let eff = f.name.split("_")[^1].dup(removeSuffix(".h5")).parseFloat
    var dfFlat = flatScale(@[f], 1e5, verbose, dropCounts = false)
    dfFlat["ε/√B"] = dfFlat["Counts", float].map_inline:
      if x > 0.0:
        eff / sqrt(x)
      else:
        0.0
    dfFlat["Eff"] = constantColumn(eff, dfFlat.len)
    df.add dfFlat
  ggplot(df.filter(f{`Energy` < 12.0}).arrange("Energy"), aes(ECol, "ε/√B", color = factor("Eff"))) +
    geom_line(aes = aes(shape = factor("Eff"))) +
    geom_point() +
    ggtitle("Comparison of signal efficiency ε / √Background") +
    ggsave(&"{outpath}/signal_efficiency_vs_sqrt_background.pdf",
           width = 1280, height = 800)
  ggplot(df.filter(f{`Energy` < 12.0}).arrange("Energy"), aes(ECol, "ε/√B", color = factor("Eff"))) +
    facet_wrap("Eff") +
    geom_line(aes = aes(shape = factor("Eff"))) +
    geom_point() +
    ggtitle("Comparison of signal efficiency ε / √Background") +
    ggsave(&"{outpath}/signal_efficiency_vs_sqrt_background_facet.pdf",
           width = 1280, height = 800)

import orgtables
proc printBackgroundRates(df: DataFrame, factor: float, energyMin: float, rateTable: string) =
  type Line = tuple[Classifier: string, ε_eff: float, Scinti, FADC, Septem, Line: bool, ε_total, Rate: float]
  proc toLine(df: DataFrame, backRate: float): Line =
    result = (Classifier : df["Classifier"].unique.item(string),
              ε_eff      : df["ε_eff"].unique.item(float),
              Scinti     : df["Scinti"].unique.item(bool),
              FADC       : df["FADC"].unique.item(bool),
              Septem     : df["Septem"].unique.item(bool),
              Line       : df["Line"].unique.item(bool),
              ε_total    : df["ε_total"].unique.item(float),
              Rate       : backRate)

  proc intBackRate(df: DataFrame, factor: float, energyRange: Slice[float]): seq[Line] {.discardable.} =
    result = newSeq[Line]()
    for tup, subDf in groups(df.group_by("Dataset")): # for each `name` argument separately!
      let f = tup[0][1].toStr
      let intBackRate = calcIntegratedBackgroundRate(subDf, factor, energyRange)
      let size = energyRange.b - energyRange.a
      log(true): # always print
        &"Dataset: {f}"
      log(true): # for some reason cannot have in one, shrug
        &"\t Integrated background rate in range: {energyRange}: {intBackRate:.4e} cm⁻²·s⁻¹"
      log(true):
        &"\t Integrated background rate/keV in range: {energyRange}: {intBackRate / size:.4e} keV⁻¹·cm⁻²·s⁻¹"

      result.add toLine(subDf, intBackRate / size)

  if rateTable.len == 0:
    intBackRate(df, factor, max(energyMin, 0.0) .. 12.0)
    intBackRate(df, factor, max(energyMin, 0.5) .. 2.5)
    intBackRate(df, factor, max(energyMin, 0.5) .. 5.0)
    intBackRate(df, factor, max(energyMin, 0.0) .. 2.5)
    intBackRate(df, factor, max(energyMin, 4.0) .. 8.0)
    intBackRate(df, factor, max(energyMin, 2.0) .. 8.0)
  let tab = intBackRate(df, factor, max(energyMin, 0.0) .. 8.0)
  if rateTable.len > 0:
    writeFile(rateTable, tab.toOrgTable(precision = 3))
  echo tab.toOrgTable(precision = 3)

proc main(files: seq[string], log = false, title = "",
          show2014 = false,
          separateFiles = false,
          suffix = "",
          # names are the names associated to each file. If len > 0 use that name instead of something derived from
          # filename. Makes for more natural way to separate combine things! Order needs to be same as `files`!
          names: seq[string] = @[],
          compareEfficiencies = false,
          plotMedianBools = false,
          combName = "", # name for combined data (i.e. `separateFiles == false`)
          combYear = 0, # year for combined data (i.e. `separateFiles == false`)
          totalTime = -1, # total time to use in hours # not supported with `names` or `separateFiles`!
          hidePoints = false, ## disables the `geom_point` call
          hideErrors = false, ## disables the `geom_errorbar` call
          fill = false, ## If true, will `fill` with alpha 0.5 for case of multiple datasets.
                        ## Else will color outline
          region: ChipRegion = crAll, # use either all data or cut to given region
          energyDset = "energyFromCharge",
          centerChip = 3,
          toaCutoff: seq[float] = @[],
          readToA = false,
          topMargin = 1.0,
          yMax = -1.0,
          energyMin = 0.0, energyMax = 12.0,
          useTeX = false,
          showPreliminary = false,
          showNumClusters = false,
          showTotalTime = false,
          genTikZ = false,
          filterNoisyPixels = false,
          logPlot = false,
          outpath = "plots",
          outfile = "",
          rateTable = "", # path to a file in which a table of rates will be printed
          applyEfficiencyNormalization = false,
          quiet = false,
          noPlot = false # skips the plot, good for background rate table
         ) =
  discard existsOrCreateDir(outpath)
  if readToA:
    DsetNames.add "toaLength"

  if filterNoisyPixels:
    addNoisyPixels()

  let verbose = not quiet

  let logLFiles = readFiles(files, names, region, energyDset, centerChip, energyMin, energyMax, readToA, toaCutoff, applyEfficiencyNormalization, verbose)
  let fnameSuffix = logLFiles.mapIt($it.year).join("_") & "_show2014_" & $show2014 & "_separate_" & $separateFiles & "_" & suffix.replace(" ", "_")

  let factor = if log: 1.0 else: 1e5

  if plotMedianBools:
    var df: DataFrame
    var totalTime: float
    for f in logLFiles:
      df.add f.df.clone() ## need to clone otherwise modify?!
      totalTime += f.totalTime
    plotMedianBools(computeMedianBools(df).flatScaleMedian(factor = factor,
                                                           totalTime = totalTime,
                                                           verbose = verbose),
                    fnameSuffix, title = "Comparison of background rate by different variables: clusters < and > than median",
                    suffix = suffix,
                    outpath = outpath,
                    verbose = verbose)

  if compareEfficiencies:
    plotEfficiencyComparison(logLFiles, outpath, verbose)
  else:
    var df = newDataFrame()
    if names.len > 0:
      doAssert names.len == files.len, "Need one name for each input file!"
      #for logL in logLFiles:
      df.add flatScale(logLFiles, factor, verbose)
    elif separateFiles:
      df = flatScale(logLFiles, factor, verbose)
    elif logLFiles.len > 0:
      # add up all input files
      var logL: LogLFile
      for l in logLFiles:
        log(true):
          &"Total time: {l.totalTime} of file: {l.name}"
        if l.totalTime == 0.0:
          logL = l
        else:
          logL.totalTime += l.totalTime
          logL.df.add l.df.clone()
      if totalTime.Hour > 0.0.Hour:
        logL.totalTime = totalTime.Hour.to(Second).float
      log(true):
        &"Total total time: {logL.totalTime}"
      logL.name = combName
      logL.year = combYear
      logL.df["File"] = combName
      if logL.name.len == 0:
        raise newException(ValueError, "separateFiles == true requires a `combName`!")
      if logL.year == 0:
        raise newException(ValueError, "separateFiles == true requires a `combYear`!")
      df = flatScale(@[logL], factor, verbose)
    # else no input files. Only 2014 data to plot
    #if separateFiles:
    #  for logL in logLFiles:
    #    df.add flatScale(@[logL], factor)
    #else:
    #  df = flatScale(logLFiles, factor)
    if show2014:
      if Ccol in df:
        df.drop(Ccol)
      let df2014 = read2014Data(log)
      df.add df2014

    if df.len == 0 and not show2014:
      echo "[INFO]: Neither any input files given, nor 2014 data to plot. Please provide either."
      return
    elif df.len > 0:
      printBackgroundRates(df, factor, energyMin, rateTable)

    if not noPlot:
      plotBackgroundRate(
        df.filter(f{idx(Ecol) <= energyMax}), # filter histogram bins to target energy
        fnameSuffix, title,
        outpath, outfile,
        show2014, suffix,
        hidePoints = hidePoints, hideErrors = hideErrors, fill = fill,
        useTeX = useTeX, showPreliminary = showPreliminary, genTikZ = genTikZ,
        showNumClusters = showNumClusters, showTotalTime = showTotalTime,
        topMargin = topMargin, yMax = yMax, energyMin = energyMin, energyMax = energyMax,
        logPlot = logPlot,
        applyEfficiencyNormalization = applyEfficiencyNormalization
      )

when isMainModule:
  dispatch main, help = {
    "files" : "Input files to plot backround rate from",
    "log" : "Create a log plot. TODO: merge with logPlot!",
    "title" : "Set the title of the resulting plot",
    "show2014" : "If true show the background rate of 2014/15 CAST run",
    "separateFiles" : "If set will plot background rates for each input file separatetely.",
    "suffix" : "Suffix will added to title and filename.",
    "names" : """
Names are the names associated to each file. If len > 0 use that name instead of something derived from
filename. Makes for more natural way to separate combine things! Order needs to be same as `files`!
""",
    "compareEfficiencies" : """Given multiple input files corresponding to different efficiencies
of the logL cut, creates a comparison plot.""",
    "plotMedianBools" : "Plot different datasets that are smaller than the median of that dataset.",
    "combName" : "Name for combined data (i.e. `separateFiles == false`).",
    "combYear" : "Year for combined data (i.e. `separateFiles == false`).",
    "totalTime" : "Total time to use in hours, not supported with `names` or `separateFiles`.",
    "hidePoints" : "If set do not show points of data (only histogram).",
    "hideErrors" : "If set disables error bars.",
    "fill" : "If true, will `fill` with alpha 0.5 for case of multiple datasets. Else will color outline.",
    "region" : "The chip region to cut to.",
    "energyDset" : "The energy dataset to base the x axis on, {energyFromCharge, energyFromPixel}.",
    "centerChip" : "The center chip of the detector stored in the files. 0 for single chip detectors.",
    "toaCutoff" : "For detectors with ToA data, cut off all clusters with larger ToA lengths than this.",
    "readToA" : "If set will also read ToA related datasets. Required for `toaCutoff.",
    "topMargin" : "Margin at the top of the plot for long titles.",
    "yMax" : "If any given, limit the y axis to this value.",
    "energyMax" : "If any given, limit the energy & x axis to this maximum value.",
    "energyMin" : "If any given, limit the energy to this minimum value (x axis remains unchanged).",
    "useTeX" : "Generate a plot using TeX",
    "showPreliminary" : "If set shows a big 'Preliminary' message in the center.",
    "showNumClusters" : "If set adds number of input clusters to title.",
    "showTotalTime" : "If set adds the total time of background data to title.",
    "genTikZ" : "If set only generate TikZ instead of compiling a TeX file to PDF.",
    "filterNoisyPixels" : "If set removes the (currently hardcoded) noisy pixels.",
    "logPlot" : "Alternate setting to activate log10 plot. TO BE REMOVED",
    "outpath" : "The path in which the plot is placed, by default './plots'",
    "outfile" : "The name of the generated plot, by default constructed from other parameters"}
