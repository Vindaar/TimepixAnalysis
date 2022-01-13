import ggplotnim, seqmath, sequtils, os, sugar, strscans, strformat, strutils, sugar
import ingrid / [tos_helpers]
from arraymancer import tensor

import nimhdf5, numericalnim

import cligen

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
echo DsetNames

proc scaleDset(data: Column, totalTime, factor: float): Column =
  ## scales the data in `data` according to the area of the gold region,
  ## total time and bin width. The data has to be pre binned of course!
  # XXX: make sure we never use wrong area if input data makes use of `--chipRegion` feature!
  let area = pow(0.95 - 0.45, 2) # area of gold region!
  const bin_width = 0.2 # 0.392
  const shutter_open = 1.0 # Does not play any role, because totalDuration takes
                           # this into account already!
  let scale = factor / (totalTime * shutter_open * area * bin_width) #* 1e5
  result = toColumn data.toTensor(float).map_inline(x * scale)

proc readTime(h5f: H5FileObj): float =
  let lhGrp = h5f["/likelihood".grp_str]
  result = lhGrp.attrs["totalDuration", float]

proc scaleDset(h5f: H5FileObj, data: Column, factor: float): Column =
  result = scaleDset(data, h5f.readTime(), factor)

proc extractYear(f: string): int =
  ## We assume the filenames
  ## - always contain a year
  ## - always separated by underscores
  let fname = f.extractFilename
  var dummy: string
  if fname.scanf("$*_$i_$*", dummy, result, dummy):
    echo "found a year ", result

proc readFiles(files: seq[string], names: seq[string], region: ChipRegion): seq[LogLFile] =
  ## reads all H5 files given and stacks them to a single
  ## DF. An additional column is added, which corresponds to
  ## the filename. That way one can later differentiate which
  ## events belong to what and decide if data is supposed to
  ## be accumulated or not
  ##
  ## After reading we will cut to the given `region`.
  doAssert names.len == 0 or names.len == files.len, "Need one name for each input file!"
  for idx, file in files:
    let h5f = H5open(file, "r")
    var df = h5f.readDsets(likelihoodBase(), some((3, DsetNames)))
      .rename(f{Ecol <- "energyFromCharge"})
    let fname = if names.len > 0: names[idx]
                else: file.extractFilename
    df["File"] = constantColumn(fname, df.len)
    echo "[INFO]: File: ", fname
    echo "[INFO]: Elements before cut: ", df.len
    df = df.filter(f{float: inRegion(`centerX`, `centerY`, region) == true})
    echo "[INFO]: Elements after cut: ", df.len
    result.add LogLFile(name: fname,
                        totalTime: readTime(h5f),
                        df: df,
                        year: extractYear(fname))
    discard h5f.close()

proc histogram(df: DataFrame): DataFrame =
  ## Calculates the histogam of the energy data in the `df` and returns
  ## a histogram of the binned data
  ## TODO: allow to do this by combining different `File` values
  let (hist, bins) = histogram(df[Ecol].toTensor(float).toRawSeq,
                               range = (0.0, 20.0), bins = 100)
  result = seqsToDf({ Ecol : bins, Ccol : concat(hist, @[0]) })

template sumIt(s: seq[typed], body: untyped): untyped =
  var res: float
  for it {.inject.} in s:
    res += body
  res

proc flatScale(files: seq[LogLFile], factor: float, dropCounts = true): DataFrame =
  #var count {.global.} = 0
  var df: DataFrame
  for f in files:
    if not f.df.isNil:
      df.add f.df
  result = newDataFrame()

  for tup, subDf in groups(df.group_by("File")):
  #for tup, subDf in groups(df.group_by("L<L_median")):
    var dfLoc = newDataFrame()
    dfLoc = subDf.histogram()
    dfLoc = dfLoc.mutate(f{float: "CountErr" ~ sqrt(`Counts`)})
    dfLoc[Rcol] = dfLoc[Ccol].scaleDset(files.sumIt(it.totalTime), factor)
    dfLoc["RateErr"] = dfLoc["CountErr"].scaleDset(files.sumIt(it.totalTime), factor)
    dfLoc["Dataset"] = constantColumn(tup[0][1].toStr, dfLoc.len) #"2017/18_" & $count, dfLoc.len)
    dfLoc = dfLoc.mutate(f{"yMin" ~ `Rate` - `RateErr`}, f{"yMax" ~ `Rate` + `RateErr`})
    dfLoc["yMin"] = dfLoc["yMin", float].map_inline:
      if x >= 0.0: x
      else: 0.0
    if dropCounts:
      dfLoc.drop("Counts")
    result.add dfLoc
  echo result.pretty(-1)

proc flatScaleMedian(df: DataFrame, factor: float, totalTime: float): DataFrame =
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
  echo result.pretty(-1)

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
  result = trapz(rate.toRawSeq, energies.toRawSeq) / factor

proc computeMedianBools(df: DataFrame): DataFrame =
  var df = df
  var medianNames = newSeq[string]()
  proc computeMedianBool(dset: string) =
    let medianVal = df[dset, float].toRawSeq.median(50)
    let nameStr = $dset & "<" & $dset & "_median"
    medianNames.add nameStr
    df = df.mutate(f{float -> bool: nameStr ~ df[dset][idx] < medianVal})
  for dset in DsetNames:
    if dset != "energyFromCharge":
      computeMedianBool(dset)

  result = df.gather(medianNames, key = "Variable", value = "Value")

proc plotMedianBools(df: DataFrame, fnameSuffix, title: string,
                     suffix: string) =
  let transparent = color(0, 0, 0, 0)
  let fname = &"plots/background_rate_median_bools_{fnameSuffix}.pdf"
  echo "INFO: storing plot in ", fname
  #let df = df.select([Ecol, Rcol, "Variable", "Value"])
  echo df.pretty(-1)

  ggplot(df, aes(Ecol, Rcol, fill = "Value")) +
    facet_wrap("Variable", scales = "free") +
    geom_histogram(stat = "identity", position = "identity", alpha = some(0.5),
                   color = some(transparent),
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
                        show2014: bool, suffix: string,
                        useTeX: bool) =
  var df = df.filter(f{c"Energy" < 12.0})
  let titleSuff = if title.len > 0: title
                  elif show2014: " compared to 2014/15"
                  else: ""
  let transparent = color(0, 0, 0, 0)
  let fname = &"plots/background_rate_{fnameSuffix}.pdf"
  echo "INFO: storing plot in ", fname
  echo df

  ggplot(df, aes(Ecol, Rcol, fill = "Dataset")) +
    geom_histogram(stat = "identity", position = "identity", alpha = some(0.5),
                   color = some(transparent),
                   hdKind = hdOutline) +
    geom_point(binPosition = "center", position = "identity") +
    geom_errorbar(binPosition = "center",
                  aes = aes(yMin = "yMin", yMax = "yMax"),
                  errorBarKind = ebLines) +
    scale_y_continuous() +
    xlim(0, 12.0)
  if useTeX:
    plt + annotate("GridPix preliminary",
             x = 4.3, y = 2.0,
             rotate = -30.0,
             font = font(16.0, color = color(0.92, 0.92, 0.92)),
             backgroundColor = transparent) +
    xlab(r"Energy [\si{keV}]") +
    ylab(r"Rate [\SI{1e-5}{keV⁻¹ cm⁻² s⁻¹}]", margin = 1.1) +
    #minorGridLines() +
    ggtitle(title) +
    theme_latex() +
    ggsave(fname, width = 800, height = 480, useTex = true, standalone = true)
  else:
    plt + xlab("Energy [keV]") +
    ylab("Rate [10⁻⁵ keV⁻¹ cm⁻² s⁻¹]") +
    ggtitle(&"Background rate of Run 2 & 3 (2017/18){titleSuff}{suffix}") +
    ggsave(fname, width = 800, height = 480)

  df = df.drop(["Dataset", "yMin", "yMax"])
  df.writeCsv("/tmp/background_rate_data.csv")

proc plotEfficiencyComparison(files: seq[LogLFile]) =
  ## creates a plot comparing different logL cut efficiencies. Each filename should contain the
  ## efficiency in the format
  ## `_eff_<eff>.h5`
  ## - for each file:
  ##   - for each bin:
  ##     - eff / sqrt(counts)
  ## plot geom_line + geom_point of these values
  var df = newDataFrame()
  for f in files:
    doAssert "_eff_" in f.name
    # ugly for now
    let eff = f.name.split("_")[^1].dup(removeSuffix(".h5")).parseFloat
    var dfFlat = flatScale(@[f], 1e5, dropCounts = false)
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
    ggsave("plots/signal_efficiency_vs_sqrt_background.pdf",
           width = 1280, height = 800)
  ggplot(df.filter(f{`Energy` < 12.0}).arrange("Energy"), aes(ECol, "ε/√B", color = factor("Eff"))) +
    facet_wrap("Eff") +
    geom_line(aes = aes(shape = factor("Eff"))) +
    geom_point() +
    ggtitle("Comparison of signal efficiency ε / √Background") +
    ggsave("plots/signal_efficiency_vs_sqrt_background_facet.pdf",
           width = 1280, height = 800)


proc main(files: seq[string], log = false, title = "", show2014 = false,
          separateFiles = false,
          suffix = "",
          # names are the names associated to each file. If len > 0 use that name instead of something derived from
          # filename. Makes for more natural way to separate combine things! Order needs to be same as `files`!
          names: seq[string] = @[],
          compareEfficiencies = false,
          plotMedianBools = false,
          region: ChipRegion = crAll, # use either all data or cut to given region
          useTeX = false
         ) =
  discard existsOrCreateDir("plots")
  let logLFiles = readFiles(files, names, region)
  let fnameSuffix = logLFiles.mapIt($it.year).join("_") & "_show2014_" & $show2014 & "_separate_" & $separateFiles & suffix

  let factor = if log: 1.0 else: 1e5

  if plotMedianBools:
    var df: DataFrame
    var totalTime: float
    for f in logLFiles:
      df.add f.df.clone() ## need to clone otherwise modify?!
      totalTime += f.totalTime
    plotMedianBools(computeMedianBools(df).flatScaleMedian(factor = factor,
                                                           totalTime = totalTime),
                    fnameSuffix, title = "Comparison of background rate by different variables: clusters < and > than median",
                    suffix = suffix)

  if compareEfficiencies:
    plotEfficiencyComparison(logLFiles)
  else:
    var df = newDataFrame()
    if names.len > 0:
      doAssert names.len == files.len, "Need one name for each input file!"
      #for logL in logLFiles:
      df.add flatScale(logLFiles, factor)
    #if separateFiles:
    #  for logL in logLFiles:
    #    df.add flatScale(@[logL], factor)
    else:
      df = flatScale(logLFiles, factor)
    echo df

    ## NOTE: this has to be calculated before we add 2014 data if we do, of course,
    ## because otherwise we have everything duplicated!
    let intBackRate = calcIntegratedBackgroundRate(df, factor)
    echo &"Integrated background rate: {intBackRate:.4e} cm⁻² s⁻¹"
    echo &"Integrated background rate/keV: {intBackRate / 10.0:.4e} keV⁻¹ cm⁻² s⁻¹"
    if show2014:
      df.drop(Ccol)
      var df2014 = toDf(readCsv(Data2014, sep = ' ', header = "#"))
        .rename(f{Rcol <- "Rate[/keV/cm²/s]"}, f{"yMax" <- "dRateUp"},
                f{"yMin" <- "dRateDown"}, f{Ecol <- "E[keV]"})
      let lastLine = seqsToDf({ Ecol : @[10.1], Rcol : @[0.0], "yMin" : @[0.0], "yMax" : @[0.0] })
      df2014.add lastLine
      if not log:
        df2014 = df2014.mutate(f{Rcol ~ 1e5 * `Rate`}, f{"yMin" ~ 1e5 * `yMin`},
                               f{"yMax" ~ 1e5 * `yMax`})
      df2014 = df2014.mutate(f{"yMin" ~ `Rate` - `yMin`}, f{"yMax" ~ `Rate` + `yMax`},
                             f{Ecol ~ `Energy` - 0.1})
      df2014["Dataset"] = constantColumn("2014/15", df2014.len)
      echo df2014
      df.add df2014
    plotBackgroundRate(df, fnameSuffix, title, show2014, suffix, useTeX = useTeX)

when isMainModule:
  dispatch main
