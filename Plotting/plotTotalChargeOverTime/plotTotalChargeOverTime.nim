import nimhdf5, ggplotnim, os, strformat, strutils, sequtils, algorithm, times, os, stats
import ingrid / [tos_helpers, ingrid_types, calibration]

import cligen

## See status and progress section
## ``Calculate total charge of background data over time``
## Just calculates the total charge binned over a specific time interval
## and plots it

type
  SupportedRead = SomeFloat | SomeInteger | string | bool | Value
  ReadProc = proc(h5f: H5File, applyRegionCut: bool, chip: int): DataFrame

let RotAngle = getEnv("ROT_ANGLE", "-20.0").parseFloat
let Width = getEnv("WIDTH", "1200").parseFloat
let Height = getEnv("HEIGHT", "800").parseFloat
let FacetMargin = getEnv("FACET_MARGIN", "0.5").parseFloat
let Top = getEnv("TOP", "1.5").parseFloat
let Left = getEnv("LEFT", "2.0").parseFloat
let Bottom = getEnv("BOTTOM", "1.0").parseFloat
let Right = getEnv("RIGHT", "2.0").parseFloat
let LegendLeft = getEnv("LEGEND_LEFT", "0.5").parseFloat
let LegendBottom = getEnv("LEGEND_BOTTOM", "0.175").parseFloat
let UseTex = getEnv("USE_TEX", "false").parseBool

const CommonDsets = { "energyFromCharge" : "energy",
                      "fractionInTransverseRms" : "fracRmsTrans",
                      "lengthDivRmsTrans" : "L_div_RMS_trans",
                      "eccentricity" : "eccentricity" }

let readKeys = concat(CommonDsets.mapIt(it[0]), @["timestamp"])
let addnKeys = ["Mean", "Median", "Variance", "Skewness", "Kurtosis"]
var allKeys: seq[string]
for k in readKeys:
  if k != "timestamp":
    for adn in addnKeys:
      allKeys.add k & adn
  else:
    allKeys.add k

proc genNames(): seq[string] {.compileTime.} =
  result = @["totalCharge", "eventNumber", "hits", "centerX", "centerY", "rmsTransverse"]
  for (d, v) in CommonDsets:
    result.add d

proc toPeriod(v: float): string =
  result = v.int.fromUnix.format("dd/MM/YYYY")

proc fromPeriod(s: string): float =
  result = s.parseTime("dd/MM/YYYY", utc()).toUnixFloat

proc uniquePeriods(df: DataFrame): seq[string] =
  ## Returns all run periods sorted in increasing order as strings
  result = df["runPeriods"]
    .unique.toTensor(string).toSeq1D
    .mapIt(fromPeriod(it)).sorted(SortOrder.Ascending)
    .mapIt(toPeriod(it))

proc closestIdx(t: float, ts: seq[float]): int =
  ## returns the index of the time from `ts` which is closest to `t`
  result = ts.mapIt(abs(it - t)).minIndex

proc splitDf(df: DataFrame,
             periods: seq[string] = @[]): DataFrame =
  var periodSeq = newSeq[string](df.len)
  let times = df["timestamp"].toTensor(int)
  var lastTimestamp = times[0]
  var period = if periods.len == 0: toPeriod(lastTimestamp.float) else: periods[0]
  if periods.len == 0:
    for i in 0 ..< times.size:
      let t = times[i]
      if t - lastTimestamp > 10 * 86400:
        period = toPeriod t.float
      periodSeq[i] = period
      lastTimestamp = t
  else:
    var idx = 0
    var nextPeriodTime = fromPeriod periods[idx+1]
    let periodTimes = periods.mapIt(it.fromPeriod)
    for i in 0 ..< times.size:
      let t = times[i]
      periodSeq[i] = periods[t.float.closestIdx(periodTimes)]
  result = df
  result["RunPeriod"] = periodSeq

proc readTstampDf(h5f: H5File, applyRegionCut: bool,
                  chip = 3): DataFrame =
  const evNames = genNames()
  const allNames = ["timestamp", "eventNumber"]

  template readIt(names, baseName, df: untyped): untyped =
    for name in names:
      let dsetName = baseName / name
      let dsetH5 = h5f[dsetName.dset_str]
      withDset(dsetH5):
        when type(dset) is seq[SupportedRead]:
          df[name] = dset

  for run, grp in runs(h5f, recoBase()):
    var dfEv = newDataFrame()
    var dfAll = newDataFrame()
    let group = h5f[grp.grp_str]
    readIt(evNames, grp / &"chip_{chip}", dfEv)
    if applyRegionCut:
      # if we do not apply the region cut the slicing does not make any sense, because
      # the indices for the slices will be wrong!
      dfEv = dfEv.applyGasGainCut

      ## For event, read gas gain slices
      var sliceNum = newSeq[int](dfEv.len)
      # Note: we only read the slice data to match the input *binning*. If the data
      # was *not* binned (full run gas gain), but we want to plot against 90 min bins,
      # using the gas gain time slice 90 min bins is _still a good idea_ as it gives us
      # the same data points as in the 90 min case, just with energies that are worse
      # calibrated.
      ## XXX: do not try to read `90`, but rather try to read the version of the input time!
      let gsname = if grp / &"chip_{chip}" / "gasGainSlices90" in h5f: "gasGainSlices90"
                   else: "gasGainSlices"
      for (idx, slice) in iterGainSlicesIdx(h5f, grp / &"chip_{chip}", gsname):
        sliceNum[slice] = repeat(idx, slice.len)
      dfEv["sliceNum"] = sliceNum

    readIt(allNames, grp, dfAll)
    var dfJoined = innerJoin(dfEv, dfAll, "eventNumber")
    dfJoined["runNumber"] = constantColumn(run, dfJoined.len)

    result.add dfJoined

proc readPhotoEscapeDf(h5f: H5File, applyRegionCut: bool, chip: int): DataFrame =
  ## These are the x positions (in charge) of the spectrum
  #const kalphaCharge = 4
  #const escapeCharge = 1
  const kalphaCharge = 3
  const escapeCharge = 0
  const parPrefix = "p"
  const dateStr = "yyyy-MM-dd'.'HH:mm:ss" # example: 2017-12-04.13:39:45
  var
    photoSeq: seq[float]
    escapeSeq: seq[float]
    dates: seq[float]
  for run, grp in runs(h5f, recoBase()):
    let group = h5f[(recoBase & $run).grp_str]
    let chpGrpName = group.name / &"chip_{chip}"
    let dset = h5f[(chpGrpName / "FeSpectrumCharge").dset_str]
    photoSeq.add dset.attrs[parPrefix & $kalphaCharge, float]
    escapeSeq.add dset.attrs[parPrefix & $escapeCharge, float]
    dates.add parseTime(group.attrs["dateTime", string],
                        dateStr,
                        utc()).toUnix.float
  result = toDf({ "timestamp" : dates, "photo" : photoSeq,
                      "escape" : escapeSeq })

proc readCdl(): DataFrame =
  const path = "/mnt/1TB/CAST/CDL_2019/calibration-cdl-2018.h5"
  result = newDataFrame()
  if fileExists(path):
    const group = "/calibration-cdl-feb2019-Mn-Cr-12kV"
    const dsets = ["CdlSpectrum", "CdlSpectrumEvents", "CdlSpectrumCharge"]
    var h5f = H5open(path, "r")
    for name in dsets:
      let dsetName = group / name
      let dsetH5 = h5f[dsetName.dset_str]
      withDset(dsetH5):
        when type(dset) is seq[SupportedRead]:
          result[name] = dset
    result["runNumber"] = constantColumn(0, result.len)
    result["runType"] = constantColumn("CDL", result.len)
    result = result.rename(f{"totalCharge" <- "CdlSpectrumCharge"},
                           f{"hits" <- "CdlSpectrum"},
                           f{"eventNumber" <- "CdlSpectrumEvents"})
    discard h5f.close()

proc readFiles(files: seq[string], runType: string,
               applyRegionCut: bool,
               readAllChips: bool,
               readProc: ReadProc = readTstampDf): DataFrame =
  ## reads all H5 files given and stacks them to a single
  ## DF. An additional column is added, which corresponds to
  ## the filename. That way one can later differentiate which
  ## events belong to what and decide if data is supposed to
  ## be accumulated or not
  result = newDataFrame()
  for file in files:
    let h5f = H5open(file, "r")
    if readAllChips:
      for chip in 0 ..< 6: ## XXX: assert correct number of chips!
        var df = readProc(h5f, applyRegionCut = applyRegionCut, chip = chip)
        df["chip"] = chip
        result.add df
    else:
      result.add readProc(h5f, applyRegionCut = applyRegionCut, chip = 3)
    discard h5f.close()
  if result.len > 0:
    result = result.arrange("timestamp")
    result["runType"] = constantColumn(runType, result.len)

template len[T](t: Tensor[T]): int = t.size.int

proc mapToRunPeriods[T](tmeanStamps: T,
                        periodTab: OrderedTable[int, string]): seq[string] =
  const splitPeriod = 86400 * 5 # 5 days?

  let periods = toSeq(keys(periodTab))
  var lastPeriod = if periods.len == 0: toPeriod tmeanStamps[0]
                   else: periodTab[periods[0]]
  result.add lastPeriod
  ## walk data again to split by splitPeriods
  for i in 1 ..< tmeanStamps.len:
    if tmeanStamps[i] - tmeanStamps[i - 1] >= splitPeriod:
      if periods.len == 0:
        lastPeriod = toPeriod tmeanStamps[i]
      else:
        # find the largest run period smaller in time than curTstamp
        let diffs = periods.mapIt(abs(it - tmeanStamps[i].int))
        let idx = diffs.argMin
        if diffs[idx] <= splitPeriod * 10:
          lastPeriod = periodTab[periods[idx]]
        else:
          lastPeriod = toPeriod tmeanStamps[i]
    result.add lastPeriod

proc calcMeanTimestamps(t: Tensor[float], interval: float): seq[float] =
  let numIntervals = ((t[t.size - 1] - t[0]) / interval).ceil.int
  result = newSeqOfCap[float](numIntervals * 2)
  var tStart = t[0]
  for i in 0 ..< t.size:
    if (t[i] - tStart) >= interval:
      # run over bin range
      result.add (t[i] + tStart) / 2.0
      tStart = t[i]

proc calcMean(df: DataFrame,
              tstamps: Tensor[float],
              outLen: int,
              n: string,
              interval: float,
              normBy: static string = "",
              toCount: static bool = false,
              useMedian: static bool = false): seq[float] =
  ## If `toCount` is used it only counts the number of entries within a time interval.
  ## If `normBy` uses the given dataset to normalize by (normalization by the
  ## *sum* of that dataset in the interval)
  ## If `useMedian` computes the median instead of the mean value.
  when not toCount: ## do not need it when we count entries
    let data = df[n].toTensor(float)
  const toNorm = normBy.len > 0
  when toNorm:
    let dataNorm = df[normBy].toTensor(float)
  result = newSeq[float](outLen)

  var
    current = 0.0
    norm = 0.0
    tStart = tstamps[0]
    numCluster = 0
    j = 0
    tStartIdx = 0
    length = 0.0
  var timesToPlot = newSeq[float]()
  when useMedian:
    var vec = newSeq[float]()
  for i in 0 ..< tstamps.size:
    let inInterval = (tstamps[i] - tStart) < interval
    if not inInterval:
      #if numCluster < 400:
      #  echo "Jump from last to this: ", tstamps[i-1], " to ", tstamps[i], " = ", tstamps[i] - tstamps[i-1], " and numCluster= ", numCluster
      #  echo "total time difference was: ", tstamps[i-1] - tStart, "\n"
      # run over bin range
      when toNorm:
        result[j] = current / norm
        norm = 0.0
      elif useMedian:
        result[j] = vec.median()
        vec = newSeq[float]()
      elif toCount:
        if length > 3600: ## Want to exclude intervals that are too short!
          result[j] = numCluster.float / length # normalize number of counts by length of this interval
      else:
        result[j] = current / numCluster.float
      numCluster = 0
      current = 0.0
      tStart = tstamps[i]
      tStartIdx = i
      timesToPlot.add tStart
      inc j
    when toNorm:
      norm += dataNorm[i]
    when useMedian:
      vec.add data[i]
    # distinct from `when` above! else here needed for both above.
    when not toCount:
      current += data[i]
    length = tstamps[i] - tStart
    inc numCluster
  let tt = timesToPlot.mapIt(it.int)
  let tstr = tt.mapIt($it.fromUnix)
  let df = toDf({"time" : tt, "tstr" : tstr, "idx" : toSeq(0 ..< timesToPlot.len)})
  ggplot(df, aes(idx, time)) + geom_point() + ggsave("/tmp/testtime.pdf")
  df.write_csv("/tmp/testtime.csv")

proc computeStats(df: DataFrame): DataFrame =
  result = newDataFrame()
  var res = initTable[string, seq[float]]() # reduced values from this, not a DF to grow better
  for k in allKeys:
    res[k] = newSeqOfCap[float](1000)
  var periods = newSeqOfCap[string](1000)
  var timesToPlot = newSeq[float]()

  for (tup, subDf) in groups(df.group_by(["runNumber", "sliceNum"])):
    ## NOTE: if we include slices with less than 100 events, we include some where the
    ## energy calibration is rather bad!
    #if subDf.len < 100: continue
    for k in readKeys:
      let data = subDf[k, float].toRawSeq
      if k != "timestamp":
        var stat: RunningStat
        stat.push(data)
        res[k & "Mean"].add stat.mean()
        res[k & "Median"].add data.percentile(50)
        res[k & "Variance"].add stat.variance()
        res[k & "Skewness"].add stat.skewness()
        res[k & "Kurtosis"].add stat.kurtosis()
      else:
        res["timestamp"].add((data[^1] + data[0]).float / 2.0)
    periods.add subDf["RunPeriod", string][0]
  for k, v in res:
    result[k] = toColumn v
  result["runType"] = constantColumn(df["runType", string][0], result.len)
  result["runPeriods"] = periods
  result = result.arrange("timestamp")

  let df2 = result.mutate(f{float -> int: "timestamp" ~ `timestamp`.int})
  var num {.global.} = 0
  df2.writeCsv(&"/tmp/nov2017_{num}.csv")
  inc num

proc calculateMeanDf(df: DataFrame, interval: float,
                     periodTab: OrderedTable[int, string] = initOrderedTable[int, string](),
                     meanColumn = "meanCharge",
                     useMedian = false,
                     toCount = false): DataFrame =
  ## interval: time in minutes
  # we're going to do a rather crude approach. Given that we have to bin
  # a continuous variable, we're not well equipped using ggplotnim's dataframe
  # I fear. So just do it manually...
  let tstamps = df["timestamp", float]
  let tmeanstamps = calcMeanTimestamps(tstamps, interval)
  let runPeriods = mapToRunPeriods(tmeanStamps, periodTab)

  let outLen = tmeanstamps.len
  let sums = df.calcMean(tstamps, outLen, "hits", interval)
  var means: seq[float]
  if useMedian:
    means = df.calcMean(tstamps, outLen, "totalCharge", interval, useMedian = true)
  elif toCount:
    means = df.calcMean(tstamps, outLen, "", interval, toCount = true)
  else:
    means = df.calcMean(tstamps, outLen, "totalCharge", interval, normBy = "hits")

  result = toDf({ "timestamp" : tmeanStamps, "sumCharge" : sums,
                  meanColumn : means, "runPeriods" : runPeriods })

  template calcAndAdd(name, dfName: untyped): untyped =
    let tmp = df.calcMean(tstamps, outLen, name, interval, useMedian = true)
    result[dfName] = tmp
  for (k, v) in CommonDsets:
    calcAndAdd(k, v)

  let filterPeriods = result.group_by("runPeriods").count("runPeriods")
    .filter(f{`n` == 1})["runPeriods"].toTensor(string).toRawSeq
  echo filterPeriods
  result = result.filter(f{string -> bool: `runPeriods` notin filterPeriods})
  result["runType"] = constantColumn(df["runType", 0].toStr, result.len)
  if "chip" in df:
    echo df
    let chip = df["chip"].unique.item(int) ## The input should only have a single chip number!
    result["chip"] = chip
  echo "Number of run periods with more than 1 entry: ", result["runPeriods"].unique

  let df2 = result.mutate(f{float -> int: "timestamp" ~ `timestamp`.int})
  var num2 {.global.} = 0
  df2.writeCsv(&"/tmp/nov2017_mean_{num2}.csv")
  inc num2

proc commonPlotFields(plt: GgPlot, periods: seq[Value]): GgPlot =
  ## Adds all those geoms / scales etc. that are common between multiple plots in this script.
  result = plt +
    facet_wrap("runPeriods", scales = "free", order = periods) +
    geom_point(alpha = some(0.75)) +
    #scale_x_continuous(labels = toPeriod) +
    scale_x_date(name = "Date", isTimestamp = true,
                 dateSpacing = initDuration(weeks = 2),
                 formatString = "dd/MM/YYYY", dateAlgo = dtaAddDuration) +
                        #dateSpacing = initDuration(weeks = 26), dateAlgo = dtaAddDuration) +
    facetMargin(FacetMargin, ukCentimeter) +
    margin(top = Top, bottom = Bottom, right = Right, left = Left) +
    legendPosition(LegendLeft, LegendBottom) +
    xlab("Date", rotate = RotAngle, alignTo = "right", margin = 0.0)
  if UseTex:
    proc th(): Theme =
      result = singlePlot()
      result.tickLabelFont = some(font(7.0))
    result = result +  #theme_scale(1.2) +
      ylab("Median value", margin = 3.0) +
      themeLatex(fWidth = 1.0, width = 1200, height = 800, baseTheme = th)
      #facetHeaderText(font = font(14.0)) # readd otherwise overwritten

proc plotOverTime(df: DataFrame, interval: float, titleSuff: string,
                  useLog = true,
                  applyRegionCut = false,
                  useMedian = false,
                  toCount = false,
                  normalizeMedian = true,
                  colorBy = "runType",
                  outpath = "out",
                  cutoff = 0.0,
                  ylabel = "", title = "") =
  let interval = interval / 60.0 # convert back to minutes for title
  createDir(outpath)
  var nameSuff = if titleSuff == "all data": "all" else: "filtered"
  if applyRegionCut:
    nameSuff.add "_crSilver"
  ## XXX: adjust prefixes!
  let meanPrefix = if useMedian: "Median" else: "Mean"
  let normSuffix = if useMedian and normalizeMedian: " each run type normalized to 1"
                   else: ""
  let title = if title.len > 0: title
              else: &"{meanPrefix} total charge within {interval:.1f} min, {normSuffix}{titleSuff}"

  let periods = uniquePeriods(df).mapIt(%~ ("runPeriods", it))

  var pltSum = ggplot(df, aes("timestamp", "sumCharge", color = factor(colorBy)))
    .commonPlotFields(periods) +
    ggtitle(&"Sum of total charge within {interval:.1f} min{titleSuff}")
  let yScale = if useMedian and normalizeMedian:
                 f{float: `meanCharge` / max(col("meanCharge"))}
               else:
                 f{"meanCharge"}
  let ylabel = if ylabel.len > 0: ylabel
               elif useMedian: "Median value"
               else: "Mean value"

  var df = df
  if toCount:
    df = df.filter(f{float: `meanCharge` >= cutoff})

  var pltMean = ggplot(df, aes("timestamp",
                               yScale,
                               color = factor(colorBy)))
    .commonPlotFields(periods) +
    ylab(ylabel) +
    ggtitle(title)
  if useLog:
    pltSum = pltSum + scale_y_log10()
    pltMean = pltMean + scale_y_log10()
  pltSum + ggsave(&"{outpath}/background_sum_charge_binned_{interval:.1f}_min_{nameSuff}.pdf",
                  width = Width, height = Height,
                  useTeX = UseTex, standalone = UseTex)
  pltMean + ggsave(&"{outpath}/background_{meanPrefix.toLowerAscii()}_charge_binned_{interval:.1f}_min_{nameSuff}.pdf",
                   width = Width, height = Height,
                   useTeX = UseTex, standalone = UseTex)



proc plotHistos(df: DataFrame, interval: float, titleSuff: string,
                useLog = true,
                applyRegionCut = false,
                outpath = "out") =
  let interval = interval / 60.0 # convert back to minutes for title
  createDir(outpath)
  var nameSuff = if titleSuff == "all data": "all" else: "filtered"
  if applyRegionCut:
    nameSuff.add "_crSilver"

  # Get unique periods and turn into form needed for ggplotnim
  let periods = uniquePeriods(df).mapIt(%~ ("runPeriods", it))

  for key in readKeys:
    if key == "timestamp": continue
  #for (k, name) in CommonDsets:
    for adn in addnKeys:
      let name = key & adn
      #let key = k
      #let adn = "median"
      #let df = df.filter(f{float -> bool: `lengthDivRmsTrans` < 1e13 and
      #                     `eccentricity` < 1e13})
      #  .filter(f{string -> bool: `runPeriods` == "30/10/2017"})
      #  .mutate(f{float -> float: "energy" ~ `energy` * 1e6})
      #let df = df.mutate(f{float -> float: "energy" ~ `energy` * 1e6})
      var pltTmp = ggplot(df, aes("timestamp", name, color = "runType"))
        .commonPlotFields(periods) +
        ylim(2.0, 6.5) +
        ggtitle(&"{adn} of cluster {key} within {interval:.1f} min{titleSuff}")
      var pltTmpHisto = ggplot(df, aes(name, fill = "runType")) +
        facet_wrap("runPeriods", scales = "free") +
        geom_histogram(bins = 300, density = true, position = "identity") +
        ggtitle(&"Histogram of {adn} cluster {key} within {interval:.1f} min{titleSuff}")
      if useLog:
        pltTmp = pltTmp + scale_y_log10()
        pltTmpHisto = pltTmpHisto + scale_y_log10()
      pltTmp + ggsave(&"{outpath}/background_{adn.normalize}_{key}_{interval:.1f}_min_{nameSuff}.pdf",
                      width = Width, height = Height,
                      useTeX = UseTex, standalone = UseTex)
                       #width = 800, height = 480)
      pltTmpHisto + ggsave(&"{outpath}/background_histogram_{adn.normalize}_{key}_{interval:.1f}_min_{nameSuff}.pdf",
                           width = Width, height = Height,
                           useTeX = UseTex, standalone = UseTex)
      echo "created ", adn.normalize, " ", name
                            #width = 800, height = 480)

  # write the DF to the output path (convert timestamp to int to get full accuracy w/ default precision)
  df.mutate(f{float -> int: "timestamp" ~ `timestamp`.int}).writeCsv(&"{outpath}/data_{interval:.1f}_min_{nameSuff}.csv")

  when false:
    ggplot(df, aes("meanCharge", color = "runType")) +
      geom_histogram(bins = 100, position = "identity", alpha = some(0.5)) +
      margin(top = 2) +
      ggtitle(&"Histogram of binned mean charge, within {interval:.1f} min, {titleSuff}") +
      ggsave(&"{outpath}/background_histo_mean_binned_{interval:.1f}_min_{nameSuff}.pdf",
              width = 800, height = 480)

proc plotPhotoDivEscape(df, dfTime: DataFrame, periods: OrderedTable[int, string],
                        interval: float,
                        outpath = "out",
                        titleSuff = "") =
  ## dfTime contains the whole time series data, df is just the data containing
  ## Fe spectra photo + escape peaks
  var df = df
  let times = toSeq(periods.keys()).sorted
  df["runPeriods"] = df["timestamp"].toTensor(float).mapToRunPeriods(periods)
  df = df.mutate(f{"photoEsc" ~ `photo` / `escape`})
    .group_by("runPeriods")
    .mutate(f{"val" ~ `photoEsc` / max(`photoEsc`)})
  let maxPhotoEsc = df["photoEsc"].toTensor(float).max
  let maxEnergy = dfTime["energy"].toTensor(float).max
  var dfTime = dfTime.group_by("runPeriods")
    .mutate(f{"val" ~ `energy` / max(`energy`)})
  df = df.select(["timestamp", "val", "runPeriods", "runType"])
  dfTime = dfTime.select(["timestamp", "val", "runPeriods", "runType"])
  let dfPlot = bind_rows([("", df), ("", dfTime)])
  let periods = uniquePeriods(df).mapIt(%~ ("runPeriods", it))
  ggplot(dfPlot, aes("timestamp", "val", color = "runType"))
    .commonPlotFields(periods) +
    ggtitle("Normalized cmp of photo/escape peak in charge & median energy " &
            &"{interval/60.0:.1f} min, max photo/esc {maxPhotoEsc:.2f}, max " &
            &"median energy {maxEnergy:.2f} {titleSuff}") +
    ggsave(&"{outpath}/photo_div_escape_vs_time.pdf", width = Width, height = Height,
           useTeX = UseTex, standalone = UseTex)

proc plotFeSpectra(df: DataFrame,
                   outpath = "out") =
  proc plt(dfHalf: DataFrame, suff: string) =
    ggplot(dfHalf, aes("hits")) +
      facet_wrap("runNumber", scales = "free") +
      geom_histogram(binWidth = 1.0) +
      ggsave(&"out/spectra_{suff}.pdf",
             width = 2000, height = 1500)
  proc pltCharge(dfHalf: DataFrame, suff: string) =
    ggplot(dfHalf, aes("totalCharge")) +
      facet_wrap("runNumber", scales = "free") +
      geom_histogram(breaks = arange(0.0, 1.5e6, 3e3).toRawSeq) +
      ggsave(&"out/spectra_charge_{suff}.pdf",
             width = 2000, height = 1500)

  let df2017 = df.filter(f{int: `runNumber` <= 187 or `runNumber` == 0})
  let df2018 = df.filter(f{int: `runNumber` > 187 or `runNumber` == 0})
  echo df2017.unique("runNumber").len
  echo df2018.unique("runNumber").len
  plt(df2017, "2017_to_run_187")
  plt(df2018, "2018_from_run_187")
  plt(df, "all")

  pltCharge(df2017, "2017_to_run_187")
  pltCharge(df2018, "2018_from_run_187")
  pltCharge(df, "all")

proc concat(dfs: varargs[DataFrame]): DataFrame =
  for df in dfs:
    result.add df

proc getPeriods(df: DataFrame): OrderedTable[int, string] =
  let meanT = df.group_by("runPeriods").summarize(f{float -> float: "t" << mean(df["timestamp"])})["t"]
    .toTensor(int).toRawSeq.sorted
  let p = df.unique("runPeriods")["runPeriods"].toTensor(string).toRawSeq
    .mapIt(it.parseTime("dd/MM/YYYY", utc()).toUnix.int).sorted
    .mapIt(it.fromUnix.format("dd/MM/YYYY"))
  for (k, v) in zip(meanT, p):
    result[k] = v
  echo result

proc main(files: seq[string],
          interval, cutoffHits: float,
          calibFiles: seq[string] = @[],
          cutoffCharge: float = 0.0,
          createSpectra: bool = false,
          timeSeries: bool = false,
          photoDivEscape: bool = false,
          applyRegionCut = false,
          useLog = false,
          useMedian = false,
          toCount = false,
          readAllChips = false, ## Whether to only read center or all Septemboard chips
          normalizeMedian = false,
          outpath = "out",
          ylabel = "", title = "", titleSuff = "",
          countCutoff = 0.0
         ) =
  ## Input should be both H5 `DataRuns*_reco.h5` data files
  ## `interval` is the time to average per bin in minutes
  ## `cutoffCharge` is a filter on an amount of charge each cluster must
  ## have to take part
  let interval = interval * 60.0 # convert to seconds
  let dfBack = readFiles(files, "background", applyRegionCut = applyRegionCut, readAllChips = readAllChips)
  let dfCalib = readFiles(calibFiles, "calibration", applyRegionCut = applyRegionCut, readAllChips = readAllChips)
  ## check if there are additional files in the toml file

  if timeSeries:
    template all(args: varargs[DataFrame]): untyped =
      #let all1 = arg1.splitDf.computeStats()
      #let all2 = arg2.splitDf(all1.uniquePeriods).computeStats()
      var df = calculateMeanDf(args[0], interval, useMedian = useMedian, toCount = toCount)
      let periods = df.getPeriods
      for i in 1 ..< args.len:
        df.add calculateMeanDf(args[i], interval, periods, useMedian = useMedian, toCount = toCount)
      df
    proc filt(arg: DataFrame): DataFrame =
      arg.filter(f{float -> bool: `totalCharge` > cutoffCharge},
                 f{float -> bool: `hits` < cutoffHits})
    proc filtConc(args: varargs[DataFrame]): DataFrame =
      var df = filt(args[0]).splitDf().computeStats()
      let periods = df.uniquePeriods
      for i in 1 ..< args.len:
        df.add filt(args[i]).splitDf(periods).computeStats()
      #let all1 = filt(arg1).calculateMeanDf(interval)
      #let all2 = filt(arg2).calculateMeanDf(interval, all1.getPeriods)
      df

    let dfAll = if readAllChips:
                  doAssert calibFiles.len == 0, "Currently not supported to use calibration files when reading all chips"
                  ## Little hacky to get a sequence of dataframes
                  var dfs = newSeq[DataFrame]()
                  for (t, subDf) in groups(dfBack.group_by("chip")):
                    var mdf = subDf
                    mdf["chip"] = t[0][1].toInt
                    dfs.add mdf
                  all(dfs)
                else:
                  all(dfBack, dfCalib)
    echo dfAll.len
    #echo dfFilter
    let regionCut = if applyRegionCut: " cut to crSilver, 0.1 < rmsTrans < 1.5"
                    else: "cut to "

    let colorBy = if readAllChips: "chip" else: "runType"
    echo dfAll
    let titleSuff = if titleSuff.len > 0: titleSuff
                    else: &", {regionCut}, charge > {cutoffCharge}, hits < {cutoffHits.int}"
    plotOverTime(dfAll, interval,
                 titleSuff = titleSuff,
                 applyRegionCut = applyRegionCut,
                 useLog = useLog,
                 useMedian = useMedian,
                 toCount = toCount,
                 normalizeMedian = normalizeMedian,
                 colorBy = colorBy,
                 outpath = outpath,
                 ylabel = ylabel, title = title,
                 cutoff = countCutoff)
    if not toCount:
      let dfFilter = filtConc(dfBack, dfCalib)
      plotHistos(dfFilter, interval,
                 titleSuff = titleSuff,
                 applyRegionCut = applyRegionCut,
                 useLog = false,
                 outpath = outpath)

  if photoDivEscape:
    let dfBackMean = calculateMeanDf(dfBack, interval)
    let periods = dfBackMean.getPeriods
    let dfCalibMean = calculateMeanDf(dfCalib, interval, periods)
    let dfPhoto = readFiles(calibFiles, "Photo/Escape", applyRegionCut = applyRegionCut,
                            readAllChips = readAllChips,
                            readProc = readPhotoEscapeDf)
    echo dfPhoto

    plotPhotoDivEscape(dfPhoto, dfCalibMean, periods, interval, outpath = outpath)

  if createSpectra:
    var dfComb: DataFrame
    dfComb.add dfCalib.select(["totalCharge", "hits", "runType", "runNumber", "eventNumber"])
    dfComb.add readCdl()
    let dfFilter = dfComb.filter(f{float -> bool: `totalCharge` > cutoffCharge},
                                  f{float -> bool: `hits` < cutoffHits})
    plotFeSpectra(dfFilter, outpath = outpath)


when isMainModule:
  dispatch main
