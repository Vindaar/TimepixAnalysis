import nimhdf5, ggplotnim, os, strformat, strutils, sequtils, algorithm, times, os, stats
import ingrid / [tos_helpers, ingrid_types, calibration]

import cligen

## See status and progress section
## ``Calculate total charge of background data over time``
## Just calculates the total charge binned over a specific time interval
## and plots it

type
  SupportedRead = SomeFloat | SomeInteger | string | bool | Value
  ReadProc = proc(h5f: H5File, applyRegionCut: bool): DataFrame

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

proc readTstampDf(h5f: H5File, applyRegionCut: bool): DataFrame =
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
    readIt(evNames, grp / "chip_3", dfEv)
    if applyRegionCut:
      # if we do not apply the region cut the slicing does not make any sense, because
      # the indices for the slices will be wrong!
      dfEv = dfEv.applyGasGainCut

      ## For event, read gas gain slices
      var sliceNum = newSeq[int](dfEv.len)
      let gsname = if grp / "chip_3" / "gasGainSlices90" in h5f: "gasGainSlices90"
                   else: "gasGainSlices"
      for (idx, slice) in iterGainSlicesIdx(h5f, grp / "chip_3", gsname):
        sliceNum[slice] = repeat(idx, slice.len)
      dfEv["sliceNum"] = sliceNum

    readIt(allNames, grp, dfAll)
    var dfJoined = inner_join(dfEv, dfAll, "eventNumber")
    dfJoined["runNumber"] = constantColumn(run, dfJoined.len)

    result.add dfJoined

proc readPhotoEscapeDf(h5f: H5File, applyRegionCut: bool): DataFrame =
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
    let chpGrpName = group.name / "chip_3"
    let dset = h5f[(chpGrpName / "FeSpectrumCharge").dset_str]
    photoSeq.add dset.attrs[parPrefix & $kalphaCharge, float]
    escapeSeq.add dset.attrs[parPrefix & $escapeCharge, float]
    dates.add parseTime(group.attrs["dateTime", string],
                        dateStr,
                        utc()).toUnix.float
  result = seqsToDf({ "timestamp" : dates, "photo" : photoSeq,
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
               readProc: ReadProc = readTstampDf): DataFrame =
  ## reads all H5 files given and stacks them to a single
  ## DF. An additional column is added, which corresponds to
  ## the filename. That way one can later differentiate which
  ## events belong to what and decide if data is supposed to
  ## be accumulated or not
  for file in files:
    let h5f = H5open(file, "r")
    result.add readProc(h5f, applyRegionCut = applyRegionCut)
    discard h5f.close()
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
              useMedian: static bool = false): seq[float] =
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
  var timesToPlot = newSeq[float]()
  when useMedian:
    var vec = newSeq[float]()
  for i in 0 ..< tstamps.size:
    let inInterval = (tstamps[i] - tStart) < interval
    if not inInterval:
      # run over bin range
      when toNorm:
        result[j] = current / norm
        norm = 0.0
      elif useMedian:
        result[j] = vec.median(50)
        vec = newSeq[float]()
      else:
        result[j] = current / numCluster.float
      numCluster = 0
      current = 0.0
      tStart = tstamps[i]
      timesToPlot.add tStart
      inc j
    when toNorm:
      norm += dataNorm[i]
    when useMedian:
      vec.add data[i]
    current += data[i]
    inc numCluster
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
                     periodTab: OrderedTable[int, string] = initOrderedTable[int, string]()): DataFrame =
  ## interval: time in minutes
  # we're going to do a rather crude approach. Given that we have to bin
  # a continuous variable, we're not well equipped using ggplotnim's dataframe
  # I fear. So just do it manually...
  let tstamps = df["timestamp", float]
  let tmeanstamps = calcMeanTimestamps(tstamps, interval)
  let runPeriods = mapToRunPeriods(tmeanStamps, periodTab)

  let outLen = tmeanstamps.len
  let sums = df.calcMean(tstamps, outLen, "hits", interval)
  let means = df.calcMean(tstamps, outLen, "totalCharge", interval, "hits")

  result = seqsToDf({"timestamp" : tmeanStamps, "sumCharge" : sums,
                      "meanCharge" : means, "runPeriods" : runPeriods })

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
  echo "Number of run periods with more than 1 entry: ", result["runPeriods"].unique

proc formatTime(v: float): string =
  result = v.int.fromUnix.format("dd/MM/YYYY")

proc plotDf(df: DataFrame, interval: float, titleSuff: string,
            useLog = true,
            applyRegionCut = false,
            outpath = "out") =
  let interval = interval / 60.0 # convert back to minutes for title
  createDir(outpath)
  var nameSuff = if titleSuff == "all data": "all" else: "filtered"
  if applyRegionCut:
    nameSuff.add "_crSilver"
  when false:
    var pltSum = ggplot(df, aes("timestamp", "sumCharge", color = "runType")) +
      facet_wrap("runPeriods", scales = "free") +
      geom_point(alpha = some(0.5)) +
      scale_x_continuous(labels = toPeriod) +
      xlab(rotate = -45, alignTo = "right") +
      ggtitle(&"Sum of total charge within {interval:.1f} min, {titleSuff}")
    var pltMean = ggplot(df, aes("timestamp", "meanCharge", color = "runType")) +
      facet_wrap("runPeriods", scales = "free") +
      scale_x_continuous(labels = toPeriod) +
      xlab(rotate = -45, alignTo = "right") +
      geom_point(alpha = some(0.5)) +
      ggtitle(&"Mean of total charge within {interval:.1f} min, {titleSuff}")
    if useLog:
      pltSum = pltSum + scale_y_log10()
      pltMean = pltMean + scale_y_log10()
    pltSum + ggsave(&"{outpath}/background_sum_charge_binned_{interval:.1f}_min_{nameSuff}.pdf",
                     width = 1920, height = 1080)
    pltMean + ggsave(&"{outpath}/background_mean_charge_binned_{interval:.1f}_min_{nameSuff}.pdf",
                      width = 1920, height = 1080)

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
      var pltTmp = ggplot(df, aes("timestamp", name, color = "runType")) +
        facet_wrap("runPeriods", scales = "free") +
        scale_x_continuous(labels = toPeriod) +
        xlab("timestamp", margin = 2.5, rotate = -45, alignTo = "right") +
        geom_point(alpha = some(0.5)) +
        ylim(2, 6.5) +
        margin(top = 1.75, bottom = 3) +
        ggtitle(&"{adn} of cluster {key} within {interval:.1f} min, {titleSuff}")
      var pltTmpHisto = ggplot(df, aes(name, fill = "runType")) +
        facet_wrap("runPeriods", scales = "free") +
        geom_histogram(bins = 300, density = true, position = "identity") +
        ggtitle(&"Histogram of {adn} cluster {key} within {interval:.1f} min, {titleSuff}")
      if useLog:
        pltTmp = pltTmp + scale_y_log10()
        pltTmpHisto = pltTmpHisto + scale_y_log10()
      pltTmp + ggsave(&"{outpath}/background_{adn.normalize}_{key}_{interval:.1f}_min_{nameSuff}.pdf",
                       width = 1920, height = 1080)
                       #width = 800, height = 480)
      pltTmpHisto + ggsave(&"{outpath}/background_histogram_{adn.normalize}_{key}_{interval:.1f}_min_{nameSuff}.pdf",
                            width = 1920, height = 1080)
      echo "created ", adn.normalize, " ", name
                            #width = 800, height = 480)

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
  ggplot(dfPlot, aes("timestamp", "val", color = "runType")) +
    facet_wrap("runPeriods", scales = "freeX") +
    geom_point(alpha = some(0.8)) +
    scale_x_continuous(labels = formatTime) +
    xlab(rotate = -45, alignTo = "right") +
    margin(top = 1.5) +
    ggtitle("Normalized cmp of photo/escape peak in charge & median energy " &
            &"{interval/60.0:.1f} min, max photo/esc {maxPhotoEsc:.2f}, max " &
            &"median energy {maxEnergy:.2f} {titleSuff}") +
    ggsave(&"{outpath}/photo_div_escape_vs_time.pdf", width = 1920, height = 1080)

proc plotFeSpectra(df: DataFrame,
                   outpath = "out") =
  proc plt(dfHalf: DataFrame, suff: string) =
    ggplot(dfHalf, aes("hits")) +
      facet_wrap("runNumber", scales = "free") +
      geom_histogram(binWidth = 1.0) +
      ggsave(&"out/spectra_{suff}.pdf",
             width = 2000, height = 1500)
  proc pltCharge(dfHalf: DataFrame, suff: string) =
    ggplot(dfHalf, aes("totalCharge"), numXTicks = 6) +
      facet_wrap("runNumber", scales = "free") +
      geom_histogram(breaks = arange(0.0, 1.5e6, 3e3).toRawSeq) +
      ggsave(&"out/spectra_charge_{suff}.pdf",
             width = 2000, height = 1500)

  let df2017 = df.filter(f{`runNumber` <= 187 or `runNumber` == 0})
  let df2018 = df.filter(f{`runNumber` > 187 or `runNumber` == 0})
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
          calibFiles: seq[string],
          interval, cutoffHits: float,
          cutoffCharge: float = 0.0,
          createSpectra: bool = false,
          timeSeries: bool = true,
          photoDivEscape: bool = false,
          applyRegionCut = false) =
  ## Input should be both H5 `DataRuns*_reco.h5` data files
  ## `interval` is the time to average per bin in minutes
  ## `cutoffCharge` is a filter on an amount of charge each cluster must
  ## have to take part
  let interval = interval * 60.0 # convert to seconds
  let dfBack = readFiles(files, "background", applyRegionCut = applyRegionCut)
  let dfCalib = readFiles(calibFiles, "acalibration", applyRegionCut = applyRegionCut)
  ## check if there are additional files in the toml file

  if timeSeries:
    template all(arg1, arg2: DataFrame): untyped =
      let all1 = calculateMeanDf(arg1, interval)
      let all2 = calculateMeanDf(arg2, interval, all1.getPeriods)
      concat(all1, all2)
    template filt(arg: DataFrame): untyped =
      arg.filter(f{float -> bool: `totalCharge` > cutoffCharge},
                 f{float -> bool: `hits` < cutoffHits})
    template filtConc(arg1, arg2: DataFrame): untyped =
      let all1 = filt(arg1).calculateMeanDf(interval)
      let all2 = filt(arg2).calculateMeanDf(interval, all1.getPeriods)
      concat(all1, all2)

    let dfAll = all(dfBack, dfCalib)
    let dfFilter = filtConc(dfBack, dfCalib)
    echo dfAll
    echo dfFilter
    let regionCut = if applyRegionCut: " cut to crSilver, 0.1 < rmsTrans < 1.5 "
                    else: ""
    plotDf(dfAll, interval, titleSuff = &"all data{regionCut}",
           applyRegionCut = applyRegionCut)

    plotDf(dfFilter, interval,
           titleSuff = &"{regionCut}charge > {cutoffCharge}, hits < {cutoffHits:.0f} filtered out",
           applyRegionCut = applyRegionCut,
           useLog = false)

  if photoDivEscape:
    let dfBackMean = calculateMeanDf(dfBack, interval)
    let periods = dfBackMean.getPeriods
    let dfCalibMean = calculateMeanDf(dfCalib, interval, periods)
    let dfPhoto = readFiles(calibFiles, "Photo/Escape", applyRegionCut = applyRegionCut,
                            readProc = readPhotoEscapeDf)
    echo dfPhoto

    plotPhotoDivEscape(dfPhoto, dfCalibMean, periods, interval)

  if createSpectra:
    var dfComb: DataFrame
    dfComb.add dfCalib.select(["totalCharge", "hits", "runType", "runNumber", "eventNumber"])
    dfComb.add readCdl()
    let dfFilter = dfComb.filter(f{float -> bool: `totalCharge` > cutoffCharge},
                                  f{float -> bool: `hits` < cutoffHits})
    plotFeSpectra(dfFilter)


when isMainModule:
  dispatch main
