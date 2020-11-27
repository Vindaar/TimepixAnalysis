import nimhdf5, ggplotnim, os, strformat, strutils, sequtils, algorithm
import ingrid / tos_helpers, times, os

import cligen

## See status and progress section
## ``Calculate total charge of background data over time``
## Just calculates the total charge binned over a specific time interval
## and plots it

type SupportedRead = SomeFloat | SomeInteger | string | bool | Value

const CommonDsets = { "energyFromCharge" : "energy",
                      "fractionInTransverseRms" : "fracRmsTrans",
                      "lengthDivRmsTrans" : "L_div_RMS_trans",
                      "eccentricity" : "eccentricity" }

proc genNames(): seq[string] {.compileTime.} =
  result = @["totalCharge", "eventNumber", "hits"]
  for (d, v) in CommonDsets:
    result.add d

proc readTstampDf(h5f: H5FileObj): DataFrame =
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
    readIt(allNames, grp, dfAll)
    var dfJoined = inner_join(dfEv, dfAll, "eventNumber")
    dfJoined["runNumber"] = constantColumn(run.parseInt, dfJoined.len)
    result.add dfJoined

proc readCdl(): DataFrame =
  const path = "/mnt/1TB/CAST/CDL_2019/calibration-cdl-2018.h5"
  result = newDataFrame()
  if existsFile(path):
    const group = "/calibration-cdl-feb2019-Mn-Cr-12kV"
    const dsets = ["CdlSpectrum", "CdlSpectrumEvents", "CdlSpectrumCharge"]
    var h5f = H5file(path, "r")
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

proc readFiles(files: seq[string], runType: string): DataFrame =
  ## reads all H5 files given and stacks them to a single
  ## DF. An additional column is added, which corresponds to
  ## the filename. That way one can later differentiate which
  ## events belong to what and decide if data is supposed to
  ## be accumulated or not
  for file in files:
    let h5f = H5File(file, "r")
    result.add readTstampDf(h5f)
    discard h5f.close()
  result = result.arrange("timestamp")
  result["runType"] = constantColumn(runType, result.len)

proc calculateMeanDf(df: DataFrame, interval: float,
                     periodTab: OrderedTable[int, string] = initOrderedTable[int, string]()): DataFrame =
  ## interval: time in minutes
  # we're going to do a rather crude approach. Given that we have to bin
  # a continuous variable, we're not well equipped using ggplotnim's dataframe
  # I fear. So just do it manually...
  let tstamps = df["timestamp"].toTensor(float)
  let periods = toSeq(keys(periodTab))
  const splitPeriod = 86400 * 5 # 5 days?
  proc toPeriod(t: float): string =
    ## convert t to unix timestamp and then format date
    result = fromUnix(t.int).format("dd/MM/YYYY")

  let numIntervals = ((tstamps[tstamps.size - 1] - tstamps[0]) / interval).ceil.int
  var tmeanStamps = newSeqOfCap[float](numIntervals * 2)
  var tStart = tstamps[0]
  for i in 0 ..< tstamps.size:
    if (tstamps[i] - tStart) >= interval:
      # run over bin range
      tmeanStamps.add (tstamps[i] + tStart) / 2.0
      tStart = tstamps[i]

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
        inc j
      when toNorm:
        #echo "adding ", dataNorm[i]
        norm += dataNorm[i]
      when useMedian:
        vec.add data[i]
      current += data[i]
      inc numCluster

  var lastPeriod = if periods.len == 0: toPeriod tmeanStamps[0]
                   else: periodTab[periods[0]]
  var runPeriods = newSeq[string]()
  runPeriods.add lastPeriod
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
    runPeriods.add lastPeriod

  let outLen = tmeanStamps.len
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
            outpath = "out") =
  let interval = interval / 60.0 # convert back to minutes for title
  createDir(outpath)
  let nameSuff = if titleSuff == "all data": "all" else: "filtered"
  var pltSum = ggplot(df, aes("timestamp", "sumCharge", color = "runType")) +
    facet_wrap("runPeriods", scales = "free") +
    geom_point(alpha = some(0.5)) +
    scale_x_continuous(labels = formatTime) +
    xlab(rotate = -45, alignTo = "right") +
    ggtitle(&"Sum of total charge within {interval:.1f} min, {titleSuff}")
  var pltMean = ggplot(df, aes("timestamp", "meanCharge", color = "runType")) +
    facet_wrap("runPeriods", scales = "free") +
    scale_x_continuous(labels = formatTime) +
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

  for (k, name) in CommonDsets:
    let df = df.filter(f{float -> bool: `L_div_RMS_trans` < 1e13 and
                         `eccentricity` < 1e13 and
                         df[name][idx] > 0.0})
    if name == "eccentricity":
      echo "NOW ", titleSuff
    var pltTmp = ggplot(df, aes("timestamp", name, color = "runType")) +
      facet_wrap("runPeriods", scales = "free") +
      scale_x_continuous(labels = formatTime) +
      xlab(rotate = -45, alignTo = "right") +
      geom_point(alpha = some(0.5)) +
      ggtitle(&"Median of cluster {name} within {interval:.1f} min, {titleSuff}")
    if useLog:
      pltTmp = pltTmp + scale_y_log10()
    pltTmp + ggsave(&"{outpath}/background_mean_{name}_binned_{interval:.1f}_min_{nameSuff}.pdf",
                     width = 1920, height = 1080)

  ggplot(df, aes("meanCharge", color = "runType")) +
    geom_histogram(bins = 100, position = "identity", alpha = some(0.5)) +
    margin(top = 2) +
    ggtitle(&"Histogram of binned mean charge, within {interval:.1f} min, {titleSuff}") +
    ggsave(&"{outpath}/background_histo_mean_binned_{interval:.1f}_min_{nameSuff}.pdf",
            width = 800, height = 480)

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
          interval, cutoffCharge, cutoffHits: float,
          createSpectra: bool = false) =
  ## Input should be both H5 `DataRuns*_reco.h5` data files
  ## `interval` is the time to average per bin in minutes
  ## `cutoffCharge` is a filter on an amount of charge each cluster must
  ## have to take part
  let interval = interval * 60.0 # convert to seconds
  let dfBack = readFiles(files, "background")
  let dfCalib = readFiles(calibFiles, "acalibration")

  block TimeSeriesPlots:
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
    plotDf(dfAll, interval, "all data")
    plotDf(dfFilter, interval,
           &"charge > {cutoffCharge}, hits < {cutoffHits:.0f} filtered out",
           useLog = false)

  if createSpectra:
    var dfComb: DataFrame
    dfComb.add dfCalib.select(["totalCharge", "hits", "runType", "runNumber", "eventNumber"])
    dfComb.add readCdl()
    let dfFilter = dfComb.filter(f{float -> bool: `totalCharge` > cutoffCharge},
                                  f{float -> bool: `hits` < cutoffHits})
    plotFeSpectra(dfFilter)


when isMainModule:
  dispatch main
