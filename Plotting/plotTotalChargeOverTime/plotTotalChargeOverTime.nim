import nimhdf5, ggplotnim, os, strformat
import ingrid / tos_helpers, times

import cligen

## See status and progress section
## ``Calculate total charge of background data over time``
## Just calculates the total charge binned over a specific time interval
## and plots it

type SupportedRead = SomeFloat | SomeInteger | string | bool | Value

proc readTstampDf(h5f: H5FileObj): DataFrame =
  const evNames = ["totalCharge", "eventNumber", "hits"]
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
    let dfJoined = inner_join(dfEv, dfAll, "eventNumber")
    result.add dfJoined

proc readFiles(files: seq[string]): DataFrame =
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

proc caclulateMeanDf(df: DataFrame, interval: float): DataFrame =
  ## interval: time in minutes
  # we're going to do a rather crude approach. Given that we have to bin
  # a continuous variable, we're not well equipped using ggplotnim's dataframe
  # I fear. So just do it manually...
  let tstamps = df["timestamp"].toTensor(float)
  let totCh = df["totalCharge"].toTensor(float)
  proc toPeriod(t: float): string =
    ## convert t to unix timestamp and then format date
    result = fromUnix(t.int).format("dd/MM/YYYY")

  var
    curSum = 0.0
    curNum = 0
    tStart = tstamps[0]
    ## NOTE: inefficient, but calculating the real number is rather annoying, so gonna live
    ## with inefficiency
    sums = newSeq[float]()
    means = newSeq[float]()
    tmeanStamps = newSeq[float]()
    runPeriods = newSeq[string]()
    lastPeriod = toPeriod tStart
  const splitPeriod = 86400 * 5 # 5 days?
  for i in 0 ..< tstamps.size:
    let inInterval = (tstamps[i] - tStart) < (interval * 60.0)
    if not inInterval:
      # run over bin range
      sums.add curSum
      means.add (curSum / curNum.float)
      let curTstamp = (tstamps[i] + tStart) / 2.0
      if curTstamp - tmeanStamps[^1] >= splitPeriod:
        lastPeriod = toPeriod curTstamp
      tmeanStamps.add curTstamp
      curSum = 0.0
      curNum = 0
      tStart = tstamps[i]
      runPeriods.add lastPeriod
    curSum += totCh[i]
    inc curNum
  result = seqsToDf({"timestamp" : tmeanStamps, "sumCharge" : sums,
                      "meanCharge" : means, "runPeriods" : runPeriods })
  let filterPeriods = result.group_by("runPeriods").count("runPeriods")
    .filter(f{`n` == 1})["runPeriods"].toTensor(string).toRawSeq
  result = result.filter(f{string -> bool: `runPeriods` notin filterPeriods})
  echo "Number of run periods with more than 1 entry: ", result["runPeriods"].unique

proc formatTime(v: float): string =
  result = v.int.fromUnix.format("dd/MM/YYYY")

proc plotDf(df: DataFrame, interval: float, titleSuff: string,
            useLog = true) =
  const outpath = "out"
  createDir(outpath)
  let nameSuff = if titleSuff == "all data": "all" else: "filtered"
  var pltSum = ggplot(df, aes("timestamp", "sumCharge")) +
    facet_wrap("runPeriods", scales = "free") +
    geom_point() +
    scale_x_continuous(labels = formatTime) +
    xlab(rotate = -45, alignTo = "right") +
    ggtitle(&"Sum of total charge within {interval:.1f} min, {titleSuff}")
  var pltMean = ggplot(df, aes("timestamp", "meanCharge")) +
    facet_wrap("runPeriods", scales = "free") +
    scale_x_continuous(labels = formatTime) +
    xlab(rotate = -45, alignTo = "right") +
    geom_point() +
    ggtitle(&"Mean of total charge within {interval:.1f} min, {titleSuff}")
  if useLog:
    pltSum = pltSum + scale_y_log10()
    pltMean = pltSum + scale_y_log10()
  pltSum + ggsave(&"{outpath}/background_sum_charge_binned_{interval:.1f}_min_{nameSuff}.pdf",
                   width = 1920, height = 1080)
  pltMean + ggsave(&"{outpath}/background_mean_charge_binned_{interval:.1f}_min_{nameSuff}.pdf",
                    width = 1920, height = 1080)

  ggplot(df, aes("meanCharge")) +
    geom_histogram() +
    ggtitle(&"Histogram of binned mean charge, within {interval:.1f} min, {titleSuff}") +
    ggsave(&"{outpath}/background_histo_mean_binned_{interval:.1f}_min_{nameSuff}.pdf",
            width = 800, height = 480)

proc main(files: seq[string], interval, cutoffCharge, cutoffHits: float) =
  ## Input should be both H5 `DataRuns*_reco.h5` data files
  ## `interval` is the time to average per bin in minutes
  ## `cutoffCharge` is a filter on an amount of charge each cluster must
  ## have to take part
  var df = readFiles(files)

  let dfAll = caclulateMeanDf(df, interval)
  let dfFilter = df.filter(f{float -> bool: `totalCharge` > cutoffCharge},
                           f{float -> bool: `hits` < cutoffHits})
    .caclulateMeanDf(interval)

  echo dfFilter
  echo dfAll
  plotDf(dfAll, interval, "all data")
  plotDf(dfFilter, interval,
         &"charge > {cutoffCharge}, hits < {cutoffHits:.0f} filtered out",
         useLog = false)


when isMainModule:
  dispatch main
