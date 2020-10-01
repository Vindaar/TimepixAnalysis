import nimhdf5, ggplotnim, os, strformat
import ingrid / tos_helpers

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

proc caclulateMeanDf(df: DataFrame, interval: float): DataFrame =
  ## interval: time in minutes
  # we're going to do a rather crude approach. Given that we have to bin
  # a continuous variable, we're not well equipped using ggplotnim's dataframe
  # I fear. So just do it manually...
  let tstamps = df["timestamp"].toTensor(float)
  let totCh = df["totalCharge"].toTensor(float)
  var
    curSum = 0.0
    curNum = 0
    tStart = tstamps[0]
    ## NOTE: inefficient, but calculating the real number is rather annoying, so gonna live
    ## with inefficiency
    sums = newSeq[float]()
    means = newSeq[float]()
    tmeanStamps = newSeq[float]()
  for i in 0 ..< tstamps.size:
    if (tstamps[i] - tStart) >= (interval * 60.0):
      # run over bin range
      sums.add curSum
      means.add (curSum / curNum.float)
      tmeanStamps.add (tstamps[i] + tStart) / 2.0
      curSum = 0.0
      curNum = 0
      tStart = tstamps[i]
    curSum += totCh[i]
    inc curNum
  result = seqsToDf({"timestamp" : tmeanStamps, "sumCharge" : sums,
                      "meanCharge" : means })

proc plotDf(df: DataFrame, interval: float, titleSuff: string) =
  const outpath = "out"
  createDir(outpath)
  let nameSuff = if titleSuff == "all data": "all" else: "filtered"
  ggplot(df, aes("timestamp", "sumCharge")) +
    geom_point() +
    scale_y_log10() +
    ggtitle(&"Sum of total charge within {interval:.1f} min, {titleSuff}") +
    ggsave(&"{outpath}/background_sum_charge_binned_{interval:.1f}_min_{nameSuff}.pdf",
            width = 800, height = 480)
  ggplot(df, aes("timestamp", "meanCharge")) +
    geom_point() +
    scale_y_log10() +
    ggtitle(&"Mean of total charge within {interval:.1f} min, {titleSuff}") +
    ggsave(&"{outpath}/background_mean_charge_binned_{interval:.1f}_min_{nameSuff}.pdf",
            width = 800, height = 480)

  ggplot(df, aes("meanCharge")) +
    geom_histogram() +
    ggtitle(&"Histogram of binned mean of charge, within {interval:.1f} min") +
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

  plotDf(dfAll, interval, "all data")
  plotDf(dfFilter, interval, &"charge > {cutoffCharge}, hits < {cutoffHits:.0f} filtered out")


when isMainModule:
  dispatch main
