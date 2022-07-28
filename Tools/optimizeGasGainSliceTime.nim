import shell, parsetoml, strutils, nimhdf5, os, times, strformat
import ggplotnim

import ingrid / [tos_helpers, calibration, ingrid_types]

#[
Need to optimize...

1. compute gas gain slices
2. they are stored in `gasGainSlices` dataset.
  -> need a way to keep all slices in the file so that we can just read any
     slicing using the `interval` from the `config.toml` file?
3. calibrate the energy based on a chosen `interval`
4. compute the median energies for each time slice
5. create histogram for each run period
6. perform fitting of a gaussian

We could store the raw median cluster data as a CSV file to quickly have
all data available later on? That way we only have to do steps 1-4 first.

Things to take note:
- computation of gas gain intervals by far the slowest part.
- computation of energy clusters fast
- computation of median energy fast

So we can just compute let's say 6 different gas gain time slicings and
compute the energies and median cluster energies. Then store these as
CSV files. Once that is done we perform the gaussian fits on all inputs
and try to find some measure for the "goodness" that matches our understanding.
From there we could in principle optimize for such a value. But doing
the optimization without having a bunch of different histograms is too
problematic.

On a practical side:
1. We have to read the `config.toml` file each time and replace it by the version
containing the interval we wish to compute.

2. Then run the reconstruction with `--only_gas_gain`, which creates a dataset
`gasGainSlices<IntervalLength>`.

3. Recompute the energy using `--only_energy_from_e`, which overrides the existing
energies. This means we have to start with computing the median values and write to
CSV for the existing 30 min stored in the files first.

4. Compute the median cluster energies according to `plotTotalChargeVsTime`.

5. Store DF as CSV.
]#

const timeIntervals = [45, 60, 90, 120, 180, 300]
const Data2017 = "/mnt/1TB/CAST/2017/DataRuns2017_Reco.h5"
const Data2018 = "/mnt/1TB/CAST/2018_2/DataRuns2018_Reco.h5"

proc percentile[T](t: Tensor[T], perc: float): float =
  let dataSorted = t.sorted
  let perIdx = min((t.size.float * perc).round.int, t.size - 1)
  result = dataSorted[perIdx]

proc rewriteToml(path: string, interval: int) =
  ## rewrites the given TOML file in the `path` to use the `interval`
  ## instead of the existing value
  var data = readFile(path).splitLines
  for l in mitems(data):
    if l.startsWith("gasGainInterval"):
      l = "gasGainInterval = " & $interval
  writeFile(path, data.join("\n"))

proc formatTime(v: float): string =
  result = v.int.fromUnix.format("dd/MM/YYYY")

proc splitDf(df: DataFrame): DataFrame =
  var periodSeq = newSeq[string](df.len)
  let times = df["timestamp"].toTensor(int)
  var lastTimestamp = times[0]
  var period = formatTime(lastTimestamp.float)
  for i in 0 ..< times.size:
    let t = times[i]
    if t - lastTimestamp > 10 * 86400:
      period = formatTime t.float
    periodSeq[i] = period
    lastTimestamp = t
  result = df
  result["RunPeriod"] = periodSeq

proc computeMedianEnergy(f: string, interval: int): DataFrame =
  result = newDataFrame()
  var h5f = H5open(f, "r")
  for run, grp in runs(h5f):
    echo "Reading run ", run
    let df = h5f.readGasGainDf(grp / "chip_3",
                               @["rmsTransverse", "centerX", "centerY",
                                 "hits", "energyFromCharge"])
      .applyGasGainCut()
    let dsetName = if interval == 30: "gasGainSlices" else: "gasGainSlices" & $(interval)
    let gainSlices = h5f[grp / "chip_3" / dsetName, GasGainIntervalResult]
    var energy = newSeq[float](gainSlices.len)
    var t = newSeq[float](gainSlices.len)
    var i = 0
    for subDf in iterGainSlicesDf(df, gainSlices):
      energy[i] = subDf["energyFromCharge", float].percentile(0.5)
      t[i] = subDf["timestamp", float].percentile(0.5)
      inc i
    result.add toDf({"timestamp" : t, "medianEnergy" : energy})
  result = result.arrange("timestamp")
  discard h5f.close()

proc computeAndWriteFinalDf(interval: int) =
  var df = newDataFrame()
  for f in [Data2017, Data2018]:
    df.add computeMedianEnergy(f, interval)
  df = df.splitDf()
  df.write_csv(&"out/median_energies_interval_{interval}.csv")

#computeAndWriteFinalDf(30)
when false:
  ggplot(df, aes("timestamp", "medianEnergy")) +
    facet_wrap("RunPeriod") +
    geom_point() +
    ggsave("/tmp/test_median.pdf")

proc computeNewGainAndEnergies(arg: tuple[f: string, interval: int]) {.thread.} =
  let (f, interval) = arg
  echo "from thread ", f, " to handle ", interval

  let (res, err, code) = shellVerboseErr:
    #reconstruction ($f) "--only_gas_gain"
    reconstruction ($f) "--only_energy_from_e"
  if code != 0:
    raise newException(Exception, "Error running the given command for " & $f)

proc generateCsv() =
  createDir("out")
  var thr: array[2, Thread[tuple[f: string, interval: int]]]
  for i, t in timeIntervals:
    #if i == 0:
    #  # need to handle currently existing interval
    #  computeAndWriteFinalDf(30)
    rewriteToml("/home/basti/CastData/ExternCode/TimepixAnalysis/Analysis/ingrid/config.toml", t)
    createThread(thr[0], computeNewGainAndEnergies, (Data2017, t))
    createThread(thr[1], computeNewGainAndEnergies, (Data2018, t))
    joinThreads(thr)

    computeAndWriteFinalDf(t)

proc fromFile(f: string): DataFrame =
  result = toDf(readCsv(f))
  let interval = f.extractFilename.split(".")[0].split("_")[^1].parseInt
  result["Interval"] = constantColumn(interval, result.len)

proc plot() =
  var df = newDataFrame()
  for (pcKind, f) in walkDir("out"):
    case pcKind
    of pcFile: df.add fromFile(f)
    else: echo "Skipping ", f, " of kind ", pcKind
  echo df
  df = df.mutate(fn {"medianEnergy" ~ `medianEnergy` * 1e6})
    .filter(fn {`medianEnergy` < 5.0})
  const bins = 150

  ggplot(df, aes("medianEnergy", color = "Interval")) +
    facet_wrap("RunPeriod", scales = "free") +
    geom_histogram(bins = bins,
                   lineWidth = some(1.25),
                   alpha = some(0.0),
                   hdKind = hdOutline,
                   # density = true
                   position = "identity") +
    xlim(2, 5) +
    ggtitle("Different interval lengths for the gas gain computation in minutes") +
    ggsave("/tmp/medianEnergy_intervals.pdf", width = 1200, height = 800)
  for (tup, subDf) in groups(df.group_by("RunPeriod")):
    let period = tup[0][1].toStr.replace("/", "_")
    ggplot(subDf, aes("medianEnergy", fill = "Interval")) +
      ggridges("Interval", overlap = 1.0) +
      geom_histogram(bins = bins,
                     hdKind = hdOutline,
                     position = "identity") +
      ggtitle("Different interval lengths for the gas gain computation in minutes") +
      ggsave(&"/tmp/medianEnergy_ridges_{period}.pdf", width = 1200, height = 800)

  for (tup, subDf) in groups(df.group_by("Interval")):
    let suff = tup[0][1].toStr
    ggplot(subDf, aes("timestamp", "medianEnergy")) +
      facet_wrap("RunPeriod", scales = "free") +
      geom_point(alpha = some(0.7)) +
      ggtitle(&"Median energy vs time {suff}") +
      ggsave(&"/tmp/medianEnergy_vs_time_{suff}.pdf", width = 1200, height = 800)


proc main(genCsv: bool = false, plot: bool = false) =
  if genCsv:
    generateCsv()
  if plot:
    plot()

import cligen
when isMainModule:
  dispatch main
