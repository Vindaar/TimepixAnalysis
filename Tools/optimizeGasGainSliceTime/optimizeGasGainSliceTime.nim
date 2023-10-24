import shell, strutils, nimhdf5, os, times, strformat, tables, algorithm, sugar
import ggplotnim, stats

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
const Data2017 = "DataRuns2017_Reco.h5"
const Data2018 = "DataRuns2018_Reco.h5"

proc percentile[T: SomeNumber](t: Tensor[T], perc: float): float =
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

proc findFile(p: string, f: string): string =
  ## Returns the path to the file `f` if it exists in one of two allowed places, else
  ## raises.
  if existsFile(p / f):
    result = p / f
  elif existsFile(p / "2017" / f):
    result = p / "2017" / f
  elif existsFile(p / "2018_2" / f):
    result = p / "2018_2" / f
  else:
    raise newException(IOError, "The given path `" & $p & "` does not contain " & $f)

proc computeAndWriteFinalDf(path: string, interval: int) =
  var df = newDataFrame()
  for f in [Data2017, Data2018]:
    df.add computeMedianEnergy(path.findFile(f), interval)
  df = df.splitDf()
  df.mutate(fn {float -> int: "timestamp" ~ `timestamp`.int})
    .write_csv(&"out/median_energies_interval_{interval}.csv")

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

proc generateCsv(path: string) =
  createDir("out")
  var thr: array[2, Thread[tuple[f: string, interval: int]]]
  for i, t in timeIntervals:
    #if i == 0:
    #  # need to handle currently existing interval
    #  computeAndWriteFinalDf(30)
    ## XXX: replace hardcoded path here! Use CT TPA introspection module I need to write!
    ## Alternative would be a custom config file that we provide the `reconstruction` call with!
    rewriteToml("/home/basti/CastData/ExternCode/TimepixAnalysis/Analysis/ingrid/config.toml", t)
    createThread(thr[0], computeNewGainAndEnergies, (path.findFile(Data2017), t))
    createThread(thr[1], computeNewGainAndEnergies, (path.findFile(Data2018), t))
    joinThreads(thr)

    path.computeAndWriteFinalDf(t)

proc fromFile(f: string): DataFrame =
  result = readCsv(f)
  let interval = f.extractFilename.split(".")[0].split("_")[^1].parseInt
  result["Interval"] = constantColumn(interval, result.len)


import fitl / gof
import arraymancer / stats / kde

let UseTeX = getEnv("USE_TEX", "false").parseBool
let Width = getEnv("WIDTH", "1000").parseFloat
let Height = getEnv("HEIGHT", "600").parseFloat
let LineWidth = getEnv("LINE_WIDTH", "2.0").parseFloat
proc plot(path, outpath: string) =
  var df = newDataFrame()
  for f in walkFiles(path / "*.csv"):
    df.add fromFile(f)
  echo df
  df = df.mutate(fn {"medianEnergy" ~ `medianEnergy` * 1e6})
    .filter(fn {`medianEnergy` < 5.0})
  const bins = 150

  let labOrd = collect:
    for i, x in timeIntervals.sorted(SortOrder.Descending):
      { %~ x : i }
  echo labOrd
  let xlab = "Median energy [keV]"
  df = df.rename(f{xlab <- "medianEnergy"})

  # start with ridgeline and individual histograms of each period
  for (tup, subDf) in groups(df.group_by("RunPeriod")):
    let period = tup[0][1].toStr.replace("/", "_")
    ## XXX: fix code to work with different label order!!
    #let subDf = subDf.filter(fn {float: idx("medianEnergy") >= percentile(col("medianEnergy"), 0.01) and
    #                             idx("medianEnergy") <= percentile(col("medianEnergy"), 0.99)})
    let ridgeBins = 100
    ggplot(subDf, aes(xlab, fill = factor("Interval"))) +
      ggridges("Interval", overlap = 2.0) + # , labelOrder = labOrd) +
      geom_histogram(bins = ridgeBins,
                     hdKind = hdOutline,
                     position = "identity",
                     density = true,
                     color = "black", lineWidth = LineWidth) +
      #yMargin(0.2) +
      ggtitle(&"Different interval lengths for the gas gain computation in minutes, data starting {tup[0][1].toStr}") +
      ggsave(&"{outpath}/medianEnergy_ridges_{period}.pdf", width = Width, height = Height,
             useTeX = UseTeX, standalone = UseTeX)

    ggplot(subDf, aes(xlab, color = factor("Interval"))) +
      geom_histogram(bins = 40,
                     lineWidth = LineWidth,
                     alpha = some(0.0),
                     hdKind = hdOutline,
                     density = true,
                     position = "identity") +
      #xlim(2, 5) +
      ggtitle("Different interval lengths for the gas gain computation in minutes") +
      ggsave(&"{outpath}/medianEnergy_intervals_{period}.pdf", width = Width, height = Height,
             useTeX = UseTeX, standalone = UseTeX)

  # and the individual plots
  for (tup, subDf) in groups(df.group_by("Interval")):
    let suff = tup[0][1].toStr
    ggplot(subDf, aes("timestamp", xlab)) +
      facet_wrap("RunPeriod", scales = "free") +
      geom_point(alpha = some(0.7)) +
      ggtitle(&"Median energy vs time {suff}") +
      ggsave(&"{outpath}/medianEnergy_vs_time_{suff}.pdf", width = Width, height = Height,
             useTeX = UseTeX, standalone = UseTeX)

  # next plot the data as histograms in a facet
  let dfF = df.group_by("Interval")
    .filter(fn {float: idx(xlab) >= percentile(col(xlab), 0.01) and
                idx(xlab) <= percentile(col(xlab), 0.99)})
  ggplot(dfF, aes(xlab, color = factor("Interval"))) +
    facet_wrap("RunPeriod", scales = "free") +
    geom_histogram(bins = 30,
                   lineWidth = LineWidth,
                   alpha = 0.0,
                   hdKind = hdOutline,
                   density = true,
                   position = "identity") +
    xlim(2, 5) +
    ggtitle("Different interval lengths for the gas gain computation in minutes") +
    ggsave(&"{outpath}/medianEnergy_intervals.pdf", width = Width, height = Height,
           useTeX = UseTeX, standalone = UseTeX)

  # and finally using a KDE approach
  ## XXX: replace by new `geom_density` in ggplotnim!
  var dfK = newDataFrame()
  let samples = linspace(2.0, 4.0, 1000).toTensor
  for (tup, subDf) in groups(df.group_by(["Interval", "RunPeriod"])):
    let data = subDf[xlab, float]
    let kde = kde(data, samples = samples)
    var dfLoc = toDf({"Median energy [keV]" : samples, "Density" : kde, "RunPeriod" : tup[1][1].toStr, "Interval" : tup[0][1].toInt})
    dfK.add dfLoc

  ggplot(dfK, aes(xlab, "Density", color = factor("Interval"))) +
    facet_wrap("RunPeriod", scales = "free") +
    geom_line() +
    ggtitle("Different interval lengths for the gas gain computation in minutes") +
    ggsave(&"{outpath}/medianEnergy_kde_intervals.pdf", width = Width, height = Height,
           useTeX = UseTeX, standalone = UseTeX)

  for (tup, subDf) in groups(dfK.group_by("RunPeriod")):
    let period = tup[0][1].toStr.replace("/", "_")
    ggplot(subDf, aes(xlab, "Density", fill = factor("Interval"))) +
      ggridges("Interval", overlap = 2.0) + # , labelOrder = labOrd) +
      geom_line(color = "black", size = LineWidth) +
      ggtitle(&"Different interval lengths for the gas gain computation in minutes, data starting {tup[0][1].toStr}") +
      ggsave(&"{outpath}/medianEnergy_kde_ridges_{period}.pdf", width = Width, height = Height,
             useTeX = UseTeX, standalone = UseTeX)

  # As a bonus compute GoF tests for the data in each case and create a comparison plot.
  block Tests:
    let allMean = df[xlab, float].mean
    let allVar = df[xlab, float].variance
    var tests = newSeq[string]()
    var gofs = newSeq[float]()
    var ints = newSeq[int]()
    var pers = newSeq[string]()
    for test in GoFTest:
      for (tup, subDf) in groups(df.group_by(["Interval", "RunPeriod"])):
        let dat = subDf[xlab, float].toSeq1D
        let gof = gofStat(dat, (allMean, allVar), test)
        echo tup, ", ", test, " = ", gof
        tests.add $test
        gofs.add gof
        ints.add tup[0][1].toInt
        pers.add tup[1][1].toStr
    let dfT = toDf({"Goodness of fit test" : tests, "GoF value" : gofs, "Interval" : ints, "RunPeriod" : pers})
    ggplot(dfT, aes("Goodness of fit test", "GoF value", color = "Interval", shape = "RunPeriod")) +
      geom_point() +
      scale_y_log10() +
      ggtitle("Goodness of fit tests for median energy data by period & interval length") +
      ggsave(&"{outpath}/gofs_for_different_binnings.pdf", width = 600, height = 360,
             useTeX = UseTeX, standalone = UseTeX)

proc main(path: string, genCsv: bool = false, plot: bool = false,
          outpath = "/tmp/") =
  ## The given `path` is either a path to a place that contains both `DataRuns` H5 files
  ## (potentially in 2017, 2018_2 subdirectories) or a path to a directory containing CSV
  ## files that are generated from the `genCsv` call. Note that in the latter case the
  ## directory *must not* contain any other CSV files!
  if genCsv:
    generateCsv(path)
  if plot:
    plot(path, outpath)

import cligen
when isMainModule:
  dispatch main
