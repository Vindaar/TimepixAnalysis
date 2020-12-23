import nimhdf5, ggplotnim, os, strformat, strutils, sequtils, algorithm, sets
import ingrid / [tos_helpers, calibration, ingrid_types], times, os

import cligen

proc percentile[T](t: Tensor[T], perc: float): float =
  let dataSorted = t.sorted
  let perIdx = min((t.size.float * perc).round.int, t.size - 1)
  result = dataSorted[perIdx]

proc readVlen(h5f: H5File,
              path: string,
              dsetName: string,
              dtype: typedesc,
              idx: seq[int]): seq[seq[dtype]] =
  ## reads variable length data `dsetName` and returns it
  ## In contrast to `read` this proc does *not* convert the data.
  let vlenDtype = special_type(dtype)
  let dset = h5f[(path / dsetName).dset_str]
  result = dset[vlenDType, dtype, idx]

proc readNorm(h5f: H5File,
              path: string,
              dsetName: string,
              dtype: typedesc,
              idx: seq[int]): seq[dtype] =
  ## reads variable length data `dsetName` and returns it
  ## In contrast to `read` this proc does *not* convert the data.
  let dset = h5f[(path / dsetName).dset_str]
  result = dset.read_hyperslab(dtype, @[idx[0], 0], @[idx.len, 1])

proc plotOccupancySlice(h5f: H5File, run, chip, idx: int, slice: Slice[int], path: string) =
  ## plots an occupancy of a single run slice
  let xD = h5f.readVlen(path, "x", uint8, toSeq(slice))
  let yD = h5f.readVlen(path, "y", uint8, toSeq(slice))
  let chD = h5f.readVlen(path, "charge", float, toSeq(slice))
  let evNums = h5f.readNorm(path, "eventNumber", int, toSeq(slice))
  var occ = newTensor[float]([NPix, NPix])
  var occCounts = newTensor[int]([NPix, NPix])
  # TODO: add option to use real event number instead of indices
  for i in 0 ..< xD.len:
    let xEv = xD[i]
    let yEv = yD[i]
    let chEv = chD[i]
    let evNum = evNums[i]
    for j in 0 ..< xEv.len:
      let x = xEv[j].int
      let y = yEv[j].int
      let ch = chEv[j]
      occ[x, y] += ch
      occCounts[x, y] += 1
    if run == 164 and idx in {57, 58}:
      let dfEv = seqsToDf(xEv, yEv, chEv)
      ggplot(dfEv, aes("xEv", "yEv", color = "chEv")) +
        geom_point() +
        xlim(0, 256) + ylim(0, 256) +
        margin(top = 1.5) +
        ggtitle(&"Run {run}, chip {chip}, slice {idx}, eventNum {evNum} eventIdx {i}, hits {dfEv.len}") +
        ggsave(&"out/event_{i}_run_{run}_chip_{chip}_slice_{idx}.pdf")

  const NPix = 256
  var
    xT = newTensorUninit[int](NPix * NPix)
    yT = newTensorUninit[int](NPix * NPix)
    zT = newTensorUninit[float](NPix * NPix)
    zCountT = newTensorUninit[int](NPix * NPix)
  var i = 0
  for idx, val in occ:
    let x = idx[0].int
    let y = idx[1].int
    xT[i] = x
    yT[i] = y
    zT[i] = val
    zCountT[i] = occCounts[x, y]
    inc i
  let df = seqsToDf(xT, yT, zT, zCountT)
  let perc = zT.percentile(0.99)
  ggplot(df, aes("xT", "yT", fill = "zT")) +
    geom_raster() +
    scale_fill_continuous(scale = (low: 0.0, high: perc)) +
    xlim(0, NPix) + ylim(0, NPix) +
    ggtitle(&"occupancy of pixel charges, run {run} chip {chip} slice {idx}") +
    ggsave(&"out/occupancy_charge_run_{run}_chip_{chip}_slice_{idx}.pdf", width = 1200,
            height = 1200)

  ggplot(df, aes("xT", "yT", fill = "zCountT")) +
    geom_raster() +
    scale_fill_continuous() +
    xlim(0, NPix) + ylim(0, NPix) +
    ggtitle(&"occupancy of pixel counts, run {run} chip {chip} slice {idx}") +
    ggsave(&"out/occupancy_counts_run_{run}_chip_{chip}_slice_{idx}.pdf", width = 1200,
            height = 1200)

proc plotSeptemEvents(df: DataFrame, run, slice: int) =
  var eventIdx = 0
  for tup, dfEv in groups(df.group_by("eventNumber")):
    let evNum = tup[0][1].toInt
    #if dfEv.len > 2000:
    let dfSeptem = dfToSeptemEvent(dfEv)
    let perc = df["charge"].toTensor(float).percentile(0.96)
    echo "Perc: ", perc
    ggplot(dfSeptem, aes("x", "y", color = "charge")) +
      geom_point(alpha = some(0.8)) +
      xlim(0, 3*256) + ylim(0, 3*256) +
      scale_color_continuous(scale = (low: 0.0, high: perc)) +
      ggtitle(&"evNumber {evNum}, run {run}, hits (all chips) {dfSeptem.len}") +
      ggsave(&"out/septem_events_run{run}_event_{evNum}_eventIdx_{eventIdx}_slice_{slice}.pdf")
    inc eventIdx

proc readGasGains(f: string): DataFrame =
  var h5f = H5open(f, "r")
  result = newDataFrame()
  const RunOfInterest = -1
  for runNumber, grp in runs(h5f, recoBase()):
    #if runNumber != RunOfInterest: continue
    let dset = h5f[(grp / "chip_3/charge").dset_str]
    let dfData = h5f.readGasGainDf(dset.name.parentDir, 3, @["rmsTransverse", "centerX", "centerY", "hits"])
    var gains: seq[float]
    var gainsFit: seq[float]
    var gainsMeanFit: seq[float]
    var times: seq[int]
    var sliceNum = 0
    let dfDataFilter = dfData.applyGasGainCut()
    #let dfRun = getSeptemDataFrame(h5f, runNumber, allowedChips = @[3])
    ## get the ``last`` gas gain slice computation (add `$(interval)` if desired for specific)
    let gasGainSlices = h5f[grp / "chip_3" / "gasGainSlices", GasGainIntervalResult]
    for g in gasGainSlices:
      # read the eventNumbers of this slice, then filter fullRun DF based on these
      # eventNumbers
      #let evNumbers = h5f.readNorm(dset.name.parentDir, "eventNumber", int,
      #                             toSeq(g.sliceStart .. g.sliceStop)).toSet
      #let evNumsFiltered = dfDataFilter
      #  .filter(f{int -> bool: `eventNumber` in evNumbers})["eventNumber"].toTensor(int)
      #  .toRawSeq.toSet
      #let dfEvs = dfRun.filter(f{int -> bool: `eventNumber` in evNumsFiltered})
      #if runNumber == RunOfInterest and g.idx in {8, 9, 10}:
      #  if g.idx == 9:
      #    echo "Event numbers: ", evNumbers.toSeq.sorted
      #    plotSeptemEvents(dfEvs, runNumber, g.idx)
      #h5f.plotOccupancySlice(runNumber, 3, sliceNum, slice, dset.name.parentDir)

      echo "Run ", runNumber, " ", g
      gainsFit.add g.G_fit
      gains.add g.G
      gainsMeanFit.add g.G_fitMean
      times.add g.tStart
      inc sliceNum
    var df = seqsToDf({ "Gain" : gains, "GainFit" : gainsFit,
                        "GainFitMean" : gainsMeanFit, "timestamp" : times })
    df["Run"] = constantColumn(runNumber, df.len)
    df["SliceIdx"] = toSeq(0 ..< sliceNum)
    result.add df
  if result.len > 0:
    result = result.arrange("timestamp")

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

proc main(files: seq[string]) =
  var df = newDataFrame()
  for f in files:
    df.add f.readGasGains
  block IndividualPlots:
    ggplot(df, aes("timestamp", "Gain")) +
      geom_point() +
      ggsave("/tmp/test.pdf")
    ggplot(df, aes("timestamp", "GainFit")) +
      geom_point() +
      ggsave("/tmp/testFit.pdf")
    ggplot(df, aes("timestamp", "GainFitMean")) +
      geom_point() +
      ggsave("/tmp/testMeanFit.pdf")

  df = df.splitDf
    .gather(@["Gain", "GainFit", "GainFitMean"], "CalcType", "GasGain")

  createDir("out")
  ggplot(df, aes("timestamp", "GasGain", color = "CalcType")) +
    facet_wrap("RunPeriod", scales = "free_x") +
    geom_point(alpha = some(0.5)) +
    #geom_text(data = df.filter(f{float -> bool: `GasGain` > perc95}),
    #          aes = aes(x = "timestamp", #f{`timestamp`},
    #                    y = "GasGain",
    #                    text = "Run")) +
    scale_x_continuous(labels = formatTime) +
    #ylim(0, 8000.0) +
    xlab(rotate = -45, alignTo = "right") +
    ggsave("out/gas_gain_slices_vs_time_30_minutes.pdf", width = 1920, height = 1080)

  let perc95 = df["GasGain"].toTensor(float).max * 0.80
  echo "Runs with worst gas gain above: ", perc95, " stored in out/bad_run_slices.csv"
  let dfBad = df.filter(f{float -> bool: `GasGain` > perc95})
  dfBad.write_csv("out/bad_run_slices.csv")

when isMainModule:
  dispatch main
