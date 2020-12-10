import nimhdf5, ggplotnim, os, strformat, strutils, sequtils, algorithm
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
  var occ = newTensor[float]([NPix, NPix])
  var occCounts = newTensor[int]([NPix, NPix])
  # TODO: add option to use real event number instead of indices
  for i in 0 ..< xD.len:
    let xEv = xD[i]
    let yEv = yD[i]
    let chEv = chD[i]
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
        ggtitle(&"Run {run}, chip {chip}, slice {idx}, event {i}, hits {dfEv.len}") +
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
  let perc75 = zT.percentile(0.75)
  ggplot(df, aes("xT", "yT", fill = "zT")) +
    geom_raster() +
    scale_fill_continuous(scale = (low: 0.0, high: perc75)) +
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

proc readGasGains(f: string): DataFrame =
  var h5f = H5open(f, "r")
  for run, grp in runs(h5f, recoBase()):
    if run.parseInt != 164: continue
    let dset = h5f[(grp / "chip_3/charge").dset_str]
    let dfData = h5f.readGasGainDf(dset.name.parentDir, 3, @[])
    var gains: seq[float]
    var gainsFit: seq[float]
    var gainsMeanFit: seq[float]
    var times: seq[int]
    var sliceNum = 0
    for (g, slice) in iterGainSlicesFromAttrs(dset, dfData, 30.0):
      h5f.plotOccupancySlice(run.parseInt, 3, sliceNum, slice, dset.name.parentDir)

      echo "Run ", run, " ", g
      let gidx = g.toAttrPrefix()
      gainsFit.add dset.attrs[gidx & "G_fit", float]
      gains.add dset.attrs[gidx & "G", float]
      gainsMeanFit.add dset.attrs[gidx & "G_fitmean", float]
      times.add dset.attrs[g.toSliceStartAttr(), int]
      inc sliceNum
    var df = seqsToDf({ "Gain" : gains, "GainFit" : gainsFit,
                        "GainFitMean" : gainsMeanFit, "timestamp" : times })
    df["Run"] = constantColumn(run, df.len)
    df["SliceIdx"] = toSeq(0 ..< sliceNum)
    result.add df
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
  #echo perc95
  let dfBad = df.filter(f{float -> bool: `GasGain` > perc95})
  dfBad.write_csv("out/bad_run_slices.csv")

when isMainModule:
  dispatch main
