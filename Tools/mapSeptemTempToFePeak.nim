import os except FileInfo
import ggplotnim, seqmath, sequtils, tables, options, times
import nimhdf5, cligen

import ingrid / [tos_helpers, ingrid_types]

type
  FeFileKind = enum
    fePixel, feCharge, feFadc

const TempFile = "../resources/cast_2017_2018_temperatures.csv"
const OrgFormat = "'<'yyyy-MM-dd ddd H:mm'>'"
var dfTemp = toDf(readCsv(TempFile))
  .filter(f{c"Temp / °" != "-"})
dfTemp["Timestamp"] = dfTemp["Date"].toTensor(string).map_inline(parseTime(x, OrgFormat, utc()).toUnix)

const Peak = "μ"
const PeakNorm = "μ/μ_max"
const TempPeak = "(μ/T) / max"

proc readFePeaks(files: seq[string], feKind: FeFileKind = fePixel): DataFrame =
  const kalphaPix = 10
  const kalphaCharge = 4
  const parPrefix = "p"
  const dateStr = "yyyy-MM-dd'.'HH:mm:ss" # example: 2017-12-04.13:39:45
  var dset: string
  var kalphaIdx: int
  case feKind
  of fePixel:
    kalphaIdx = kalphaPix
    dset = "FeSpectrum"
  of feCharge:
    kalphaIdx = kalphaCharge
    dset = "FeSpectrumCharge"
  of feFadc:
    kalphaIdx = kalphaCharge
    dset = "FeSpectrumFadcPlot" # raw dataset is `minVal` instead of `FeSpectrumFadc`

  var h5files = files.mapIt(H5open(it, "r"))
  var fileInfos = newSeq[FileInfo]()
  for h5f in mitems(h5files):
    let fi = h5f.getFileInfo()
    fileInfos.add fi
  var
    peakSeq = newSeq[float]()
    dateSeq = newSeq[float]()
  for (h5f, fi) in zip(h5files, fileInfos):
    for r in fi.runs:
      let group = h5f[(recoBase() & $r).grp_str]
      let chpGrpName = if feKind in {fePixel, feCharge}: group.name / "chip_3"
                       else: group.name / "fadc"
      peakSeq.add h5f[(chpGrpName / dset).dset_str].attrs[
        parPrefix & $kalphaIdx, float
      ]
      dateSeq.add parseTime(group.attrs["dateTime", string],
                            dateStr,
                            utc()).toUnix.float
  result = toDf({ Peak : peakSeq,
                  "Timestamp" : dateSeq })
    .arrange("Timestamp", SortOrder.Ascending)
    .mutate(f{float: PeakNorm ~ idx(Peak) / max(col(Peak))},
            f{"Type" <- $feKind})

proc findTemp(d: int, dates: seq[int], temps: Tensor[float]): float =
  let idx = dates.lowerBound(d)
  let idxMin = max(0, idx - 1)
  if idx == idxMin:
    result = temps[idx]
  else:
    let
      idxMoreDiff = dates[idx] - d
      idxLessDiff = dates[idxMin] - d
      # interpolate temp between the two temps corresponding both idx
      t0 = temps[idxMin]
      t1 = temps[idx]
      totDiff = dates[idx] - dates[idxMin]
    echo "tot diff ", totDiff, "   ", fromUnix(dates[idx]), "    ", fromUnix(dates[idxMin]), "   ", fromUnix(d)
    result = (t0 * (abs(totDiff - idxLessDiff).float) +
              t1 * abs(totDiff - idxMoreDiff).float) /
      totDiff.float

proc main(calibFiles: seq[string], inputs: set[FeFileKind]) =
  var dfPeaks = newDataFrame()
  for input in inputs:
    var df = readFePeaks(calibFiles, input)
    ggplot(df, aes("Timestamp", Peak)) +
      geom_point() +
      ggsave("/tmp/time_vs_peak_pos_$#.pdf" % $input)
    dfPeaks.add df
  if inputs.card > 1:
    ggplot(dfPeaks, aes("Timestamp", PeakNorm, color = "Type")) +
      geom_point() +
      ggsave("/tmp/time_vs_peak_pos.pdf")

  let runs = dfTemp["Run number"].toTensor(int)
  let temps = dfTemp["Temp / °"].toTensor(float)
  let dates = dfTemp["Timestamp"].toTensor(int).toRawSeq
  let feDates = dfPeaks["Timestamp"].toTensor(int)
  let peaks = dfPeaks[Peak].toTensor(float)
  var peakNorm = newSeq[float](peaks.size)

  var idx = 0
  for date in feDates:
    let temp = findTemp(date, dates, temps)
    peakNorm[idx] = peaks[idx] / (temp + 273.15)
    echo "Temp ", temp, " peak ", peaks[idx], " norm ", peakNorm[idx]
    inc idx
  var df = toDf({ TempPeak : peakNorm, "Timestamp" : feDates })
  ggplot(df, aes(Timestamp, TempPeak)) +
    geom_point() +
    ggsave("/tmp/time_vs_peak_norm_by_temp.pdf")

  # finally combine the twe into one plot
  dfPeaks = dfPeaks.mutate(f{Peak ~ df[Peak][idx] / max(df[Peak])})
  df = df.mutate(f{TempPeak ~ df[TempPeak][idx] / max(df[TempPeak])})
  df = innerJoin(df, dfPeaks, by = "Timestamp").gather([TempPeak, Peak], "Normalization", "peak")
  ggplot(df, aes(Timestamp, peak, color = "Normalization")) +
    geom_point() +
    ylab("Normalized ⁵⁵Fe peak position") +
    ggtitle("Variation of peak position of ⁵⁵Fe peak over time, normalized by temperature") +
    # ggsave("/tmp/time_vs_peak_temp_normed_comparison.pdf")
    ggsave("/tmp/time_vs_peak_temp_normed_comparison.pdf", width = 800, height = 480)


when isMainModule:
  dispatch(main, echoResult = false, noAutoEcho = true)
