import nimhdf5, ggplotnim, strformat, sequtils, os
import ingrid / [tos_helpers, ingrid_types]

import cligen

proc toDf[T: object](data: seq[T]): DataFrame =
  ## Converts a seq of objects that (may only contain scalar fields) to a DF
  result = newDataFrame()
  for i, d in data:
    for field, val in fieldPairs(d):
      if field notin result:
        result[field] = newColumn(toColKind(type(val)), data.len)
      result[field, i] = val

proc main(files: seq[string], parameter = "theta") =

  var df = newDataFrame()
  for f in files:
    let h5f = H5file(f, "r")
    let fInfo = h5f.getFileInfo()
    for r in fInfo.runs:
      for c in fInfo.chips:
        let group = recoDataChipBase(r) & $c
        var gainSlicesDf = h5f[group & "/gasGainSlices90", GasGainIntervalResult].toDf
        gainSlicesDf["Chip"] = c
        gainSlicesDf["Run"] = r
        gainSlicesDf["File"] = f
        df.add gainSlicesDf
    discard h5f.close()
  df.showBrowser()
  df.writeCsv("/t/gas_gain_fit_parameters.csv")
  ggplot(df, aes(parameter, fill = factor("Chip"))) +
    geom_histogram(position = "identity", bins = 60,
                   hdKind = hdOutline, alpha = some(0.5)) +
    ggtitle("Input files: " & $files.mapIt(it.extractFilename).join(", ")) +
    ggsave(&"/tmp/{parameter}_values.pdf", width = 1000, height = 600)

  if files.len > 1:
    ggplot(df, aes(parameter, fill = factor("Chip"))) +
      facet_wrap("File", scales = "free") +
      geom_histogram(position = "identity", bins = 60,
                     hdKind = hdOutline, alpha = some(0.5)) +
      ggtitle("Input files: " & $files.mapIt(it.extractFilename).join(", ")) +
      ggsave(&"/tmp/{parameter}_values_by_files.pdf", width = 1200, height = 800)

    ggplot(df, aes("theta", "G_fit", fill = factor("Chip"))) +
      facet_wrap("File", scales = "free") +
      geom_point(size = some(2.0), alpha = some(0.9)) +
      ggtitle("Input files: " & $files.mapIt(it.extractFilename).join(", ")) +
      ggsave(&"/tmp/theta_vs_G_fit_values_by_files.pdf", width = 1200, height = 800)
  else:
    ggplot(df, aes(parameter, fill = factor("Chip"))) +
      geom_histogram(position = "identity", bins = 60,
                     hdKind = hdOutline, alpha = some(0.5)) +
      ggtitle("Input files: " & $files.mapIt(it.extractFilename).join(", ")) +
      ggsave(&"/tmp/{parameter}_values_by_files.pdf", width = 1200, height = 800)

    ggplot(df, aes(parameter, "G_fit", fill = factor("Chip"))) +
      geom_point(size = some(2.0), alpha = some(0.9)) +
      ggtitle("Input files: " & $files.mapIt(it.extractFilename).join(", ")) +
      ggsave(&"/tmp/{parameter}_vs_G_fit_values_by_files.pdf", width = 1200, height = 800)


  #let dfG = df.gather(["N", "G_fit", "theta", "redChiSq"], key = "Parameter", value = "Value")
  #ggplot(dfG, aes("Value", fill = factor("Chip"))) +
  #  facet_wrap("Parameter", scales = "free") +
  #  geom_histogram(position = "identity", bins = 60,
  #                 hdKind = hdOutline, alpha = some(0.5)) +
  #  ggtitle("Input files: " & $files.mapIt(it.extractFilename).join(", ")) +
  #  ggsave(&"/tmp/{parameter}_values.pdf", width = 1200, height = 800)


  #  N*: float # amplitude (?) of the polya fit, parameter 0
  #  G_fit*: float # gas gain as determined from fit, parameter 1
  #  G_fitmean*: float # gas gain as determined from mean of fit "distribution"
  #  G*: float # fit parameter as determined from mean of data histogram
  #  theta*: float # theta parameter of polya fit, parameter 2
  #  redChiSq*: float # reduced χ² of the polya fit of this slice


when isMainModule:
  dispatch main
