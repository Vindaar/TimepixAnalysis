import nimhdf5
import ggplotnim
#import .. / karaPlot / [plotData, common_types]
import ingrid / [tos_helpers, ingrid_types]
import os, strutils, strformat
import sequtils, options
import chroma
import seqmath
import times
import cligen

#from ../karaPlot / plotData import readEventSparse

proc readDsets(files: seq[string],
               names: seq[string],
               dsets: varargs[string],
               xrayCuts = true,
               tfKind: TargetFilterKind = tfMnCr12,
               run = -1,
               pathBase = likelihoodBase(),
              ): DataFrame =
  ## reads all likelihood data in the given `h5f` file as well as the
  ## corresponding energies. Flattened to a 1D seq.
  ## This proc is for TPA generated H5 files! (i.e. containing run_* groups, ...)
  # iterate over all groups, read all likelihood and energy dsets
  result = newDataFrame()
  let dsets = some((3, @dsets))
  for i, file in files:
    let h5f = H5open(file, "r")
    var df = h5f.readDsets(pathBase, chipDsets = dsets, run = run)
    if xrayCuts:
      df = df.cutXrayCleaning(tfKind, 5.0, 7.0)
    let fname = file.extractFilename
    df["File"] = constantColumn(fname, df.len)
    if names.len > 0:
      df["Type"] = names[i]
    result.add df

proc readFeVsTime(h5f: H5FileObj): DataFrame =
  ## reads all likelihood data in the given `h5f` file as well as the
  ## corresponding energies. Flattened to a 1D seq.
  ## This proc is for TPA generated H5 files! (i.e. containing run_* groups, ...)
  # iterate over all groups, read all likelihood and energy dsets
  result = newDataFrame()
  var data = newColumn()
  var means: seq[float]
  var dates: seq[float]
  const dateStr = "yyyy-MM-dd'.'HH:mm:ss" # example: 2017-12-04.13:39:45
  for run, grp in runs(h5f):
    echo "Accesing ", grp
    let group = h5f[grp.grp_str]
    let centerChip = "chip_" & $group.attrs["centerChip", int]
    let fedset = h5f[(group.name / centerChip / "FeSpectrum").dset_str]
    means.add fedset.attrs["p10", float]
    dates.add parseTime(group.attrs["dateTime", string], dateStr, utc()).toUnix.float

  result["means"] = means
  result["dates"] = dates
  ggplot(result, aes("dates", "means")) +
    geom_point() +
    ggsave("/tmp/fe_vs_time.pdf")

proc plotHisto(df: DataFrame, dsets: seq[string]) =

  ## split by high and low
  echo df
  var df = df
  df.drop("eventNumber")

  df = df.mutate(f{float -> bool: "8<x<10" ~ `energyFromCharge` >= 8.0 and `energyFromCharge` <= 10.0})
    .gather(dsets, key = "Dset", value = "vals")
  echo df
  ggplot(df, aes("vals", fill = "8<x<10")) +
    facet_wrap("Dset", scales = "free") +
    scale_x_continuous() +
    geom_histogram(bins = 100, binBy = "subset", position = "identity", alpha = some(0.5),
                   hdKind = hdOutline) +
    ggsave(&"out/dsets_facet_lhood.pdf", height = 1080, width = 1920)

proc plotKdeRidges(df: DataFrame, dsets: seq[string],
                   plotPath, prefix, suffix: string) =
  var df = df
  for dset in dsets:
    let d = dset
    df = df.mutate(f{float: d ~ idx(d) / col(d).max})
  df = df.gather(dsets, key = "Dset", value = "Value")
  echo "Plotting ridgeline"

  let by = if "Type" in df: "Type" else: "File"
  ggplot(df, aes("Value", fill = factor(by))) +
    ggridges("Dset", overlap = 4.0) +
    geom_density(normalize = true, alpha = 0.6, color = "black", size = 0.5) +
    ggtitle(&"Ridgeline plot for all InGrid properties in a.u.: x/max(x) {suffix}") +
    margin(left = 4) +
    ggsave(&"{plotPath}/{prefix}_ridgeline_kde_by_run.pdf", width = 900, height = 600)

proc main(files: seq[string],
          names: seq[string],
          eventDisplay: bool = false,
          xrayCuts = true,
          cutTfKind: TargetFilterKind = tfMnCr12,
          run = -1,
          cutLow = 0.0, cutHigh = Inf,
          chip = 3, dset = "",
          prefix = "ingrid_properties",
          plotPath = "/tmp/",
          suffix = "",
          feVsTime = false) =
  if names.len > 0 and files.len != names.len:
    raise newException(ValueError, "If any `names` given, there must be one name for each file. " &
      "Given files: " & $files & ", but only the following names: " & $names)
  var dsets = newSeq[string]()
  for dkKind in InGridDsetKind:
    if dkKind notin { igNumClusters, igFractionInHalfRadius, igRadiusDivRmsTrans,
                      igRadius, igBalance, igLengthDivRadius, igInvalid, igEventNumber,
                      igDiffusion, igGasGain, igEnergyFromPixel, igEnergyFromCdlFit,
                      igLikelihood }:
      dsets.add dkKind.toDset(fkTpa)

  var df = readDsets(files,
                     names = names,
                     dsets = dsets,
                     pathBase = recoBase(),
                     run = run,
                     tfKind = cutTfKind,
                     xrayCuts = xrayCuts)
  #df.plotHisto(dsets)

  df.plotKdeRidges(dsets, plotPath, prefix, suffix)

  #if feVsTime:
  #  let dftime = readFeVsTime(h5f)
  #  let df = h5f.readDsets(dset)
  #  ggplot(df,#df.filter(f{float: df[dset] > binLow and
  #         #                   df[dset] < binHigh}),
  #                   aes(dset)) +
  #    geom_histogram(bins = 50) +
  #    ggtitle(&"File: {h5file}, all runs, dsets: {dset}") +
  #    ggsave("/tmp/" & dset & ".pdf")
  #if eventDisplay:
  #  for i in 0 .. 100:
  #    let df = readEventSparse(h5f, 305, 3, i)
  #    ggplot(df, aes("x", "y", color = "ch")) +
  #      geom_point() + #aes(width = 1.0, height = 1.0)) +
  #      xlim(0, 256) +
  #      ylim(0, 256) +
  #      ggsave("event_" & $i & ".pdf")



when isMainModule:
  dispatch main
