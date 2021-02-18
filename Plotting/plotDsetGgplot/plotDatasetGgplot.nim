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
               names: varargs[string],
               pathBase = likelihoodBase()): DataFrame =
  ## reads all likelihood data in the given `h5f` file as well as the
  ## corresponding energies. Flattened to a 1D seq.
  ## This proc is for TPA generated H5 files! (i.e. containing run_* groups, ...)
  # iterate over all groups, read all likelihood and energy dsets
  result = newDataFrame()
  let dsets = some((3, @names))
  for file in files:
    let h5f = H5open(file, "r")
    var df = h5f.readDsets(pathBase, chipDsets = dsets)
    let fname = file.extractFilename
    df["File"] = constantColumn(fname, df.len)
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

proc plotHisto(df: DataFrame, names: seq[string]) =

  ## split by high and low
  echo df
  var df = df
  df.drop("eventNumber")

  df = df.mutate(f{float -> bool: "8<x<10" ~ `energyFromCharge` >= 8.0 and `energyFromCharge` <= 10.0})
    .gather(names, key = "Dset", value = "vals")
  echo df
  ggplot(df, aes("vals", fill = "8<x<10")) +
    facet_wrap("Dset", scales = "free") +
    scale_x_continuous() +
    geom_histogram(bins = 100, binBy = "subset", position = "identity", alpha = some(0.5),
                   hdKind = hdOutline) +
    ggsave(&"out/dsets_facet_lhood.pdf", height = 1080, width = 1920)

proc main(files: seq[string],
          eventDisplay: bool = false,
          cutLow = 0.0, cutHigh = Inf,
          chip = 3, dset = "",
          feVsTime = false) =
  var names = newSeq[string]()
  for dkKind in InGridDsetKind:
    if dkKind notin {igNumClusters, igFractionInHalfRadius, igRadiusDivRmsTrans,
                      igRadius, igBalance, igLengthDivRadius, igInvalid, igHits, igTotalCharge, igEventNumber}:
      names.add dkKind.toDset(fkTpa)

  var df = readDsets(files, names = names)
  df.plotHisto(names)


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
