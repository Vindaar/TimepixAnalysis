import nimhdf5
import ggplotnim
#import .. / karaPlot / [plotData, common_types]
import ingrid / [tos_helpers, ingrid_types]
import os, strutils, strformat
import sequtils
import docopt
import chroma
import seqmath
import times

from ../karaPlot / plotData import readEventSparse



const doc = """
A tool to create a plot of a InGrid / FADC dataset of a single or
two files.

Usage:
  plotDatasets <H5file> --dset=NAME --chip=NUMBER --eventDisplay [--file2 <H5File2>] [options]

Options:
  --file2 <H5File2>   Optional second file to compare with
  --dset=NAME         The name of the dataset to plot.
  --chip=NUMBER       The number of the chip we want to plot data from.
  --binLow=LOW        Low range of binning
  --binHigh=HIGH      High range of binning
  -h, --help          Show this help.
  --version           Show the version number.
"""

proc readDsets(h5f: H5FileObj, names: varargs[string]): DataFrame =
  ## reads all likelihood data in the given `h5f` file as well as the
  ## corresponding energies. Flattened to a 1D seq.
  ## This proc is for TPA generated H5 files! (i.e. containing run_* groups, ...)
  # iterate over all groups, read all likelihood and energy dsets
  result = newDataFrame()
  for name in names:
    var data = newColumn()
    for run, grp in runs(h5f):
      let group = h5f[grp.grp_str]
      let centerChip = "chip_" & $group.attrs["centerChip", int]
      doAssert grp / centerChip / name in h5f[(group.name / centerChip).grp_str]
      data = data.add toColumn(h5f[grp / centerChip / name, float])
    result[name] = data

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

proc main() =

  let args = docopt(doc)
  echo args

  let h5file = $args["<H5file>"]
  let h5file2 = $args["--file2"]
  let dset = $args["--dset"]
  let chip = ($args["--chip"]).parseInt
  var outfile = $dset
  let evDisplay = $args["--eventDisplay"]

  let binLowS = $args["--binLow"]
  let binHighS = $args["--binHigh"]
  var
    binLow: float
    binHigh: float
  try:
    binLow = binLowS.parseFloat
  except ValueError: discard
  try:
    binHigh = binHighS.parseFloat
  except ValueError: discard

  var h5f = H5open(h5file, "r")
  if evDisplay == "nil":
    let dftime = readFeVsTime(h5f)
    let df = h5f.readDsets(dset)
    ggplot(df,#df.filter(f{float: df[dset] > binLow and
           #                   df[dset] < binHigh}),
                     aes(dset)) +
      geom_histogram(bins = 50) +
      ggtitle(&"File: {h5file}, all runs, dsets: {dset}") +
      ggsave("/tmp/" & dset & ".pdf")
  else:
    for i in 0 .. 100:
      let df = readEventSparse(h5f, 305, 3, i)
      ggplot(df, aes("x", "y", color = "ch")) +
        geom_point() + #aes(width = 1.0, height = 1.0)) +
        xlim(0, 256) +
        ylim(0, 256) +
        ggsave("event_" & $i & ".pdf")



when isMainModule:
  main()
