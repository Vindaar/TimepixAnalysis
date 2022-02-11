import nimhdf5
import ggplotnim
#import .. / karaPlot / [plotData, common_types]
import ingrid / [tos_helpers, ingrid_types]
import os, strutils, strformat
import sequtils
import docopt
import chroma
import seqmath

const doc = """
A tool to create a plot of a InGrid / FADC dataset of a single or
two files.

Usage:
  plotDatasets <H5file> --dset=NAME --chip=NUMBER [--file2 <H5File2>] [options]

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
  result = initDataFrame()
  for name in names:
    var data = newSeq[float]()
    for run, grp in runs(h5f):
      let group = h5f[grp.grp_str]
      let centerChip = "chip_" & $group.attrs["centerChip", int]
      doAssert grp / centerChip / name in h5f[(group.name / centerChip).grp_str]
      data.add h5f[grp / centerChip / data, float64]
    result[name] = toColumn data



proc main() =

  let args = docopt(doc)
  echo args

  let h5file = $args["<H5file>"]
  let h5file2 = $args["--file2"]
  let dset = $args["--dset"]
  let chip = ($args["--chip"]).parseInt
  var outfile = $dset
  let evDisplay = $args["--eventDisplay"].toBool

  let binLowS = $args["--binLow"]
  let binHighS = $args["--binHigh"]
  var
    binLow: float
    binHigh: float
  if binLowS != "nil":
    binLow = binLowS.parseFloat
    binHigh = binHighS.parseFloat

  var h5f = H5open(h5file, "r")
  let df = h5f.readDsets(dset)
  ggplot(df, aes(dset)) +
    geom_histogram() +
    ggtitle(&"File: {h5file}, all runs, dsets: {dset}") +
    ggsave("/tmp/" & dset & ".pdf")



when isMainModule:
  main()
