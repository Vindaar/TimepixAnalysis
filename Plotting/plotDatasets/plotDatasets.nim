import plotly, nimhdf5
import .. / karaPlot / [plotData, common_types]
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

const GridColor = color(0.8, 0.8, 0.8, 0.8)
const LegendBg = color(1.0, 1.0, 1.0, 0.5)
const LegendBorder = color(0.6, 0.6, 0.6, 0.6)
const Color1 = color(1.0, 0.0, 102.0 / 256.0)
const Color2 = color(0.0, 153.0 / 256.0, 204 / 256.0)

proc makePlot(fname, dset: string, chip: int,
              binRangeCustom = (0.0, 0.0)): Plot[float] =
  var h5f = H5file(fname, "r")
  var (binSize, binRange) = getBinSizeAndBinRange(dset)
  if binRangeCustom[0] != binRangeCustom[1]:
    binRange = binRangeCustom
  let fInfo = getFileInfo(h5f)
  let pd = PlotDescriptor(runType: rtNone,
                          name: dset,
                          runs: fInfo.runs,
                          plotKind: pkInGridDset,
                          chip: chip,
                          range: (low: -Inf, high: Inf, name: "All"),
                          binSize: binSize,
                          binRange: binRange)
  let (name, plt) = createPlot(h5f, fInfo, pd)
  result = plt.plPlot

proc main() =

  let args = docopt(doc)
  echo args

  let h5file = $args["<H5file>"]
  let h5file2 = $args["--file2"]
  let dset = $args["--dset"]
  let chip = ($args["--chip"]).parseInt
  var outfile = $dset

  let binLowS = $args["--binLow"]
  let binHighS = $args["--binHigh"]
  var
    binLow: float
    binHigh: float
  if binLowS != "nil":
    binLow = binLowS.parseFloat
    binHigh = binHighS.parseFloat

  BKind = bPlotly
  var plt = makePlot(h5file, dset, chip, (binLow, binHigh))
  if h5file2 != "nil":
    let plt2 = makePlot(h5file2, dset, chip, (binLow, binHigh))
    plt = plt.addTrace(plt2.traces[0])
      .legendLocation(0.85, 1.05)
      .legendBgColor(LegendBg)
      .legendBorderColor(LegendBorder)
      .legendBorderWidth(1)
      .gridColor(GridColor)
      .name("Signal", idx = 0)
      .name("Background", idx = 1)
      .markerColor(@[Color1], idx = 0)
      .markerColor(@[Color2], idx = 1)
      .width(800)
      .height(500)
      .title("Comparison of signal / background of " & $dset)

    plt.layout.barMode = BarMode.Overlay
    plt.traces[1].opacity = 0.5

    let font = Font(size: 18)
    #plt.layout.xaxis.font = font
    #plt.layout.yaxis.font = font
    #plt.layout.legend.font = font
    plt.layout.font = font
    outfile = outfile & "_merged"

  plt.layout.showLegend = true
  plt.show(outfile & ".svg")


when isMainModule:
  main()
