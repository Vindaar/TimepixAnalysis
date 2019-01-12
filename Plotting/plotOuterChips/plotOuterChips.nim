import .. / karaPlot / plotData
import docopt, strutils
import nimhdf5, plotly, chroma
import ingrid / ingrid_types

when defined(linux):
  const commitHash = staticExec("git rev-parse --short HEAD")
else:
  const commitHash = ""

# get date using `CompileDate` magic
const currentDate = CompileDate & " at " & CompileTime

const docTmpl = """
Version: $# built on: $#
A tool to plot a histogram of events on outer chips

Usage:
  plotOuterChips <H5file> <H5file2> [options]

Options:
  -h, --help             Show this help
  --version              Show the version number
"""
const doc = docTmpl % [commitHash, currentDate]

proc makeOuterChipPlot(fname: string, runType: RunTypeKind): PlotV =
  var h5f = H5file(fname, "r")
  let fInfo = getFileInfo(h5f)
  let pds = createOuterChipHistograms(h5f,
                                      runType = runType,
                                      fileInfo = fInfo)
  let (name, plt) = createPlot(h5f, fInfo, pds[0])
  result = plt

proc main =

  # somewhat ugly hack
  plotData.BKind = bPlotly

  let args = docopt(doc)
  let backFname = $args["<H5file>"]
  let calibFname = $args["<H5file2>"]

  let pltBack = makeOuterChipPlot(backFname, rtBackground)
  let pltCalib = makeOuterChipPlot(calibFname, rtCalibration)

  let gridColor = color(0.8, 0.8, 0.8, 0.8)

  var plt = pltBack.plPlot
  plt = plt.addTrace(pltCalib.plPlot.traces[0])
    .title("# pixels on outer chips for blob on center")
    .xlabel("Outer chips # pixels hit")
    .gridColor(gridColor)
    .markerColor(@[color(1.0, 0.0, 102.0 / 256.0)], idx = 0)
    .markerColor(@[color(0.0, 153.0 / 256.0, 204 / 256.0)], idx = 1)

  plt.traces[0].opacity = 1.0
  plt.traces[1].opacity = 0.5
  plt.layout.barMode = BarMode.Overlay

  plt.traces[0].histNorm = HistNorm.None
  plt.traces[1].histNorm = HistNorm.None
  plt.show("outerHits_blobCenter.svg")

  plt.traces[0].histNorm = HistNorm.ProbabilityDensity
  plt.traces[1].histNorm = HistNorm.ProbabilityDensity
  plt.layout.yaxis.title = "Probability density"

  plt.show("outerHits_blobCenter_normalizedPDF.svg")

when isMainModule:
  main()
