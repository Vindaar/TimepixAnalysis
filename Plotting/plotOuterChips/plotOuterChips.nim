import .. / karaPlot / [plotData, common_types]
import docopt, strutils, os, sets, json
import nimhdf5, plotly, chroma
import ingrid / [ingrid_types, tos_helpers]


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
  plotOuterChips <H5file> (<H5file2> | --lhood <file>) [options]

Options:
  --lhood <file>         If given uses the file to extract events, which pass
                         the likelihood cuts instead of all events.
  -h, --help             Show this help
  --version              Show the version number
"""
const doc = docTmpl % [commitHash, currentDate]
const GridColor = color(0.8, 0.8, 0.8, 0.8)
const LegendBg = color(1.0, 1.0, 1.0, 0.5)
const LegendBorder = color(0.6, 0.6, 0.6, 0.6)
const Color1 = color(1.0, 0.0, 102.0 / 256.0)
const Color2 = color(0.0, 153.0 / 256.0, 204 / 256.0)

proc makeOuterChipPlot(fname: string, runType: RunTypeKind,
                       cutRange: CutRange): PlotV =
  var h5f = H5file(fname, "r")
  let fInfo = getFileInfo(h5f)
  let pds = createOuterChipHistograms(h5f,
                                      runType = runType,
                                      fileInfo = fInfo,
                                      cutRange = cutRange)
  let (name, plt) = createPlot(h5f, fInfo, pds[0])
  result = plt

proc makeOuterChipPlotLhood(fname, lhoodFname: string) =
  ## creates a similar plot to the one in `makeOuterChipPlot`, but performs
  ## a cut on the events passing the likelihood cut
  var h5f = H5file(fname, "r")
  let fInfo = getFileInfo(h5f)
  var h5L = H5file(lhoodFname, "r")
  let fLhoodInfo = getFileInfo(h5L, likelihoodGroupGrpStr())
  var data: seq[int]
  for r in fLhoodInfo.runs:
    # iterate all runs, extract events and
    let centerChip = fInfo.centerChip
    let dsetName = (likelihoodBase() & $r) / "chip_" &
      $centerChip / "eventNumber"
    if dsetName notin h5L:
      # this run has no events left, continue
      continue
    let evNumLhood = toSet(h5L[dsetName, int])
    for c in fInfo.chips:
      if c != centerChip:
        let
          evNum = h5f[recoDataChipBase(r) & $c / "eventNumber", int]
          hits = h5f[recoDataChipBase(r) & $c / "hits", int]
        for i, ev in evNum:
          if ev in evNumLhood:
            data.add hits[i]

  histPlot(data)
    .binRange(0.0, 400.0)
    .binSize(10.0)
    #.title("# pixels on outer chips for events passing LogL cut for Run 2 (2017/18)")
    .title("# pixels on outer chips for events passing LogL cut for Run 3 (Oct-Dec 2018)")
    .name("Background data, outer chips # pix") # "Outer chips # pix hit")
    .xlabel("Hits / #")
    .ylabel("#")
    .gridColor(GridColor)
    .markerColor(@[Color2], idx = 0)
    .legendLocation(0.55, 0.95)
    .legendBgColor(LegendBg)
    .legendBorderColor(LegendBorder)
    .legendBorderWidth(1)
    .width(800)
    .height(500)
    .show("outerHits_passingLogLCuts.svg")

proc main =
  # somewhat ugly hack
  plotData.BKind = bPlotly

  let args = docopt(doc)
  let backFname = $args["<H5file>"]
  let calibFname = $args["<H5file2>"]
  let lHoodFname = $args["--lhood"]

  if calibFname != "nil":
    let cutRangeNil = (low: 0.0, high: Inf, name: "noCut")
    let cutRangePhoto = (low: 5.5, high: 6.3, name: "photoPeak")
    let pltBack = makeOuterChipPlot(backFname, rtBackground, cutRangeNil)
    let pltCalib = makeOuterChipPlot(calibFname, rtCalibration, cutRangePhoto)

    var plt = pltCalib.plPlot
    plt = plt.addTrace(pltBack.plPlot.traces[0])
      .title("# pixels on outer chips for 'X-ray like' on center for Run 2 (2017/18)")
      #.title("# pixels on outer chips for 'X-ray like' on center for Run 3 (Oct-Dec 2018)")
      .xlabel("Hits / #")
      .ylabel("#")
      .name("Calibration data, outer chips # pix", idx = 0)
      .name("Background data, outer chips # pix", idx = 1)
      .gridColor(GridColor)
      .legendLocation(0.55, 0.95)
      .markerSize(10)
      .legendBgColor(LegendBg)
      .legendBorderColor(LegendBorder)
      .legendBorderWidth(1)
      .markerColor(@[Color1], idx = 0)
      .markerColor(@[Color2], idx = 1)
      .width(800)
      .height(500)

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
  else:
    makeOuterChipPlotLhood(backFname, lHoodFname)

when isMainModule:
  main()
