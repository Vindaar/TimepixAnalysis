import plotly
import os, strutils, strformat, times, sequtils, math, macros, algorithm, sets
import options, logging, typeinfo, json
# import websocket, asynchttpserver, asyncnet, asyncdispatch

import shell
import arraymancer
import zero_functional
import nimhdf5
import seqmath

import nimpy
import docopt
import chroma
import parsetoml

import dataset_helpers
import ingrid/ingrid_types
import ingrid/tos_helpers
import ingrid/calibration
import helpers/utils
import ingridDatabase / [databaseDefinitions, databaseRead]

type
  # for clarity define a type for the Docopt argument table
  DocoptTab = Table[string, Value]

template canImport(x: untyped): bool =
  compiles:
    import x

when canImport(ingridDatabase):
  # only in this case give option to use plot from database
  import ingridDatabase

when defined(linux):
  const commitHash = staticExec("git rev-parse --short HEAD")
else:
  const commitHash = ""

# get date using `CompileDate` magic
const currentDate = CompileDate & " at " & CompileTime

const docTmpl = """
Version: $# built on: $#
A tool to plot data from H5 files

Usage:
  plotData <H5file> [--runType=<type>] [--backend=<val>] [--no_fadc] [--no_ingrid] [options]

Options:
  --runType=<type>    Select run type (Calib | Back | Xray)
                      The following are parsed case insensetive:
                      Calib = {"calib", "calibration", "c"}
                      Back = {"back", "background", "b"}
                      Xray = {"xray", "xrayfinger", "x"}
  --backend=<val>     Select the plotting backend to be chosen.
                      The followoing backends are available:
                      Python / Matplotlib: {"python", "py", "matplotlib", "mpl"}
                      Nim / Plotly: {"nim", "plotly"}
  --no_fadc           If set no FADC plots will be created.
  --no_ingrid         If set no InGrid plots will be created.
  --no_occupancy      If set no occupancy plots will be created.
  --no_polya          If set no polya plots will be created.
  --no_fe_spec        If set no Fe spectrum will be created.
  -h, --help          Show this help
  --version           Show the version number
"""
const doc = docTmpl % [commitHash, currentDate]

const
  GoldenMean = (sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
  FigWidth = 1200.0                     # width in inches
  FigHeight = FigWidth * GoldenMean     # height in inches

type
  ConfigFlagKind = enum
    cfNone, cfNoFadc, cfNoInGrid, cfNoOccupancy, cfNoPolya, cfNoFeSpectrum

  # enum listing all available `plot types` we can produce
  PlotKind = enum
    pkInGridDset, pkFadcDset, pkPolya, pkCombPolya, pkOccupancy, pkOccCluster,
    pkFeSpec, pkEnergyCalib, pkFeChargeSpec, pkFeVsTime, pkCalibRandom,
    pkAnyScatter, pkMultiDset, pkInGridCluster

  BackendKind = enum
    bNone, bMpl, bPlotly

  Scatter[T] = object
    case kind: BackendKind
    of bMpl:
      x: PyObject
      y: PyObject
      xerr: PyObject
      yerr: PyObject
    of bPlotly:
      # for plotly we store everything in a trace already
      tr: Trace[T]
    else: discard

  # variant object for the layout combining both
  # TODO: make generic or always use float?
  PlotV = object
    case kind: BackendKind
    of bMpl:
      # what needs to go here?
      plt: PyObject
      fig: PyObject
      ax: PyObject
    of bPlotly:
      plLayout: Layout
      plPlot: Plot[float]
    else: discard

  DataKind = enum
    dkInGrid, dkFadc

  ShapeKind = enum
    Rectangle, Square

  ClampKind = enum
    ckFullRange, ckAbsolute, ckQuantile

  # a simple object storing the runs, chips etc. from a given
  # H5 file
  FileInfo = object
    runs: seq[int]
    chips: seq[int]
    runType: RunTypeKind
    rfKind: RunFolderKind
    centerChip: int
    centerChipName: string
    hasFadc: bool # reads if FADC group available
    # NOTE: add other flags for other optional plots?
    # if e.g. FeSpec not available yet, we can just call the
    # procedure to create it for us

# global variable which stores the backend the user selected
var BKind: BackendKind = bNone
# global OrderedSet to store all files we save to later
# combine to single PDF
var imageSet = initOrderedSet[string]()
var plotlyJson = newJObject()

var ShowPlots = false
const PlotlySaveSvg = false

# let server = newAsyncHttpServer()

# create directories, if not exist
if not dirExists("logs"):
  createDir("logs")
if not dirExists("figs"):
  createDir("figs")

# set up the logger
var L = newConsoleLogger()
var fL = newFileLogger("logs/plotData.log", fmtStr = verboseFmtStr)
addHandler(L)
addHandler(fL)


const InGridFnameTemplate = "$1_run$2_chip$3_$4"
const InGridTitleTemplate = "Dataset: $1 for run $2, chip $3 in range: $4"
const FadcFnameTemplate = "fadc_$1_run$2_$3"
const FadcTitleTemplate = "Dataset: $1 for run $2, fadc in range: $3"
const PolyaFnameTemplate = "polya_run$1_chip$2"
const PolyaTitleTemplate = "Polya for run $1 of chip $2"
const CombPolyaFnameTemplate = "combined_polya_run$1"
const CombPolyaTitleTemplate = "Polyas for all chips of run $1"
const OccupancyFnameTemplate = "occupancy_run$1_chip$2_$3"
const OccupancyTitleTemplate = "Occupancy of chip $1 for run $2, $3"
const OccClusterFnameTemplate = "occupancy_clusters_run$1_chip$2_$3"
const OccClusterTitleTemplate = "Occupancy of cluster centers for run $1, chip $2, $3"
const FeSpecFnameTemplate = @["fe_spectrum_run$runNumber",
                              "fe_energy_calib_run$runNumber"]
type
  CutRange = tuple[low, high: float, name: string]

  PlotDescriptor = object
    runType: RunTypeKind
    name: string
    runs: seq[int]
    chip: int
    case plotKind: PlotKind
    of pkInGridDset, pkFadcDset:
      range: CutRange
    of pkOccupancy, pkOccCluster:
      case clampKind: ClampKind
      of ckAbsolute:
        # absolute clamp tp `clampA`
        clampA: float
      of ckQuantile:
        # clamp to `clampQ` quantile
        clampQ: float
      of ckFullRange:
        # no field for ckFullRange
        discard
    of pkCombPolya:
      chipsCP: seq[int]
    else:
      discard

# PlotDescriptor object
# storing all needed information to create a specific plot
func `%`*(pd: PlotDescriptor): JsonNode =
  ## serialize `PlotDescriptor` to Json
  result = newJObject()
  result["runType"] = % pd.runType
  result["name"] = % pd.name
  result["runs"] = % pd.runs
  result["chip"] = % pd.chip
  result["plotKind"] = % pd.plotKind
  case pd.plotKind
  of pkInGridDset, pkFadcDset:
    result["range"] = %* { "low": % pd.range[0],
                           "high": % pd.range[1],
                           "name": % pd.range[2] }
  of pkOccupancy, pkOccCluster:
    result["clampKind"] = % pd.clampKind
    case pd.clampKind
    of ckAbsolute:
      result["clampA"] = % pd.clampA
    of ckQuantile:
      result["clampQ"] = % pd.clampQ
    else: discard
  of pkCombPolya:
    result["chipsCP"] = % pd.chipsCP
  else: discard

func `$`*(pd: PlotDescriptor): string =
  ## use JSON representation as string
  result = (% pd).pretty

func parsePd*(pd: JsonNode): PlotDescriptor =
  ## parse a PlotDescriptor stored as JsonNode back to an object
  result.runType = parseEnum[RunTypeKind](pd["runType"].getStr, rtNone)
  result.name = pd["name"].getStr
  result.runs = to(pd["runs"], seq[int])
  result.chip = pd["chip"].getInt
  result.plotKind = parseEnum[PlotKind](pd["plotKind"].getStr)
  case result.plotKind
  of pkInGridDset, pkFadcDset:
    result.range = (low: pd["range"]["low"].getFloat,
                    high: pd["range"]["high"].getFloat,
                    name: pd["range"]["name"].getStr)
  of pkOccupancy, pkOccCluster:
    result.clampKind = parseEnum[ClampKind](pd["clampKind"].getStr)
    case result.clampKind
    of ckAbsolute:
      result.clampA = pd["clampA"].getFloat
    of ckQuantile:
      result.clampQ = pd["clampQ"].getFloat
    else: discard
  of pkCombPolya:
    result.chipsCP = to(pd["chipsCP"], seq[int])
  else: discard

# karax client
# - static table of
#   Table[PlotDescriptor, string]
#   where the string is an SVG
#   Karax waits to receive plots via sockets after they're
#   created as they are now, read and sent
# - "dynamic" table of
#   Table[PlotDescriptor, JsonNode]
#   where the JsonNode is the output of plotly before pasting it
#   into the HTML template
# Upon startup server sends Karax client the FileInfo, which should
# be enough to have basic information to request plots
# - can provide fields to enter information needed for PlotDescriptor
#   where all information is not from a text field, but rather a dropdown
#   menu? With FileInfo we have everything needed (in addition to some
#   definitions like dataset names, which we will store in a TOML file

proc jsonPlotly(pltV: PlotV): JsonNode =
  ## returns JsonNode from the PlotV
  result = newJObject()
  case pltV.kind
  of bPlotly:
    let trSeq = pltV.plPlot.traces.mapIt(% it)
    result["Traces"] = % trSeq
    result["Layout"] = % pltV.plPlot.layout
  else:
    warn "Unsuitable type for Json " & $(pltV.kind)

proc savePlot(p: PlotV, outfile: string, fullPath = false) =
  var fname = outfile
  if not fullPath:
    fname = "figs" / outfile & ".svg"
  info &"Saving file: {fname}"
  case BKind
  of bPlotly:
    if not PlotlySaveSvg:
      plotlyJson[outfile] = jsonPlotly(p)
    else:
      p.plPlot.saveImage(fname)
      imageSet.incl(fname)
  of bMpl:
    discard p.plt.savefig(fname)
    imageSet.incl(fname)
    if ShowPlots:
      discard callMethod(p.plt, "show")
    else:
      discard callMethod(p.plt, "close", "all")
  else: discard

proc importPyplot(): PyObject =
  let mpl = pyImport("matplotlib")
  discard mpl.use("TKagg")
  result = pyImport("matplotlib.pyplot")

proc initPlotV(title: string, xlabel: string, ylabel: string, shape = ShapeKind.Rectangle): PlotV =
  ## returns the plotly layout of the given shape
  var width: int
  case shape
  of ShapeKind.Rectangle:
    width = FigWidth.int
  of ShapeKind.Square:
    # for square take height as width
    width = FigHeight.int
  case BKind
  of bPlotly:
    let plLayout = Layout(title: title,
                          width: width, height: FigHeight.int,
                          xaxis: Axis(title: xlabel),
                          yaxis: Axis(title: ylabel),
                          autosize: false)
    result = PlotV(kind: bPlotly,
                   plLayout: plLayout)
  of bMpl:
    let plt = importPyplot()
    let (fig, ax) = plt.subplots(1, 1).to((PyObject, PyObject))
    let dpi = fig.get_dpi()
    discard fig.set_size_inches(width.float / dpi.to(float), FigHeight.float / dpi.to(float))
    discard ax.set_title(title)
    discard ax.set_xlabel(xlabel)
    discard ax.set_ylabel(ylabel)
    result = PlotV(kind: bMpl,
                   plt: plt,
                   fig: fig,
                   ax: ax)
  else: discard

proc read(h5f: var H5FileObj,
          runNumber: int,
          dsetName: string,
          chipNumber = 0,
          isFadc = false,
          dtype: typedesc = float): seq[dtype] =
  ## reads the given dataset from the H5 file and returns it
  ## (potentially) converted to ``dtype``
  var dset: H5DataSet
  if isFadc:
    dset = h5f[(fadcRunPath(runNumber) / dsetName).dset_str]
  else:
    dset = h5f[(recoPath(runNumber, chipNumber).string / dsetName).dset_str]
  when dtype is SomeNumber:
    let convert = dset.convertType(dtype)
    result = dset.convert
  else:
    type subtype = utils.getInnerType(dtype)
    # NOTE: only support 2D seqs
    let convert = dset.convertType(subtype)
    when subType isnot SomeNumber:
      raise newException(Exception, "Cannot convert N-D sequence to type: " &
        subtype.name)
    else:
      # convert to subtype and reshape to dsets shape
      result = dset.convert.reshape2D(dset.shape)

proc getFileInfo(h5f: var H5FileObj): FileInfo =
  ## returns a set of all run numbers in the given file
  # virist file
  h5f.visitFile()
  var readAux = false
  # get reconstruction group
  let reco = h5f[recoGroupGrpStr()]
  result.runType = parseEnum[RunTypeKind](reco.attrs["runType", string], rtNone)
  result.rfKind = parseEnum[RunFolderKind](reco.attrs["runFolderKind", string],
                                           rfUnknown)
  result.centerChip = reco.attrs["centerChip", int]
  result.centerChipName = reco.attrs["centerChipName", string]

  # get allowed chips from toml
  let tomlConfig = parseToml.parseFile("config.toml")
  let allowedChips = tomlConfig["General"]["allowedChips"].getElems
  var chipSet: set[uint16]
  for el in allowedChips:
    chipSet.incl el.getInt.uint16

  for runNumber, group in runs(h5f):
    result.runs.add runNumber.parseInt
    if not readAux:
      let grp = h5f[group.grp_str]
      let nChips = grp.attrs["numChips", int]
      result.chips = toSeq(0 ..< nChips)#.mapIt(it)
      readAux = true

  result.chips = result.chips.filterIt(it.uint16 in chipSet)
  # sort the run numbers
  result.runs.sort
  echo result

proc plotHist[T](xIn: seq[seq[T]], title, dset, outfile: string) =
  ## plots the data in `x` as a histogram
  let xs = xIn.mapIt(it.mapIt(it.float))
  let binRangeO = getBinRangeForDset(dset)
  let nbinsO = getNumBinsForDset(dset)
  var binRange: tuple[start, stop: float]
  var nbins: int
  if binRangeO.isSome:
    binRange = get(binRangeO)
  else:
    binRange = (xs[0].percentile(5), xs[0].percentile(95))
  if nBinsO.isSome:
    nBins = get(nbinsO)
  else:
    nBins = 100
  info &"Bin range {binRange} for dset: {dset}"
  var pltV = initPlotV(title, dset, "#")
  case BKind
  of bPlotly:
    var traces: seq[Trace[float]]
    for x in xs:
      let binSize = (binRange[1] - binRange[0]) / nbins.float
      traces.add Trace[float](`type`: PlotType.Histogram,
                              bins: binRange,
                              binSize: binSize,
                              #nbins: nBins,
                              xs: x,
                              name: dset)
    pltV.plPlot = Plot[float](layout: pltV.plLayout, traces: traces)
    pltV.savePlot(outfile, fullPath = true)
  of bMpl:
    for x in xs:
      discard pltV.ax.hist(x,
                            bins = nbins,
                            range = binRange)
    pltV.savePlot(outfile, fullPath = true)
  else: discard

proc plotBar[T](binsIn, countsIn: seq[seq[T]], title: string,
                xlabel: string, dsets: seq[string],
                outfile: string, drawPlots = false): PlotV =
  ## plots the data in `x` as a histogram
  let bins = binsIn.mapIt(it.mapIt(it.float))
  let counts = countsIn.mapIt(it.mapIt(it.float))
  result = initPlotV(title, xlabel, "#")
  var traces: seq[Trace[float]]
  for i in 0 .. bins.high:
    case BKind
    of bPlotly:
      let
        bin = bins[i]
        count = counts[i]
        dset = dsets[i]
      traces.add Trace[float](`type`: PlotType.Bar,
                              xs: bin,
                              ys: count,
                              name: dset)
      result.plPlot = Plot[float](layout: result.plLayout, traces: traces)
      if drawPlots:
        result.plPlot.show()
    of bMpl:
      let
        # cut of last bin edge in matplotlib case. Expects same shape for edges
        # and counts
        bin = bins[i][0 .. ^2]
        count = counts[i]
        dset = dsets[i]
      let width = toSeq(0 ..< bins[i].high).mapIt(bins[i][it + 1] - bins[i][it])
      discard result.ax.bar(bin,
                            count,
                            align = "edge",
                            width = width,
                            label = dset)
      # we cannot call `show` directly, because the Nim compiler tries to
      # call the Nim `show` procedure and fails
      if drawPlots:
        discard callMethod(result.plt, "show")

    else: discard

iterator chips(group: var H5Group): (int, H5Group) =
  ## returns all chip groups within the given Run group
  doAssert group.name.parentDir == "/reconstruction", "Bad group : " & group.name
  doAssert "run_" in group.name, "Bad group : " & group.name
  for grp in groups(group):
    if "chip_" in grp.name:
      let chipNum = grp.attrs["chipNumber", int]
      yield (chipNum, grp)

proc histograms(h5f: var H5FileObj, runType: RunTypeKind,
                fileInfo: FileInfo,
                flags: set[ConfigFlagKind]): seq[PlotDescriptor] =
  const dsets = ["length", "width", "skewnessLongitudinal", "skewnessTransverse",
                 "kurtosisLongitudinal", "kurtosisTransverse", "rotationAngle",
                 "eccentricity", "fractionInTransverseRms", "lengthDivRmsTrans"]
  const fadcDsets = ["minvals", "fallTime", "riseTime"]
  # TODO: perform cut on photo peak and escape peak, done by getting the energy
  # and performing a cut around peak +- 300 eV maybe
  # NOTE: need photo peak and escape peak plenty of times here. Just write wrapper
  # which takes energy and performs cuts, *then* creates plots
  var ranges: seq[CutRange]
  case runType
  of rtCalibration:
    # creates 3 plots for the given dataset.
    # - 1 around photo peak
    # - 1 around escape peak
    # - 1 without cut
    ranges = @[(low: 5.5, high: 6.2, name: "Photopeak"),
               (low: 2.7, high: 3.2, name: "Escapepeak"),
               (low: -Inf, high: Inf, name: "All")]
  else:
    # else just take all
    ranges = @[(low: -Inf, high: Inf, name: "All")]

  for r in ranges:
    if cfNoInGrid notin flags:
      for ch in fileInfo.chips:
        for dset in dsets:
          result.add PlotDescriptor(runType: runType,
                                    name: dset,
                                    runs: fileInfo.runs,
                                    plotKind: pkInGridDset,
                                    chip: ch,
                                    range: r)
    if cfNoFadc notin flags:
      for dset in fadcDsets:
        result.add PlotDescriptor(runType: runType,
                                  name: dset,
                                  runs: fileInfo.runs,
                                  plotKind: pkFadcDset,
                                  range: r)

proc calcOccupancy[T](x, y: seq[T]): Tensor[float] =
  ## calculates the occupancy of the given x and y datasets
  ## Either for a `seq[seq[T: SomeInteger]]` in which case we're calculating
  ## the occupancy of a raw clusters or `seq[T: SomeFloat]` in which case
  ## we're dealing with center positions of clusters
  result = newTensor[float]([256, 256])
  # iterate over events
  for i in 0 .. x.high:
    let
      xEv = x[i]
      yEv = y[i]
    when T is seq:
      for j in 0 .. xEv.high:
        result[xEv[j].int, yEv[j].int] += 1.0
    elif T is SomeFloat:
      result[xEv.round.int, yEv.round.int] += 1.0

proc plotHist2D(data: Tensor[float], title, outfile: string) = #descr: string) =
  ## creates a 2D histogram plot (basically an image) of the given 2D
  ## Tensor
  doAssert data.rank == 2
  var pltV = initPlotV(title, "pixel x", "pixel y", ShapeKind.Square)
  case BKind
  of bPlotly:
    let tr = Trace[float](`type`: PlotType.HeatMap,
                      colormap: ColorMap.Viridis,
                      zs: data.toRawSeq.reshape2D([256, 256]))
    pltV.plPlot = Plot[float](layout: pltV.plLayout, traces: @[tr])
    pltV.savePlot(outfile, fullPath = true)
  of bMpl:
    discard pltV.plt.imshow(data.toRawSeq.reshape([256, 256]),
                       cmap = "viridis")
    discard pltV.plt.colorbar()
    pltV.savePlot(outfile, fullPath = true)
  else:
    discard

proc plotScatter(pltV: PlotV, x, y: seq[float], name, outfile: string) =
  case BKind
  of bPlotly:
    # in this case plot is defined
    let trFit = Trace[float](`type`: PlotType.Scatter,
                             xs: x,
                             ys: y,
                             name: name)
    pltV.plPlot.traces.add trFit
    pltV.savePlot(outfile, fullPath = true)
  of bMpl:
    # in this case `ax` is defined
    discard pltV.ax.plot(x,
                         y,
                         label = name,
                         color = "r")
    pltV.savePlot(outfile, fullPath = true)
  else:
    warn &"Unsupported backend kind: {BKind}"

proc plotScatter(x, y: seq[float], title, name, outfile: string) =
  ## wrapper around the above if no `PlotV` yet defined
  var pltV = initPlotV(title, "x", "y", ShapeKind.Rectangle)
  pltV.plotScatter(x, y, name, outfile)

# TODO: also plot occupancies without full frames (>4095 hits)?
proc occupancies(h5f: var H5FileObj, runType: RunTypeKind,
                 fileInfo: FileInfo,
                 flags: set[ConfigFlagKind]): seq[PlotDescriptor] =
  ## creates occupancy plots for the given HDF5 file, iterating over
  ## all runs and chips
  const
    quantile = 95
    clamp3 = 10.0
  for ch in fileInfo.chips:
    let basePd = PlotDescriptor(runType: runType,
                                name: "occupancy",
                                runs: fileInfo.runs,
                                plotKind: pkOccupancy,
                                chip: ch)
    let fullPd = replace(basePd):
      clampKind = ckFullRange
    let clampAPd = replace(basePd):
      clampKind = ckAbsolute
      clampA = clamp3
    let clampQPd = replace(basePd):
      clampKind = ckQuantile
      clampQ = quantile
    let clusterPd = replace(basePd):
      plotKind = pkOccCluster
      # only consider full range for cluster centers
      clampKind = ckFullRange
    let clusterClampPd = replace(basePd):
      plotKind = pkOccCluster
      # only consider full range for cluster centers
      clampKind = ckQuantile
      clampQ = 80
    result.add @[fullPd, clampAPd, clampQPd, clusterPd, clusterClampPd]
  echo result

proc plotPolyas(h5f: var H5FileObj, group: H5Group,
                chipNum: int, runNumber: string): seq[Trace[float]] =
  ## perform the plots of the polya distribution again (including the fits)
  ## and return the polya data as to allow a combined plot afterwards
  let
    polya = h5f.read(runNumber.parseInt, "polya", chipNum,
                     dtype = seq[float])
    polyaFit = h5f.read(runNumber.parseInt, "polyaFit", chipNum,
                        dtype = seq[float])
  let nbins = polya.shape[0] - 1
  var
    bins = newSeq[float](nbins + 1)
    binCenterFit = newSeq[float](nbins)
    counts = newSeq[float](nbins)
    countsFit = newSeq[float](nbins)
  for i in 0 .. polya.high:
    bins[i] = polya[i][0]
    if i != nbins:
      # do not take last element of polya datasets, since it's just 0
      counts[i] = polya[i][1]
      binCenterFit[i] = polyaFit[i][0]
      countsFit[i] = polyaFit[i][1]
  let names = @[&"Polya of chip {chipNum} for run {runNumber}"]
  let title = &"Polya distribution of chip {chipNum} for run {runNumber}"
  let xlabel = "Number of electrons"
  let outfile = &"polya_run{runNumber}_chip{chipNum}"
  var pltV = plotBar(@[bins], @[counts], title, xlabel, names, outfile)
  # now add scatter plot for fit result
  let nameFit = &"Polya fit of chip {chipNum} for run {runNumber}"
  plotScatter(binCenterFit, countsFit, title, nameFit, outfile)

proc polya(h5f: var H5FileObj, runType: RunTypeKind,
           fileInfo: FileInfo,
           flags: set[ConfigFlagKind]): seq[PlotDescriptor] =
  ## creates the plots for the polya distribution of the datasets
  for ch in fileInfo.chips:
    result.add PlotDescriptor(runType: runType,
                              name: "polya",
                              runs: fileInfo.runs,
                              chip: ch,
                              plotKind: pkPolya)
  if fileInfo.chips.len > 1:
    result.add PlotDescriptor(runType: runType,
                              name: "polya",
                              runs: fileInfo.runs,
                              plotKind: pkCombPolya,
                              chipsCP: fileInfo.chips)

  #for runNumber, grpName in runs(h5f):
  #  var group = h5f[grpName.grp_str]
  #  var traces: seq[Trace[float]]
  #  for chipNum, chpgrp in chips(group):
  #    traces.add plotPolyas(h5f, chpgrp, chipNum, runNumber)
  #  case BKind
  #  of bPlotly:
  #    let title = &"Polyas for all chips of run {runNumber}"
  #    let xlabel = "Number of electrons"
  #    var pltV = initPlotV(title = title, xlabel = xlabel, ylabel = "#")
  #    pltV.plPlot = Plot[float](layout: pltV.plLayout,
  #                              traces: traces)
  #    let outfile = &"combined_polya_{runType}_run{runNumber}"
  #    pltV.savePlot(outfile)
  #  else:
  #    warn "Combined polya only available for Plotly backend!"

proc plotDates[T, U](x: seq[U], y: seq[T],
                  title, xlabel, dsets, outfile: string,
                  ylabel = "",
                  drawPlots = false) =
  var pltV = initPlotV(title, xlabel, ylabel, ShapeKind.Rectangle)
  case BKind
  of bPlotly:
    let tr = Trace[T](mode: PlotMode.LinesMarkers,
                      `type`: PlotType.Scatter,
                      xs: x,
                      ys: y,
                      marker: Marker[float](size: @[12.0]))
    pltV.plPlot = Plot[T](layout: pltV.plLayout, traces: @[tr])
    pltV.savePlot(outfile)
  of bMpl:
    let dtm = pyImport("datetime")
    let xP = x.mapIt(dtm.datetime.utcfromtimestamp(int(it)))
    discard pltV.ax.plot_date(xP, y, label = title)
    pltV.savePlot(outfile)
  else:
    discard

proc feSpectrum(h5f: var H5FileObj, flags: set[ConfigFlagKind]) =
  ## creates the plot of the Fe spectrum
  const kalphaPix = 10
  const kalphaCharge = 4
  const parPrefix = "p"
  const dateStr = "yyyy-MM-dd'.'HH:mm:ss" # example: 2017-12-04.13:39:45
  var
    pixSeq: seq[float]
    chSeq: seq[float]
    dates: seq[float] #string]#Time]
  for runNumber, grpName in runs(h5f):
    var group = h5f[grpName.grp_str]
    let centerChip = h5f[group.parent.grp_str].attrs["centerChip", int]
    # call `importPyplot` because it sets the appropriate backend for us
    discard importPyplot()
    let path = "figs/"
    var outfiles = @[&"fe_spectrum_run{runNumber}.svg",
                     &"fe_energy_calib_run{runNumber}.svg"]
    let outfilesCh = outfiles.mapIt(path / ("charge_" & it))
    outfiles = outfiles.mapIt(path / it)
    for o in outfiles:
      imageSet.incl o
    for o in outfilesCh:
      imageSet.incl o
    h5f.fitToFeSpectrum(runNumber.parseInt, centerChip,
                        fittingOnly = not ShowPlots,
                        outfiles = outfiles,
                        writeToFile = false)
    # extract fit parameters from center chip group
    let chpGrpName = group.name / "chip_" & $centerChip
    #let feSpec =
    pixSeq.add h5f[(chpGrpName / "FeSpectrum").dset_str].attrs[
      parPrefix & $kalphaPix, float
    ]
    chSeq.add h5f[(chpGrpName / "FeSpectrumCharge").dset_str].attrs[
      parPrefix & $kalphaCharge, float
    ]
    dates.add parseTime(group.attrs["dateTime", string],
                        dateStr,
                        utc()).toUnix.float

  # now plot
  # calculate ratio and convert to string to workaround plotly limitation of
  # only one type for Trace
  let ratio = zip(pixSeq, chSeq) --> map(it[0] / it[1])
  plotDates(dates, ratio,
            title = "Photopeak pix / charge vs time",
            xlabel = "Date",
            dsets = "a",
            outfile = "photopeak_vs_time")
            #ylabel = "# pix / charge in e^-")

func cKindStr(pd: PlotDescriptor, sep: string): string =
  case pd.clampKind
  of ckAbsolute:
    result = &"{pd.clampKind}{sep}{pd.clampA}"
  of ckQuantile:
    result = &"{pd.clampKind}{sep}{pd.clampQ}"
  of ckFullRange:
    result = &"{pd.clampKind}"

proc buildOutfile(pd: PlotDescriptor): string =
  var name = ""
  let runsStr = pd.runs.foldl($a & "_" & $b, "").strip(chars = {'_'})
  case pd.plotKind
  of pkInGridDset:
    name = InGridFnameTemplate % [pd.name,
                                  runsStr,
                                  $pd.chip,
                                  $pd.range[2]]
  of pkFadcDset:
    name = FadcFnameTemplate % [pd.name,
                                runsStr,
                                $pd.range[2]]
  of pkOccupancy:
    let clampStr = cKindStr(pd, "_")
    name = OccupancyFnameTemplate % [runsStr,
                                     $pd.chip,
                                     clampStr]
  of pkOccCluster:
    let clampStr = cKindStr(pd, "_")
    name = OccClusterFnameTemplate % [runsStr,
                                      $pd.chip,
                                      clampStr]
  of pkPolya:
    name = PolyaFnameTemplate % [runsStr,
                                 $pd.chip]
  of pkCombPolya:
    name = CombPolyaFnameTemplate % [runsStr]
  else:
    discard
  result = "figs" / (name & ".svg")

proc buildTitle(pd: PlotDescriptor): string =
  let runsStr = pd.runs.foldl($a & " " & $b, "").strip(chars = {' '})
  case pd.plotKind
  of pkInGridDset:
    result = InGridTitleTemplate % [pd.name,
                                    runsStr,
                                    $pd.chip,
                                    $pd.range[2]]
  of pkFadcDset:
    result = FadcTitleTemplate % [pd.name,
                                  runsStr,
                                  $pd.range[2]]
  of pkOccupancy:
    let clampStr = cKindStr(pd, "@")
    result = OccupancyTitleTemplate % [runsStr,
                                       $pd.chip,
                                       clampStr]
  of pkOccCluster:
    let clampStr = cKindStr(pd, "@")
    result = OccClusterTitleTemplate % [runsStr,
                                       $pd.chip,
                                       clampStr]
  of pkPolya:
    result = PolyaTitleTemplate % [runsStr,
                                   $pd.chip]
  of pkCombPolya:
    result = CombPolyaTitleTemplate % [runsStr]
  else:
    discard

func clampedOccupancy[T](x, y: seq[T], pd: PlotDescriptor): Tensor[float] =
  ## calculates the occupancy given `x`, `y`, which may be seqs of clusters
  ## or seqs of center positions
  result = calcOccupancy(x, y)
  case pd.clampKind
  of ckAbsolute:
    result = result.clamp(0.0, pd.clampA)
  of ckQuantile:
    let quant = result.toRawSeq.percentile(pd.clampQ.round.int)
    result = result.clamp(0.0, quant)
  else: discard

proc readPolya(h5f: var H5FileObj, pd: PlotDescriptor):
     (seq[float], seq[float], seq[float], seq[float]) =
  ## reads the `polya` and `polyaFit` datasets for all runs in
  ## `pd`, stacks it and returns bins, counts for both
  doAssert (pd.plotKind == pkPolya or pd.plotKind == pkCombPolya),
     &"Only supported for {pkPolya} and {pkCombPolya}. This plotKind " &
     &"is {pd.plotKind}"
  var
    lastBins: seq[float]
    lastBCFit: seq[float]
    countsFull: seq[float]
    countsFitFull: seq[float]
  for r in pd.runs:
    # TODO: replace this after rewriting `calcGasGain` in calibration.nim
    # There: disentangle data reading from the fitting data processing etc.
    # Then we can read raw data and hand that to the fitting proc here
    let
      polya = h5f.read(r, pd.name, pd.chip,
                       dtype = seq[float])
      polyaFit = h5f.read(r, pd.name & "Fit", pd.chip,
                          dtype = seq[float])
    let nbins = polya.shape[0] - 1
    var
      bins = newSeq[float](nbins + 1)
      binCenterFit = newSeq[float](nbins)
      counts = newSeq[float](nbins)
      countsFit = newSeq[float](nbins)
    for i in 0 .. polya.high:
      bins[i] = polya[i][0]
      if i != nbins:
        # do not take last element of polya datasets, since it's just 0
        counts[i] = polya[i][1]
        binCenterFit[i] = polyaFit[i][0]
        countsFit[i] = polyaFit[i][1]
    let diffR = zip(lastBins, bins) --> map(it[0] - it[1])
    if lastBins.len > 0:
      doAssert lastBins == bins, "The ToT calibration changed between the " &
        "last two runs!"
      doAssert lastBCFit == binCenterFit
    lastBins = bins
    lastBCFit = binCenterFit
    if countsFull.len > 0:
      doAssert countsFull.len == counts.len
      countsFull = zip(countsFull, counts) --> map(it[0] + it[1])
      countsFitFull = zip(countsFitFull, countsFit) --> map(it[0] + it[1])
    else:
      countsFull = counts
      countsFitFull = countsFit
  result = (lastBins, countsFull, lastBCFit, countsFitFull)

proc createPlot(h5f: var H5FileObj,
                fileInfo: FileInfo,
                pd: PlotDescriptor): string =
  ## creates a plot of kind `plotKind` for the data from all runs in `runs`
  ## for chip `chip`
  case pd.plotKind
  of pkInGridDset:
    let ranges = @[pd.range]
    var allData: seq[float]
    for r in pd.runs:
      let data = h5f.read(r, pd.name, pd.chip, dtype = float)
      # perform cut on range
      let group = h5f[recoPath(r, pd.chip)]
      let idx = cutOnProperties(h5f, group,
                      ("energyFromCharge", pd.range[0], pd.range[1]))
      allData.add idx.mapIt(data[it])
    result = buildOutfile(pd)
    let title = buildTitle(pd)
    plotHist(@[allData], title, pd.name, result)
  of pkFadcDset:
    # get the center chip group
    var allData: seq[float]
    for r in pd.runs:
      let group = h5f[recoPath(r, fileInfo.centerChip)]
      let idx = cutOnProperties(h5f, group,
                    ("energyFromCharge", pd.range[0], pd.range[1]))
      let evNumInGrid = h5f.read(r, "eventNumber", fileInfo.centerChip, dtype = int)
      # filter out correct indices passing cuts
      var inGridSet = initSet[int]()
      for i in idx:
        inGridSet.incl evNumInGrid[i]
      let evNumFadc = h5f.read(r, "eventNumber", isFadc = true, dtype = int) # [group.name / "eventNumber", int]
      let idxFadc = (toSeq(0 .. evNumFadc.high)) --> filter(evNumFadc[it] in inGridSet)
      let data = h5f.read(r, pd.name, pd.chip, isFadc = true, dtype = float)
      allData.add idxFadc --> map(data[it])
    result = buildOutfile(pd)
    let title = buildTitle(pd)
    plotHist(@[allData], title, pd.name, result)
  of pkOccupancy:
    # get x and y datasets, stack and get occupancies
    let vlenDtype = special_type(uint8)
    var occFull = newTensor[float]([256, 256])
    for r in pd.runs:
      let
        group = h5f[recoPath(r, pd.chip)]
        xGroup = h5f[(group.name / "x").dset_str]
        yGroup = h5f[(group.name / "y").dset_str]
        x = xGroup[vlenDtype, uint8]
        y = yGroup[vlenDtype, uint8]
      let occ = clampedOccupancy(x, y, pd)
      # stack this run onto the full data tensor
      occFull = occFull .+ occ
    let title = buildTitle(pd)
    result = buildOutfile(pd)
    plotHist2D(occFull, title, result)
  of pkOccCluster:
    # plot center positions
    var occFull = newTensor[float]([256, 256])
    for r in pd.runs:
      let
        group = h5f[recoPath(r, pd.chip)]
        # get centers and rescale to 256 max value
        centerX = h5f[(group.name / "centerX"), float].mapIt(it * 256.0 / 14.0)
        centerY = h5f[(group.name / "centerY"), float].mapIt(it * 256.0 / 14.0)
      let occ = clampedOccupancy(centerX, centerY, pd)
      # stack this run onto the full data tensor
      occFull = occFull .+ occ
    let title = buildTitle(pd)
    result = buildOutfile(pd)
    plotHist2D(occFull, title, result)
  of pkPolya:
    let (bins, counts, binsFit, countsFit) = h5f.readPolya(pd)
    let xlabel = "Number of electrons"
    let title = buildTitle(pd)
    result = buildOutfile(pd)
    var pltV = plotBar(@[bins], @[counts], title, xlabel, @[title], result)
    # now add fit to the existing plot
    let nameFit = &"Polya fit of chip {pd.chip}"
    pltV.plotScatter(binsFit, countsFit, nameFit, result)
  of pkCombPolya:
    var
      binsSeq: seq[seq[float]]
      countsSeq: seq[seq[float]]
      dsets: seq[string]
    let title = buildTitle(pd)
    result = buildOutfile(pd)
    for ch in pd.chipsCP:
      # get a local PlotDescriptor, which has this chip number
      let localPd = block:
        var tmp = pd
        tmp.chip = ch
        tmp
      let (bins, counts, binsFit, countsFit) = h5f.readPolya(localPd)
      binsSeq.add bins
      countsSeq.add counts
      dsets.add "Chip " & $ch
    let xlabel = "Number of electrons"
    var pltV = plotBar(binsSeq, countsSeq, title, xlabel, dsets, result)
    # now add fit to the existing plot
    pltV.savePlot(result, fullPath = true)
  else:
    discard

proc createOrg(outfile: string) =
  ## creates a simple org file consisting of headings and images
  ## SVGs are implemented using raw inline SVG
  const header = """
* $1

"""
  const tmpl = """
** $1

#+BEGIN_EXPORT html
$2
#+END_EXPORT

"""
  # now build pdf
  var orgStr = header % [outfile]
  for im in imageSet:
    info &"Adding image {im}"
    let (dir, name, ext) = im.splitFile()
    let data = readFile(im)
    orgStr = orgStr & tmpl % [name, data]
  var f = open(outfile, fmWrite)
  f.write(orgStr)
  f.close()

  shell:
    emacs `$outfile` "--batch -f org-html-export-to-html --kill"

proc jsonDump(outfile: string) =
  ## reads either the created images from the `imageSet` and dumps them
  ## contained in a JsonNode to a file.
  ## Additionally stores the plotly plots

  var jdump = newJObject()
  jdump["svg"] = newJObject()
  var svgJ = newJObject()
  for im in imageSet:
    svgJ[im] = % readFile(im)

  jdump["svg"] = svgJ
  #echo jdump.pretty
  #echo "\n\n\n\n"
  #echo "plt ", plotlyJson.pretty
  jdump["plotly"] = plotlyJson
  echo jdump.len
  var f = open(outfile, fmWrite)
  var outstr = ""
  outstr.toUgly(jdump)
  f.write(outstr)
  f.close()

proc createCalibrationPlots(h5file: string,
                            bKind: BackendKind,
                            runType: RunTypeKind,
                            flags: set[ConfigFlagKind]) =
  ## creates QA plots for calibration runs
  var h5f = H5file(h5file, "r")
  let fileInfo = getFileInfo(h5f)
  # var imageSet = initOrderedSet[string]()
  var pds: seq[PlotDescriptor]

  const length = "length"
  if cfNoOccupancy notin flags:
    pds.add occupancies(h5f, runType, fileInfo, flags) # plus center only
  if cfNoPolya notin flags:
    pds.add polya(h5f, runType, fileInfo, flags)
  if cfNoFeSpectrum notin flags:
    feSpectrum(h5f, flags)
  # energyCalib(h5f) # ???? plot of gas gain vs charge?!
  pds.add histograms(h5f, runType, fileInfo, flags) # including fadc

  for p in pds:
    let fileF = h5f.createPlot(fileInfo, p)
    case BKind
    of bPlotly:
      if PlotlySaveSvg:
        imageSet.incl fileF
    else: discard
  echo "Image set is ", imageSet.card

  # likelihoodHistograms(h5f) # need to cut on photo peak and esc peak
  # neighborPixels(h5f)
  discard h5f.close()
  var outfile = "calibration"
  for fl in flags:
    outfile &= "_" & $fl
  #outfile &= ".org"
  outfile &= ".json"
  jsonDump(outfile)
  #createOrg(outfile)

proc createBackgroundPlots(h5file: string,
                           bKind: BackendKind,
                           runType: RunTypeKind,
                           flags: set[ConfigFlagKind]) =
  ## creates QA plots for calibration runs
  var h5f = H5file(h5file, "r")
  let fileInfo = getFileInfo(h5f)
  var pds: seq[PlotDescriptor]
  const length = "length"
  if cfNoOccupancy notin flags:
    #occupancies(h5f, flags) # plus center only
    pds.add occupancies(h5f, runType, fileInfo, flags) # plus center only
  if cfNoPolya notin flags:
    pds.add polya(h5f, runType, fileInfo, flags)
  # energyCalib(h5f) # ???? plot of gas gain vs charge?!
  pds.add histograms(h5f, runType, fileInfo, flags) # including fadc
  echo pds
  # likelihoodHistograms(h5f) # need to cut on photo peak and esc peak
  # neighborPixels(h5f)
  discard h5f.close()
  var outfile = "background"
  for fl in flags:
    outfile &= "_" & $fl
  #outfile &= ".org"
  outfile &= ".json"
  jsonDump(outfile)
  #createOrg(outfile)

proc createXrayFingerPlots(bKind: BackendKind, flags: set[ConfigFlagKind]) =
  discard

proc parseBackendType(backend: string): BackendKind =
  ## given a string describing a run type, return the correct
  ## `RunTypeKind`
  if backend.normalize in ["python", "py", "matplotlib", "mpl"]:
    result = bMpl
  elif backend.normalize in ["nim", "plotly"]:
    result = bPlotly
  else:
    result = bNone

proc plotData*() =
  ## the main workhorse of the server end
  let args = docopt(doc)
  info &"Received arguments:\n  {args}"
  let h5file = $args["<H5file>"]
  let runTypeStr = $args["--runType"]
  let backendStr = $args["--backend"]
  var flags: set[ConfigFlagKind]
  if $args["--no_fadc"] == "true":
    flags.incl cfNoFadc
  if $args["--no_ingrid"] == "true":
    flags.incl cfNoIngrid
  if $args["--no_occupancy"] == "true":
    flags.incl cfNoOccupancy
  if $args["--no_polya"] == "true":
    flags.incl cfNoPolya
  if $args["--no_fe_spec"] == "true":
    flags.incl cfNoFeSpectrum
  info "Flags are:\n  {flags}"

  var runType: RunTypeKind
  var bKind: BackendKind
  if runTypeStr != "nil":
    runType = parseRunType(runTypeStr)
  if backendStr != "nil":
    BKind = parseBackendType(backendStr)

  case runType
  of rtCalibration:
    createCalibrationPlots(h5file, bKind, runType, flags)
  of rtBackground:
    createBackgroundPlots(h5file, bKind, runType, flags)
  of rtXrayFinger:
    createXrayFingerPlots(bKind, flags)
  else:
    discard

# proc cb(req: Request) {.async.} =
#   echo "cb"

# proc main =
#   var thr: Thread[void]
#   thr.createThread(plotData)

#   waitFor server.serve(Port(8080), )

when isMainModule:
  plotData()
