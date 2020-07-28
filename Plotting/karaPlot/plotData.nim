import plotly, ggplotnim
import os except FileInfo
import strutils, strformat, times, sequtils, math, macros, algorithm, sets, stats, base64
import options, logging, typeinfo, json
import websocket, asynchttpserver, asyncnet, asyncdispatch

import shell
import arraymancer
import zero_functional
import nimhdf5
import seqmath

import nimpy
import docopt
import chroma
import parsetoml
# import nimdata

import dataset_helpers
import ingrid/ingrid_types
import ingrid/tos_helpers
import ingrid / calibration
import ingrid / calibration / calib_fitting
import helpers/utils
import ingridDatabase / [databaseDefinitions, databaseRead]
import protocol, common_types

type
  # for clarity define a type for the Docopt argument table
  DocoptTab = Table[string, docopt.Value]

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
  plotData <H5file> [--runType=<type>] [--eventDisplay=<run>] [--server] [--backend=<val>] [options]

Options:
  --runType=<type>       Select run type (Calib | Back | Xray)
                         The following are parsed case insensetive:
                         Calib = {"calib", "calibration", "c"}
                         Back = {"back", "background", "b"}
                         Xray = {"xray", "xrayfinger", "x"}
  --backend=<val>        Select the plotting backend to be chosen.
                         The followoing backends are available:
                         Python / Matplotlib: {"python", "py", "matplotlib", "mpl"}
                         Nim / Plotly: {"nim", "plotly"}
  --eventDisplay=<run>   If given will show event displays of the given run.
  --server               If flag given, will launch client and send plots individually,
                         instead of creating all plots and dumping them.
  --no_fadc              If set no FADC plots will be created.
  --no_ingrid            If set no InGrid plots will be created.
  --no_occupancy         If set no occupancy plots will be created.
  --no_polya             If set no polya plots will be created.
  --no_fe_spec           If set no Fe spectrum will be created.
  -h, --help             Show this help
  --version              Show the version number
"""
const doc = docTmpl % [commitHash, currentDate]

const
  GoldenMean = (sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
  FigWidth = 800.0                      # width in pixels
  FigHeight = FigWidth * GoldenMean     # height in pixels

  InGridDsets = ["length", "width", "skewnessLongitudinal", "skewnessTransverse",
                 "kurtosisLongitudinal", "kurtosisTransverse", "rotationAngle",
                 "eccentricity", "fractionInTransverseRms", "lengthDivRmsTrans",
                 "rmsLongitudinal", "rmsTransverse", "hits"]
  FadcDsets = ["minvals", "fallTime", "riseTime"]
  AllFadcDsets = ["argMinval", "baseline", "eventNumber", "fallStop", "fallTime",
                  "minvals", "noisy", "riseStart", "riseTime"]

type
  ConfigFlagKind = enum
    cfNone, cfNoFadc, cfNoInGrid, cfNoOccupancy, cfNoPolya, cfNoFeSpectrum, cfProvideServer

  BackendKind* = enum
    bNone, bMpl, bPlotly, bGgPlot

  # variant object for the layout combining both
  # TODO: make generic or always use float?
  PlotV* = object
    annotations*: seq[string]
    case kind*: BackendKind
    of bMpl:
      # what needs to go here?
      plt*: PyObject
      fig*: PyObject
      ax*: PyObject
    of bPlotly:
      plLayout*: Layout
      plPlot*: Plot[float]
      plPlotJson*: PlotJson
    of bGgPlot:
      pltGg*: GgPlot[DataFrame]
      width*: float
      height*: float
      theme*: Theme
    else: discard

  ShapeKind = enum
    Rectangle, Square

  OutputFiletypeKind = enum
    ofUnknown,
    ofJson = "json"
    ofOrg = "org"


# global variable which stores the backend the user selected
var BKind*: BackendKind = bNone
# global OrderedSet to store all files we save to later
# combine to single PDF
var imageSet = initOrderedSet[string]()
var plotlyJson = newJObject()

const NPix = 256

var channel: Channel[(JsonNode, PacketKind)]
var stopChannel: Channel[bool]

# channel to send DataPackets received from client to work thread
var dpChannel: Channel[DataPacket]
var serveNewClientCh: Channel[bool]

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

proc genPlotDirname(h5name, outdir: string): string =
  ## generates unique name for the given input file based on its name and
  ## the current time
  let (_, name, _) = splitFile(h5name)
  let timeStr = format(now(), "yyyy-MM-dd'_'HH-mm-ss")
  result = outdir / name & "_" & timeStr

## is the directory for the current input file to store the plots
var fileDir: string
var fileType: string

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
  # TODO: move elsewhere!!!
  var
    plotlySaveSvg = false
    mplShowPlots = false
  if existsFile("config.toml"):
    let tomlConfig = parseToml.parseFile("config.toml")
    plotlySaveSvg = tomlConfig["General"]["plotlySaveSvg"].getBool
    mplShowPlots = tomlConfig["General"]["mplShowPlots"].getBool
  var fname = outfile
  if not fullPath:
    fname = "figs" / outfile & ".svg"
  info &"Saving file: {fname}"
  case BKind
  of bPlotly:
    if not plotlySaveSvg:
      plotlyJson[outfile] = jsonPlotly(p)
    else:
      if not p.plPlot.isNil:
        p.plPlot.saveImage(fname)
        imageSet.incl(fname)
  of bMpl:
    if not p.plt.isNil:
      discard p.plt.savefig(fname)
      imageSet.incl(fname)
    if mplShowPlots:
      discard callMethod(p.plt, "show")
    else:
      discard callMethod(p.plt, "close", "all")
  of bGgPlot:
    if p.kind == bGgPlot: # if plot was not initialized
      p.pltGg.ggsave(fname, p.width, p.height)
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
  of bGgPlot:
    result = PlotV(kind: bGgPlot,
                   theme: Theme(title: some(title),
                                xlabel: some(xlabel),
                                ylabel: some(ylabel)),
                   width: width.float,
                   height: FigHeight)
  else: discard

proc read*(h5f: var H5FileObj,
           runNumber: int,
           dsetName: string,
           chipNumber = 0,
           isFadc = false,
           dtype: typedesc = float,
           idx: seq[int] = @[]): seq[dtype] =
  ## reads the given dataset from the H5 file and returns it
  ## (potentially) converted to ``dtype``
  var dset: H5DataSet
  if isFadc:
    dset = h5f[(fadcRecoPath(runNumber) / dsetName).dset_str]
  else:
    dset = h5f[(recoPath(runNumber, chipNumber).string / dsetName).dset_str]
  when dtype is SomeNumber:
    if idx.len == 0:
      result = dset.readAs(dtype)
    else:
      # manual conversion required
      result = dset.readAs(idx, dtype)
  else:
    type subtype = getInnerType(dtype)
    # NOTE: only support 2D seqs
    when subType isnot SomeNumber:
      raise newException(Exception, "Cannot convert N-D sequence to type: " &
        subtype.name)
    else:
      if idx.len == 0:
        # convert to subtype and reshape to dsets shape
        result = dset.readAs(subtype).reshape2D(dset.shape)
      else:
        # manual conversion required
        result = dset[idx, subtype].reshape2D(dset.shape)

proc readVlen(h5f: var H5FileObj,
              runNumber: int,
              dsetName: string,
              chipNumber = 0,
              dtype: typedesc = float,
              idx: seq[int] = @[]): seq[seq[dtype]] =
  ## reads variable length data `dsetName` and returns it
  ## In contrast to `read` this proc does *not* convert the data.
  let vlenDtype = special_type(dtype)
  let dset = h5f[(recoPath(runNumber, chipNumber).string / dsetName).dset_str]
  if idx.len > 0:
    result = dset[vlenDType, dtype, idx]
  else:
    result = dset[vlenDtype, dtype]

template applyFilter(key: string, tomlConfig, theSet, fileInfo, field: untyped): untyped =
  if tomlConfig["General"].hasKey(key):
    let allowed = tomlConfig["General"][key].getElems
    for el in allowed:
      theSet.incl el.getInt.uint16
  if theSet.card > 0:
    fileInfo.field = fileInfo.field.filterIt(it.uint16 in theSet)


proc applyConfig(fileInfo: var FileInfo) =
  ## reads the  `config.toml` of the project and applies certain
  ## filters / settings to the fileInfo object. This allows us to
  ## disable certain plots / configurations
  # get allowed chips from toml
  let tomlConfig = parseToml.parseFile("config.toml")
  var chipSet: set[uint16]
  var runSet: set[uint16]
  # TODO: move elsewhere
  fileInfo.plotlySaveSvg = tomlConfig["General"]["plotlySaveSvg"].getBool
  applyFilter("allowedChips", tomlConfig, chipSet, fileInfo, chips)
  applyFilter("allowedRuns", tomlConfig, runSet, fileInfo, runs)

proc appliedConfig(fileInfo: FileInfo): FileInfo =
  ## returns a copy of the given `FileInfo` with the config.toml applied
  result = fileInfo
  result.applyConfig()

proc plotHist[T](xIn: seq[T], title, dset, outfile: string,
                 binS: float, binR: (float, float)): PlotV =
  ## plots the data in `x` as a histogram
  let xs = xIn.mapIt(it.float)

  var binRange = binR
  var binSize = binS
  if binRange[0] == binRange[1]:
    binRange = (xs.percentile(3), xs.percentile(97))
  if binSize == 0.0:
    # use 100 bins as a backup
    binSize = (binRange[1] - binRange[0]) / 100

  info &"Bin range {binRange} for dset: {dset}"
  result = initPlotV(title, dset, "#")
  case BKind
  of bPlotly:
    var traces: seq[Trace[float]]
    traces.add Trace[float](`type`: PlotType.Histogram,
                            bins: binRange,
                            binSize: binSize,
                            #nbins: nBins,
                            xs: xs,
                            name: dset)
    result.plPlot = Plot[float](layout: result.plLayout, traces: traces)
  of bMpl:
    let nbins = ((binRange[1] - binRange[0]) / binSize).round.int
    discard result.ax.hist(xs,
                         bins = nbins,
                         range = binRange)
  of bGgPlot:
    let df = seqsToDf(xs).filter(fn {float: `xs` >= binRange[0] and
                                            `xs` <= binRange[1]})
    result.pltGg = ggplot(df, aes("xs")) +
        geom_histogram(binWidth = binSize) +
        scale_x_continuous() +
        result.theme # just add the theme directly
  else: discard

proc plotBar[T](binsIn, countsIn: seq[seq[T]], title: string,
                xlabel: string, dsets: seq[string],
                outfile: string, drawPlots = false): PlotV =
  ## plots the data in `x` as a histogram
  let bins = binsIn.mapIt(it.mapIt(it.float))
  let counts = countsIn.mapIt(it.mapIt(it.float))
  result = initPlotV(title, xlabel, "#")
  var traces: seq[Trace[float]]
  if BKind in {bPlotly, bMpl}:
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
                                name: dset,
                                autoWidth: true)
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
      else: doAssert false
  else:
    var df = newDataFrame()
    # stack all data frames
    for i in 0 .. bins.high:
      let ldf = seqsToDf({ "bins" : binsIn[i],
                           "counts" : countsIn[i],
                           "dset" : constantColumn(dsets[i], binsIn[i].len) })
      df.add ldf
    if bins.len > 1:
      result.pltGg = ggplot(df, aes("bins", "counts")) +
          geom_histogram(aes(fill = "dset"), stat = "identity") +
          result.theme # just add the theme directly
    else:
      # don't classify by `dset`
      result.pltGg = ggplot(df, aes("bins", "counts")) +
          geom_histogram(stat = "identity") +
          result.theme # just add the theme directly

iterator chips(group: var H5Group): (int, H5Group) =
  ## returns all chip groups within the given Run group
  doAssert group.name.parentDir == "/reconstruction", "Bad group : " & group.name
  doAssert "run_" in group.name, "Bad group : " & group.name
  for grp in groups(group):
    if "chip_" in grp.name:
      let chipNum = grp.attrs["chipNumber", int]
      yield (chipNum, grp)

proc createCustomPlots(fileInfo: FileInfo): seq[PlotDescriptor] =
  ## define any plot that doesn't fit any of the other descriptions here
  block:
    let range = (low: -Inf, high: Inf, name: "All")
    result.add PlotDescriptor(runType: fileInfo.runType,
                              name: "riseTime",
                              runs: fileInfo.runs,
                              plotKind: pkFadcDset,
                              range: range,
                              binSize: 1.0,
                              binRange: (2.0, 502.0))
  block:
    let range = (low: -Inf, high: Inf, name: "All")
    result.add PlotDescriptor(runType: fileInfo.runType,
                              name: "fallTime",
                              runs: fileInfo.runs,
                              plotKind: pkFadcDset,
                              range: range,
                              binSize: 1.0,
                              binRange: (7.0, 707.0))

proc getBinSizeAndBinRange*(dset: string): (float, (float, float)) =
  ## accesses the helper procs to unwrap the options from the `dataset_helper`
  ## procs to return bin size and a bin range
  let binRangeO = getBinRangeForDset(dset)
  let binSizeO = getBinSizeForDset(dset)
  var binRange: tuple[low, high: float]
  var binSize: float
  var nbins: int
  if binRangeO.isSome:
    binRange = get(binRangeO)
  else:
    binRange = (0.0, 0.0)
  if binSizeO.isSome:
    # all good
    binSize = get(binSizeO)
  else:
    # else use number of bins for this dataset
    let nbinsO = getNumBinsForDset(dset)
    if nBinsO.isSome:
      nBins = get(nbinsO)
    else:
      nBins = 100
    # then calculate bin size
    if binRange[1] != binRange[0]:
      binSize = (binRange[1] - binRange[0]) / nBins.float
  result = (binSize, binRange)

proc histograms(h5f: var H5FileObj, runType: RunTypeKind,
                fileInfo: FileInfo,
                flags: set[ConfigFlagKind]): seq[PlotDescriptor] =
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
        for dset in InGridDsets:
          let (binSize, binRange) = getBinSizeAndBinRange(dset)
          result.add PlotDescriptor(runType: runType,
                                    name: dset,
                                    runs: fileInfo.runs,
                                    plotKind: pkInGridDset,
                                    chip: ch,
                                    range: r,
                                    binSize: binSize,
                                    binRange: binRange)
    if cfNoFadc notin flags:
      for dset in FadcDsets:
        let (binSize, binRange) = getBinSizeAndBinRange(dset)
        result.add PlotDescriptor(runType: runType,
                                  name: dset,
                                  runs: fileInfo.runs,
                                  plotKind: pkFadcDset,
                                  range: r,
                                  binSize: binSize,
                                  binRange: binRange)

proc calcOccupancy[T](x, y: seq[T]): Tensor[float] =
  ## calculates the occupancy of the given x and y datasets
  ## Either for a `seq[seq[T: SomeInteger]]` in which case we're calculating
  ## the occupancy of a raw clusters or `seq[T: SomeFloat]` in which case
  ## we're dealing with center positions of clusters
  result = newTensor[float]([NPix, NPix])
  # iterate over events
  for i in 0 .. x.high:
    let
      xEv = x[i]
      yEv = y[i]
    when T is seq:
      for j in 0 .. xEv.high:
        result[xEv[j].int, yEv[j].int] += 1.0
    elif T is SomeFloat:
      let xIdx = min(xEv.round.int, 255)
      let yIdx = min(yEv.round.int, 255)
      result[xIdx, yIdx] += 1.0

proc plotHist2D(data: Tensor[float], title, outfile: string): PlotV =
  ## creates a 2D histogram plot (basically an image) of the given 2D
  ## Tensor
  doAssert data.rank == 2
  result = initPlotV(title, "pixel x", "pixel y", ShapeKind.Square)
  case BKind
  of bPlotly:
    let tr = Trace[float](`type`: PlotType.HeatMap,
                          colormap: ColorMap.Viridis,
                          zs: data.toRawSeq.reshape2D([NPix, NPix]))
    result.plPlot = Plot[float](layout: result.plLayout, traces: @[tr])
  of bMpl:
    discard result.plt.imshow(data.toRawSeq.reshape2D([NPix, NPix]),
                       cmap = "viridis")
    discard result.plt.colorbar()
  of bGgPlot:
    # inefficient so far!
    # build 3 cols, x, y, z
    var
      x = newTensorUninit[int](data.shape[0])
      y = newTensorUninit[int](data.shape[0])
      z = newTensorUninit[float](data.shape[0])
    var i = 0
    for idx, val in data:
      x[i] = idx[0]
      y[i] = idx[1]
      z[i] = val
    let df = seqsToDf(x, y, z)
    result.pltGg = ggplot(df, aes("x", "y", fill = "z")) +
        geom_tile() +
        result.theme # just add the theme directly
  else:
    discard

proc plotScatter(pltV: var PlotV, x, y: seq[float], name, outfile: string,
                 isNew = false) =
  case BKind
  of bPlotly:
    # in this case plot is defined
    let trFit = Trace[float](`type`: PlotType.Scatter,
                             xs: x,
                             ys: y,
                             name: name)
    pltV.plPlot.traces.add @[trFit]
  of bMpl:
    # in this case `ax` is defined
    discard pltV.ax.plot(x,
                         y,
                         label = name,
                         color = "r")
  of bGgPlot:
    let df = seqsToDf(x, y)
    if not isNew:
      # in this case a plot already exists!
      pltV.pltGg = pltV.pltGg +
        geom_line(data = df, aes = aes("x", "y"),
                  color = some(parseHex("FF00FF"))) +
        pltV.theme # just add the theme directly
    else:
      # in this case a plot already exists!
      pltV.pltGg = ggplot(df, aes("x", "y")) +
        geom_line() +
        pltV.theme # just add the theme directly
  else:
    warn &"Unsupported backend kind: {BKind}"

proc plotScatter(x, y: seq[float], title, name, outfile: string): PlotV =
  ## wrapper around the above if no `PlotV` yet defined
  result = initPlotV(title, "x", "y", ShapeKind.Rectangle)
  case BKind
  of bPlotly:
    result.plPlot = Plot[float](layout: result.plLayout, traces: @[])
  else: discard
  result.plotScatter(x, y, name, outfile, isNew = true)

proc makeSubplot(pd: PlotDescriptor, plts: seq[PlotV]): PlotV =
  result = initPlotV(pd.title, pd.xlabel, pd.ylabel)
  let baseLayout = plts[0].plLayout
  baseLayout.width = 1600
  baseLayout.height = 800
  # move all potential annotations from plts[> 0] to baseLayout
  for i in 1 .. plts.high:
    baseLayout.annotations.add plts[i].plPlot.layout.annotations
  case plts.len
  of 2:
    result.plPlotJson = subplots:
      baseLayout: baseLayout
      plot:
        plts[0].plPlot
        pd.domain[0]
      plot:
        plts[1].plPlot
        pd.domain[1]
  else: echo "Subplots currently only supported for 2 plots!"

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

#proc plotPolyas(h5f: var H5FileObj, group: H5Group,
#                chipNum: int, runNumber: string): seq[Trace[float]] =
#  ## perform the plots of the polya distribution again (including the fits)
#  ## and return the polya data as to allow a combined plot afterwards
#  let
#    polya = h5f.read(runNumber.parseInt, "polya", chipNum,
#                     dtype = seq[float])
#    polyaFit = h5f.read(runNumber.parseInt, "polyaFit", chipNum,
#                        dtype = seq[float])
#  let nbins = polya.shape[0] - 1
#  var
#    bins = newSeq[float](nbins + 1)
#    binCenterFit = newSeq[float](nbins)
#    counts = newSeq[float](nbins)
#    countsFit = newSeq[float](nbins)
#  for i in 0 .. polya.high:
#    bins[i] = polya[i][0]
#    if i != nbins:
#      # do not take last element of polya datasets, since it's just 0
#      counts[i] = polya[i][1]
#      binCenterFit[i] = polyaFit[i][0]
#      countsFit[i] = polyaFit[i][1]
#  let names = @[&"Polya of chip {chipNum} for run {runNumber}"]
#  let title = &"Polya distribution of chip {chipNum} for run {runNumber}"
#  let xlabel = "Number of electrons"
#  let outfile = &"polya_run{runNumber}_chip{chipNum}"
#  var pltV = plotBar(@[bins], @[counts], title, xlabel, names, outfile)
#  # now add scatter plot for fit result
#  let nameFit = &"Polya fit of chip {chipNum} for run {runNumber}"
#  plotScatter(binCenterFit, countsFit, title, nameFit, outfile)

proc polya(h5f: var H5FileObj, runType: RunTypeKind,
           fileInfo: FileInfo,
           flags: set[ConfigFlagKind]): seq[PlotDescriptor] =
  ## creates the plots for the polya distribution of the datasets
  for ch in fileInfo.chips:
    result.add PlotDescriptor(runType: runType,
                              name: "polya",
                              xlabel: "Number of electrons",
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
                     title, xlabel, outfile: string,
                     ylabel = ""): PlotV =
  result = initPlotV(title, xlabel, ylabel, ShapeKind.Rectangle)
  case BKind
  of bPlotly:
    let tr = Trace[T](mode: PlotMode.Markers,
                      `type`: PlotType.Scatter,
                      xs: x,
                      ys: y,
                      marker: Marker[float](size: @[12.0]))
    result.plPlot = Plot[T](layout: result.plLayout, traces: @[tr])
  of bMpl:
    let dtm = pyImport("datetime")
    let xP = x.mapIt(dtm.datetime.utcfromtimestamp(int(it)))
    discard result.ax.plot_date(xP, y, label = title)
  of bGgPlot:
    let df = seqsToDf(x, y)
    result.pltGg = ggplot(df, aes("x", "y")) +
        geom_point() +
        result.theme # just add the theme directly
  else:
    discard

proc feSpectrum(h5f: var H5FileObj, runType: RunTypeKind,
                fileInfo: FileInfo,
                flags: set[ConfigFlagKind]): seq[PlotDescriptor] =
  ## creates the plot of the Fe spectrum
  for r in fileInfo.runs:
    let basePd = PlotDescriptor(runType: runType,
                                name: "FeSpectrumPlot",
                                xlabel: "# pixels",
                                runs: @[r],
                                chip: fileInfo.centerChip,
                                plotKind: pkFeSpec)
    let energyPd = replace(basePd):
      name = "EnergyCalib"
      xlabel = "# pixels"
      ylabel = "Energy / keV"
      plotKind = pkEnergyCalib
    let feChargePd = replace(basePd):
      name = "FeSpectrumChargePlot"
      xlabel = "Charge / 10^3 electrons"
      plotKind = pkFeSpecCharge
    let energyChargePd = replace(basePd):
      name = "EnergyCalibCharge"
      xlabel = "Charge / 10^3 electrons"
      ylabel = "Energy / keV"
      plotKind = pkEnergyCalibCharge
    result.add @[basePd, energyPd, feChargePd, energyChargePd]



  let photoVsTime = PlotDescriptor(runType: runType,
                                   name: "PhotoPeakVsTime",
                                   xlabel: "Time / unix",
                                   runs: fileInfo.runs,
                                   chip: fileInfo.centerChip,
                                   plotKind: pkFeVsTime)
  let photoChVsTime = PlotDescriptor(runType: runType,
                                     name: "PhotoPeakChargeVsTime",
                                     xlabel: "Time / unix",
                                     runs: fileInfo.runs,
                                     chip: fileInfo.centerChip,
                                     plotKind: pkFeChVsTime)
  let phPixDivChVsTime = PlotDescriptor(runType: runType,
                                   name: "PhotoPixDivChVsTime",
                                   xlabel: "Time / unix",
                                   runs: fileInfo.runs,
                                   chip: fileInfo.centerChip,
                                   plotKind: pkFePixDivChVsTime)
  let photoVsTimeHalfH = PlotDescriptor(runType: runType,
                                        name: "PhotoPeakVsTimeHalfHour",
                                        xlabel: "Time / unix",
                                        runs: fileInfo.runs,
                                        chip: fileInfo.centerChip,
                                        plotKind: pkFeVsTime,
                                        splitBySec: 1800,
                                        lastSliceError: 0.2,
                                        dropLastSlice: false)
  let photoChVsTimeHalfH = PlotDescriptor(runType: runType,
                                          name: "PhotoPeakChargeVsTimeHalfHour",
                                          xlabel: "Time / unix",
                                          runs: fileInfo.runs,
                                          chip: fileInfo.centerChip,
                                          plotKind: pkFeChVsTime,
                                          splitBySec: 1800,
                                          lastSliceError: 0.2,
                                          dropLastSlice: false)
  let phPixDivChVsTimeHalfH = PlotDescriptor(runType: runType,
                                             name: "PhotoPixDivChVsTimeHalfHour",
                                             xlabel: "Time / unix",
                                             runs: fileInfo.runs,
                                             chip: fileInfo.centerChip,
                                             plotKind: pkFePixDivChVsTime,
                                             splitBySec: 1800,
                                             lastSliceError: 0.2,
                                             dropLastSlice: false)
  result.add @[photoVsTime, phPixDivChVsTime,
               photoChVsTime, photoChVsTimeHalfH,
               photoVsTimeHalfH, phPixDivChVsTimeHalfH]
  #for runNumber, grpName in runs(h5f):
  #  var group = h5f[grpName.grp_str]
  #  let centerChip = h5f[group.parent.grp_str].attrs["centerChip", int]
  #  # call `importPyplot` because it sets the appropriate backend for us
  #  discard importPyplot()
  #  let path = "figs/"
  #  var outfiles = @[&"fe_spectrum_run{runNumber}.svg",
  #                   &"fe_energy_calib_run{runNumber}.svg"]
  #  let outfilesCh = outfiles.mapIt(path / ("charge_" & it))
  #  outfiles = outfiles.mapIt(path / it)
  #  for o in outfiles:
  #    imageSet.incl o
  #  for o in outfilesCh:
  #    imageSet.incl o
  #  h5f.fitToFeSpectrum(runNumber.parseInt, centerChip,
  #                      fittingOnly = not ShowPlots,
  #                      outfiles = outfiles,
  #                      writeToFile = false)
  #  # extract fit parameters from center chip group
  #  let chpGrpName = group.name / "chip_" & $centerChip
  #  #let feSpec =
  #  pixSeq.add h5f[(chpGrpName / "FeSpectrum").dset_str].attrs[
  #    parPrefix & $kalphaPix, float
  #  ]
  #  chSeq.add h5f[(chpGrpName / "FeSpectrumCharge").dset_str].attrs[
  #    parPrefix & $kalphaCharge, float
  #  ]
  #  dates.add parseTime(group.attrs["dateTime", string],
  #                      dateStr,
  #                      utc()).toUnix.float
  #
  ## now plot
  ## calculate ratio and convert to string to workaround plotly limitation of
  ## only one type for Trace
  #let ratio = zip(pixSeq, chSeq) --> map(it[0] / it[1])
  #plotDates(dates, ratio,
  #          title = "Photopeak pix / charge vs time",
  #          xlabel = "Date",
  #          dsets = "a",
  #          outfile = "photopeak_vs_time")
  #          #ylabel = "# pix / charge in e^-")

proc buildEvents[T, U](x, y: seq[seq[T]], ch: seq[seq[U]],
                       events: OrderedSet[int]): Tensor[float] =
  ## takes all events via their *indices* (not real event numbers) and returns
  ## a (256, 256) tensor for each event
  let nEvents = events.card
  result = newTensor[float]([nEvents, NPix, NPix])
  # TODO: add option to use real event number instead of indices
  for i in 0 ..< nEvents:
    if i in events:
      let
        xi = x[i]
        yi = y[i]
        chi = ch[i]
      for (ix, iy, ich) in zipEm(xi, yi, chi):
        result[i, iy.int, ix.int] = ich.float

proc readEvent*(h5f: var H5FileObj, run, chip: int, idx: seq[int]): Tensor[float] =
  ## helper proc to read data for given indices of events `idx` and
  ## builds a tensor of `idx.len` events.
  let
    x = h5f.readVlen(run, "x",
                     chipNumber = chip,
                     dtype = uint8, idx = idx)
    y = h5f.readVlen(run, "y",
                     chipNumber = chip,
                     dtype = uint8, idx = idx)
    ch = h5f.readVlen(run, "ToT",
                      chipNumber = chip,
                      dtype = uint16, idx = idx)
  result = buildEvents(x, y, ch, toOrderedSet(@[0]))

proc readEventSparse*(h5f: var H5FileObj, run, chip: int, idx: int): DataFrame =
  ## helper proc to read data for given indices of events `idx` and
  ## builds a tensor of `idx.len` events.
  let
    x = h5f.readVlen(run, "x",
                     chipNumber = chip,
                     dtype = uint8, idx = @[idx])
    y = h5f.readVlen(run, "y",
                     chipNumber = chip,
                     dtype = uint8, idx = @[idx])
    ch = h5f.readVlen(run, "ToT",
                      chipNumber = chip,
                      dtype = uint16, idx = @[idx])
  result = seqsToDf({"x" : x[0], "y" : y[0], "ch" : ch[0]})


proc createEventDisplayPlots(h5f: var H5FileObj,
                             run: int,
                             runType: RunTypeKind,
                             fileInfo: FileInfo,
                             events: OrderedSet[int]): seq[PlotDescriptor] =
  for ch in fileInfo.chips:
    for ev in events:
      result.add PlotDescriptor(runType: runType,
                                name: "EventDisplay",
                                xlabel: "x",
                                ylabel: "y",
                                runs: @[run],
                                chip: ch,
                                plotKind: pkInGridEvent,
                                event: ev)

proc createFadcPlots(h5f: var H5FileObj,
                     run: int,
                     runType: RunTypeKind,
                     fileInfo: FileInfo,
                     events: OrderedSet[int]): seq[PlotDescriptor] =
  for ev in events:
    result.add PlotDescriptor(runType: runType,
                              name: "EventDisplay FADC",
                              xlabel: "Clock cycles FADC",
                              ylabel: "U / V",
                              runs: @[run],
                              chip: fileInfo.centerChip,
                              plotKind: pkFadcEvent,
                              event: ev)

proc createOuterChipHistograms*(h5f: var H5FileObj,
                               runType: RunTypeKind,
                               fileInfo: FileInfo,
                               cutRange: CutRange): seq[PlotDescriptor] =
  var chips = fileInfo.chips
  let oldLen = chips.len
  chips.delete(chips.find(fileInfo.centerChip))
  doAssert chips.len + 1 == oldLen
  result.add PlotDescriptor(runType: runType,
                            name: &"Outer chips # pix hit {runType}",
                            xlabel: "# pixels on outer chips",
                            ylabel: "#",
                            runs: fileInfo.runs,
                            outerChips: chips,
                            rangeCenter: cutRange,
                            plotKind: pkOuterChips)

proc createInGridFadcEvDisplay(h5f: var H5FileObj,
                               run: int,
                               runType: RunTypeKind,
                               fileInfo: FileInfo,
                               events: OrderedSet[int]): seq[PlotDescriptor] =
  var finfo = fileInfo
  finfo.chips = @[finfo.centerChip]
  let ingridPds = createEventDisplayPlots(h5f, run, runType, fInfo, events)
  let ingridDomain = (left: 0.0, bottom: 0.0, width: 0.45, height: 1.0)
  let fadcPds = createFadcPlots(h5f, run, runType, fInfo, events)
  let fadcDomain = (left: 0.525, bottom: 0.05, width: 0.575, height: 0.6)
  for tup in zip(ingridPds, fadcPds):
    let
      ingrid = tup[0]
      fadc = tup[1]
    result.add PlotDescriptor(runType: runType,
                              name: "Event display InGrid + FADC",
                              runs: @[run],
                              chip: fileInfo.centerChip,
                              plotKind: pkSubPlots,
                              plots: @[ingrid, fadc],
                              domain: @[ingridDomain, fadcDomain])

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

proc readPlotFit(h5f: var H5FileObj, pd: PlotDescriptor):
     (seq[float], seq[float], seq[float], seq[float]) =
  ## reads the `plot` and `Fit` datasets for all runs in
  ## `pd`, stacks it and returns bins, counts for both
  ## This is currently used for:
  ## - ploya + polyaFit
  ## - FeSpectrum + FeSpectrumFit
  ## - FeSpectrumCharge + FeSpectrumChargeFit
  let plotSet = {pkPolya, pkCombPolya, pkFeSpec, pkFeSpecCharge}
  doAssert (pd.plotKind in plotSet),
     &"Only supported for {plotSet}. This plotKind is {pd.plotKind}"
  var
    lastBins: seq[float]
    lastBinsFit: seq[float]
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
    let nbins = polya.shape[0]
    let (bins, counts) = polya.split(SplitSeq.Seq2Col)
    let (binsFit, countsFit) = polyaFit.split(SplitSeq.Seq2Col)
    if lastBins.len > 0:
      doAssert lastBins == bins, "The ToT calibration changed between the " &
        "last two runs!"
      doAssert lastBinsFit == binsFit
    lastBins = bins
    lastBinsFit = binsFit
    if countsFull.len > 0:
      doAssert countsFull.len == counts.len
      countsFull = zip(countsFull, counts) --> map(it[0] + it[1])
      countsFitFull = zip(countsFitFull, countsFit) --> map(it[0] + it[1])
    else:
      countsFull = counts
      countsFitFull = countsFit
  result = (lastBins, countsFull, lastBinsFit, countsFitFull)

template createFeVsTimeDataFrame(h5f: var H5FileObj,
                                 group: H5Group,
                                 pd: PlotDescriptor): untyped {.dirty.} =
  ## template to create the HDF5 data frame to read FeVsTime data
  ## This can neither be a proc, because of the static type defined by
  ## feSchema, evSchema as well as the joined data frame
  ## Nor can it be a normal template due to a `cannot instantiate DynamicStackArray`
  ## from the joinTheta call for the 3rd argument.
  ## Put everything in a block to manually create some hygiene
  when false:
    const feSchema = [
      intCol("FeSpectrum"),
      intCol("FeSpectrumEvents")
    ]
    const evSchema = [
      intCol("eventNumber"),
      intCol("timestamp")
    ]
    type feType = schemaType(feSchema)
    type evType = schemaType(evSchema)
    let dfFe = fromHDF5[feType](DF, h5f, h5f.name,
                                group.name / "chip_" & $pd.chip).cache()
    let dfFeRenamed = dfFe.map(x => (
      FeSpectrum: x.FeSpectrum,
      eventNumber: x.FeSpectrumEvents
    ))
    let dfEv = fromHDF5[evType](DF, h5f, h5f.name,
                                group.name).cache()
    let joined = joinTheta(
      dfFeRenamed,
      dfEv,
      (a, b) => a.eventNumber == b.eventNumber,
      (a, b) => mergeTuple(a, b, ["eventNumber"])
    )
    joined
  else:
    # implement ggplotnim DF usage
    var path = group.name
    let evNum = h5f[path / "eventNumber", int]
    let tstamp = h5f[path / "timestamp", int]
    let dfEv = seqsToDf({ "eventNumber" : evNum, "timestamp" : tstamp })
    path &= "/chip_" & $pd.chip
    let feSpec = h5f[path / "FeSpectrum", int]
    let feSpecCh = h5f[path / "FeSpectrumCharge", float]
    let feSpecEvNum = h5f[path / "FeSpectrumEvents", int]
    let dfFeSpec = seqsToDf({ "eventNumber" : feSpecEvNum,
                              "FeSpectrum" : feSpec,
                              "FeSpectrumCharge" : feSpecCh})
    # return the joined df from template
    inner_join(dfEv, dfFeSpec, by = "eventNumber")

proc determineStartStopFeVsTime(df: DataFrame): (BiggestInt, BiggestInt) =
  ## returns the start and stop times of a dataframe with timestamps
  # now get timeslice from dataframe
  when false:
    let tStopDf = joined.sort(x => x.timestamp, SortOrder.Descending)
      .take(1)
      .collect()
    let tStop = tStopDf[0].timestamp

    let tStartDf = joined.sort(x => x.timestamp, SortOrder.Ascending)
      .take(1)
      .collect()
    let tStart = tStartDf[0].timestamp
  else:
    #let cc = df.collect()
    let tStart = df["timestamp"][0, int]
    let tStop = df["timestamp"][df.high, int]
  result = (tStart.BiggestInt, tStop.BiggestInt)

proc determineNumBatchesFeVsTime(length: int, pd: PlotDescriptor): int =
  result = length div pd.splitBySec
  var useLastBatch = false
  if not pd.dropLastSlice:
    inc result
    useLastBatch = true
  else:
    # get size of last batch
    let lastBatch = length mod pd.splitBySec
    if lastBatch > (length.float * (1.0 - pd.lastSliceError)).round.int:
      # then also take it
      inc result
      useLastBatch = true

proc handleInGridDset(h5f: var H5FileObj,
                      fileInfo: FileInfo,
                      pd: PlotDescriptor): (string, PlotV) =
  let ranges = @[pd.range]
  var allData: seq[float]
  for r in pd.runs:
    let data = h5f.read(r, pd.name, pd.chip, dtype = float)
    # perform cut on range
    let group = h5f[recoPath(r, pd.chip)]
    if pd.range[0] != -Inf and pd.range[1] != Inf:
      let idx = cutOnProperties(h5f, group,
                    ("energyFromCharge", pd.range[0], pd.range[1]))
      allData.add idx.mapIt(data[it])
    else:
      allData.add data
  result[0] = buildOutfile(pd, fileDir, fileType)
  let title = buildTitle(pd)
  result[1] = plotHist(allData, title, pd.name, result[0], pd.binSize, pd.binRange)

proc handleFadcDset(h5f: var H5FileObj,
                    fileInfo: FileInfo,
                    pd: PlotDescriptor): (string, PlotV) =
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
  result[0] = buildOutfile(pd, fileDir, fileType)
  let title = buildTitle(pd)
  result[1] = plotHist(allData, title, pd.name, result[0], pd.binSize, pd.binRange)

proc handleOccupancy(h5f: var H5FileObj,
                     fileInfo: FileInfo,
                     pd: PlotDescriptor): (string, PlotV) =
  # get x and y datasets, stack and get occupancies
  let vlenDtype = special_type(uint8)
  var occFull = newTensor[float]([NPix, NPix])
  for r in pd.runs:
    let
      x = h5f.readVlen(r, "x", pd.chip, dtype = uint8)
      y = h5f.readVlen(r, "y", pd.chip, dtype = uint8)
    let occ = clampedOccupancy(x, y, pd)
    # stack this run onto the full data tensor
    occFull = occFull .+ occ
  let title = buildTitle(pd)
  result[0] = buildOutfile(pd, fileDir, fileType)
  result[1] = plotHist2D(occFull, title, result[0])

proc handleOccCluster(h5f: var H5FileObj,
                      fileInfo: FileInfo,
                      pd: PlotDescriptor): (string, PlotV) =
  # plot center positions
  var occFull = newTensor[float]([NPix, NPix])
  for r in pd.runs:
    let
      group = h5f[recoPath(r, pd.chip)]
      # get centers and rescale to 256 max value
      centerX = h5f[(group.name / "centerX"), float].mapIt(it * NPix.float / 14.0)
      centerY = h5f[(group.name / "centerY"), float].mapIt(it * NPix.float / 14.0)
    let occ = clampedOccupancy(centerX, centerY, pd)
    # stack this run onto the full data tensor
    occFull = occFull .+ occ
  let title = buildTitle(pd)
  result[0] = buildOutfile(pd, fileDir, fileType)
  result[1] = plotHist2D(occFull, title, result[0])

proc handleBarScatter(h5f: var H5FileObj,
                      fileInfo: FileInfo,
                      pd: PlotDescriptor): (string, PlotV) =
  result[0] = buildOutfile(pd, fileDir, fileType)
  let xlabel = pd.xlabel
  let title = buildTitle(pd)
  let (bins, counts, binsFit, countsFit) = h5f.readPlotFit(pd)
  result[1] = plotBar(@[bins], @[counts], title, xlabel, @[title], result[0])
  # now add fit to the existing plot
  let nameFit = &"Fit of chip {pd.chip}"
  result[1].plotScatter(binsFit, countsFit, nameFit, result[0])

proc handleCombPolya(h5f: var H5FileObj,
                     fileInfo: FileInfo,
                     pd: PlotDescriptor): (string, PlotV) =
  var
    binsSeq: seq[seq[float]]
    countsSeq: seq[seq[float]]
    dsets: seq[string]
  let title = buildTitle(pd)
  result[0] = buildOutfile(pd, fileDir, fileType)
  for ch in pd.chipsCP:
    # get a local PlotDescriptor, which has this chip number
    let localPd = block:
      var tmp = pd
      tmp.chip = ch
      tmp
    let (bins, counts, binsFit, countsFit) = h5f.readPlotFit(localPd)
    binsSeq.add bins
    countsSeq.add counts
    dsets.add "Chip " & $ch
  let xlabel = "Number of electrons"
  result[1] = plotBar(binsSeq, countsSeq, title, xlabel, dsets, result[0])

proc handleFeVsTime(h5f: var H5FileObj,
                    fileInfo: FileInfo,
                    pd: PlotDescriptor): (string, PlotV) =
  const kalphaPix = 10
  const kalphaCharge = 4
  var kalphaIdx: int
  var dset: string
  case pd.plotKind
  of pkFeVsTime:
    kalphaIdx = kalphaPix
    dset = "FeSpectrum"
  of pkFeChVsTime:
    kalphaIdx = kalphaCharge
    dset = "FeSpectrumCharge"
  else: doAssert false, "Invalid for handle fe vs time!"
  const parPrefix = "p"
  const dateStr = "yyyy-MM-dd'.'HH:mm:ss" # example: 2017-12-04.13:39:45
  var
    pixSeq: seq[float]
    dates: seq[float] #string]#Time]

  let title = buildTitle(pd)
  result[0] = buildOutfile(pd, fileDir, fileType)
  for r in pd.runs:
    let group = h5f[(recoBase & $r).grp_str]
    let chpGrpName = group.name / "chip_" & $pd.chip
    if pd.splitBySec == 0:
      pixSeq.add h5f[(chpGrpName / dset).dset_str].attrs[
        parPrefix & $kalphaIdx, float
      ]
      dates.add parseTime(group.attrs["dateTime", string],
                          dateStr,
                          utc()).toUnix.float
    else:
      # split `FeSpectrum` hits by `splitBySec` and perform fit
      #let pyFitFe = pyImport("ingrid.fit_fe_spectrum")
      let joined = createFeVsTimeDataFrame(h5f, group, pd)

      let (tStart, tStop) = determineStartStopFeVsTime(joined)
      let length = (tStop - tStart).int
      let nBatches = determineNumBatchesFeVsTime(length, pd)

      for i in 0 ..< nBatches:
        # extract the correct data
        let slStart = i * pd.splitBySec + tStart
        let slStop = if (i + 1) == nBatches:
                       i * pd.splitBySec + length mod pd.splitBySec + tStart
                     else:
                       (i + 1) * pd.splitBySec + tStart
        if slStop - slStart < 300:
          echo "Skipping last batch of less than 5 min length: ",
            (slStop - slStart), " seconds."
          continue
        echo "Run: ", r, " in batch ", i
        echo "Starting from ", slStart
        echo "Stopping at ", slStop

        when false:
          let hits = joined.map(x => x.projectTo(FeSpectrum, timestamp))
            .filter(x => (x.timestamp >= slStart and x.timestamp < slStop))
            .map(x => x.FeSpectrum.int)
            .collect()
        else:
          let hitsTensor = joined.filter(fn {int: `timestamp` >= slStart and
                                                  `timestamp` < slStop})[dset]
            .toTensor(int)
          var hits: seq[int]
          if hitsTensor.size > 0:
            hits = hitsTensor.clone.toRawSeq
          else:
            # skip this batch
            echo "WARNING: Skipping empty batch!"
            continue

        var res: FeSpecFitData
        case pd.plotKind
        of pkFeVsTime: res = fitFeSpectrum(hits)
        of pkFeChVsTime: res = fitFeSpectrumCharge(hits)
        else: doAssert false, "not possible"
        #let res0 = res[0].mapIt(it.float)
        #let res1 = res[1].mapIt(it.float)

        #let kalphaLoc = res[0].popt.toNimSeq(float)[kalphaPix]

        #let ppp = barPlot(res0, res1)
        #ppp.traces[0].autowidth = true
        #
        ## now create fit from fit results, add to plot
        #let popt = res[2]
        #let countFit = res0.mapIt(feSpectrumFunc(popt, it))
        #let p2 = scatterPlot(res0, countFit).mode(PlotMode.Lines)
        #ppp.traces.add p2.traces[0]
        #ppp.show()

        #let pyRes = pyFitFe.fitAndPlotFeSpectrum([hits], "", ".", r,
        #                                         true)
        let kalphaLoc = res.kalpha
        let tstamp = (slStop + slStart).float / 2.0
        pixSeq.add kalphaLoc
        dates.add tstamp

  # calculate mean and STD of peak location
  let
    stdPeak = standardDeviation(pixSeq)
    xloc = max(dates)
    yloc = max(pixSeq)
    # calculate mean of peak location
  let meanPeak = mean(pixSeq)
  case BKind
  of bPlotly:
    # now plot
    result[1] = plotDates(dates, pixSeq,
                          title = title,
                          xlabel = pd.xlabel,
                          outfile = result[0])
    let stdAnn = plotly.Annotation(x: xloc,
                                   y: yloc,
                                   text: &"Peak variation σ: {stdPeak:.2f}")
    let meanAnn = replace(stdAnn):
      y = yloc - (yloc - min(pixSeq)) * 0.05
      text = &"Mean location µ: {meanPeak:.2f}"
    result[1].plPlot.layout.annotations.add [stdAnn, meanAnn]
  of bGgPlot:
    let stdAnn = ggplotnim.Annotation(x: some(xloc),
                                      y: some(yloc),
                                      text: &"Peak variation σ: {stdPeak:.2f}")
    let meanAnn = replace(stdAnn):
      x = some(yloc - (yloc - min(pixSeq)) * 0.05)
      text = &"Mean location µ: {meanPeak:.2f}"
    # now plot
    result[1] = plotDates(dates, pixSeq,
                          title = title,
                          xlabel = pd.xlabel,
                          outfile = result[0])
    result[1].pltGg.annotations.add @[stdAnn, meanAnn]
  else: echo "WARNING: annotations unsupported on backend: " & $BKind


proc handleFePixDivChVsTime(h5f: var H5FileObj,
                            fileInfo: FileInfo,
                            pd: PlotDescriptor): (string, PlotV) =
  const kalphaPix = 10
  const kalphaCharge = 4
  const parPrefix = "p"
  const dateStr = "yyyy-MM-dd'.'HH:mm:ss" # example: 2017-12-04.13:39:45
  var
    pixSeq: seq[float]
    chSeq: seq[float]
    dates: seq[float] #string]#Time]
  let title = buildTitle(pd)
  result[0] = buildOutfile(pd, fileDir, fileType)
  for r in pd.runs:
    let group = h5f[(recoBase & $r).grp_str]
    let chpGrpName = group.name / "chip_" & $pd.chip
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
  result[1] = plotDates(dates, ratio,
                        title = title,
                        xlabel = pd.xlabel,
                        outfile = result[0])
                        #ylabel = "# pix / charge in e^-")

proc handleOuterChips(h5f: var H5FileObj,
                      fileInfo: FileInfo,
                      pd: PlotDescriptor): (string, PlotV) =
  var data: seq[int]
  const
    rmsTransHigh = 1.2
    eccHigh = 1.3

  for r in pd.runs:
    # cut on `centerChip` for rms and eccentricity
    let
      grp = h5f[recoPath(r, fileInfo.centerChip)]
      idx = cutOnProperties(h5f, grp,
                            ("eccentricity", -Inf, eccHigh),
                            ("rmsTransverse", -Inf, rmsTransHigh),
                            ("energyFromCharge", pd.rangeCenter.low, pd.rangeCenter.high))
      evNumCenterAll = h5f.read(r, "eventNumber", fileInfo.centerChip, dtype = int)
      evNumCenter = toSet(idx.mapIt(evNumCenterAll[it]))
    for c in pd.outerChips:
      let
        hitsAll = h5f.read(r, "hits", c, dtype = int)
        evNumAll = h5f.read(r, "eventNumber", c, dtype = int)
      # now add all hits passing
      for i, ev in evNumAll:
        if ev in evNumCenter and hitsAll[i] > 3:
          data.add hitsAll[i]
  result[0] = pd.title
  const binSize = 3.0
  const binRange = (0.0, 400.0)
  let outfile = buildOutfile(pd, fileDir, fileType)
  result[1] = plotHist(data, pd.title, pd.name, pd.title, binSize, binRange)

iterator ingridEventIter(h5f: var H5FileObj,
                         fileInfo: FileInfo,
                         pds: seq[PlotDescriptor]): (string, PlotV) =
  var events {.global.}: seq[int]
  # all PDs are guaranteed to be from  the same run!
  let run = pds[0].runs[0]
  # TODO: check if all PDs necessarily have same chip. Will be the case for
  # InGrid = FADC, but not necessarily only chip?
  let chip = pds[0].chip

  var lastRun {.global.} = 0
  var evNums {.global.}: seq[int]
  var evTab {.global.}: Table[int, int]
  if run != lastRun:
    evNums = h5f.read(run, "eventNumber", dtype = int,
                      chipNumber = chip)
    evTab = initTable[int, int]()
    for i, ev in evNums:
      evTab[i] = ev
      events.add ev
    lastRun = run

  events.setLen(0)
  for pd in pds:
    events.add evTab[pd.event]

  let ev = readEvent(h5f, run, chip, events)
  for i, pd in pds:
    var texts: seq[string]
    let event = evTab[pd.event]
    for d in InGridDsets:
      let val = h5f.read(run, d,
                         chipNumber = chip,
                         dtype = float, idx = @[event])[0]
      let s = &"{d:>20}: {val:6.4f}"
      texts.add s

    let title = buildTitle(pd)
    var outfile = buildOutfile(pd, fileDir, fileType)
    var pltV = plotHist2D(ev[i,_,_].squeeze.clone, title, outfile)
    case BKind
    of bPlotly:
      for i, a in texts:
        pltV.plPlot.layout.annotations.add plotly.Annotation(x: 0.0,
                                                             y: 0.0,
                                                             xshift: 800.0,
                                                             yshift: 600.0 - (i.float * 20.0),
                                                             text: texts[i])
    of bGgPlot:
      for i, a in texts:
        pltV.pltGg.annotations.add ggplotnim.Annotation(left: some(0.1),
                                                        bottom: some(0.9 - i.float * 0.05),
                                                        text: texts[i])
    else:
      echo "InGrid Event property annotations not supported on " &
        "Matplotlib backend yet!"
    pltV.annotations = texts
    yield (outfile, pltV)

iterator fadcEventIter(h5f: var H5FileObj,
                       fileInfo: FileInfo,
                       pds: seq[PlotDescriptor]): (string, PlotV) =
  var events {.global.}: seq[int]
  # all PDs are guaranteed to be from  the same run!
  let run = pds[0].runs[0]

  var lastRun {.global.} = 0
  var evNums {.global.}: seq[int]
  var evTab {.global.}: Table[int, int]
  if run != lastRun:
    evNums = h5f.read(run, "eventNumber", dtype = int,
                      isFadc = true)
    evTab = initTable[int, int]()
    for i, ev in evNums:
      evTab[ev] = i
    lastRun = run
  events.setLen(0)
  for pd in pds:
    events.add evTab[pd.event]

  let
    dset = h5f[(fadcRecoPath(run) / "fadc_data").dset_str]
    fadc = dset.read_hyperslab(float64, @[events[0], 0], @[events.len, 2560])
    fShape = [fadc.len div 2560, 2560]
    fTensor = fadc.toTensor.reshape(fShape)

  for i, pd in pds:
    var texts: seq[string]
    let event = evTab[pd.event]
    for d in AllFadcDsets:
      let val = h5f.read(run, d, isFadc = true,
                         dtype = float, idx = @[event])[0]
      let s = &"{d:15}: {val:6.4f}"
      texts.add s

    let xFadc = toSeq(0 ..< 2560).mapIt(it.float)
    let title = buildTitle(pd)
    var outfile = buildOutfile(pd, fileDir, fileType)
    let yFadc = fTensor[i,_].squeeze.clone.toRawSeq
    var pltV = plotScatter(xFadc, yFadc, title, title, outfile)
    case BKind
    of bPlotly:
      pltV.plPlot = pltV.plPlot.mode(PlotMode.Lines).lineWidth(1)
      for i, a in texts:
        pltV.plPlot.layout.annotations.add plotly.Annotation(x: 0.0,
                                                             y: 0.0,
                                                             xshift: 1100.0,
                                                             yshift: 600.0 - (i.float * 20.0),
                                                             text: texts[i])
      pltV.annotations = texts
    of bGgPlot:
      for i, a in texts:
        pltV.pltGg.annotations.add ggplotnim.Annotation(left: some(0.1),
                                                        bottom: some(0.9 - (i.float * 0.05)),
                                                        text: texts[i])
      pltV.annotations = texts
    else:
      echo "FADC property annotations not supported on Matplotlib backend yet!"
    yield (outfile, pltV)

proc handleIngridEvent(h5f: var H5FileObj,
                       fileInfo: FileInfo,
                       pd: PlotDescriptor): (string, PlotV) =
  doAssert pd.plotKind == pkInGridEvent
  for outfile, pltV in ingridEventIter(h5f, fileInfo, @[pd]):
    # only a single pd
    result = (outfile, pltV)

proc handleFadcEvent(h5f: var H5FileObj,
                     fileInfo: FileInfo,
                     pd: PlotDescriptor): (string, PlotV) =
  doAssert pd.plotKind == pkFadcEvent
  for outfile, pltV in fadcEventIter(h5f, fileInfo, @[pd]):
    # only a single pd
    result = (outfile, pltV)

proc createPlot*(h5f: var H5FileObj, fileInfo: FileInfo,
                 pd: PlotDescriptor): (string, PlotV)

proc handleSubPlots(h5f: var H5FileObj,
                    fileInfo: FileInfo,
                    pd: PlotDescriptor): (string, PlotV) =
#    result.add PlotDescriptor(runType: runType,
#                              name: "Event display InGrid + FADC",
#                              runs: @[run],
#                              chip: fileInfo.centerChip,
#                              plotKind: pkSubPlots,
#                              plots: @[ingrid, fadc],
#                              domain: @[ingridDomain, fadcDomain])
  var fnames: seq[string]
  var plts: seq[PlotV]
  for p in pd.plots:
    let (outfile, pltV) = createPlot(h5f, fileInfo, p)
    fnames.add outfile
    plts.add pltV

  # now combine plots
  let plt = makeSubplot(pd, plts)
  plt.plPlotJson.show()

proc createPlot*(h5f: var H5FileObj,
                 fileInfo: FileInfo,
                 pd: PlotDescriptor): (string, PlotV) =
  ## creates a plot of kind `plotKind` for the data from all runs in `runs`
  ## for chip `chip`
  # TODO: think: have createPlot return the `PlotV` object. Then we could
  # call this proc recursively and add other plots to the existing figure
  let test = % pd
  echo "Test is ", test.pretty
  # test reverse
  #let test2 = parsePd(test)
  #doAssert test2 == pd

  case pd.plotKind
  of pkInGridDset:
    result = handleInGridDset(h5f, fileInfo, pd)
  of pkFadcDset:
    result = handleFadcDset(h5f, fileInfo, pd)
  of pkOccupancy:
    result = handleOccupancy(h5f, fileInfo, pd)
  of pkOccCluster:
    result = handleOccCluster(h5f, fileInfo, pd)
  of pkPolya, pkFeSpec, pkFeSpecCharge:
    result = handleBarScatter(h5f, fileInfo, pd)
  of pkCombPolya:
    result = handleCombPolya(h5f, fileInfo, pd)
  of pkFeVsTime, pkFeChVsTime:
    result = handleFeVsTime(h5f, fileInfo, pd)
  of pkFePixDivChVsTime:
    result = handleFePixDivChVsTime(h5f, fileInfo, pd)
  of pkInGridEvent:
    # TODO: make this into a call to an `InGridEventIterator`
    result = handleIngridEvent(h5f, fileInfo, pd)
    #for outfile, pltV in ingridEventIter(h5f, fileInfo, pd):
    #  yield (outfile, pltV)
  of pkFadcEvent:
    result = handleFadcEvent(h5f, fileInfo, pd)
    #for outfile, pltV in fadcEventIter(h5f, fileInfo, pd):
    #  yield (outfile, pltV)
  of pkSubPlots:
    result = handleSubPlots(h5f, fileInfo, pd)
    #for outfile, pltV in subPlotsIter(h5f, fileInfo, pd):
    #  yield (outfile, pltV)
  of pkOuterChips:
    result = handleOuterChips(h5f, fileInfo, pd)
  else:
    discard

  # finally call savePlot if we actually created a plot
  #if not result[1].plPlot.isNil:
  info &"Calling savePlot for {pd.plotKind} with filename {result[0]}"
  savePlot(result[1], result[0], fullPath = true)


proc createOrg(outfile, fileType: string) =
  ## creates a simple org file consisting of headings and images
  ## SVGs are implemented using raw inline SVG
  const header = """
* =$1=

"""
  const tmpl = """
** =$1=

#+BEGIN_EXPORT html
$2
#+END_EXPORT

"""
  # now build pdf
  var orgStr = header % [outfile]
  for im in imageSet:
    if im.len > 0:
      info &"Adding image {im} and {im.repr}"
      let (dir, name, ext) = im.splitFile()
      if fileType == "svg":
        let data = readFile(im)
        orgStr = orgStr & tmpl % [name, data]
      elif fileType == "pdf":
        orgStr = orgStr & tmpl % [
          name,
          &"""<embed src="{im}" width="{FigWidth}px" height="{FigHeight * 1.1}px"/>"""
        ]
      elif fileType == "png":
        # encode via base64
        let data = readFile(im)
        orgStr = orgStr & tmpl % [
          name,
          &"""<img alt="My Image" src="data:image/png;base64,{encode(data)}" />"""
        ]

      else: doAssert false, "Invalid filetype: " & $filetype

  var f = open(outfile, fmWrite)
  f.write(orgStr)
  f.close()

  shell:
    emacs ($outfile) "--batch -f org-html-export-to-html --kill"

proc clearPlots() =
  ## clears the imageSet as well as the JsonNode for plotly
  imageSet = initOrderedSet[string]()
  plotlyJson = newJObject()

proc jsonDump(imSet: OrderedSet[string], plotly: JsonNode): JsonNode=
  ## combines the image set for SVGs and the Plotly JsonNode to a combined
  ## JsonNode, which can either be sent via socket or dumped to a file
  result = newJObject()
  result["svg"] = newJObject()
  var svgJ = newJObject()
  for im in imSet:
    svgJ[im] = % readFile(im)

  result["svg"] = svgJ
  result["plotly"] = plotly

proc jsonDump(outfile = "", clear = false): string =
  ## reads either the created images from the `imageSet` and dumps them
  ## contained in a JsonNode to a file.
  ## Additionally stores the plotly plots

  let jDump = jsonDump(imageSet, plotlyJson)

  if outfile.len > 0:
    info "Num SVG plots: ", imageSet.card
    info "Num Json plots: ", plotlyJson.len

    info "Writing JSON file: ", outfile
    var f = open(outfile, fmWrite)
    var outstr = ""
    outstr.toUgly(jdump)
    f.write(outstr)
    f.close()
  else:
    result = jdump.pretty

  if clear:
    clearPlots()

proc handleOutput(basename: string, flags: set[ConfigFlagKind]): string =
  ## handles output file creation. Reads the config.toml file
  ## and uses it to decide whether to store data as Org file or
  ## dump to Json
  # TODO: create some kind of Config object, which does this once at
  # startup of the program
  # TODO: when Org is used, we need to save all plots as SVG!
  result = basename
  let tomlConfig = parseToml.parseFile("config.toml")
  let outfileKind = parseEnum[OutputFiletypeKind](
    tomlConfig["General"]["outputFormat"].getStr,
    ofUnknown
  )
  for fl in flags:
    result &= "_" & $fl
  case outfileKind
  of ofOrg:
    result &= ".org"
    createOrg(result, fileType)
  of ofJson:
    result &= ".json"
    discard jsonDump(result)
  else:
    warn "Unsupported output file format!"

proc plotsFromPds(h5f: var H5FileObj,
                  fileInfo: FileInfo,
                  pds: seq[PlotDescriptor]) =
  ## calls `createPlot` for each PD and saves filename to
  ## `imageSet` if necessary
  for p in pds:
    let (fileF, pltV) = h5f.createPlot(fileInfo, p)
    case BKind
    of bPlotly:
      if fileInfo.plotlySaveSvg:
        imageSet.incl fileF
    of bGgPlot:
      imageSet.incl fileF
    else: discard

proc serveNewClient(): bool =
  let (newClient, _)  = serveNewClientCh.tryRecv()
  result = newClient

proc awareSend[T](ch: var Channel[T], data: JsonNode, pKind: PacketKind) =
  while not trySend(ch, (data, pKind)):
    echo "DataThread: sleeping to send ", pKind
    sleep(100)

proc serveRequests(h5f: var H5FileObj,
                   fileInfo: FileInfo,
                   pds: seq[PlotDescriptor]) =
  while true:
    # receive data packet
    # TODO: replace by a `reqChannel`?
    echo "Serving client requests!"

    let dp = dpChannel.recv()
    echo "Received kind ", dp.reqKind
    case dp.reqKind
    of rqPlot:
      let pd = dp.payload.parseJson.parsePd
      #for sp in createPlotIter(h5f, fileInfo, pd):
      discard createPlot(h5f, fileInfo, pd)
      let jData = jsonDump(imageSet, plotlyJson)
      # echo "Jdata corresponding: ", jData
      awareSend(channel, jData, PacketKind.Plots)
      info "Clear plots"
      clearPlots()
    of rqPlotDescriptors:
      let pdJson = % pds
      awareSend(channel, pdJson, PacketKind.Descriptors)
    of rqFileInfo:
      let fiJson = % fileInfo
      awareSend(channel, fiJson, PacketKind.FileInfo)
    else:
      echo "Unknown req kind in context: ", dp.reqKind
      discard

proc serve(h5f: var H5FileObj,
           fileInfo: FileInfo,
           pds: seq[PlotDescriptor],
           flags: set[ConfigFlagKind] = {}) =
  ## serves the client
  # before we do anything, send all PDs to the client

  while true:
    while not serveNewClient():
      sleep(100)

    if serveNewClient():
      continue

    for p in pds:
      #for sp in createPlotIter(h5f, fileInfo, p):
      discard createPlot(h5f,fileInfo, p)
      echo "Number of pds ", pds.len
      echo "Current plot ", p
      let jData = jsonDump(imageSet, plotlyJson)
      # echo "Jdata corresponding: ", jData
      while not trySend(channel, (jData, PacketKind.Plots)):
        echo "DataThread: sleep..."
        sleep(100)
      info "Clear plots"
      clearPlots()

    echo "Done with all plots!, sending stop!"
    stopChannel.send(true)

proc eventDisplay(h5file: string,
                  run: int,
                  runType: RunTypeKind,
                  bKind: BackendKind,
                  flags: set[ConfigFlagKind]): string =
  ## use as event display tool
  var h5f = H5file(h5file, "r")
  let fileInfo = getFileInfo(h5f)
  let fInfoConfig = fileInfo.appliedConfig()
  let events = toOrderedSet(toSeq(0 .. 50))#initOrderedSet[int]())
  let pds = createEventDisplayPlots(h5f, run, runType, fInfoConfig, events)

  if cfProvideServer in flags:
    serve(h5f, fileInfo, pds, flags)
  else:
    plotsFromPds(h5f, fInfoConfig, pds)
    let outfile = "eventDisplay"  & "_" & getRunsStr(fInfoConfig.runs)
    result = handleOutput(outfile, flags)

  discard h5f.close()

proc createCalibrationPlots(h5file: string,
                            bKind: BackendKind,
                            runType: RunTypeKind,
                            flags: set[ConfigFlagKind]): string =
  ## creates QA plots for calibration runs
  var h5f = H5file(h5file, "r")
  var fileInfo = getFileInfo(h5f)
  let fInfoConfig = fileInfo.appliedConfig()
  # var imageSet = initOrderedSet[string]()
  var pds: seq[PlotDescriptor]

  const length = "length"
  if cfNoOccupancy notin flags:
    pds.add occupancies(h5f, runType, fInfoConfig, flags) # plus center only
  if cfNoPolya notin flags:
    pds.add polya(h5f, runType, fInfoConfig, flags)
  if cfNoFeSpectrum notin flags:
    pds.add feSpectrum(h5f, runType, fInfoConfig, flags)
  # energyCalib(h5f) # ???? plot of gas gain vs charge?!
  pds.add histograms(h5f, runType, fInfoConfig, flags) # including fadc

  if cfProvideServer in flags:
    serveRequests(h5f, fileInfo, pds)
    #serve(h5f, fileInfo, pds, flags)
  else:
    plotsFromPds(h5f, fInfoConfig, pds)
    echo "Image set is ", imageSet.card

    # likelihoodHistograms(h5f) # need to cut on photo peak and esc peak
    # neighborPixels(h5f)
    discard h5f.close()
    let outfile = "calibration" & "_" & getRunsStr(fInfoConfig.runs)
    result = handleOutput(outfile, flags)

proc createBackgroundPlots(h5file: string,
                           bKind: BackendKind,
                           runType: RunTypeKind,
                           flags: set[ConfigFlagKind]): string =
  ## creates QA plots for calibration runs
  var h5f = H5file(h5file, "r")
  var fileInfo = getFileInfo(h5f)
  let fInfoConfig = fileInfo.appliedConfig()
  var pds: seq[PlotDescriptor]
  const length = "length"
  if cfNoOccupancy notin flags:
    #occupancies(h5f, flags) # plus center only
    pds.add occupancies(h5f, runType, fInfoConfig, flags) # plus center only
  if cfNoPolya notin flags:
    pds.add polya(h5f, runType, fInfoConfig, flags)
  # energyCalib(h5f) # ???? plot of gas gain vs charge?!
  pds.add histograms(h5f, runType, fInfoConfig, flags) # including fadc

  # pds.add createCustomPlots(fInfoConfig)

  if cfProvideServer in flags:
    serveRequests(h5f, fileInfo, pds)
    #serve(h5f, fileInfo, pds, flags)
  else:
    plotsFromPds(h5f, fInfoConfig, pds)
    echo "Image set is ", imageSet.card

    # likelihoodHistograms(h5f) # need to cut on photo peak and esc peak
    # neighborPixels(h5f)
    discard h5f.close()
    let outfile = "background" & "_" & getRunsStr(fInfoConfig.runs)
    result = handleOutput(outfile, flags)

proc createXrayFingerPlots(bKind: BackendKind, flags: set[ConfigFlagKind]): string =
  discard

proc handlePlotTypes(h5file: string,
                     bKind: BackendKind,
                     runType: RunTypeKind,
                     flags: set[ConfigFlagKind],
                     evDisplayRun = none[int]()) =
  ## handles dispatch of the correct data type / kind / mode to be plotted
  ## calibration, X-ray, background, event display...
  ## If not run in server mode, launch staticDisplay with output JSON file
  ## TODO: what happens for SVG? Need to check too.
  var outfile = ""
  if evDisplayRun.isSome():
    outfile = eventDisplay(h5file, evDisplayRun.get, runType, BKind, flags)
  else:
    case runType
    of rtCalibration:
      outfile = createCalibrationPlots(h5file, bKind, runType, flags)
    of rtBackground:
      outfile = createBackgroundPlots(h5file, bKind, runType, flags)
    of rtXrayFinger:
      outfile = createXrayFingerPlots(bKind, flags)
    else:
      discard
  if cfProvideServer notin flags:
    # if not using server, start client, in case the data is being stored as
    # JSON
    let filecall = &"--file:{outfile}"
    shell:
      ./karaRun "-r" ($filecall) staticClient.nim

proc parseBackendType(backend: string): BackendKind =
  ## given a string describing a run type, return the correct
  ## `RunTypeKind`
  if backend.normalize in ["python", "py", "matplotlib", "mpl"]:
    result = bMpl
  elif backend.normalize in ["plotly"]:
    result = bPlotly
  elif backend.normalize in ["nim", "ggplot"]:
    result = bGgPlot
  else:
    result = bNone

proc sendDataPacket(ws: AsyncWebSocket, data: JsonNode, kind: PacketKind) =
  # convert payload to a string
  var sendData = ""
  toUgly(sendData, data)
  # split into several parts (potentially)
  echo "Now sending ", sendData[0 .. (min(150, sendData.high))]
  let nParts = ceil(sendData.len.float / FakeFrameSize.float).int
  if nParts > 1:
    for i in 0 ..< nParts:
      let i_start = i * FakeFrameSize
      let i_stop = min((i + 1) * FakeFrameSize, sendData.len)
      let dPart = sendData[i_start ..< i_stop]
      var packet: DataPacket
      if i == 0:
        packet = initDataPacket(kind,
                                Messages.DataStart,
                                payload = dPart,
                                recvData = true)
        asyncCheck ws.sendText(packet.asData)
      elif i < nParts - 1:
        packet = initDataPacket(kind,
                                Messages.Data,
                                payload = dPart,
                                recvData = true)
        asyncCheck ws.sendText(packet.asData)
      else:
        doAssert i == nParts - 1
        packet = initDataPacket(kind,
                                Messages.DataStop,
                                payload = dPart,
                                recvData = false,
                                done = true)
        waitFor ws.sendText(packet.asData)
  else:
    # handle the single packet case differently
    let packet = initDataPacket(kind,
                                Messages.DataSingle,
                                payload = sendData,
                                recvData = false,
                                done = true)
    waitFor ws.sendText(packet.asData)


proc processClient(req: Request) {.async.} =
  ## handle a single client
  let (ws, error) = await verifyWebsocketRequest(req)
  if ws.isNil:
    echo "WS negotiation failed: ", error
    await req.respond(Http400, "Websocket negotiation failed: " & error)
    req.client.close()
    return
  else:
    # receive connection successful package
    let (opcodeConnect, dataConnect) = await ws.readData()
    if dataConnect != $Messages.Connected:
      echo "Received wrong packet, quit early"
      return
    echo "New websocket customer arrived! ", dataConnect
  # let worker know that new client connected
  serveNewClientCh.send(true)

  # TODO: send PlotDescriptors first and then go into the loop
  # Include a channel to tell the main thread to restart plots, if
  # already underway

  var dp: DataPacket
  while true:
    let (opcode, data) = await ws.readData()
    echo "(opcode: ", opcode, ", data length: ", data

    case opcode
    of Opcode.Text, Opcode.Cont:
      # parse the `DataPacket` we just received
      dp = parseDataPacket(data)
      # now just send this DataPacket and await response from worker
      case dp.kind
      of PacketKind.Request:
        case dp.reqKind
        of rqPlot, rqPlotDescriptors, rqFileInfo:
          echo "Sending via dpChannel!"
          dpChannel.send(dp)
          echo "Sent"
        else: discard
      else: discard
      # first check whether main process is done
      #let (stopAvailable, toStop) = stopChannel.tryRecv()
      let (hasData, dataTup) = channel.tryRecv()
      if hasData:
        let (jData, packetKind) = dataTup
        ws.sendDataPacket(jData, packetKind)
      sleep(50)
      #if toStop:
      #  break
    of Opcode.Close:
      asyncCheck ws.close()
      let (closeCode, reason) = extractCloseData(data)
      echo "Socket went away, close code: ", closeCode, ", reason: ", reason
      req.client.close()
      return
    else:
      echo "Unkown error: ", opcode

  echo "This client dies now!"
  asyncCheck ws.close()
  req.client.close()
  stopChannel.send(true)

proc serveClient(server: AsyncHttpServer) {.async.} =
  var
    stopAvailable = false
    stop = false
  var clientFut = server.serve(Port(8081), processClient)
  while not stopAvailable:
    # NOTE: this causes the puzzling stops of the program! `serve` sends stop via the
    # `stopChannel`, which in reality is just being catched here! That just force stops
    # the server itself...
    #(stopAvailable, stop) = stopChannel.tryRecv()
    #if stop:
    #  server.close()
    #  break
    if not clientFut.finished:
      # client still connected, continue
      poll(500)
    else:
      # this client disconnected early, so accept another one
      clientFut = server.serve(Port(8081), processClient)

proc serve() =
  var server: AsyncHttpServer
  server = newAsyncHttpServer(reuseAddr = true, reusePort = true)
  channel.open(1)
  stopChannel.open(1)
  serveNewClientCh.open(1)
  dpChannel.open(1)
  shell:
    ./karaRun "-d:client -r staticClient.nim"

  waitFor server.serveClient()

proc plotData*() =
  ## the main workhorse of the server end
  let args = docopt(doc)
  info &"Received arguments:\n  {args}"
  let h5file = $args["<H5file>"]
  fileDir = genPlotDirName(h5file, "figs")
  fileType = parseToml.parseFile("config.toml")["General"]["filetype"].getStr
  discard existsOrCreateDir(fileDir)
  let runTypeStr = $args["--runType"]
  let backendStr = $args["--backend"]
  let evDisplayStr = $args["--eventDisplay"]
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
  if $args["--server"] == "true":
    flags.incl cfProvideServer

  info &"Flags are:\n  {flags}"

  var thr: Thread[void]
  if cfProvideServer in flags:
    # set up the server and launch the client
    # create and run the websocket server
    thr.createThread(serve)

  var runType: RunTypeKind
  var bKind: BackendKind
  var evDisplayRun: Option[int] = none[int]()
  if runTypeStr != "nil":
    runType = parseRunType(runTypeStr)
  if backendStr != "nil":
    BKind = parseBackendType(backendStr)
  if evDisplayStr != "nil":
    evDisplayRun = some(parseInt(evDisplayStr))

  handlePlotTypes(h5file, bKind, runType, flags, evDisplayRun = evDisplayRun)

  if cfProvideServer in flags:
    joinThread(thr)

when isMainModule:
  plotData()
