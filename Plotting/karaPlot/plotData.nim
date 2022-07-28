## A tool to plot data from H5 files
import plotly
import ggplotnim
import os except FileInfo
import strutils, strformat, times, sequtils, math, macros, algorithm, sets, stats, base64
import options, logging, typeinfo, json
import ws, asynchttpserver, asyncnet, asyncdispatch

import shell
import arraymancer
import zero_functional
import nimhdf5
import seqmath

import nimpy
import cligen
import chroma
import parsetoml

import dataset_helpers
import ingrid/[ingrid_types, tos_helpers, calibration]
import ingrid / calibration / calib_fitting
import helpers/utils
#import ingridDatabase / [databaseDefinitions, databaseRead]
import protocol, common_types

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
const Version = commitHash & " built on: " & currentDate

const
  GoldenMean = (sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
  FigWidth = 1200.0                      # width in pixels
  FigHeight = FigWidth * GoldenMean     # height in pixels

  InGridDsets = ["length", "width", "skewnessLongitudinal", "skewnessTransverse",
                 "kurtosisLongitudinal", "kurtosisTransverse", "rotationAngle",
                 "eccentricity", "fractionInTransverseRms", "lengthDivRmsTrans",
                 "rmsLongitudinal", "rmsTransverse", "hits", "energyFromPixel",
                 "energyFromCharge", "likelihood", "centerX", "centerY"]
  ToADsets = ["toaLength", "toaRms", "toaMean", "toaMin", "toaSkewness", "toaKurtosis"]
  FadcDsets = ["minvals", "fallTime", "riseTime"]
  AllFadcDsets = ["argMinval", "baseline", "eventNumber", "fallStop", "fallTime",
                  "minvals", "noisy", "riseStart", "riseTime"]

type
  ShapeKind = enum
    Rectangle, Square

## is the directory for the current input file to store the plots
var fileDir: string
var fileType: string

proc dataPath*(fileInfo: FileInfo, run, chip: int): grp_str =
  case fileInfo.tpaFileKind
  of tpkRawData:    result = rawPath(run, chip)
  of tpkReco:       result = recoPath(run, chip)
  of tpkLogL:       result = likelihoodPath(run, chip)
  of tpkTpx3Interp: result = tpx3Path(run, chip)
  of tpkTpx3Raw:
    doAssert false, "Invalid file kind to plot data from! " & $tpkTpx3Raw

proc dataBase*(fileInfo: FileInfo): string =
  case fileInfo.tpaFileKind
  of tpkRawData:    result = rawDataBase()
  of tpkReco:       result = recoBase()
  of tpkLogL:       result = likelihoodBase()
  of tpkTpx3Interp: result = tpx3Base()
  of tpkTpx3Raw:
    doAssert false, "Invalid file kind to plot data from! " & $tpkTpx3Raw

template buildFilter(key: string, tomlConfig, theSet: untyped): untyped =
  if tomlConfig["General"].hasKey(key):
    let allowed = tomlConfig["General"][key].getElems
    for el in allowed:
      theSet.incl el.getInt.uint16

proc contains*(dset: string, cuts: seq[GenericCut]): bool =
  for c in cuts:
    if dset == c.dset:
      return true

proc getCuts*(selector: DataSelector, file: string, dset: string): seq[GenericCut] =
  # only add those cuts where the file is valid & the dset are. Dset can
  # be overwritten by `applyAll`
  for c in selector.cuts:
    if (c.applyFile.len == 0 or file in c.applyFile) and
       (selector.applyAll or c.applyDset.len == 0 or dset in c.applyDset):
      result.add c

proc getMasks*(selector: DataSelector, file: string, dset: string): seq[MaskRegion] =
  # only add those masks where the file is valid & the dset are. Dset can
  # be overwritten by `applyAll`
  for m in selector.maskRegion:
    if (m.applyFile.len == 0 or file in m.applyFile) and
       (selector.applyAll or m.applyDset.len == 0 or dset in m.applyDset):
      result.add m

proc parseCut(tab: TomlValueRef): GenericCut =
  result = GenericCut(applyFile: tab["applyFile"].getElems.mapIt(it.getStr),
                      applyDset: tab["applyDset"].getElems.mapIt(it.getStr),
                      dset: tab["dset"].getStr,
                      min: tab["min"].getFloat,
                      max: tab["max"].getFloat)

proc parseMaskRegion(tab: TomlValueRef): MaskRegion =
  let x = tab["x"].getElems
  let y = tab["y"].getElems
  result = MaskRegion(applyFile: tab["applyFile"].getElems.mapIt(it.getStr),
                      applyDset: tab["applyDset"].getElems.mapIt(it.getStr),
                      x: (min: ChipCoord(x[0].getInt), max: ChipCoord(x[1].getInt)),
                      y: (min: ChipCoord(y[0].getInt), max: ChipCoord(y[1].getInt)))

proc parseTomlCuts(cfg: TomlValueRef): seq[GenericCut] =
  doAssert "Cuts" in cfg, "Cuts table not defined in your config.toml. Please add it."
  let cutTab = cfg["Cuts"]
  let cuts = cutTab["cuts"]
  for c in cuts:
    let cut = cutTab[$(c.getInt)]
    result.add parseCut(cut)

proc parseTomlMaskRegions(cfg: TomlValueRef): seq[MaskRegion] =
  doAssert "MaskRegions" in cfg, "Mask region table not defined in your config.toml. Please add it."
  let regionTab = cfg["MaskRegions"]
  let regions = regionTab["regions"]
  for r in regions:
    let region = regionTab[$(r.getInt)]
    result.add parseMaskRegion(region)

proc initConfig*(chips: set[uint16],
                 runs: set[uint16],
                 flags: set[ConfigFlagKind],
                 cuts: seq[GenericCut],
                 maskRegion: seq[MaskRegion],
                 tomlConfig: TomlValueRef,
                 ingridDsets: set[IngridDsetKind] = {},
                 fadcDsets: seq[string] = @[],
                 head: int = 0,
                 xDset: string = "",
                 yDset: string = "",
                 zDset: string = ""
                ): Config =
  let fType = tomlConfig["General"]["filetype"].getStr
  let outputType = tomlConfig["General"]["outputFormat"].getStr
  let combinedFormat = parseEnum[OutputFiletypeKind](outputType)
  let fileFormat = parseEnum[PlotFiletypeKind](fType)

  var cuts = cuts
  cuts.add parseTomlCuts(tomlConfig)

  var maskRegion = maskRegion
  maskRegion.add parseTomlMaskRegions(tomlConfig)

  ## XXX: remove this, set global
  fileType = fType

  var customPlots: seq[CustomPlot]
  if xDset.len > 0:
    customPlots = @[CustomPlot(x: xDset, y: yDset, color: zDset)]
  elif yDset.len > 0:
    doAssert false, "If a y axis is given, also supply an x axis! Otherwise only give x"

  ## Add the runs & chips from the configuration file to the user input
  var chips = chips
  var runs = runs
  buildFilter("allowedChips", tomlConfig, chips)
  buildFilter("allowedRuns", tomlConfig, runs)
  result = Config(flags: flags,
                  chips: chips,
                  runs: runs,
                  ingridDsets: ingridDsets,
                  fadcDsets: fadcDsets,
                  cuts: cuts,
                  maskRegion: maskRegion,
                  idxs: toSeq(0 ..< head),
                  fileType: fileFormat,
                  outputType: combinedFormat,
                  customPlots: customPlots)

proc initSelector(config: Config, cuts: seq[GenericCut] = @[],
                  chipRegion: ChipRegion = crAll,
                  applyAll: Option[bool] = none[bool]()
                 ): DataSelector =
  let appAll = if applyAll.isSome: applyAll.get
               else: cfApplyAllCuts in config.flags
  result = DataSelector(region: chipRegion,
                        cuts: concat(config.cuts, cuts),
                        maskRegion: config.maskRegion, #concat(config.cuts, maskRegion),
                        idxs: config.idxs,
                        applyAll: appAll)

## #############################################
## Globals
## #############################################

# global variable which stores the backend the user selected
var BKind*: PlottingBackendKind = bNone
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

# global storing the path to the config file
var ConfigFile = currentSourcePath().parentDir / "config.toml" # default in current dir

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

proc showPlot(p: PlotV, config: Config, fname: string) =
  ## If the `show` flag is set, show it
  case BKind
  of bMpl:
    discard callMethod(p.plt, "show")
  of bGgPlot:
    echo "Opening plot: ", fname
    case config.fileType
    of ftSvg:
      shell:
        inkview ($fname)
    of ftPdf:
      shell:
        evince ($fname)
    of ftPng:
      shell:
        nomacs ($fname)
  else:
    discard # nothing to show / already showing (plotly)

proc savePlot(p: PlotV, outfile: string, config: Config, fullPath = false) =
  # TODO: move elsewhere!!!
  var fname = outfile
  if not fullPath:
    fname = "figs" / outfile & ".svg"
  info &"Saving file: {fname}"
  case BKind
  of bPlotly:
    if not config.plotlySaveSvg:
      plotlyJson[outfile] = jsonPlotly(p)
    else:
      if not p.plPlot.isNil:
        p.plPlot.saveImage(fname)
        imageSet.incl(fname)
  of bMpl:
    if not p.plt.isNil:
      discard p.plt.savefig(fname)
      imageSet.incl(fname)
      if cfShow notin config.flags:
        discard callMethod(p.plt, "close", "all")
  of bGgPlot:
    if p.kind == bGgPlot: # if plot was not initialized
      if not p.invalid:
        try:
          p.pltGg.ggsave(fname, p.width, p.height)
        except AssertionError:
          # just continue here
          echo "WARNING: raised AssertionError trying to plot " & $fname
          discard
  else: discard
  if cfShow in config.flags:
    showPlot(p, config, fname)

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
    width = (FigHeight * 1.2).int
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

import sets

proc toIdx(x: float): int = (x / 14.0 * 256.0).round.int.clamp(0, 255)
proc applyMaskRegion(h5f: H5File, selector: DataSelector, dset: string,
                     group: H5Group, idx: seq[int]): seq[int] =
  ## applies the masking region in `maskRegion` of `selector`
  # build set of masked pixels from selector
  let masks = selector.getMasks(h5f.name.extractFilename, dset)
  if masks.len == 0:
    return idx # short cut to avoid data reading

  var noisyPixels = newSeq[(int, int)]()
  for maskRegion in masks:
    for x in maskRegion.x.min .. maskRegion.x.max:
      for y in maskRegion.y.min .. maskRegion.y.max:
        noisyPixels.add (x, y)
  if noisyPixels.len == 0: return idx
  let noiseSet = noisyPixels.toHashSet
  # read data
  let
    posX = h5f.readAs(group.name / "centerX", float)
    posY = h5f.readAs(group.name / "centerY", float)
  result = newSeqOfCap[int](idx.len)
  let idxs = if idx.len == 0: toSeq(0 ..< posX.len) else: idx
  for i in idxs:
    if (posX[i].toIdx, posY[i].toIdx) notin noiseSet:
      result.add i

proc hasCuts(selector: DataSelector, file: string, dset: H5Dataset): bool =
  ## returns whether we apply cuts to this (or all) datasets. To differentiate
  ## between an empty `idx` due to cuts removing everything or not using cuts
  # if the data selector has selected a region, we already have a cut
  if selector.region != crAll:
    result = true

  for c in selector.cuts:
    if (c.applyFile.len == 0 or file in c.applyFile) and
       (selector.applyAll or c.applyDset.len == 0 or dset.name in c.applyDset):
      result = true
      break
  for m in selector.maskRegion:
    if (m.applyFile.len == 0 or file in m.applyFile) and
       (selector.applyAll or m.applyDset.len == 0 or dset.name in m.applyDset):
      result = true
      break

proc toTuple(c: GenericCut): tuple[dset: string, lower, upper: float] =
  (dset: c.dset, lower: c.min, upper: c.max)

proc toTuple(c: seq[GenericCut]): seq[tuple[dset: string, lower, upper: float]] =
  result = c.mapIt(it.toTuple)

proc applyCuts(h5f: H5File, selector: DataSelector, dset: H5Dataset, idx: seq[int]): seq[int] =
  ## Apply potential cut to this dataset
  result = idx

  let parentGrp = h5f[dset.parent.grp_str]
  var performedCut = false
  if idx.len == 0: # only apply if no indices already specified
    let cuts = selector.getCuts(h5f.name.extractFilename, dset.name)
    # now apply the correct cuts
    # get the cut indices to only read passing data
    if cuts.len > 0:
      result = cutOnProperties(h5f, parentGrp, selector.region, cuts.toTuple)
      performedCut = true
    elif cuts.len == 0 and selector.region != crAll:
      result = cutOnProperties(h5f, parentGrp, selector.region)
      performedCut = true

  if selector.idxs.len > 0:
    if result.len == 0 and not performedCut:
      let maxIdx = min(selector.idxs.len, dset.shape[0])
      result = selector.idxs.filterIt(it < maxIdx)
    elif result.len > 0:
      let maxIdx = min(selector.idxs.len, result.len)
      let idxs = selector.idxs.filterIt(it < maxIdx)
      result = idxs.mapIt(result[it])
    # else leave as is
  # now possibly apply region cut
  if performedCut and result.len > 0:
    result = applyMaskRegion(h5f, selector, dset.name, parentGrp, result)
  elif not performedCut:
    result = applyMaskRegion(h5f, selector, dset.name, parentGrp, @[])

proc readFull*(h5f: H5File,
               fileInfo: FileInfo,
               runNumber: int,
               dsetName: string,
               selector: DataSelector,
               chipNumber = 0,
               isFadc = false,
               dtype: typedesc = float,
               idx: seq[int] = @[]): (seq[int], seq[dtype]) =
  ## reads the given dataset from the H5 file and returns it
  ## (potentially) converted to ``dtype``
  var dset: H5DataSet
  if isFadc:
    dset = h5f[(fadcRecoPath(runNumber) / dsetName).dset_str]
  else:
    dset = h5f[(fileInfo.dataPath(runNumber, chipNumber).string / dsetName).dset_str]

  ## possibly set indices to a subset based on cuts
  let idx = h5f.applyCuts(selector, dset, idx)
  when dtype is SomeNumber:
    if not selector.hasCuts(h5f.name.extractFilename, dset): # read everything and idx.len == 0:
      result[1] = dset.readAs(dtype)
    elif idx.len > 0:
      # manual conversion required
      result[1] = dset.readAs(idx, dtype)
    # else nothing to read, will remain empty
  else:
    type subtype = getInnerType(dtype)
    # NOTE: only support 2D seqs
    when subType isnot SomeNumber:
      raise newException(Exception, "Cannot convert N-D sequence to type: " &
        subtype.name)
    else:
      if not selector.hasCuts(h5f.name.extractFilename, dset): # read everything idx.len == 0:
        # convert to subtype and reshape to dsets shape
        result[1] = dset.readAs(subtype).reshape2D(dset.shape)
      elif idx.len > 0:
        # manual conversion required
        result[1] = dset[idx, subtype].reshape2D(dset.shape)
      # else nothing to read, will remain empty
  result[0] = idx

proc read*(h5f: H5File,
           fileInfo: FileInfo,
           runNumber: int,
           dsetName: string,
           selector: DataSelector,
           chipNumber = 0,
           isFadc = false,
           dtype: typedesc = float,
           idx: seq[int] = @[]): seq[dtype] =
  let (idxs, res) = readFull(h5f, fileInfo, runNumber, dsetName, selector, chipNumber, isFadc, dtype, idx)
  result = res

proc readVlen(h5f: H5File,
              fileInfo: FileInfo,
              runNumber: int,
              dsetName: string,
              selector: DataSelector,
              chipNumber = 0,
              dtype: typedesc = float,
              idx: seq[int] = @[]): seq[seq[dtype]] =
  ## reads variable length data `dsetName` and returns it
  ## In contrast to `read` this proc does *not* convert the data.
  let vlenDtype = special_type(dtype)
  let dset = h5f[(fileInfo.dataPath(runNumber, chipNumber).string / dsetName).dset_str]
  let idx = h5f.applyCuts(selector, dset, idx)
  if not selector.hasCuts(h5f.name.extractFilename, dset):
    result = dset[vlenDType, dtype]
  elif idx.len > 0:
    result = dset[vlenDtype, dtype, idx]
  # else nothing to read, remain empty

template applyFilter(theSet, fileInfo, field: untyped): untyped =
  if theSet.card > 0:
    fileInfo.field = fileInfo.field.filterIt(it.uint16 in theSet)

proc applyConfig(fileInfo: var FileInfo, config: Config) =
  ## Applies the given configuration to the `fileInfo` object, i.e.
  ## restricting the allowed runs & chips from the file.
  applyFilter(config.chips, fileInfo, chips)
  applyFilter(config.runs, fileInfo, runs)

proc appliedConfig(fileInfo: FileInfo, config: Config): FileInfo =
  ## returns a copy of the given `FileInfo` with the config.toml applied
  result = fileInfo
  result.applyConfig(config)

proc plotHist[T](xIn: seq[T], title, dset, outfile: string,
                 binS: float, binR: (float, float)): PlotV =
  ## plots the data in `x` as a histogram
  let xs = xIn.mapIt(it.float).filterIt(classify(it) notin {fcInf, fcNegInf, fcNaN})

  var binRange = binR
  var binSize = binS
  if binRange[0] == binRange[1]:
    binRange = (xs.percentile(3), xs.percentile(97))
  if binSize == 0.0:
    # use 100 bins as a backup
    binSize = (binRange[1] - binRange[0]) / 100

  info &"Bin range {binRange} for dset: {dset}"
  case BKind
  of bPlotly:
    result = initPlotV(title, dset, "#")
    var traces: seq[Trace[float]]
    traces.add Trace[float](`type`: PlotType.Histogram,
                            bins: binRange,
                            binSize: binSize,
                            #nbins: nBins,
                            xs: xs,
                            name: dset)
    result.plPlot = Plot[float](layout: result.plLayout, traces: traces)
  of bMpl:
    result = initPlotV(title, dset, "#")
    let nbins = ((binRange[1] - binRange[0]) / binSize).round.int
    discard result.ax.hist(xs,
                         bins = nbins,
                         range = binRange)
  of bGgPlot:
    let df = toDf(xs).filter(fn {float: `xs` >= binRange[0] and
                                            `xs` <= binRange[1]})
    result = initPlotV(title & " # entries: " & $df.len, dset, "#")
    result.pltGg = ggplot(df, aes("xs")) +
        geom_histogram(binWidth = binSize, hdKind = hdOutline) +
        margin(top = 2) +
        scale_x_continuous() + scale_y_continuous() +
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
      let ldf = toDf({ "bins" : binsIn[i],
                           "counts" : countsIn[i],
                           "dset" : constantColumn(dsets[i], binsIn[i].len) })
      df.add ldf
    if bins.len > 1:
      result.pltGg = ggplot(df, aes("bins", "counts")) +
          geom_histogram(aes(fill = "dset"), stat = "identity", hdKind = hdOutline) +
          scale_y_continuous() +
          result.theme # just add the theme directly
    else:
      # don't classify by `dset`
      result.pltGg = ggplot(df, aes("bins", "counts")) +
          geom_histogram(stat = "identity", hdKind = hdOutline) +
          scale_y_continuous() +
          result.theme # just add the theme directly

iterator chips(group: var H5Group): (int, H5Group) =
  ## returns all chip groups within the given Run group
  doAssert group.name.parentDir == "/reconstruction", "Bad group : " & group.name
  doAssert "run_" in group.name, "Bad group : " & group.name
  for grp in groups(group):
    if "chip_" in grp.name:
      let chipNum = grp.attrs["chipNumber", int]
      yield (chipNum, grp)

proc createCustomPlots(fileInfo: FileInfo, config: Config): seq[PlotDescriptor] =
  ## define any plot that doesn't fit any of the other descriptions here
  let selector = initSelector(config)
  for plt in config.customPlots:
    if plt.x.len > 0 and plt.y.len > 0:
      let customPlot = CustomPlot(kind: cpScatter,
                                  x: plt.x,
                                  y: plt.y,
                                  color: plt.color)
      result.add PlotDescriptor(runType: fileInfo.runType,
                                name: plt.x & "_vs_" & plt.y,
                                xLabel: plt.x,
                                yLabel: plt.y,
                                title: plt.x & "/" & plt.y,
                                selector: selector,
                                isCenterChip: true,
                                chip: -1, ## will be set to center chip when reading data!
                                runs: fileInfo.runs,
                                plotKind: pkCustomPlot,
                                customPlot: customPlot)
    elif plt.x.len > 0 and plt.y.len == 0:
      let customPlot = CustomPlot(kind: cpHistogram,
                                  x: plt.x,
                                  color: plt.color)
      result.add PlotDescriptor(runType: fileInfo.runType,
                                name: plt.x,
                                xLabel: plt.x,
                                yLabel: "Counts",
                                title: "Histogram: " & plt.x,
                                selector: selector,
                                isCenterChip: true,
                                chip: -1, ## will be set to center chip when reading data!
                                runs: fileInfo.runs,
                                plotKind: pkCustomPlot,
                                customPlot: customPlot)

  when false:
    block:
      result.add PlotDescriptor(runType: fileInfo.runType,
                                name: "riseTime",
                                selector: selector,
                                runs: fileInfo.runs,
                                plotKind: pkFadcDset,
                                binSize: 1.0,
                                binRange: (2.0, 502.0))
    block:
      result.add PlotDescriptor(runType: fileInfo.runType,
                                name: "fallTime",
                                selector: selector,
                                runs: fileInfo.runs,
                                plotKind: pkFadcDset,
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

proc histograms(h5f: H5File, runType: RunTypeKind,
                fileInfo: FileInfo,
                config: Config): seq[PlotDescriptor] =
  # TODO: perform cut on photo peak and escape peak, done by getting the energy
  # and performing a cut around peak +- 300 eV maybe
  # NOTE: need photo peak and escape peak plenty of times here. Just write wrapper
  # which takes energy and performs cuts, *then* creates plots
  var selectors: seq[DataSelector]
  let energyDset = "energyFromCharge" ## XXX: make configurable!!!
  case runType
  of rtCalibration:
    # creates 3 plots for the given dataset.
    # - 1 around photo peak
    # - 1 around escape peak
    # - 1 without cut
    if cfCutFePeak in config.flags:
      let Ed = "energyFromCharge"
      let ranges = @[GenericCut(dset: Ed, applyDset: @[Ed], min: 5.5, max: 6.2),
                     GenericCut(dset: Ed, applyDset: @[Ed], min: 2.7, max: 3.2)]
      selectors = ranges.mapIt(initSelector(config, @[it], applyAll = some(false)))
      selectors.add initSelector(config)
    else:
      selectors = @[initSelector(config)]
  else:
    # else just take all
    selectors = @[initSelector(config)]

  ## TODO: make datasets and chip regions selectable!
  for selector in selectors:
    if cfInGrid in config.flags:
      ## XXX: make region selection CLI argument & config file!
      for region in [crAll, crGold]: #ChipRegion:
        var sel = selector
        sel.region = region
        for ch in fileInfo.chips:
          for dset in concat(@InGridDsets, @ToADsets):
            let (binSize, binRange) = getBinSizeAndBinRange(dset)
            result.add PlotDescriptor(runType: runType,
                                      name: dset,
                                      selector: sel,
                                      runs: fileInfo.runs,
                                      plotKind: pkInGridDset,
                                      chip: ch,
                                      isCenterChip: fileInfo.centerChip == ch,
                                      binSize: binSize,
                                      binRange: binRange)
    if cfFadc in config.flags:
      for dset in FadcDsets:
        let (binSize, binRange) = getBinSizeAndBinRange(dset)
        result.add PlotDescriptor(runType: runType,
                                  name: dset,
                                  selector: selector,
                                  runs: fileInfo.runs,
                                  plotKind: pkFadcDset,
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
      ## continue if full event.
      ## TODO: replace by solution that also works for clusters!!
      if xEv.len >= 4095: continue

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
      x = newTensorUninit[int](NPix * NPix)
      y = newTensorUninit[int](NPix * NPix)
      z = newTensorUninit[float](NPix * NPix)
    var i = 0
    for idx, val in data:
      x[i] = idx[0]
      y[i] = idx[1]
      z[i] = val
      inc i
    let df = toDf(x, y, z)
    if data.max > 0.0:
      result.pltGg = ggplot(df, aes("x", "y", fill = "z")) +
          geom_raster() +
          scale_fill_continuous(scale = (low: 0.0, high: data.max)) +
          xlim(0, NPix) + ylim(0, NPix) +
          result.theme # just add the theme directly
    else:
      result.invalid = true
  else:
    discard

proc plotSparseEvent(df: DataFrame, title, outfile: string): PlotV =
  ## creates a 2D histogram plot (basically an image) of the given 2D
  ## Tensor
  result = initPlotV(title, "pixel x", "pixel y", ShapeKind.Rectangle)
  case BKind
  of bPlotly:
    discard
  of bMpl:
    discard
  of bGgPlot:
    result.pltGg = ggplot(df, aes("x", "y", color = "ch")) +
      geom_point() +
      margin(left = 9.0, right = 3, top = 2) +
      xlim(0, NPix) + ylim(0, NPix) +
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
    let df = toDf(x, y)
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

proc plotCustomScatter(x, y: seq[float],
                       pd: PlotDescriptor,
                       #name: string,
                       title: string,
                       z: seq[float] = @[] # z is the coloring
                      ): PlotV =
  result = initPlotV(title, pd.xlabel, pd.ylabel, ShapeKind.Rectangle)
  case BKind
  of bPlotly:
    # in this case plot is defined
    let trFit = Trace[float](`type`: PlotType.Scatter,
                             xs: x,
                             ys: y,
                             name: title)
    result.plPlot.traces.add @[trFit]
  of bMpl:
    # in this case `ax` is defined
    discard result.ax.plot(x,
                         y,
                         label = title,
                         color = "r")
  of bGgPlot:
    if z.len > 0:
      let df = toDf(x, y, z)
      let pSize = if df.len > 1_000: some(1.5) else: some(3.0)
      result.pltGg = ggplot(df, aes("x", "y", color = "z")) +
        geom_point(size = pSize) +
        result.theme # just add the theme directly
    else:
      let df = toDf(x, y)
      let pSize = if df.len > 1_000: some(1.5) else: some(3.0)
      result.pltGg = ggplot(df, aes("x", "y")) +
        geom_point(size = pSize) +
        result.theme # just add the theme directly
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
proc occupancies(h5f: H5File, runType: RunTypeKind,
                 fileInfo: FileInfo,
                 config: Config): seq[PlotDescriptor] =
  ## creates occupancy plots for the given HDF5 file, iterating over
  ## all runs and chips
  const
    quantile = 95
    clamp3 = 10.0
  let selector = initSelector(config)
  for ch in fileInfo.chips:
    let basePd = PlotDescriptor(runType: runType,
                                name: "occupancy",
                                selector: selector,
                                runs: fileInfo.runs,
                                plotKind: pkOccupancy,
                                chip: ch,
                                isCenterChip: fileInfo.centerChip == ch)
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
      clampQ = 95
    result.add @[fullPd, clampAPd, clampQPd, clusterPd, clusterClampPd]

#proc plotPolyas(h5f: H5File, group: H5Group,
#                chipNum: int, runNumber: string): seq[Trace[float]] =
#  ## perform the plots of the polya distribution again (including the fits)
#  ## and return the polya data as to allow a combined plot afterwards
#  let
#    polya = h5f.read(runNumber.parseInt, fileInfo, "polya", chipNum,
#                     dtype = seq[float])
#    polyaFit = h5f.read(runNumber.parseInt, fileInfo, "polyaFit", chipNum,
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

proc polya(h5f: H5File, runType: RunTypeKind,
           fileInfo: FileInfo,
           config: Config): seq[PlotDescriptor] =
  ## creates the plots for the polya distribution of the datasets
  let selector = initSelector(config)
  for ch in fileInfo.chips:
    result.add PlotDescriptor(runType: runType,
                              name: "polya",
                              selector: selector,
                              xlabel: "Number of electrons",
                              runs: fileInfo.runs,
                              chip: ch,
                              isCenterChip: fileInfo.centerChip == ch,
                              plotKind: pkPolya)
  if fileInfo.chips.len > 1:
    result.add PlotDescriptor(runType: runType,
                              name: "polya",
                              selector: selector,
                              runs: fileInfo.runs,
                              plotKind: pkCombPolya,
                              chipsCP: fileInfo.chips)

proc totPerPixel(h5f: H5File, runType: RunTypeKind,
                 fileInfo: FileInfo,
                 config: Config): seq[PlotDescriptor] =
  ## creates the plots for the ToT per pixel distribution of the datasets
  ## i.e. the raw ToT histogram (equivalent to polya for raw ToT)
  let selector = initSelector(config)
  let (binSize, binRange) = getBinSizeAndBinRange("totPerPixel")
  for ch in fileInfo.chips:
    result.add PlotDescriptor(runType: runType,
                              name: "totPerPixel",
                              selector: selector,
                              xlabel: "ToT [clock cycles]",
                              runs: fileInfo.runs,
                              chip: ch,
                              isCenterChip: fileInfo.centerChip == ch,
                              plotKind: pkToTPerPixel,
                              binSize: binSize,
                              binRange: binRange)

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
    raise newException(Exception, "`dtm.datetime` is broken. Hence a date plot on matplotlib is currently unsupported!")
    #let xP = x.mapIt(dtm.datetime.utcfromtimestamp(int(it)))
    #discard result.ax.plot_date(xP, y, label = title)
  of bGgPlot:
    let df = toDf(x, y)
    result.pltGg = ggplot(df, aes("x", "y")) +
        geom_point() +
        result.theme # just add the theme directly
  else:
    discard

proc feSpectrum(h5f: H5File, runType: RunTypeKind,
                fileInfo: FileInfo,
                config: Config): seq[PlotDescriptor] =
  ## creates the plot of the Fe spectrum
  let selector = initSelector(config)
  for r in fileInfo.runs:
    let basePd = PlotDescriptor(runType: runType,
                                name: "FeSpectrumPlot",
                                selector: selector,
                                xlabel: "# pixels",
                                runs: @[r],
                                chip: fileInfo.centerChip,
                                isCenterChip: true,
                                plotKind: pkFeSpec)
    #let energyPd = replace(basePd):
    #  name = "EnergyCalib"
    #  xlabel = "# pixels"
    #  ylabel = "Energy / keV"
    #  plotKind = pkEnergyCalib
    #let feChargePd = replace(basePd):
    #  name = "FeSpectrumChargePlot"
    #  xlabel = "Charge / 10^3 electrons"
    #  plotKind = pkFeSpecCharge
    #let energyChargePd = replace(basePd):
    #  name = "EnergyCalibCharge"
    #  xlabel = "Charge / 10^3 electrons"
    #  ylabel = "Energy / keV"
    #  plotKind = pkEnergyCalibCharge
    result.add @[basePd] #, energyPd, feChargePd, energyChargePd]

  ## TODO: turn these into a separate thing from the `feSpectrum`
  #let photoVsTime = PlotDescriptor(runType: runType,
  #                                 name: "PhotoPeakVsTime",
  #                                 xlabel: "Time / unix",
  #                                 selector: selector,
  #                                 runs: fileInfo.runs,
  #                                 chip: fileInfo.centerChip,
  #                                 isCenterChip: true,
  #                                 plotKind: pkFeVsTime)
  #let photoChVsTime = PlotDescriptor(runType: runType,
  #                                   name: "PhotoPeakChargeVsTime",
  #                                   xlabel: "Time / unix",
  #                                   selector: selector,
  #                                   runs: fileInfo.runs,
  #                                   chip: fileInfo.centerChip,
  #                                   isCenterChip: true,
  #                                   plotKind: pkFeChVsTime)
  #let phPixDivChVsTime = PlotDescriptor(runType: runType,
  #                                 name: "PhotoPixDivChVsTime",
  #                                 xlabel: "Time / unix",
  #                                 selector: selector,
  #                                 runs: fileInfo.runs,
  #                                 chip: fileInfo.centerChip,
  #                                 isCenterChip: true,
  #                                 plotKind: pkFePixDivChVsTime)
  #let photoVsTimeHalfH = PlotDescriptor(runType: runType,
  #                                      name: "PhotoPeakVsTimeHalfHour",
  #                                      xlabel: "Time / unix",
  #                                      selector: selector,
  #                                      runs: fileInfo.runs,
  #                                      chip: fileInfo.centerChip,
  #                                      isCenterChip: true,
  #                                      plotKind: pkFeVsTime,
  #                                      splitBySec: 1800,
  #                                      lastSliceError: 0.2,
  #                                      dropLastSlice: false)
  #let photoChVsTimeHalfH = PlotDescriptor(runType: runType,
  #                                        name: "PhotoPeakChargeVsTimeHalfHour",
  #                                        xlabel: "Time / unix",
  #                                        selector: selector,
  #                                        runs: fileInfo.runs,
  #                                        chip: fileInfo.centerChip,
  #                                        isCenterChip: true,
  #                                        plotKind: pkFeChVsTime,
  #                                        splitBySec: 1800,
  #                                        lastSliceError: 0.2,
  #                                        dropLastSlice: false)
  #let phPixDivChVsTimeHalfH = PlotDescriptor(runType: runType,
  #                                           name: "PhotoPixDivChVsTimeHalfHour",
  #                                           xlabel: "Time / unix",
  #                                           selector: selector,
  #                                           runs: fileInfo.runs,
  #                                           chip: fileInfo.centerChip,
  #                                           isCenterChip: true,
  #                                           plotKind: pkFePixDivChVsTime,
  #                                           splitBySec: 1800,
  #                                           lastSliceError: 0.2,
  #                                           dropLastSlice: false)
  #
  #result.add @[photoVsTime, phPixDivChVsTime,
  #             photoChVsTime, photoChVsTimeHalfH,
  #             photoVsTimeHalfH, phPixDivChVsTimeHalfH]

proc fePhotoDivEscape(h5f: H5File, runType: RunTypeKind,
                      fileInfo: FileInfo,
                      config: Config
                     ): seq[PlotDescriptor] =
  let selector = initSelector(config)
  result = @[PlotDescriptor(runType: runType,
                            name: "Photopeak / Escape peak",
                            xlabel: "Time / unix",
                            selector: selector,
                            runs: fileInfo.runs,
                            chip: fileInfo.centerChip,
                            isCenterChip: true,
                            plotKind: pkFePhotoDivEscape)]

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

#proc readEvents*(h5f: H5File, run, chip: int, idx: seq[int],
#                config: Config): Tensor[float] =
#  ## helper proc to read data for given indices of events `idx` and
#  ## builds a tensor of `idx.len` events.
#  let
#    x = h5f.readVlen(fileInfo, run, "x", pd.
#                     chipNumber = chip,
#                     dtype = uint8, idx = idx)
#    y = h5f.readVlen(fileInfo, run, "y", config,
#                     chipNumber = chip,
#                     dtype = uint8, idx = idx)
#    ch = h5f.readVlen(fileInfo, run, "ToT", config,
#                      chipNumber = chip,
#                      dtype = uint16, idx = idx)
#  result = buildEvents(x, y, ch, toOrderedSet(@[0]))

proc readEventsSparse*(h5f: H5File, fileInfo: FileInfo, run, chip: int, #idx: int,
                       selector: DataSelector): DataFrame =
  ## helper proc to read data for given indices of events `idx` and
  ## builds a tensor of `idx.len` events.
  let
    x = h5f.readVlen(fileInfo, run, "x", selector,
                     chipNumber = chip,
                     dtype = uint8)
    y = h5f.readVlen(fileInfo, run, "y", selector,
                     chipNumber = chip,
                     dtype = uint8)
    ch = h5f.readVlen(fileInfo, run, "ToT", selector,
                      chipNumber = chip,
                      dtype = uint16)
    (events, _) = h5f.readFull(fileInfo, run, "eventNumber", selector,
                               chipNumber = chip,
                               dtype = int)
  result = newDataFrame()
  doAssert events.len == x.len, "Events are: " & $events.len & " and x.len " & $x.len
  for i in 0 ..< x.len:
    let dfLoc = toDf({ "x" : x[i], "y" : y[i], "ch" : ch[i],
                           "Index" : events[i] })
    echo "Adding df ", dfLoc
    result.add dfLoc

proc readIngridForEventDisplay*(h5f: H5File, fileInfo: FileInfo, run, chip: int, #idx: int,
                                selector: DataSelector): DataFrame =
  ## helper proc to read data for given indices of events `idx` and
  ## builds a tensor of `idx.len` events.
  result = newDataFrame()
  for dset in concat(@InGridDsets, @ToADsets):
    echo "Attpmting to read: ", dset
    try:
      result[dset] = h5f.read(fileInfo, run, dset, selector,
                              chipNumber = chip,
                              dtype = float)
    except KeyError:
      echo "Could not read dataset ", dset
      continue # drop this key

proc createEventDisplayPlots(h5f: H5File,
                             run: int,
                             runType: RunTypeKind,
                             fileInfo: FileInfo,
                             config: Config,
                             events: OrderedSet[int]): seq[PlotDescriptor] =
  let selector = initSelector(config)
  for ch in fileInfo.chips:
    #for ev in events:
    result.add PlotDescriptor(runType: runType,
                              name: "EventDisplay",
                              selector: selector,
                              xlabel: "x",
                              ylabel: "y",
                              runs: @[run],
                              chip: ch,
                              isCenterChip: fileInfo.centerChip == ch,
                              plotKind: pkInGridEvent)
                              #event: ev)

proc createFadcPlots(h5f: H5File,
                     run: int,
                     runType: RunTypeKind,
                     fileInfo: FileInfo,
                     config: Config,
                     events: OrderedSet[int]): seq[PlotDescriptor] =
  let selector = initSelector(config)
  for ev in events:
    result.add PlotDescriptor(runType: runType,
                              name: "EventDisplay FADC",
                              selector: selector,
                              xlabel: "Clock cycles FADC",
                              ylabel: "U / V",
                              runs: @[run],
                              chip: fileInfo.centerChip,
                              isCenterChip: true,
                              plotKind: pkFadcEvent,
                              event: ev)

proc createOuterChipHistograms*(h5f: H5File,
                               runType: RunTypeKind,
                               fileInfo: FileInfo,
                               config: Config,
                               cutRange: GenericCut): seq[PlotDescriptor] =
  let selector = initSelector(config, @[cutRange])
  var chips = fileInfo.chips
  let oldLen = chips.len
  chips.delete(chips.find(fileInfo.centerChip))
  doAssert chips.len + 1 == oldLen
  result.add PlotDescriptor(runType: runType,
                            name: &"Outer chips # pix hit {runType}",
                            xlabel: "# pixels on outer chips",
                            ylabel: "#",
                            selector: selector,
                            runs: fileInfo.runs,
                            outerChips: chips,
                            plotKind: pkOuterChips)

proc createInGridFadcEvDisplay(h5f: H5File,
                               run: int,
                               runType: RunTypeKind,
                               fileInfo: FileInfo,
                               config: Config,
                               events: OrderedSet[int]): seq[PlotDescriptor] =
  let selector = initSelector(config)
  var finfo = fileInfo
  finfo.chips = @[finfo.centerChip]
  let ingridPds = createEventDisplayPlots(h5f, run, runType, fInfo, config, events)
  let ingridDomain = (left: 0.0, bottom: 0.0, width: 0.45, height: 1.0)
  let fadcPds = createFadcPlots(h5f, run, runType, fInfo, config, events)
  let fadcDomain = (left: 0.525, bottom: 0.05, width: 0.575, height: 0.6)
  for tup in zip(ingridPds, fadcPds):
    let
      ingrid = tup[0]
      fadc = tup[1]
    result.add PlotDescriptor(runType: runType,
                              name: "Event display InGrid + FADC",
                              runs: @[run],
                              selector: selector,
                              chip: fileInfo.centerChip,
                              isCenterChip: true,
                              plotKind: pkSubPlots,
                              plots: @[ingrid, fadc],
                              domain: @[ingridDomain, fadcDomain])

proc percentile[T](t: Tensor[T], perc: float): float =
  let dataSorted = t.reshape(t.size.int).sorted
  let perIdx = min((t.size.float * perc).round.int, t.size - 1)
  result = dataSorted[perIdx]

proc clampedOccupancy[T](x, y: seq[T], pd: PlotDescriptor): Tensor[float] =
  ## calculates the occupancy given `x`, `y`, which may be seqs of clusters
  ## or seqs of center positions
  result = calcOccupancy(x, y)
  case pd.clampKind
  of ckAbsolute:
    result = result.clamp(0.0, pd.clampA)
  of ckQuantile:
    let quant = result.percentile(pd.clampQ / 100.0)
    result = result.clamp(0.0, quant)
  else: discard

proc readPlotFit(h5f: H5File, fileInfo: FileInfo, pd: PlotDescriptor):
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
    if pd.selector.cuts.len > 0:
      echo "WARN: Cutting on fit datasets. This may have unintended consequences!"
    let
      data = h5f.read(fileInfo, r, pd.name, pd.selector, pd.chip,
                       dtype = seq[float])
      dataFit = h5f.read(fileInfo, r, pd.name & "Fit", pd.selector, pd.chip,
                          dtype = seq[float])
    let nbins = data.shape[0]
    let (bins, counts) = data.split(SplitSeq.Seq2Col)
    let (binsFit, countsFit) = dataFit.split(SplitSeq.Seq2Col)
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

template createFeVsTimeDataFrame(h5f: H5File,
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
    let dfEv = toDf({ "eventNumber" : evNum, "timestamp" : tstamp })
    path &= "/chip_" & $pd.chip
    let feSpec = h5f[path / "FeSpectrum", int]
    let feSpecCh = h5f[path / "FeSpectrumCharge", float]
    let feSpecEvNum = h5f[path / "FeSpectrumEvents", int]
    let dfFeSpec = toDf({ "eventNumber" : feSpecEvNum,
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

proc handleInGridDset(h5f: H5File,
                      fileInfo: FileInfo,
                      pd: PlotDescriptor,
                      config: Config): (string, PlotV) =
  var allData: seq[float]
  for r in pd.runs:
    allData.add h5f.read(fileInfo, r, pd.name, pd.selector, pd.chip, dtype = float)
    # perform cut on range
  result[0] = buildOutfile(pd, fileDir, fileType)
  let title = buildTitle(pd)
  result[1] = plotHist(allData, title, pd.name, result[0], pd.binSize, pd.binRange)

proc handleFadcDset(h5f: H5File,
                    fileInfo: FileInfo,
                    pd: PlotDescriptor,
                    config: Config): (string, PlotV) =
  # get the center chip group
  var allData: seq[float]
  for r in pd.runs:
    let group = h5f[fileInfo.dataPath(r, fileInfo.centerChip)]
    let inGridSet = h5f.read(fileInfo, r, "eventNumber", pd.selector, fileInfo.centerChip, dtype = int)
      .toSet()
    let evNumFadc = h5f.read(fileInfo, r, "eventNumber", pd.selector,
                             isFadc = true, dtype = int) # [group.name / "eventNumber", int]
    let idxFadc = (toSeq(0 .. evNumFadc.high)) --> filter(evNumFadc[it] in inGridSet)
    let data = h5f.read(fileInfo, r, pd.name, pd.selector, pd.chip, isFadc = true, dtype = float)
    allData.add idxFadc --> map(data[it])
  result[0] = buildOutfile(pd, fileDir, fileType)
  let title = buildTitle(pd)
  result[1] = plotHist(allData, title, pd.name, result[0], pd.binSize, pd.binRange)

proc handleOccupancy(h5f: H5File,
                     fileInfo: FileInfo,
                     pd: PlotDescriptor,
                     config: Config): (string, PlotV) =
  # get x and y datasets, stack and get occupancies
  let vlenDtype = special_type(uint8)
  var occFull = newTensor[float]([NPix, NPix])
  for r in pd.runs:
    let
      x = h5f.readVlen(fileInfo, r, "x", pd.selector, pd.chip, dtype = uint8)
      y = h5f.readVlen(fileInfo, r, "y", pd.selector, pd.chip, dtype = uint8)
    let occ = clampedOccupancy(x, y, pd)
    # stack this run onto the full data tensor
    occFull = occFull .+ occ
  let title = buildTitle(pd)
  result[0] = buildOutfile(pd, fileDir, fileType)
  result[1] = plotHist2D(occFull, title, result[0])

proc handleOccCluster(h5f: H5File,
                      fileInfo: FileInfo,
                      pd: PlotDescriptor,
                      config: Config): (string, PlotV) =
  # plot center positions
  var occFull = newTensor[float]([NPix, NPix])
  for r in pd.runs:
    let
      group = h5f[fileInfo.dataPath(r, pd.chip)]
      # get centers and rescale to 256 max value
      centerX = h5f[(group.name / "centerX"), float].mapIt(it * NPix.float / TimepixSize)
      centerY = h5f[(group.name / "centerY"), float].mapIt(it * NPix.float / TimepixSize)
    let occ = clampedOccupancy(centerX, centerY, pd)
    # stack this run onto the full data tensor
    occFull = occFull .+ occ
  let title = buildTitle(pd)
  result[0] = buildOutfile(pd, fileDir, fileType)
  result[1] = plotHist2D(occFull, title, result[0])

proc handleBarScatter(h5f: H5File,
                      fileInfo: FileInfo,
                      pd: PlotDescriptor,
                      config: Config): (string, PlotV) =
  result[0] = buildOutfile(pd, fileDir, fileType)
  let xlabel = pd.xlabel
  let title = buildTitle(pd)
  let (bins, counts, binsFit, countsFit) = h5f.readPlotFit(fileInfo, pd)
  result[1] = plotBar(@[bins], @[counts], title, xlabel, @[title], result[0])
  # now add fit to the existing plot
  let nameFit = &"Fit of chip {pd.chip}"
  result[1].plotScatter(binsFit, countsFit, nameFit, result[0])

proc handleCombPolya(h5f: H5File,
                     fileInfo: FileInfo,
                     pd: PlotDescriptor,
                     config: Config): (string, PlotV) =
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
    let (bins, counts, binsFit, countsFit) = h5f.readPlotFit(fileInfo, localPd)
    binsSeq.add bins
    countsSeq.add counts
    dsets.add "Chip " & $ch
  let xlabel = "Number of electrons"
  result[1] = plotBar(binsSeq, countsSeq, title, xlabel, dsets, result[0])

proc handleFeVsTime(h5f: H5File,
                    fileInfo: FileInfo,
                    pd: PlotDescriptor,
                    config: Config): (string, PlotV) =
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
    let group = h5f[(fileInfo.dataBase & $r).grp_str]
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

proc handleFePixDivChVsTime(h5f: H5File,
                            fileInfo: FileInfo,
                            pd: PlotDescriptor,
                            config: Config): (string, PlotV) =
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
    let group = h5f[(fileInfo.dataBase & $r).grp_str]
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

proc handleFePhotoDivEscape(h5f: H5File,
                            fileInfo: FileInfo,
                            pd: PlotDescriptor,
                            config: Config): (string, PlotV) =
  const kalphaCharge = 3 # for the ``amplitude``! not the mean position
  const escapeCharge = 0
  const parPrefix = "p"
  const dateStr = "yyyy-MM-dd'.'HH:mm:ss" # example: 2017-12-04.13:39:45
  var
    photoSeq: seq[float]
    escapeSeq: seq[float]
    dates: seq[float] #string]#Time]
  let title = buildTitle(pd)
  result[0] = buildOutfile(pd, fileDir, fileType)
  for r in pd.runs:
    let group = h5f[(fileInfo.dataBase & $r).grp_str]
    let chpGrpName = group.name / "chip_" & $pd.chip
    let dset = h5f[(chpGrpName / "FeSpectrumCharge").dset_str]
    photoSeq.add dset.attrs[parPrefix & $kalphaCharge, float]
    escapeSeq.add dset.attrs[parPrefix & $escapeCharge, float]
    dates.add parseTime(group.attrs["dateTime", string],
                        dateStr,
                        utc()).toUnix.float
  # now plot
  # calculate ratio and convert to string to workaround plotly limitation of
  # only one type for Trace
  let ratio = zip(photoSeq, escapeSeq) --> map(it[0] / it[1])
  result[1] = plotDates(dates, ratio,
                        title = title,
                        xlabel = pd.xlabel,
                        outfile = result[0])

## XXX: Outer chip is broken, as our `DataSelector` doesn't yet allow to select data
## on *another* chip!

#proc handleOuterChips(h5f: H5File,
#                      fileInfo: FileInfo,
#                      pd: PlotDescriptor,
#                      config: Config): (string, PlotV) =
#  var data: seq[int]
#  const
#    rmsTransHigh = 1.2
#    eccHigh = 1.3
#
#  for r in pd.runs:
#    # cut on `centerChip` for rms and eccentricity
#    let
#      grp = h5f[fileInfo.dataPath(r, fileInfo.centerChip)]
#      idx = cutOnProperties(h5f, grp,
#                            ("eccentricity", -Inf, eccHigh),
#                            ("rmsTransverse", -Inf, rmsTransHigh),
#                            ("energyFromCharge", pd.rangeCenter.low, pd.rangeCenter.high))
#      evNumCenterAll = h5f.read(fileInfo, r, "eventNumber", pd.selector, fileInfo.centerChip, dtype = int)
#      evNumCenter = toSet(idx.mapIt(evNumCenterAll[it]))
#    for c in pd.outerChips:
#      let
#        hitsAll = h5f.read(fileInfo, r, "hits", pd.selector, c, dtype = int)
#        evNumAll = h5f.read(fileInfo, r, "eventNumber", pd.selector, c, dtype = int)
#      # now add all hits passing
#      for i, ev in evNumAll:
#        if ev in evNumCenter and hitsAll[i] > 3:
#          data.add hitsAll[i]
#  result[0] = pd.title
#  const binSize = 3.0
#  const binRange = (0.0, 400.0)
#  let outfile = buildOutfile(pd, fileDir, fileType)
#  result[1] = plotHist(data, pd.title, pd.name, pd.title, binSize, binRange)

proc handleToTPerPixel(h5f: H5File,
                       fileInfo: FileInfo,
                       pd: PlotDescriptor,
                       config: Config): (string, PlotV) =
  result[0] = buildOutfile(pd, fileDir, fileType)
  let xlabel = pd.xlabel
  let title = buildTitle(pd)
  var tots: seq[uint16]
  for run in pd.runs:
    tots.add flatten(h5f.readVlen(fileInfo, run, "ToT", pd.selector,
                                  chipNumber = pd.chip,
                                  dtype = uint16))
  # if tots longer than `10_000_000`, compute in batches
  const batchSize = 10_000_000
  if tots.len > batchSize:
    var hist: seq[float]
    var bins: seq[float]
    for batch in 0 ..< (tots.len.float / batchSize).ceil.int:
      let lower = batch * batchSize
      let upper = min((batch + 1) * batchSize, tots.len)
      echo "lower ", lower, " and higher ", upper, " and totslen ", tots.len
      var (h, b) = histogram(tots[lower ..< upper].mapIt(it.int),
                             bins = ((pd.binRange[1] - pd.binRange[0]) / pd.binSize).round.int,
                             range = (pd.binRange[0], pd.binRange[1]))
      if bins.len == 0:
        bins = b
        hist = newSeq[float](bins.len)
      for i, el in h:
        hist[i] += el.float
    result[1] = plotBar(@[bins], @[hist], title, xlabel, @[pd.name], result[0])
  else:
    result[1] = plotHist(tots.mapIt(it.int), title, pd.name, result[0],
                         binS = pd.binSize, binR = pd.binRange)

iterator ingridEventIter(h5f: H5File,
                         fileInfo: FileInfo,
                         #pds: seq[PlotDescriptor],
                         pd: PlotDescriptor,
                         config: Config): (string, PlotV) =
  #var events {.global.}: seq[int]
  # all PDs are guaranteed to be from  the same run!
  let run = pd.runs[0] #pds[0].runs[0]
  # TODO: check if all PDs necessarily have same chip. Will be the case for
  # InGrid = FADC, but not necessarily only chip?
  let chip = pd.chip #pds[0].chip

  when false:
    var lastRun {.global.} = 0
    var evNums {.global.}: seq[int]
    var evTab {.global.}: Table[int, int]
    if run != lastRun:
      evNums = h5f.read(fileInfo, run, "eventNumber", pd.selector, dtype = int,
                        chipNumber = chip)
      evTab = initTable[int, int]()
      for i, ev in evNums:
        evTab[i] = ev
        events.add ev
      lastRun = run

    events.setLen(0)
    for pd in pds:
      events.add evTab[pd.event]

    let ev = readEvent(h5f, run, chip, events, config)
    for i, pd in pds:
      discard
  let events = readEventsSparse(h5f, fileInfo, run, chip, pd.selector)
  let dfDsets = readIngridForEventDisplay(h5f, fileInfo, run, chip, pd.selector)
  var idx = 0
  for (tup, subDf) in groups(group_by(events, "Index")):
    let event = tup[0][1].toInt

    echo "Generating event display for ", tup, " is ", subDf
    ## XXX: HACK: change to sometingn better to generate filenames
    var pd = pd
    pd.event = event
    let title = buildTitle(pd)
    var outfile = buildOutfile(pd, fileDir, fileType)
    var pltV = plotSparseEvent(subDf, title, outfile)

    var texts: seq[string]
    if pd.selector.cuts.len > 0:
      texts.add &"{\" \":>15}Cuts used: "
      for c in pd.selector.cuts:
        texts.add &"{c.dset:>15}: [{c.min:.2f}, {c.max:.2f}]"
    if pd.selector.maskRegion.len > 0:
      texts.add &"{\" \":>15}Masks used: "
      for m in pd.selector.maskRegion:
        texts.add &"    [x: ({m.x.min}, {m.x.max}), y: ({m.y.min}, {m.y.max})]"
    for dset in concat(@InGridDsets, @ToADsets):
      try:
        let val = dfDsets[dset, idx, float]
        let s = &"{dset:>25}: {val:6.4f}"
        texts.add s
      except KeyError:
        echo "Ignoring missing key: ", dset, " for annotation!"
        continue # ignore this field
    echo "Have texts : ", texts
    case BKind
    of bPlotly:
      for i, a in texts:
        pltV.plPlot.layout.annotations.add plotly.Annotation(x: 0.0,
                                                             y: 0.0,
                                                             xshift: 800.0,
                                                             yshift: 600.0 - (i.float * 20.0),
                                                             text: texts[i])
    of bGgPlot:
      # create a font to use using the `ggplotnim.font` helper
      let font = font(10.0, family = "monospace")
      for i, a in texts:
        pltV.pltGg.annotations.add ggplotnim.Annotation(left: some(-0.3),
                                                        bottom: some(0.15 + i.float * 0.03),
                                                        font: font,
                                                        text: texts[i])
    else:
      echo "InGrid Event property annotations not supported on " &
        "Matplotlib backend yet!"
    pltV.annotations = texts
    echo "Plot has annotations: ", pltV.annotations
    yield (outfile, pltV)
    inc idx

iterator fadcEventIter(h5f: H5File,
                       fileInfo: FileInfo,
                       pds: seq[PlotDescriptor],
                       config: Config): (string, PlotV) =
  var events {.global.}: seq[int]
  # all PDs are guaranteed to be from  the same run!
  let run = pds[0].runs[0]

  var lastRun {.global.} = 0
  var evNums {.global.}: seq[int]
  var evTab {.global.}: Table[int, int]
  if run != lastRun:
    evNums = h5f.read(fileInfo, run, "eventNumber", pds[0].selector, dtype = int,
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
      let val = h5f.read(fileInfo, run, d, pd.selector, isFadc = true,
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

proc handleIngridEvent(h5f: H5File,
                       fileInfo: FileInfo,
                       pd: PlotDescriptor,
                       config: Config): (string, PlotV) =
  doAssert pd.plotKind == pkInGridEvent
  for outfile, pltV in ingridEventIter(h5f, fileInfo, pd, config):
    # result will be last event, HACK
    result = (outfile, pltV)
    # call save here
    ## XXX: HACK!!!
    info &"Calling savePlot for {pd.plotKind} with filename {result[0]}"
    try:
      savePlot(result[1], result[0], config, fullPath = true)
    except Exception as e:
      echo "Failed to generate plot with error ", e.msg
      continue

proc handleFadcEvent(h5f: H5File,
                     fileInfo: FileInfo,
                     pd: PlotDescriptor,
                     config: Config): (string, PlotV) =
  doAssert pd.plotKind == pkFadcEvent
  for outfile, pltV in fadcEventIter(h5f, fileInfo, @[pd], config):
    # only a single pd
    result = (outfile, pltV)

proc handleCustomPlot(h5f: H5File,
                      fileInfo: FileInfo,
                      pd: PlotDescriptor,
                      config: Config): (string, PlotV) =
  doAssert pd.plotKind == pkCustomPlot
  if pd.processData.isNil:
    ## Assume this is a scatter plot
    let cPlt = pd.customPlot
    case cPlt.kind
    of cpScatter:
      var allX: seq[float]
      var allY: seq[float]
      var allZ: seq[float]
      for r in pd.runs:
        allX.add h5f.read(fileInfo, r, cPlt.x, pd.selector, dtype = float,
                         chipNumber = fileInfo.centerChip)
        allY.add h5f.read(fileInfo, r, cPlt.y, pd.selector, dtype = float,
                         chipNumber = fileInfo.centerChip)
        if cPlt.color.len > 0:
          allZ.add h5f.read(fileInfo, r, cPlt.color, pd.selector, dtype = float,
                            chipNumber = fileInfo.centerChip)
      let title = buildTitle(pd)
      result[0] = buildOutfile(pd, fileDir, fileType)
      result[1] = plotCustomScatter(allX, allY, pd, title, allZ)
    of cpHistogram:
      var allX: seq[float]
      var allZ: seq[float]
      for r in pd.runs:
        allX.add h5f.read(fileInfo, r, cPlt.x, pd.selector, dtype = float,
                         chipNumber = fileInfo.centerChip)
        if cPlt.color.len > 0:
          allZ.add h5f.read(fileInfo, r, cPlt.color, pd.selector, dtype = float,
                            chipNumber = fileInfo.centerChip)
      let title = buildTitle(pd)
      result[0] = buildOutfile(pd, fileDir, fileType)
      result[1] = plotHist(allX, title, cPlt.x, result[0], binS = 1.0, binR = (0.0, 0.0))
  else:
    result = pd.processData(h5f, fileInfo, pd, config)

proc createPlot*(h5f: H5File, fileInfo: FileInfo,
                 pd: PlotDescriptor, config: Config): (string, PlotV)

proc handleSubPlots(h5f: H5File,
                    fileInfo: FileInfo,
                    pd: PlotDescriptor,
                    config: Config): (string, PlotV) =
#    result.add PlotDescriptor(runType: runType,
#                              name: "Event display InGrid + FADC",
#                              runs: @[run],
#                              chip: fileInfo.centerChip,
#                              isCenterChip: true,
#                              plotKind: pkSubPlots,
#                              plots: @[ingrid, fadc],
#                              domain: @[ingridDomain, fadcDomain])
  var fnames: seq[string]
  var plts: seq[PlotV]
  for p in pd.plots:
    let (outfile, pltV) = createPlot(h5f, fileInfo, p, config)
    fnames.add outfile
    plts.add pltV

  # now combine plots
  let plt = makeSubplot(pd, plts)
  plt.plPlotJson.show()

proc createPlot*(h5f: H5File,
                 fileInfo: FileInfo,
                 pd: PlotDescriptor,
                 config: Config
                ): (string, PlotV) =
  ## creates a plot of kind `plotKind` for the data from all runs in `runs`
  ## for chip `chip`
  # TODO: think: have createPlot return the `PlotV` object. Then we could
  # call this proc recursively and add other plots to the existing figure
  let test = % pd
  info "Generating plot for: ", test.pretty
  # test reverse
  #let test2 = parsePd(test)
  #doAssert test2 == pd
  try:
    case pd.plotKind
    of pkInGridDset:
      result = handleInGridDset(h5f, fileInfo, pd, config)
    of pkFadcDset:
      result = handleFadcDset(h5f, fileInfo, pd, config)
    of pkOccupancy:
      result = handleOccupancy(h5f, fileInfo, pd, config)
    of pkOccCluster:
      result = handleOccCluster(h5f, fileInfo, pd, config)
    of pkPolya, pkFeSpec, pkFeSpecCharge:
      result = handleBarScatter(h5f, fileInfo, pd, config)
    of pkCombPolya:
      result = handleCombPolya(h5f, fileInfo, pd, config)
    of pkFeVsTime, pkFeChVsTime:
      result = handleFeVsTime(h5f, fileInfo, pd, config)
    of pkFePixDivChVsTime:
      result = handleFePixDivChVsTime(h5f, fileInfo, pd, config)
    of pkFePhotoDivEscape:
      result = handleFePhotoDivEscape(h5f, fileInfo, pd, config)
    of pkInGridEvent:
      # TODO: make this into a call to an `InGridEventIterator`
      result = handleIngridEvent(h5f, fileInfo, pd, config)
      #for outfile, pltV in ingridEventIter(h5f, fileInfo, pd, config):
      #  yield (outfile, pltV)
    of pkFadcEvent:
      result = handleFadcEvent(h5f, fileInfo, pd, config)
      #for outfile, pltV in fadcEventIter(h5f, fileInfo, pd, config):
      #  yield (outfile, pltV)
    of pkSubPlots:
      result = handleSubPlots(h5f, fileInfo, pd, config)
      #for outfile, pltV in subPlotsIter(h5f, fileInfo, pd, config):
      #  yield (outfile, pltV)
    #of pkOuterChips:
    #  result = handleOuterChips(h5f, fileInfo, pd, config)
    of pkToTPerPixel:
      result = handleToTPerPixel(h5f, fileInfo, pd, config)
    of pkCustomPlot:
      result = handleCustomPlot(h5f, fileInfo, pd, config)
    else:
      discard

    # finally call savePlot if we actually created a plot
    #if not result[1].plPlot.isNil:
    info &"Calling savePlot for {pd.plotKind} with filename {result[0]}"
    savePlot(result[1], result[0], config, fullPath = true)
  except KeyError as e:
    echo "WARNING: Could not generate the plot: " & $pd & ". Skipping it."
    echo "Exception message: ", e.msg

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

proc handleOutput(basename: string, config: Config): string =
  ## handles output file creation. Reads the config.toml file
  ## and uses it to decide whether to store data as Org file or
  ## dump to Json
  # TODO: create some kind of Config object, which does this once at
  # startup of the program
  # TODO: when Org is used, we need to save all plots as SVG!
  result = basename
  for fl in config.flags:
    result &= "_" & $fl
  case config.outputType
  of ofOrg:
    result &= ".org"
    createOrg(result, fileType)
  of ofJson:
    result &= ".json"
    discard jsonDump(result)
  else:
    warn "Unsupported output file format!"

proc plotsFromPds(h5f: H5File,
                  fileInfo: FileInfo,
                  pds: seq[PlotDescriptor],
                  config: Config) =
  ## calls `createPlot` for each PD and saves filename to
  ## `imageSet` if necessary
  for p in pds:
    let (fileF, pltV) = h5f.createPlot(fileInfo, p, config)
    case BKind
    of bPlotly:
      if config.plotlySaveSvg:
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

proc serveRequests(h5f: H5File,
                   fileInfo: FileInfo,
                   pds: seq[PlotDescriptor],
                   config: Config) =
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
      discard createPlot(h5f, fileInfo, pd, config)
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

proc serve(h5f: H5File,
           fileInfo: FileInfo,
           pds: seq[PlotDescriptor],
           config: Config) =
  ## serves the client
  # before we do anything, send all PDs to the client

  while true:
    while not serveNewClient():
      sleep(100)

    if serveNewClient():
      continue

    for p in pds:
      #for sp in createPlotIter(h5f, fileInfo, p):
      discard createPlot(h5f,fileInfo, p, config)
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
                  runOpt: Option[int],
                  runType: RunTypeKind,
                  bKind: PlottingBackendKind,
                  config: Config): string =
  ## use as event display tool
  var h5f = H5open(h5file, "r")
  let fileInfo = getFileInfo(h5f)
  let fInfoConfig = fileInfo.appliedConfig(config)
  let events = toOrderedSet(toSeq(0 .. 50))#initOrderedSet[int]())
  ## XXX: events irrelevant
  let run = runOpt.get
  var pds: seq[PlotDescriptor]
  for r in fileInfo.runs:
    if run > 0 and r == run:
      pds.add createEventDisplayPlots(h5f, r, runType, fInfoConfig, config, events)
    elif run < 0:
      pds.add createEventDisplayPlots(h5f, r, runType, fInfoConfig, config, events)

  if cfProvideServer in config.flags:
    serve(h5f, fileInfo, pds, config)
  else:
    plotsFromPds(h5f, fInfoConfig, pds, config)
    let outfile = "eventDisplay"  & "_" & getRunsStr(fInfoConfig.runs)
    result = handleOutput(outfile, config)

  discard h5f.close()

include ./moreCustomPlots

proc genCalibrationPlotPDs(h5f: H5File,
                           runType: RunTypeKind,
                           config: Config
                          ): seq[PlotDescriptor] =
  ## Creates the PlotDescriptors for all calibration plots in the `h5file`
  var fileInfo = getFileInfo(h5f)
  let fInfoConfig = fileInfo.appliedConfig(config)
  # var imageSet = initOrderedSet[string]()

  if cfOccupancy in config.flags:
    result.add occupancies(h5f, runType, fInfoConfig, config) # plus center only
  if cfPolya in config.flags:
    result.add polya(h5f, runType, fInfoConfig, config)
  if cfToTPerPixel in config.flags:
    result.add totPerPixel(h5f, runType, fInfoConfig, config)
  if cfFeSpectrum in config.flags:
    result.add feSpectrum(h5f, runType, fInfoConfig, config)
  #result.add fePhotoDivEscape(h5f, runType, fInfoConfig, config)
  # energyCalib(h5f) # ???? plot of gas gain vs charge?!
  if cfIngrid in config.flags:
    result.add histograms(h5f, runType, fInfoConfig, config) # including fadc
  if config.customPlots.len > 0:
    result.add createCustomPlots(fInfoConfig, config)

  # now deal with custom plots
  if cfCompiledCustom in config.flags:
    const path = currentSourcePath().parentDir
    const FileAtCompileTime = readFile(path / "moreCustomPlots.nim")
    if fileExists(path / "moreCustomPlots.nim"):
      let fileAtRunTime = readFile(path / "moreCustomPlots.nim")
      if FileAtCompileTime != fileAtRunTime:
        raise newException(Exception, "The existing file `moreCustomPlots.nim` changed compared to the " &
          "version we compiled against! You gave the `--compiledCustom` option, which implies you wish to " &
          "generate these plots. Please recompile this program (`plotData.nim`).")
    else:
      info "Cannot find `moreCustomPlots.nim` source file. Using compiled version!"
    result.add moreCustom(fInfoConfig, config)

proc createCalibrationPlots(h5file: string,
                            bKind: PlottingBackendKind,
                            runType: RunTypeKind,
                            config: Config): string =
  ## creates QA plots for calibration runs
  withH5(h5file, "r"):
    var fileInfo = getFileInfo(h5f)
    let fInfoConfig = fileInfo.appliedConfig(config)

    let pds = genCalibrationPlotPDs(h5f, runType, config)
    if cfProvideServer in config.flags:
      serveRequests(h5f, fileInfo, pds, config)
      #serve(h5f, fileInfo, pds, flags)
    else:
      plotsFromPds(h5f, fInfoConfig, pds, config)
      echo "Image set is ", imageSet.card

      # likelihoodHistograms(h5f) # need to cut on photo peak and esc peak
      # neighborPixels(h5f)
      discard h5f.close()
      let outfile = "calibration" & "_" & getRunsStr(fInfoConfig.runs)
      result = handleOutput(outfile, config)

proc genBackgroundPlotPDs(h5f: H5File, runType: RunTypeKind,
                          config: Config): seq[PlotDescriptor] =
  var fileInfo = getFileInfo(h5f)
  let fInfoConfig = fileInfo.appliedConfig(config)
  if cfOccupancy in config.flags:
    result.add occupancies(h5f, runType, fInfoConfig, config) # plus center only
  if cfPolya in config.flags:
    result.add polya(h5f, runType, fInfoConfig, config)
  if cfToTPerPixel in config.flags:
    result.add totPerPixel(h5f, runType, fInfoConfig, config)
  # energyCalib(h5f) # ???? plot of gas gain vs charge?!
  if cfIngrid in config.flags:
    result.add histograms(h5f, runType, fInfoConfig, config) # including fadc
  # result.add createCustomPlots(fInfoConfig, config)
  if config.customPlots.len > 0:
    result.add createCustomPlots(fInfoConfig, config)
    # now deal with custom compiled plots
    result.add moreCustom(fInfoConfig, config)

proc createBackgroundPlots(h5file: string,
                           bKind: PlottingBackendKind,
                           runType: RunTypeKind,
                           config: Config): string =
  ## creates QA plots for calibration runs
  withH5(h5file, "r"):
    var fileInfo = getFileInfo(h5f)
    let fInfoConfig = fileInfo.appliedConfig(config)

    let pds = genBackgroundPlotPDs(h5f, runType, config)
    if cfProvideServer in config.flags:
      serveRequests(h5f, fileInfo, pds, config)
      #serve(h5f, fileInfo, pds, config.flags)
    else:
      plotsFromPds(h5f, fInfoConfig, pds, config)
      echo "Image set is ", imageSet.card

      # likelihoodHistograms(h5f) # need to cut on photo peak and esc peak
      # neighborPixels(h5f)
      discard h5f.close()
      let outfile = "background" & "_" & getRunsStr(fInfoConfig.runs)
      result = handleOutput(outfile, config)

proc createXrayFingerPlots(bKind: PlottingBackendKind,
                           config: Config): string =
  discard

proc handlePlotTypes(h5file: string,
                     bKind: PlottingBackendKind,
                     runType: RunTypeKind,
                     config: Config,
                     evDisplayRun = none[int]()) =
  ## handles dispatch of the correct data type / kind / mode to be plotted
  ## calibration, X-ray, background, event display...
  ## If not run in server mode, launch staticDisplay with output JSON file
  ## TODO: what happens for SVG? Need to check too.
  var outfile = ""
  if evDisplayRun.isSome():
    outfile = eventDisplay(h5file, evDisplayRun, runType, BKind, config)
  else:
    case runType
    of rtCalibration:
      outfile = createCalibrationPlots(h5file, bKind, runType, config)
    of rtBackground:
      outfile = createBackgroundPlots(h5file, bKind, runType, config)
    of rtXrayFinger:
      outfile = createXrayFingerPlots(bKind, config)
    else:
      discard
  if cfProvideServer in config.flags:
    # if not using server, start client, in case the data is being stored as
    # JSON
    let filecall = &"--file:{outfile}"
    shell:
      ./karaRun "-r" ($filecall) staticClient.nim

proc find(pds: seq[PlotDescriptor], toFind: PlotDescriptor): Option[PlotDescriptor] =
  ## Tries to find an equivalent PlotDescriptor as `toFind` in `pds`.
  ##
  ## This implies the same type of plot, but different input file / runs etc
  for pdc in pds:
    if pdc == toFind:
      return some(pdc)

proc add(p: var PlotV, p2: PlotV, f1, f2: string) =
  ## adds `p2` to `p`, which implies combining the two plots onto one
  ##
  ## `f1` and `f2` are the names of the `From` field used for each DF
  doAssert p.kind == p2.kind, $p.kind & " vs " & $p2.kind
  case p.kind
  of bGgPlot:
    ## XXX: handle geom_raster!
    let hasRaster = p.pltGg.geoms.anyIt(it.kind == gkRaster)
    # get the df of `p` re-emit with a different df based on a combination of two
    var dfP = p.pltGg.data
    var dfP2 = p2.pltGg.data
    if "From" notin dfP:
      dfP["From"] = f1
    dfP2["From"] = f2
    dfP.add dfP2
    p.pltGg.data = dfP
    # now update `aes` to use new `From` column
    if not hasRaster:
      p.pltGg.aes.color = some(Scale(scKind: scColor, col: f{"From"}, hasDiscreteness: true,
                                     ids: p.pltGg.aes.x.get.ids))
      p.pltGg.aes.fill = some(Scale(scKind: scFillColor, col: f{"From"}, hasDiscreteness: true,
                                     ids: p.pltGg.aes.x.get.ids))
    else:
      p.pltGg.facet = some(Facet(columns: @["From"]))
    ## add the existing geoms if they have non trivial data (otherwise they will be the same)
    for g in p2.pltGg.geoms:
      if g.data.isSome: #
        p.pltGg.geoms.add g

    for g in mitems(p.pltGg.geoms):
      # update geom to include alpha and identity position
      g.userStyle.alpha = some(0.5)
      g.position = pkIdentity
  else:
    raise newException(Exception, "Not implemented to add plots of kind " & $p.kind)

proc createComparePlots(h5file: string, h5Compare: seq[string],
                        bKind: PlottingBackendKind,
                        runType: RunTypeKind,
                        compareRunTypes: seq[RunTypeKind],
                        config: Config
                       ) =
  ## Generates plots comparing all plots that can be built from `h5file`
  ## and add equivalent plot from data in each `h5Compare`
  var h5f = H5open(h5file, "r")
  var fileInfo = getFileInfo(h5f).appliedConfig(config)
  if compareRunTypes.len > 0 and compareRunTypes.len != h5Compare.len:
    raise newException(Exception, "If `compareRunTypes` given, must be the same length " &
      "as `h5Compare` files!")
  var h5fs = h5Compare.mapIt(H5open(it, "r"))
  var fInfos = h5fs.mapIt(getFileInfo(it).appliedConfig(config))

  template genPds(file, runType, flags: untyped): untyped =
    var res: seq[PlotDescriptor]
    case runType
    of rtCalibration: res = genCalibrationPlotPDs(file, runType, config)
    of rtBackground: res = genBackgroundPlotPDs(file, runType, config)
    else: discard
    res
  let pds = genPds(h5f, runType, config.flags)
  var pdComp = newSeq[seq[PlotDescriptor]](h5Compare.len)
  for i in 0 ..< h5Compare.len:
    let rt = if compareRunTypes.len > 0: compareRunTypes[i] else: runType
    pdComp[i] = genPds(h5fs[i], rt, config.flags)

  for pd in pds:
    var (outfile, plt) = createPlot(h5f, fileInfo, pd, config)
    if outfile.len == 0: continue # apparantly failed to create plot, skipping
    var validCompare = true
    for i in 0 ..< h5Compare.len:
      let pdc = pdComp[i].find(pd)
      if pdc.isSome:
        let (_, cPlt) = createPlot(h5fs[i], fInfos[i], pdc.get, config)
        # add to `plt`
        if cPlt.kind != bNone:
          plt.add(cPlt, h5f.name.extractFilename, h5fs[i].name.extractFilename & "_" & $i)
        else:
          validCompare = false
    # now generate the plot
    if validCompare:
      let outf = outfile & "_combined.pdf"
      plt.savePlot(outf, config, fullPath = true)

proc parseBackendType(backend: string): PlottingBackendKind =
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

proc sendDataPacket(ws: WebSocket, data: JsonNode, kind: PacketKind) =
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
        asyncCheck ws.send(packet.asData)
      elif i < nParts - 1:
        packet = initDataPacket(kind,
                                Messages.Data,
                                payload = dPart,
                                recvData = true)
        asyncCheck ws.send(packet.asData)
      else:
        doAssert i == nParts - 1
        packet = initDataPacket(kind,
                                Messages.DataStop,
                                payload = dPart,
                                recvData = false,
                                done = true)
        waitFor ws.send(packet.asData)
  else:
    # handle the single packet case differently
    let packet = initDataPacket(kind,
                                Messages.DataSingle,
                                payload = sendData,
                                recvData = false,
                                done = true)
    waitFor ws.send(packet.asData)


proc processClient(req: Request) {.async.} =
  ## handle a single client
  let ws = await newWebSocket(req)
  if ws.isNil:
    echo "WS negotiation failed: ", ws.readyState
    await req.respond(Http400, "Websocket negotiation failed: " & $ws.readyState)
    req.client.close()
    return
  else:
    # receive connection successful package
    let (opcodeConnect, dataConnect) = await ws.receivePacket()
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
    let (opcode, data) = await ws.receivePacket()
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
      ws.close()
      echo "Socket went away, close code: ", $ws.readyState
      req.client.close()
      return
    else:
      echo "Unkown error: ", opcode

  echo "This client dies now!"
  ws.close()
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

proc plotData*(
  h5file: string,
  runType: RunTypeKind,
  backend: PlottingBackendKind = bGgPlot,
  eventDisplay: int = int.high,
  h5Compare: seq[string] = @[], # additional files to compare against
  compareRunTypes: seq[RunTypeKind] = @[],
  server = false, fadc = false, ingrid = false, occupancy = false,
  polya = false, totPerPixel = false, fe_spec = false,
  compiledCustom = false, config = "",
  version = false,
  chips: set[uint16] = {},
  runs: set[uint16] = {},
  show = false,
  cuts: seq[GenericCut] = @[],
  maskRegion: seq[MaskRegion] = @[],
  applyAllCuts = false,
  cutFePeak = false,
  head: int = 0,
  x: string = "",
  y: string = "",
  z: string = ""
              ) =
  ## the main workhorse of the server end
  if version:
    echo "Version: ", version
    return

  if config.len > 0:
    ConfigFile = config

  fileDir = genPlotDirName(h5file, "figs")
  let tomlConfig = parseToml.parseFile(ConfigFile)
  discard existsOrCreateDir(fileDir)
  var flags: set[ConfigFlagKind]
  if fadc:
    flags.incl cfFadc
  if ingrid:
    flags.incl cfIngrid
  if occupancy:
    flags.incl cfOccupancy
  if polya:
    flags.incl cfPolya
  if fe_spec:
    flags.incl cfFeSpectrum
  if totPerPixel:
    flags.incl cfTotPerPixel
  if server:
    flags.incl cfProvideServer
  if show:
    flags.incl cfShow
  if applyAllCuts:
    flags.incl cfApplyAllCuts
  if cutFePeak:
    flags.incl cfCutFePeak
  if compiledCustom:
    flags.incl cfCompiledCustom

  let cfg = initConfig(chips, runs, flags,
                       cuts = cuts,
                       maskRegion = maskRegion,
                       tomlConfig = tomlConfig,
                       head = head,
                       xDset = x,
                       yDset = y,
                       zDset = z
  )

  info &"Flags are:\n  {flags}"

  if cfShow in flags and cfProvideServer in flags:
    echo "Please either use `--show` or `--server`, but not both!"
    return

  var thr: Thread[void]
  if cfProvideServer in flags:
    # set up the server and launch the client
    # create and run the websocket server
    thr.createThread(serve)

  var evDisplayRun: Option[int] = none[int]()
  BKind = backend
  if eventDisplay != int.high:
    evDisplayRun = some(eventDisplay)

  if h5Compare.len == 0:
    handlePlotTypes(h5file, backend, runType, cfg, evDisplayRun = evDisplayRun)
  else:
    createComparePlots(h5file, h5Compare, backend, runType, compareRunTypes, cfg)

  if cfProvideServer in flags:
    joinThread(thr)

  echo "All plots were saved to:\n", fileDir

when isMainModule:
  import cligen/argcvt
  proc argParse(dst: var GenericCut, dfl: GenericCut,
                a: var ArgcvtParams): bool =
    var vals = a.val.strip(chars = {'(', ')'}).split(',')
    if vals.len != 3: return false
    try:
      let dset = vals[0].strip(chars = {'"'})
      dst = GenericCut(dset: dset,
                       applyDset: @[dset], ## XXX: set this?
                       min: parseFloat(vals[1].strip),
                       max: parseFloat(vals[2].strip))
      result = true
    except:
      result = false

  proc argHelp*(dfl: GenericCut; a: var ArgcvtParams): seq[string] =
    result = @[ a.argKeys, "(dset: string, lower, upper: float)", $dfl ]

  proc argParse(dst: var MaskRegion, dfl: MaskRegion,
                a: var ArgcvtParams): bool =
    proc stripAll(s: string): string = s.strip(chars = {'(', ')', ' '})
    var xy = a.val.stripAll.split(',')
    if xy.len != 4: return false
    try:
      let xR = (min: ChipCoord(xy[0].stripAll.parseInt), max: ChipCoord(xy[1].stripAll.parseInt))
      let yR = (min: ChipCoord(xy[2].stripAll.parseInt), max: ChipCoord(xy[3].stripAll.parseInt))
      dst = MaskRegion(x: xR, y: yR)
      result = true
    except:
      result = false

  proc argHelp*(dfl: MaskRegion; a: var ArgcvtParams): seq[string] =
    result = @[ a.argKeys, "(x: (min, max: int), y: (min, max: int))", $dfl ]


  dispatch(plotData, help = {
    "h5file" : "The h5 input file to use",
    "runType" : """Select run type (Calib | Back | Xray)
  The following are parsed case insensetive:
  Calib = {"calib", "calibration", "c"}
  Back = {"back", "background", "b"}
  Xray = {"xray", "xrayfinger", "x"}""",

    "backend" : """Select the plotting backend to be chosen.
  The followoing backends are available:
  Python / Matplotlib: {"python", "py", "matplotlib", "mpl"}
  Nim / Plotly: {"nim", "plotly"}""",

    "eventDisplay" : """If given will show event displays of the given run. If a negative
run number is given, will generate events for all runs in the file.""",
    "h5Compare" : "If any given, all plots will compare with same data from these files.",
    "server" : """If flag given, will launch client and send plots individually,
  instead of creating all plots and dumping them.""",
    "compareRunTypes" : """If any given will use these as the run types of the compare files.
  Interpreted in the same order as `runType`. If missing assume same type as `runType`.""",

    "chips" : """If any given overwrites the `config.toml` field of `allowedChips`.
  That means only generate plots for the given chips.""",
    "runs" : """If any given overwrites the `config.toml` field of `allowedRuns`.
  That means only generate plots for the given runs.""",

    "fadc" : "If set FADC plots will be created.",
    "ingrid" : "If set InGrid plots will be created.",
    "occupancy" : "If set occupancy plots will be created.",
    "polya" : "If set polya plots will be created.",
    "totPerPixel" : "If set totPerPixel plots will be created.",
    "fe_spec" : "If set Fe spectrum will be created.",
    "cutFePeak" : "If set will create plots cut to Photo & Escape peak for calibration data",
    "compiledCustom" : "If set will create the custom plots that require recompilation (moreCustomPlots.nim)",

    "x" : "Generate plot of dataset `x` against `y`.",
    "y" : "Generate plot of dataset `x` against `y`.",
    "z" : "Generate plot of dataset `x` against `y` colored by `z`",

    "show" : "If given will open each generated plot using inkview (svg) / evince (pdf) / nomacs (png).",
    "cuts" : "Allows to cut data used for event display. Only shows events of data passing these cuts.",
    "maskRegion" : "Allows to mask a region (in pixel coords). Any pixel in the region will be cut.",
    "applyAllCuts" : "If given will apply all given cuts to all data reads.",
    "head" : "Only process the first this many elements (mainly useful for event display).",

    "config" : "Path to the TOML config file.",
    "version" : "Show the version number"})
