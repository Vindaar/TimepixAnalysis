import plotly
import os, strutils, strformat
import sequtils
import algorithm, sets
import typeinfo
import math
import nimpy
import macros
import docopt
import chroma
import seqmath
import options
import mpfit
import arraymancer
import zero_functional
import nimhdf5
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
    cfNone, cfNoFadc, cfNoInGrid, cfNoOccupancy, cfNoPolya

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

  DataKind = enum
    dkInGrid, dkFadc

  ShapeKind = enum
    Rectangle, Square

# global variable which stores the backend the user selected
var BKind: BackendKind = bNone

proc getTrace[T](thl, hits: seq[T], voltage: string): Trace[float] =
  result = Trace[float](`type`: PlotType.Scatter)
  # filter out clock cycles larger 300 and assign to `Trace`
  result.xs = thl.asType(float)
  result.ys = hits.asType(float)
  result.name = &"Voltage {voltage}"

proc plotHist*[T](traces: seq[Trace[T]], voltages: set[int16], chip = "") =
  ## given a seq of traces (SCurves and their fits), plot
  let
    layout = Layout(title: &"SCurve of Chip {chip} for voltage {voltages}",
                    width: FigWidth.int, height: FigHeight.int,
                    xaxis: Axis(title: "Threshold value"),
                    yaxis: Axis(title: "# hits$"),
                    autosize: false)
    p = Plot[float](layout: layout, traces: traces)
  let filename = &"out/scurves_{chip}.svg"
  p.saveImage(filename)

proc getLayout(title: string, xlabel: string, ylabel: string, shape = ShapeKind.Rectangle): Layout =
  ## returns the plotly layout of the given shape
  var width: int
  case shape
  of ShapeKind.Rectangle:
    width = FigWidth.int
  of ShapeKind.Square:
    # for square take height as width
    width = FigHeight.int
  result = Layout(title: title,
                  width: width, height: FigHeight.int,
                  xaxis: Axis(title: xlabel),
                  yaxis: Axis(title: ylabel),
                  autosize: false)

proc plotHist[T](xIn: seq[T], title, dset: string) =
  ## plots the data in `x` as a histogram
  let x = xIn.mapIt(it.float)
  let binRangeO = getBinRangeForDset(dset)
  let nbinsO = getNumBinsForDset(dset)
  var binRange: tuple[start, stop: float]
  var nbins: int
  if binRangeO.isSome:
    binRange = get(binRangeO)
  else:
    binRange = (0.0, x.max.float)
  if nBinsO.isSome:
    nBins = get(nbinsO)
  else:
    nBins = 50
  case BKind
  of bPlotly:
    let lyout = getLayout(title, dset, "#")
    let tr = Trace[float](`type`: PlotType.Histogram,
                          bins: binRange,
                          nbins: nbins,
                          xs: x)
    let p = Plot[float](layout: lyout, traces: @[tr])
    p.show()
  of bMpl:
    let mpl = pyImport("matplotlib")
    discard mpl.use("TKagg")
    let plt = pyImport("matplotlib.pyplot")
    let np = pyImport("numpy")
    let (fig, ax) = plt.subplots(1, 1).to((PyObject, PyObject))
    discard ax.hist(x,
                    bins = nbins,
                    range = binRange)
    # we cannot call `show` directly, because the Nim compiler tries to
    # call the Nim `show` procedure and fails
    discard callMethod(plt, "show")

  else: discard


#plotData(data, nbins, range,
#         "plots/{0}_run_{1}_chip_{2}".format(dset, run_number, chip),
#         "{0} of run {1} for chip {2}".format(dset, run_number, chip),
#         dset,
#         "# hits",
#         fitting_only = True)

iterator chips(group: var H5Group): (int, H5Group) =
  ## returns all chip groups within the given Run group
  doAssert group.name.parentDir == "/reconstruction", "Bad group : " & group.name
  doAssert "run_" in group.name, "Bad group : " & group.name
  for grp in groups(group):
    if "chip_" in grp.name:
      let chipNum = grp.attrs["chipNumber", int]
      yield (chipNum, grp)

proc plotHistsFeSpec(h5f: var H5FileObj,
                     group: H5Group,
                     dsetName: string,
                     dKind: DataKind,
                     chipNum = 0) =
  ## creates 3 plots for the given dataset.
  ## - 1 around photo peak
  ## - 1 around escape peak
  ## - 1 without cut
  const ranges = [(5.5, 6.2, "Photopeak"), (2.7, 3.2, "Escapepeak"), (-Inf, Inf, "All")]
  for r in ranges:
    case dKind
    of dkInGrid:
      let idx = cutOnProperties(h5f, group,
                      ("energyFromCharge", r[0], r[1]))
      # extract dataset and correct indices
      var ldata = h5f[group.name / dsetName, float]
      # extract indices corresponding to peak
      ldata = idx.mapIt(ldata[it])
      let title = &"Dataset: {dsetName} for chip {chipNum} in range: {r[2]}"
      plotHist(ldata, title, dsetName)
    of dkFadc:
      # get the center chip group
      let centerChip = h5f[recoGroupGrpStr()].attrs["centerChip", int]
      let centerGroup = h5f[(group.name.parentDir / &"chip_{centerChip}").grp_str]
      let idx = cutOnProperties(h5f, centerGroup,
                      ("energyFromCharge", r[0], r[1]))
      # given these idx, get event numbers
      let evNumInGrid = h5f[centerGroup.name / "eventNumber", int]
      # filter out correct indices passing cuts
      var inGridSet = initSet[int]()
      for i in idx:
        inGridSet.incl evNumInGrid[i]
      let evNumFadc = h5f[group.name / "eventNumber", int]
      #let evNumFilter = evNumFadc --> filter(it in inGridSet)
      let idxFadc = (toSeq(0 .. evNumFadc.high)) --> filter(evNumFadc[it] in inGridSet)
      let dset = h5f[(group.name / dsetName).dset_str]
      case dset.dtypeAnyKind
      of akUint16: # fallTime
        var ldata = dset[uint16]
        ldata = idxFadc --> map(ldata[it])
        let title = &"Dataset: {dsetName} for fadc in range: {r[2]}"
        plotHist(ldata, title, dsetName)
      of akFloat64: # minVals
        var ldata = dset[float]
        ldata = idxFadc --> map(ldata[it])
        let title = &"Dataset: {dsetName} for fadc in range: {r[2]}"
        plotHist(ldata, title, dsetName)
      else:
        echo "Unsupported type ", dset.dtypeAnyKind
        discard
    else:
      echo "No data kind selected: ", dKind

proc histograms(h5f: var H5FileObj, flags: set[ConfigFlagKind]) =
  const dsets = ["length", "width", "skewnessLongitudinal", "skewnessTransverse",
                 "kurtosisLongitudinal", "kurtosisTransverse", "rotationAngle",
                 "eccentricity", "fractionInTransverseRms", "lengthDivRmsTrans"]
  const fadcDsets = ["minvals", "fallTime"]
  # TODO: perform cut on photo peak and escape peak, done by getting the energy
  # and performing a cut around peak +- 300 eV maybe
  # NOTE: need photo peak and escape peak plenty of times here. Just write wrapper
  # which takes energy and performs cuts, *then* creates plots


  for runNumber, group in runs(h5f):
    var group = h5f[group.grp_str]
    if cfNoIngrid notin flags:
      for chipNum, chpgrp in chips(group):
        for dset in dsets:
          plotHistsFeSpec(h5f, chpgrp, dset, dkInGrid, chipNum)

    # fadc plots
    if cfNoFadc notin flags:
      var fadcGrp = h5f[(group.name / "fadc").grp_str]
      for dsetName in fadcDsets:
        plotHistsFeSpec(h5f, fadcGrp, dsetName, dkFadc)

proc calcOccupancy(x, y: seq[seq[uint8]]): Tensor[int] =
  ## calculates the occupancy of the given x and y datasets
  result = newTensor[int]([256, 256])
  # iterate over events
  for i in 0 .. x.high:
    let
      xEv = x[i]
      yEv = y[i]
    for j in 0 .. xEv.high:
      result[xEv[j].int, yEv[j].int] += 1

proc plotHist2D[T](data: Tensor[T], title, descr: string) =
  ## creates a 2D histogram plot (basically an image) of the given 2D
  ## Tensor
  doAssert data.rank == 2
  case BKind
  of bPlotly:
    let lyout = getLayout(title, "pixel x", "pixel y", ShapeKind.Square)
    let tr = Trace[T](`type`: PlotType.HeatMap,
                      colormap: ColorMap.Viridis,
                      zs: data.toRawSeq.reshape2D([256, 256]))
    let p = Plot[T](layout: lyout, traces: @[tr])
    p.show()
  of bMpl:
    discard
  else:
    discard

proc plotOccupancies(h5f: var H5FileObj,
                     group: H5Group,
                     chipNum: int,
                     runNumber: string) =
  ## creates 3 occupancy plots for the given run & chip for different ranges
  ## - 1 full range
  ## - 1 clamped to .75 percentile
  ## - 1 clamped to 10 hits
  const
    quantile = 95
    clamp3 = 10
  # get x and y datasets, stack and get occupancies
  let vlenDtype = special_type(uint8)
  let
    xGroup = h5f[(group.name / "x").dset_str]
    yGroup = h5f[(group.name / "y").dset_str]
    x = xGroup[vlenDtype, uint8]
    y = yGroup[vlenDtype, uint8]
  let occ = calcOccupancy(x, y)
  # calculate 3 clamp values based on occ
  let quant = occ.toRawSeq.percentile(quantile)
  let clamps = [int.high, quant.int, clamp3]
  for cl in clamps:
    var title = &"Occupancy of chip {chipNum} for run {runNumber}"
    if cl == int.high:
      title = title & " in full range"
      plotHist2D(occ, title, "Occupancy")
    else:
      let occClamped = occ.clamp(0, cl)
      title = title & &" clamped to: {cl}"
      plotHist2D(occClamped, title, "Occupancy")

proc occupancies(h5f: var H5FileObj, flags: set[ConfigFlagKind]) =
  ## creates occupancy plots for the given HDF5 file, iterating over
  ## all runs and chips
  for runNumber, grpName in runs(h5f):
    var group = h5f[grpName.grp_str]
    # now build occupancy
    for chipNum, chpgrp in chips(group):
      plotOccupancies(h5f, chpgrp, chipNum, runNumber)

proc polya(h5f: var H5FileObj, flags: set[ConfigFlagKind]) =
  ## creates the plots for the polya distribution of the datasets
  # for now try to read polya datasets?!
  discard

proc createCalibrationPlots(h5file: string, bKind: BackendKind,
                            flags: set[ConfigFlagKind]) =
  ## creates QA plots for calibration runs
  var h5f = H5file(h5file, "r")
  const length = "length"
  if cfNoOccupancy notin flags:
    occupancies(h5f, flags) # plus center only
  if cfNoPolya notin flags:
    polya(h5f, flags)
  # feSpectrum(h5f)
  # centerFePerTime(h5f)
  # energyCalib(h5f) # ????
  histograms(h5f, flags) # including fadc
  # likelihoodHistograms(h5f) # need to cut on photo peak and esc peak
  # neighborPixels(h5f)

  discard h5f.close()


proc createBackgroundPlots(bKind: BackendKind, flags: set[ConfigFlagKind]) =
  discard

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

proc main() =
  let args = docopt(doc)
  echo args
  let h5file = $args["<H5file>"]
  let runTypeStr = $args["--runType"]
  let backendStr = $args["--backend"]
  var flags: set[ConfigFlagKind]
  if $args["--no_fadc"] == "true":
    flags.incl cfNoFadc
  if $args["--no_ingrid"] == "true":
    flags.incl cfNoIngrid

  var runType: RunTypeKind
  var bKind: BackendKind
  if runTypeStr != "nil":
    runType = parseRunType(runTypeStr)
  if backendStr != "nil":
    BKind = parseBackendType(backendStr)

  case runType
  of rtCalibration:
    createCalibrationPlots(h5file, bKind, flags)
  of rtBackground:
    createBackgroundPlots(bKind, flags)
  of rtXrayFinger:
    createXrayFingerPlots(bKind, flags)
  else:
    discard


when isMainModule:
  main()
