import plotly
import os, strutils, strformat
import sequtils
import algorithm
import math
import docopt
import chroma
import seqmath
import mpfit
import zero_functional
import ingrid/ingrid_types
import ingrid/tos_helpers
import ingrid/calibration
import helpers/utils
import ingridDatabase

type
  # for clarity define a type for the Docopt argument table
  DocoptTab = Table[string, Value]

template canImport(x: untyped): bool =
  compiles:
    import x

when canImport(ingridDatabase):
  # only in this case give option to use plot from database  
  import ingridDatabase

const doc = """
A simple tool to plot SCurves or ToT calibrations.

Usage:
  plotCalibration (--scurve | --tot) (--db=chipNumber | --file=FILE | --folder=FOLDER) [options]

Options:
  --scurve         If set, perform SCurve analysis
  --tot            If set, perform ToT calibration analysis
  --db=chipNumber  If given will read information from InGrid database, if 
                   available
  --file=FILE      If given will read from a single file
  --folder=FOLDER  If given will read all voltage files from the given folder
  --chip=NUMBER    The number of this chip
  --startFit=FIT   Start the TOT fit from this measurement. If differs from 
                   StartToT constant in source code (or the --startTot value),
                   the other datapoints will still be plotted.
  --startTot=TOT   Read the ToT file from this pulse height
  -h, --help       Show this help
  --version        Show the version number
"""

const
  StopToT = 450.0

const
  GoldenMean = (sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
  FigWidth = 1200.0                     # width in inches
  FigHeight = FigWidth * GoldenMean     # height in inches

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
  let filename = &"scurves_{chip}.svg"
  p.saveImage(filename)

proc plotThlCalib*(thlCalib: FitResult, charge, thl, thlErr: seq[float], chip = "") =
  let
    data = Trace[float](mode: PlotMode.Markers, `type`: PlotType.Scatter)
    fit = Trace[float](mode: PlotMode.Lines, `type`: PlotType.Scatter)
  # flip the plot, i.e. show THL on x instead of y as done for the fit to
  # include the errors
  data.ys = charge
  data.xs = thl
  data.xs_err = newErrorBar(thlErr, color = Color(r: 0.5, g: 0.5, b: 0.5, a: 1.0))
  data.name = "THL calibration"

  # flip the plot
  fit.ys = thlCalib.x
  fit.xs = thlCalib.y
  fit.name = "THL calibration fit"

  let
    layout = Layout(title: &"THL calibration of Chip {chip}",
                    width: FigWidth.int, height: FigHeight.int,
                    yaxis: Axis(title: "Charge / e-"),
                    xaxis: Axis(title: "THL"),
                    autosize: false)
    p = Plot[float](layout: layout, traces: @[data, fit])
  let filename = &"thl_calib_{chip}.svg"
  p.saveImage(filename)

proc plotToTCalib*(totCalib: FitResult, tot: Tot, chip = 0) =
  let
    data = Trace[float](mode: PlotMode.Markers, `type`: PlotType.Scatter)
    fit = Trace[float](mode: PlotMode.Lines, `type`: PlotType.Scatter)
  # flip the plot, i.e. show THL on x instead of y as done for the fit to
  # include the errors
  data.xs = tot.pulses.mapIt(it.float)
  data.ys = tot.mean
  let stdValid = tot.std.mapIt(if classify(it) == fcNaN: 0.0 else: it)
  echo stdValid
  data.ys_err = newErrorBar(stdValid, color = Color(r: 0.5, g: 0.5, b: 0.5, a: 1.0))
  data.name = "ToT calibration"

  # flip the plot
  fit.xs = totCalib.x
  fit.ys = totCalib.y
  fit.name = "ToT calibration fit"

  let
    layout = Layout(title: &"ToT calibration of Chip {chip}",
                    width: FigWidth.int, height: FigHeight.int,
                    xaxis: Axis(title: "U_inj / mV"),
                    yaxis: Axis(title: "ToT / Clock cycles"),
                    autosize: false)
    p = Plot[float](layout: layout, traces: @[data, fit])
  let filename = &"tot_calib_{chip}.svg"
  p.saveImage(filename)

iterator sCurves(args: DocoptTab): SCurve =
  ## yields the traces of the correct argument given to
  ## the program for SCurves
  let file = $args["--file"]
  let folder = $args["--folder"]
  let db = $args["--db"]
  if db != "nil":
    when declared(ingridDatabase):
      let chipName = getSeptemHChip(db.parseInt)
      let scurves = getScurveSeq(chipName)
      for curve in scurves.curves:
        yield curve
    else:
      discard
  elif file != "nil":
    let curve = readScurveVoltageFile(file)
    yield curve
  elif folder != "nil":
    for f in walkFiles(folder.expandTilde & "/*.txt"):
      let curve = readScurveVoltageFile(f)
      yield curve

proc sCurve(args: DocoptTab, chip: string) =
  ## perform plotting and fitting of SCurves
  var
    voltages: set[int16]
    traces: seq[Trace[float]] = @[]
    curve: SCurve

    # calibration seqs
    charge: seq[float] = @[]
    thlMean: seq[float] = @[]
    thlErr: seq[float] = @[]

  for curve in sCurves(args):
    let trace = getTrace(curve.thl, curve.hits, $curve.voltage)
    traces.add trace
    voltages.incl int16(curve.voltage)
    if curve.voltage > 0:
      # now fit the normal cdf to the SCurve and add its data
      let fitRes = fitSCurve(curve)
      let
        traceFit = getTrace(fitRes.x,
                            fitRes.y,
                            &"fit {curve.voltage} / chiSq {fitRes.redChiSq}")
      traces.add traceFit

      # given `fitRes`, add fit parameters to calibration seqs
      # multiply by 50 to take into account capacitance of test pulse injection
      # to get number of electrons
      charge.add (curve.voltage * 50).float
      thlMean.add (fitRes.pRes[1])
      thlErr.add (fitRes.pErr[1])

  if traces.len > 0:
    plotHist(traces, voltages, chip)

  # now fit the calibration function
  # however, the voltage data is in a wrong order. need to sort it, before
  # we can fit
  let thlTup = zip(thlMean, thlErr)
  var chThl = zip(charge, thlTup)
  let sortedChThl = chThl.sortedByIt(it[0])
  let
    chSort = sortedChThl.mapIt(it[0])
    thlSort = sortedChThl.mapIt(it[1][0])
    # increase errors artifically by factor 100... Mostly for visualization,
    # but also as a rough guess for typical deviation visible
    thlErrSort = sortedChThl.mapIt(it[1][1] * 100)

  if chSort.len > 0:
    let thlCalib = fitThlCalib(chSort, thlSort, thlErrSort)
    plotThlCalib(thlCalib, chSort, thlSort, thlErrSort, chip)

proc parseTotInput(args: DocoptTab, startTot = 0.0): (int, Tot) =
  let file = $args["--file"]
  let folder = $args["--folder"]
  let db = $args["--db"]

  var chip = 0
  var tot: Tot

  if db != "nil":
    when declared(ingridDatabase):
      chip = db.parseInt
      let chipName = getSeptemHChip(chip)
      tot = getTotCalib(chipName)
  elif file != "nil":
    (chip, tot) = readToTFile(file, startTot)
 
  else:
    for f in walkFiles(folder.expandTilde & "/*.txt"):
      # TODO: implement multiple in same plot?
      (chip, tot) = readToTFile(f, startTot)
    
  result = (chip, tot)

proc totCalib(args: DocoptTab, startFit = 0.0, startTot = 0.0) =
  ## perform plotting and analysis of ToT calibration
  var
    bins: seq[float]
    hist: seq[float]

  let (chip, tot) = parseTotInput(args, startTot)
  let totCalib = fitToTCalib(tot, startFit)
  plotToTCalib(totCalib, tot, chip)

proc main() =

  let args = docopt(doc)
  echo args
  let db = $args["--db"]
  let file = $args["--file"]
  let folder = $args["--folder"]
  let chip = $args["--chip"]
  let startFitStr = $args["--startFit"]
  let startTotStr = $args["--startTot"]

  if db != "nil":
    when not declared(ingridDatabase):
      quit("Cannot import InGrid database. --db option not supported.")

  let
    scurve = ($args["--scurve"]).parseBool
    tot = ($args["--tot"]).parseBool
  if scurve == true:
    sCurve(args, chip)
  elif tot == true:
    var startFit = 0.0
    var startTot = 0.0    
    if startFitStr != "nil":
      startFit = startFitStr.parseFloat
    if startTotStr != "nil":
      startTot = startTotStr.parseFloat
    totCalib(args, startFit, startTot)


when isMainModule:
  main()
