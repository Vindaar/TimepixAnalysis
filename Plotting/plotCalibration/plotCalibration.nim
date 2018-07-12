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
import ingrid/tos_helpers
import ingrid/calibration

const doc = """
A simple tool to plot SCurves or ToT calibrations.

Usage:
  plotCalibration (--scurve | --tot) (--file=FILE | --folder=FOLDER) [options]

Options:
  --scurve         If set, perform SCurve analysis
  --tot            If set, perform ToT calibration analysis
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

type
  FitResult = object
    x: seq[float]
    y: seq[float]
    pRes: seq[float]
    pErr: seq[float]
    redChiSq: float

const
  NTestPulses = 1000.0
  ScalePulses = 1.01
  CurveHalfWidth = 15
  StopToT = 450.0

const
  GoldenMean = (sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
  FigWidth = 1200.0                     # width in inches
  FigHeight = FigWidth * GoldenMean     # height in inches

proc getTrace[T](bins, hist: seq[T], voltage: string): Trace[T] =
  result = Trace[T](`type`: PlotType.Scatter)
  # filter out clock cycles larger 300 and assign to `Trace`
  result.xs = bins
  result.ys = hist
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
  p.show()

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
  p.show()

proc plotToTCalib*(totCalib: FitResult, pulses, mean, std: seq[float], chip = 0) =
  let
    data = Trace[float](mode: PlotMode.Markers, `type`: PlotType.Scatter)
    fit = Trace[float](mode: PlotMode.Lines, `type`: PlotType.Scatter)
  # flip the plot, i.e. show THL on x instead of y as done for the fit to
  # include the errors
  data.xs = pulses
  data.ys = mean
  let stdValid = std.mapIt(if classify(it) == fcNaN: 0.0 else: it)
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
  p.show() 

proc readVoltageFile(filename: string): (string, seq[float], seq[float]) =
  let file = filename.expandTilde
  # - read file as string
  # - split all lines after header at \n
  # - filter lines with no content
  # - create tuple of (THL, Counts) for each line
  let dataTuple = readFile(file).splitLines[2..^1].filterIt(it.len > 0).mapIt(
    ((it.split('\t')[0].parseFloat, it.split('\t')[1].parseFloat))
  )
  result[0] = file.extractFilename.strip(chars = {'a'..'z', '_', '.'})
  result[1] = dataTuple.mapIt(it[0])
  result[2] = dataTuple.mapIt(it[1])

proc sCurve(file, folder, chip: string) =
  ## perform plotting and fitting of SCurves
  var
    voltages: set[int16]
    v = ""
    bins: seq[float] = @[]
    hist: seq[float] = @[]
    traces: seq[Trace[float]] = @[]

    # calibration seqs
    charge: seq[float] = @[]
    thlMean: seq[float] = @[]
    thlErr: seq[float] = @[]

  if file != "nil":
    (v, bins, hist) = readScurveVoltageFile(file)
    voltages.incl int16(v.parseInt)
    let trace = getTrace(bins, hist, v)
    traces.add trace
  else:
    echo "folder is ", folder
    for f in walkFiles(folder.expandTilde & "/*.txt"):
      echo f
      (v, bins, hist) = readScurveVoltageFile(f)
      voltages.incl int16(v.parseInt)
      let trace = getTrace(bins, hist, v)
      traces.add trace

      if v.parseInt > 0:
        # now fit the normal cdf to the SCurve and add its data
        let fitRes = fitSCurve(bins, hist, v.parseInt)
        let
          traceFit = getTrace(fitRes.x, fitRes.y, &"fit {v} / chiSq {fitRes.redChiSq}")
        traces.add traceFit

        # given `fitRes`, add fit parameters to calibration seqs
        # why do we multiply voltage in mV by 50?
        charge.add (v.parseFloat * 50)
        thlMean.add (fitRes.pRes[1])
        thlErr.add (fitRes.pErr[1])

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

  let thlCalib = fitThlCalib(chSort, thlSort, thlErrSort)
  plotThlCalib(thlCalib, chSort, thlSort, thlErrSort, chip)

proc totCalib(file, folder, chip: string, startFit = 0.0, startTot = 0.0) =
  ## perform plotting and analysis of ToT calibration
  var
    bins: seq[float]
    hist: seq[float]

  if file != "nil":
    let (chip, pulses, mean, std) = readToTFile(file, startTot)
    let totCalib = fitToTCalib(pulses, mean, std, startFit)
    plotToTCalib(totCalib, pulses, mean, std, chip)
    
    
    #voltages.incl int16(v.parseInt)
    #let trace = getTrace(bins, hist, v)
    #traces.add trace
  else:
    echo "folder is ", folder
    for f in walkFiles(folder.expandTilde & "/*.txt"):
      let (chip, pulses, mean, std) = readToTFile(f)
      #voltages.incl int16(v.parseInt)
      #let trace = getTrace(bins, hist, v)
      #traces.add trace


proc main() =

  let args = docopt(doc)
  let file = $args["--file"]
  let folder = $args["--folder"]
  let chip = $args["--chip"]
  let startFitStr = $args["--startFit"]
  let startTotStr = $args["--startTot"]  

  let
    scurve = ($args["--scurve"]).parseBool
    tot = ($args["--tot"]).parseBool
  if scurve == true:
    sCurve(file, folder, chip)
  elif tot == true:
    var startFit = 0.0
    var startTot = 0.0    
    if startFitStr != "nil":
      startFit = startFitStr.parseFloat
    if startTotStr != "nil":
      startTot = startTotStr.parseFloat
    totCalib(file, folder, chip, startFit, startTot)


when isMainModule:
  main()
