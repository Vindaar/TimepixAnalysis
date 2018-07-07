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

func findDrop(thl, count: seq[float]): (float, int, int) =
  ## search for the drop of the SCurve, i.e.
  ##
  ##    ___/\__
  ## __/       \__
  ## 0123456789ABC
  ## i.e. the location of feature A (the center of it)
  ## returns the thl value at the center of the drop
  ## and the ranges to be used for the fit (min and max
  ## as indices for the thl and count seqs)
  # zip thl and count
  let
    thlZipped = zip(toSeq(0 .. thl.high), thl)
    thlCount = zip(thlZipped, count)
    # helper, which is the scaled bound to check whether we need
    # to tighten the view on the valid THLs
    scaleBound = NTestPulses * ScalePulses

  # extract all elements  of thlCount where count larger than 500, return
  # the largest element
  let
    drop = thlCount.filterIt(it[1] > (NTestPulses / 2.0))[^1]
    (pCenterInd, pCenter) = drop[0]

  var minIndex = max(pCenterInd - CurveHalfWidth, 0)
  # test the min index
  if count[minIndex.int] > scaleBound:
    # in that case get the first index larger than minIndex, which
    # is smaller than the scaled test pulses
    let
      thlView = thlCount[minIndex.int .. thlCount.high]
      thlSmallerBound = thlView.filterIt(it[1] < scaleBound)
    # from all elements smaller than  the bound, extract the index,
    # [0][0][0]:
    # - [0]: want element at position 0, contains smallest THL values
    # - [0]: above results in ( (`int`, `float`), `float` ), get first tupl
    # - [0]: get `int` from above, which corresponds to indices of THL
    minIndex = thlSmallerBound[0][0][0]
  let maxIndex = min(pCenterInd + CurveHalfWidth, thl.high)

  # index is thl of ind
  result = (pCenter, minIndex, maxIndex)

func sCurveFunc(p: seq[float], x: float): float =
  ## we fit the complement of a cumulative distribution function
  ## of the normal distribution
  # parameter p[2] == sigma
  # parameter p[1] == x0
  # parameter p[0] == scale factor
  result = normalCdfC(x, p[2], p[1]) * p[0]

func thlCalibFunc(p: seq[float], x: float): float =
  ## we fit a linear function to the charges and mean thl values
  ## of the SCurves
  result = p[0] + x * p[1]

func totCalibFunc(p: seq[float], x: float): float =
  ## we fit a combination of a linear and a 1 / x function
  #debugecho "Call with ", p
  result = p[0] * x + p[1] - p[2] / (x - p[3])

proc fitSCurve[T](thl, count: seq[T], voltage: int): FitResult =
  ## performs the fit of the `sCurveFunc` to the given `thl` and `count`
  ## seqs. Returns a `FitScurve` object of the fit result
  const pSigma = 5.0
  let
    (pCenter, minIndex, maxIndex) = findDrop(thl, count)
    err = thl.mapIt(1.0)
    p = @[NTestPulses, pCenter.float, pSigma]

  let
    thlCut = thl[minIndex .. maxIndex].mapIt(it.float)
    countCut = count[minIndex .. maxIndex].mapIt(it.float)

    (pRes, res) = fit(sCurveFunc,
                      p,
                      thlCut,
                      countCut,
                      err)
  echoResult(pRes, res = res)
  # create data to plot fit as result
  result.x = linspace(thl[minIndex].float, thl[maxIndex].float, 1000)
  result.y = result.x.mapIt(sCurveFunc(pRes, it))
  result.pRes = pRes
  result.pErr = res.error
  result.redChiSq = res.reducedChiSq

proc fitThlCalib(charge, thl, thlErr: seq[float]): FitResult =

  # determine start parameters
  let p = @[0.0, (thl[1] - thl[0]) / (charge[1] - charge[0])]

  echo "Fitting ", charge, " ", thl, " ", thlErr

  let (pRes, res) = fit(thlCalibFunc, p, charge, thl, thlErr)
  # echo parameters
  echoResult(pRes, res = res)

  result.x = linspace(charge[0], charge[^1], 100)
  result.y = result.x.mapIt(thlCalibFunc(pRes, it))
  result.pRes = pRes
  result.pErr = res.error
  result.redChiSq = res.reducedChiSq

proc fitToTCalib(pulses, mean, std: seq[float], startFit = 0.0): FitResult =
  var
    # local mutable variables to potentially remove unwanted data for fit
    mPulses = pulses
    mMean = mean
    mStd = std
  
  # define the start parameters. Use a good guess...
  let p = [0.149194, 23.5359, 205.735, -100.0]
  var pLimitBare: mp_par
  var pLimit: mp_par
  pLimit.limited = [1.cint, 1]
  pLimit.limits = [-100.0, 0.0]

  if startFit > 0:
    # in this case cut away the undesired parameters
    let ind = mPulses.find(startFit) - 1
    if ind > 0:
      mPulses.delete(0, ind)
      mMean.delete(0, ind)
      mStd.delete(0, ind)      
  
  #let p = [0.549194, 23.5359, 50.735, -1.0]
  let (pRes, res) = fit(totCalibFunc,
                        p,
                        mPulses,
                        mMean,
                        mStd,
                        bounds = @[pLimitBare, pLimitBare, pLimitBare, pLimit])
  echoResult(pRes, res = res)

  # the plot of the fit is performed to the whole pulses range anyways, even if 
  result.x = linspace(pulses[0], pulses[^1], 100)
  result.y = result.x.mapIt(totCalibFunc(pRes, it))
  result.pRes = pRes
  result.pErr = res.error
  result.redChiSq = res.reducedChiSq  
  
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
