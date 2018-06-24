import plotly
import os, strutils, strformat
import sequtils
import math
import docopt
import seqmath
import mpfit

const doc = """
A simple tool to plot an SCurve
  
Usage:
  plotSCurve (--file=FILE | --folder=FOLDER) [options]

Options:
  --file=FILE      If given will read from a single file
  --folder=FOLDER  If given will read all voltage files from the given folder 
  --chip=NUMBER    The number of this chip
  -h, --help       Show this help
  --version        Show the version number
"""

type
  FitScurve = object
    thl: seq[float]
    count: seq[float]
    pRes: seq[float]
    pErr: seq[float]
    redChiSq: float

const
  NTestPulses = 1000.0
  ScalePulses = 1.01
  CurveHalfWidth = 15

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

proc fitSCurve[T](thl, count: seq[T], voltage: int): FitScurve =
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
  result.thl = linspace(thl[minIndex].float, thl[maxIndex].float, 1000)
  result.count = result.thl.mapIt(sCurveFunc(pRes, it))
  result.pRes = pRes
  result.pErr = res.error
  result.redChisq = res.reducedChisq

proc getTrace[T](bins, hist: seq[T], voltage: string): Trace[T] =
  result = Trace[T](`type`: PlotType.Scatter)
  # filter out clock cycles larger 300 and assign to `Trace`
  result.xs = bins
  result.ys = hist
  result.name = &"Voltage {voltage}"

proc plotHist*[T](traces: seq[Trace[T]], voltages: set[int16], chip = "") =
  ## given a seq of scintillator counts, plot them as a histogram
  let 
    goldenMean = (sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
    figWidth = 1200.0                     # width in inches
    figHeight = figWidth * goldenMean     # height in inches

  # take the first trace of the seq as the one to extract the min and max
  # ranges. They should all be the same
  let
    binMin = traces[0].xs[0].float
    binMax = traces[0].xs[^1].float

  let 
    layout = Layout(title: &"SCurve of Chip {chip} for voltage {voltages}",
                    width: figWidth.int, height: figHeight.int,
                    xaxis: Axis(title: "Threshold value"),
                    yaxis: Axis(title: "# hits$"),
                    autosize: false)
    p = Plot[float](layout: layout, traces: traces)
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

proc main() =

  let args = docopt(doc)
  let file = $args["--file"]
  let folder = $args["--folder"]  
  let chip = $args["--chip"]

  var
    voltages: set[int16]
    v = ""
    bins: seq[float] = @[]
    hist: seq[float] = @[]
    traces: seq[Trace[float]] = @[]
  if file != "nil":
    (v, bins, hist) = readVoltageFile(file)
    voltages.incl int16(v.parseInt)
    let trace = getTrace(bins, hist, v)
    traces.add trace
  else:
    echo "folder is ", folder
    for f in walkFiles(folder.expandTilde & "/*.txt"):
      echo f
      (v, bins, hist) = readVoltageFile(f)
      voltages.incl int16(v.parseInt)
      let trace = getTrace(bins, hist, v)
      traces.add trace

      if v.parseInt > 0:
        # now fit the normal cdf to the SCurve and add its data
        let fitRes = fitSCurve(bins, hist, v.parseInt)
        let
          traceFit = getTrace(fitRes.thl, fitRes.count, &"fit {v} / chiSq {fitRes.redChiSq}")
        traces.add traceFit

        # given `fitRes`, add fit parameters to calibration seqs
        
  plotHist(traces, voltages, chip)

when isMainModule:
  main()
