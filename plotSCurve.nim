import plotly
import os, strutils, strformat
import sequtils
import math
import docopt


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

proc getTrace[T](hist, bins: seq[T], voltage: string): Trace[T] =

  result = Trace[int](`type`: PlotType.Scatter)
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
    p = Plot[int](layout: layout, traces: traces)
  p.show()

proc readVoltageFile(filename: string): (string, seq[int], seq[int]) =
  let file = filename.expandTilde
  let dataTuple = readFile(file).splitLines[2..^1].filterIt(it.len > 0).mapIt(
    ((it.split('\t')[0].parseInt, it.split('\t')[1].parseInt))
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
    bins: seq[int] = @[]
    hist: seq[int] = @[]
    traces: seq[Trace[int]] = @[]
  if file != "nil":
    (v, bins, hist) = readVoltageFile(file)
    voltages.incl int16(v.parseInt)
    let trace = getTrace(hist, bins, v)
    traces.add trace
  else:
    echo "folder is ", folder
    for f in walkFiles(folder.expandTilde & "/*.txt"):
      echo f
      (v, bins, hist) = readVoltageFile(f)
      voltages.incl int16(v.parseInt)
      let trace = getTrace(hist, bins, v)
      traces.add trace
      
  plotHist(traces, voltages, chip)

when isMainModule:
  main()
