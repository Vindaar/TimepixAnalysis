import nimsvg
import docopt
import ingrid / cdl_spectrum_creation
import helpers / utils
import parseutils, strutils, os, xmlparser, xmltree, sequtils, strtabs

const docStr = """
Usage:
  embedFitParameters <svgPath> <fitParamsFile> [options]
  embedFitParameters -h | --help
  embedFitParameters --version

Options:
  -h, --help   Show this help
  --version    Show the version number
"""
const doc = withDocopt(docStr)

## This is a simple helper tool, which can be used to embed the fit paramters obtained
## for the CDL fits onto the SVG files created by plotly.
## We both need the path to the SVG files and the corresponding `fitparameters_<timestamp>.txt`
## file.

type
  Params = object
    tfKind: TargetFilterKind
    svg: string
    params: string # string containing the lines to dump on SVG

proc parseFitParams(fname: string): seq[Params] =
  ## parses the `fitParameter` file given to the program and stores the
  ## result in `Params` objects
  var param: Params
  var data = readFile(fname).splitLines
  var lineCnt = 0
  echo data, " from fname ", fname
  func incLine(line: var string,
               lineCnt: var int) =
    line = data[lineCnt]
    inc lineCnt
  var line: string
  incLine(line, lineCnt)
  var lineName: string
  let maxLen = data.mapIt(it.len).max
  while lineCnt < data.len:
    param = Params(svg: line.replace("svg: ", ""))
    incLine(line, lineCnt)
    param.tfKind = parseEnum[TargetFilterKind](line.replace("tfKind: ", ""))
    incLine(line, lineCnt)
    var pStr = ""
    while lineCnt < data.len and not line.startsWith("svg:"):
      if '=' notin line:
        pStr = pStr & "\n" & line.alignLeft(maxLen)
      else:
        pStr = pStr & "\n" & line
      incLine(line, lineCnt)
    param.params = pStr.strip
    result.add param
  echo "Result is ", result

proc dumpOnSvg(svgPath: string, param: Params) =
  ## reads the SVG file as XML node, builds the text to dump on using nimSVG
  ## and then combines the two for a new SVG.
  let svgFname = svgPath / param.svg
  for p in param.params.splitLines:
    echo "p ", p

  # create a svg
  let nodes = buildSvg:
    svg(width=100, height=100, xmlns="http://www.w3.org/2000/svg", version="1.1"):
      text(x=50,
           y=0,
           `text-anchor`="middle",
           `dominant-baseline`="central",
           style="font-family: 'Monospace'"):
           #transform="rotate(-90 50 50)"):
        for p in param.params.splitLines:
          tspan(x = 0, `white-space`="pre"):
            t: p.replace("\\pm", "±")

  # convert SVG to XML tree
  var fitDumpNodes = nodes.render.parseXml
  # load the SVG as XML tree
  let xmlPlot = loadXml svgFname

  # create a transform attribute, which will place the fit parameter dump
  # (`fitDumpNodes`) onto the `xmlPlot`
  let att = {"transform" : "translate(850, 150)"}.toXmlAttributes
  # create new tree combining `att`, `fitDumpNodes` into a `g` XML element
  var fitDumpTree = newXmlTree("g", fitDumpNodes.mapIt(it), att)
  # create new `g` element, into which we will put the old SVG plot
  var nnPNew = newElement("g")
  for x in xmlPlot:
    nnPNew.add x
  # delete the 0th entry on all `xmlPlot` children
  for i in 0 ..< xmlPlot.len:
    xmlPlot.delete(0)
  # finally readd everything via the `g` elements
  xmlPlot.add nnPNew
  xmlPlot.add fitDumpTree
  # write to file
  var f = open("fitDump_" & param.svg, fmWrite)
  f.write(xmlPlot)
  f.close()

proc main =
  let args = docopt(doc)
  echo args
  let svgPath = $args["<svgPath>"]
  let fitParamsFile = $args["<fitParamsFile>"]
  let params = parseFitParams(fitParamsFile)

  # given the params, iterate params, open SVG file and dump on it
  for p in params:
    dumpOnSvg(svgPath, p)


when isMainModule:
  main()
