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

          #tspan(x=50, dy="1.2em"):
          #  t: "Hello World"
          #tspan(x=50, dy="1.2em"):
          #  t: "Ok!"
  echo nodes

  var xmlNodes = nodes.render.parseXml

  #echo xmlNodes
  #echo xmlNodes.len

  let xmlPlot = loadXml svgFname
  echo xmlPlot.len

  let att = {"transform" : "translate(850, 150)"}.toXmlAttributes

  #for x in mitems(xmlNodes):
  #  var tups: StringTableRef
  #  tups = x.attrs
  #  tups["xml:space"] = "preserved"
  #  x.attrs = tups.XmlAttributes

  var nnNew = newXmlTree("g", xmlNodes.mapIt(it), att)
  #for x in mitems(nnNew):
  #  x.attr("xml:space") = "preserve"
    #echo "X is ", x.attrs
  #for x in xmlNodes:
  #  nnNew.add x
  var nnPNew = newElement("g")
  for x in xmlPLot:
    nnPNew.add x

  for i in 0 ..< xmlPlot.len:
    xmlPlot.delete(0)

  xmlPlot.add nnPNew
  xmlPlot.add nnNew
  echo xmlPlot.len
  #echo xmlPlot

  var f = open("test_" & param.svg, fmWrite)
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

#buildSvgFile("examples/text1.svg"):
#  let s = "b"
#  svg(width=100, height=100, xmlns="http://www.w3.org/2000/svg", version="1.1", baseProfile="full"):
#    text(x=50,
#         y=50,
#         `text-anchor`="middle",
#         `dominant-baseline`="central",
#         style="font-family: 'Open Sans'"):
#         #transform="rotate(-90 50 50)"):
#      tspan(x=50, dy="1.2em"):
#        t: "Hello World"
#      tspan(x=50, dy="1.2em"):
#        t: "Ok!"
