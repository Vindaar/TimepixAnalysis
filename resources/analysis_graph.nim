import nimgraphviz, strutils, sequtils, macros

#[
This little piece of code generates a GraphViz flow chart representing
the whole TPA analysis chain.
]#

proc sanitize(s: string): string =
  result = "\"" & s & "\""

proc connectNodes(g: var Graph, s: seq[string]) =
  for i in 0 ..< s.high:
    g.addEdge(s[i].sanitize, s[i + 1].sanitize)

template `<-`(lhs, rhs: untyped): untyped {.dirty.} =
  g.addEdge(rhs.sanitize, lhs.sanitize)
  lhs

func `<-`(g: var Graph, lhs, rhs: string) =
  g.addEdge(rhs.sanitize, lhs.sanitize)

func `<-`(g: var Graph, h: Graph) =
  g.addEdge(h.name.sanitize, g.name.sanitize,
            attrs = [("ltail", h.clusterName),
                     ("lhead", g.clusterName)])

macro graph(stmts: untyped): untyped =
  result = newStmtList()
  for stmt in stmts:
    echo stmt.treerepr
    expectKind(stmt, nnkCall)
    expectKind(stmt[1], nnkStmtList)
    let arg = stmt[1][0]
    expectKind(arg, nnkInfix)
    case arg[0].strVal
    of "<-":
      result.add nnkCall.newTree(bindSym"addEdge", stmt[0],
                                 nnkCall.newTree(bindSym"sanitize", arg[1]),
                                 nnkCall.newTree(bindSym"sanitize", arg[2]))
    of "->":
      result.add nnkCall.newTree(bindSym"addEdge", stmt[0],
                                 nnkCall.newTree(bindSym"sanitize", arg[2]),
                                 nnkCall.newTree(bindSym"sanitize", arg[1]))
    else:
      error("Invalid identifier: " & arg[0].repr & " for graph macro.")

var g = newGraph(directed = true)

# set some attributes of the graph:
g.graphAttr["fontsize"] = "32"
g.graphAttr["label"] = "Timepix Analysis pipeline"
g.graphAttr["compound"] = "true"

var sData = g.newSubgraph("readout")
let ascii = "ASCII files"
block Readout:
  let dataParts = @["GridPix", "scintillator", "FADC"]
  for d in dataParts:
    graph:
      sData: "FPGA" -> d
  let nodes = @["FPGA", "TOS", ascii]
  sData.connectNodes(nodes)

var types = g.newSubgraph("datasets")
block TypesOfDatasets:
  let dataTypeNodes = @["calibration (⁵⁵Fe)",
                        "background+tracking", #
                        "x-ray finger",
                        "CDL"]

  for d in dataTypeNodes:
    graph:
      types: d -> "datasets"

  types <- sData

let rawStr = "raw data manipulation"
var raw = g.newSubgraph(rawStr)
block RawDataManipulation:
  #graph:
  #  raw: "raw data manipulation" <- "data"
  raw <- types
  let nodes = @[ascii, rawStr, "HDF5 files"]
  raw.connectNodes(nodes)

let eCalib = "energy calibration"
block Reconstruction:
  # these are the nodes general to all data
  let reconstruction = @["cluster finding" <- "50 pix radius",
                       "eccentricity optimization",
                       "calc of other geometric properties",
                       "charge calibration" <- "ToT calibration",
                       "gas gain calculation" <- "Polya fit",
                       eCalib <- "Fit to all ⁵⁵Fe spectra vs. gas gain"]

  var reco = g.newSubgraph("reconstruction")
  reco.connectNodes(reconstruction)

  reco <- raw


let calibrationPath = @["generalNodes",
                        "⁵⁵Fe pixel fit",
                        "⁵⁵Fe charge fit",
                        "linear fit to photo+escape peak"]

block LikelihoodLimit:
  let recoCdl = "Reco'd CDL data"
  let CDLdata = "CDL data"
  let softwareSteps = @[ "likelihood" <- (recoCDL <- CDLdata) <- ("CDL Reference spectra" <- CDLdata),
                        "limit calculation"]
  g.connectNodes(softwareSteps)

  discard recoCdl <- eCalib



# Export graph as PNG:
g.exportImage("/tmp/test_G.svg")
