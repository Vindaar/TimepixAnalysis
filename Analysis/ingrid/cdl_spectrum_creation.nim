import parsecsv, os, streams, strutils, strformat, nimhdf5, tables, sequtils, macros
import plotly, mpfit
import ingrid / [ingrid_types, tos_helpers, calibration]
import docopt
import helpers / utils

const docStr = """
Usage:
  cdl_spectrum_creation <h5file> [options]
  cdl_spectrum_creation -h | --help
  cdl_spectrum_creation --version

Options:
  -h, --help   Show this help
  --version    Show the version number
"""
const doc = withDocopt(docStr)

type
  TargetKind = enum
    tEmpty = ""
    tCu = "Cu"
    tMn = "Mn"
    tTi = "Ti"
    tAg = "Ag"
    tAl = "Al"
    tC = "C"

  FilterKind = enum
    fEmpty = ""
    fEpic = "EPIC"
    fCr = "Cr"
    fNi = "Ni"
    fAg = "Ag"
    fAl = "Al"
    fTi = "Ti"

  TargetFilterKind = enum
    tfCuNi15 = "Cu-Ni-15kV"
    tfMnCr12 = "Mn-Cr-12kV"
    tfTiTi9 = "Ti-Ti-9kV"
    tfAgAg6 = "Ag-Ag-6kV"
    tfAlAl4 = "Al-Al-4kV"
    tfCuEpic2 = "Cu-EPIC-2kV"
    tfCuEpic0_9 = "Cu-EPIC-0.9kV"
    tfCEpic0_6 =  "C-EPIC-0.6kV"

  CdlRun = object
    number: int
    runType: RunTypeKind
    hasFadc: bool
    target: TargetKind
    filter: FilterKind
    hv: float

  CdlFitFunc = proc(p_ar: seq[float], x: float): float

  CdlFit = object
    target: TargetKind
    filter: FilterKind
    hv: float
    fit: CdlFitFunc

  ##dont need this anymore?!
    #cutsB = object
    #  target: TargetKind
    #  filter: FilterKind
    #  hv: float
    #  ck: int
    #  length: float
    #  rms_min: float
    #  rms_max: float
    #  eccentricity: float
    #
    #cutsC = object
    #  target: TargetKind
    #  filter: FilterKind
    #  hv: float
    #  ck: int
    #  charge_min: float
    #  charge_max: float
    #  length: float
    #  rms_min: float
    #  rms_max: float

  FitFuncKind = enum
    ffConst, ffPol1, ffPol2, ffGauss, ffExpGauss

  FitFuncArgs = object
    name: string
    case kind: FitFuncKind
    of ffConst:
      c: float
    of ffPol1:
      cp: float
    of ffPol2:
      cpp: float
    of ffGauss:
      gN: float
      gmu: float
      gs: float
    of ffExpGauss:
      ea: float
      eb: float
      eN: float
      emu: float
      es: float

func getLines(hist, binning: seq[float], tfKind: TargetFilterKind): seq[FitFuncArgs] =
  ## this is a runtime generator for the correct fitting function prototype,
  ## i.e. it returns a seq of parts, which need to be combined to the complete
  ## function at runtime
  let muIdx = argmax(hist)
  case tfKind
  of tfCuNi15: #need to add for rest combinatons
    result.add FitFuncArgs(name: "Cu-esc",
                           kind: ffExpGauss,
                           ea: hist[muIdx],
                           eb: hist[muIdx],
                           eN: hist[muIdx] / 3.0,
                           emu: binning[muIdx] * 2.0,
                           es: hist[muIdx] / 30.0)
    result.add FitFuncArgs(name: "Cu-Kalpha",
                           kind: ffExpGauss,
                           ea: hist[muIdx],
                           eb: hist[muIdx],
                           eN: hist[muIdx],
                           emu: binning[muIdx],
                           es: hist[muIdx] / 30.0)
  of tfMnCr12:
    result.add FitFuncArgs(name: "Mn-esc",
                           kind: ffExpGauss,
                           ea: hist[muIdx],
                           eb: hist[muIdx],
                           eN: hist[muIdx] / 3.0,
                           emu: binning[muIdx] / 2.0,
                           es: hist[muIdx] / 30.0)
    result.add FitFuncArgs(name: "Mn-Kalpha",
                           kind: ffExpGauss,
                           ea: hist[muIdx],
                           eb: hist[muIdx],
                           eN: hist[muIdx],
                           emu: binning[muIdx],
                           es: hist[muIdx] / 30.0)
  of tfTiTi9:
    result.add FitFuncArgs(name: "Ti-esc-alpha",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx],
                           gs: hist[muIdx] / 30.0)
    result.add FitFuncArgs(name: "Ti-esc-beta",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx],
                           gs: hist[muIdx] / 30.0)
    result.add FitFuncArgs(name: "Ti-Kalpha",
                           kind: ffExpGauss,
                           ea: hist[muIdx],
                           eb: hist[muIdx],
                           eN: hist[muIdx] / 3.0,
                           emu: binning[muIdx] * 2.0,
                           es: hist[muIdx] / 30.0)
    result.add FitFuncArgs(name: "Ti-Kbeta",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx],
                           gs: hist[muIdx] / 30.0)

  of tfAgAg6:
   result.add FitFuncArgs(name: "Ag-Lalpha",
                          kind: ffExpGauss,
                          ea: hist[muIdx],
                          eb: hist[muIdx],
                          eN: hist[muIdx] / 3.0,
                          emu: binning[muIdx] * 2.0,
                          es: hist[muIdx] / 30.0)
   result.add FitFuncArgs(name: "Ag-Lbeta",
                          kind: ffGauss,
                          gmu: binning[muIdx],
                          gN: hist[muIdx],
                          gs: hist[muIdx] / 10.0)
  of tfAlAl4:
    result.add FitFuncArgs(name: "Al-Kalpha",
                          kind: ffExpGauss,
                          ea: hist[muIdx],
                          eb: hist[muIdx],
                          eN: hist[muIdx] / 3.0,
                          emu: binning[muIdx] * 2.0,
                          es: hist[muIdx] / 30.0)
  of tfCuEpic2:
   result.add FitFuncArgs(name: "Cu-Lalpha",
                          kind: ffGauss,
                          gmu: binning[muIdx],
                          gN: hist[muIdx],
                          gs: hist[muIdx] / 10.0)
  of tfCuEpic0_9:
   result.add FitFuncArgs(name: "",
                          kind: ffGauss,
                          gmu: binning[muIdx],
                          gN: hist[muIdx],
                          gs: hist[muIdx] / 5.0)
   result.add FitFuncArgs(name: "",
                          kind: ffGauss,
                          gmu: binning[muIdx],
                          gN: hist[muIdx],
                          gs: hist[muIdx] / 10.0)
   result.add FitFuncArgs(name: "",
                          kind: ffGauss,
                          gmu: binning[muIdx],
                          gN: hist[muIdx],
                          gs: hist[muIdx] / 15.0)
   result.add FitFuncArgs(name: "",
                          kind: ffGauss,
                          gmu: binning[muIdx],
                          gN: hist[muIdx],
                          gs: hist[muIdx] / 20.0)
   # TODO: FIX correct function
  of tfCEpic0_6:
    result.add FitFuncArgs(name: "C-Kalpha",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx],
                           gs: hist[muIdx] / 30.0)
    result.add FitFuncArgs(name: "O-Kalpha",
                           kind: ffGauss,
                           gN: hist[muIdx] / 3.0,
                           gmu: binning[muIdx] * 2.0,
                           gs: hist[muIdx] / 30.0)
  #else:
    #discard


proc genFitFuncImpl(resultNode, idx, paramsNode, xNode, pFitNode: NimNode): NimNode =
  ## the compilet time procedure that creates the implementation lines for
  ## the <target><Filter>Funcs that we create, which calls the correct functions,
  ## e.g.
  ##   result += p_ar[i] * gauss(x, p[i + 1], p[i + 2])
  ##   inc i, 3
  expectKind(pFitNode, nnkObjConstr)
  let fkind = FitFuncKind(pFitNode[2][1].intVal)
  case fKind
  of ffConst:
    result = quote do:
      `resultNode` += `paramsNode`[`idx`]
      inc `idx`
  of ffPol1:
    result = quote do:
      `resultNode` += `paramsNode`[`idx`] * `xNode`
      inc `idx`
  of ffPol2:
    result = quote do:
      `resultNode` += `paramsNode`[`idx`] * `xNode` * `xNode`
      inc `idx`
  of ffGauss:
    result = quote do:
      `resultNode` += `paramsNode`[`idx`] * gauss(`xNode`,
                                                  `paramsNode`[`idx` + 1],
                                                  `paramsNode`[`idx` + 2])
      inc `idx`, 3 # increase by number of consumed parameters
  of ffExpGauss:
    result = quote do:
      `resultNode` += expGauss(`paramsNode`[`idx` .. `idx` + 4], `xNode`)
      inc `idx`, 5
  #else: discard

macro buildFitFunc(name: untyped, parts: openArray[FitFuncArgs]): untyped =
  ## builds a CDL fit function based on the function described by
  ## the `seq[FitFuncArgs]` at compile time. Using the `FitFuncKind` of
  ## each part, it'll write the needed implementation lines for the
  ## call to the correct functions, e.g. `gauss`, `expGauss` etc.
  # define the variables needed in the implementation function and
  # for the parameters
  let
    idx = ident"i"
    paramsNode = ident"p_ar"
    xNode = ident"x"
    resultNode = ident"result"
  # define parameters and return type of the proc we create
  let
    retType = ident"float"
    retParNode = nnkIdentDefs.newTree(paramsNode,
                                      nnkBracketExpr.newTree(
                                        ident"seq",
                                        ident"float"),
                                      newEmptyNode())
    retXNode = nnkIdentDefs.newTree(xNode,
                                    ident"float",
                                    newEmptyNode())
  # create a node to hold the procedure body
  var procBody = newStmtList()
  # declare the index variable we use
  procBody.add quote do:
    var `idx` = 0
    #debugecho `paramsNode`

  for p in parts.getImpl:
    # add the lines for the function calls
    procBody.add genFitFuncImpl(resultNode, idx, paramsNode, xNode, p)

  # now define the result variable as a new proc
  result = newProc(name = name,
                   params = [retType, retParNode, retXNode],
                   body = procBody,
                   procType = nnkFuncDef)
  #echo result.repr

macro declareFitFunc(name, stmts: untyped): untyped =
  ## DSL to declare the fit functions without having to go the
  ## const of `seq[FitFuncArgs]` + buildFitFunc macro route.
  ##
  ## .. code-block::
  ##   declareFitFunc(cEpic):
  ##     ffGauss: "C-Kalpha"
  ##     ffGauss: "O-Kalpha"
  ##   # will be expanded to
  ##   var cEpicF {.compileTime.} = @[FitFuncArgs(name: "C-Kalpha", kind: ffGauss),
  ##                                  FitFuncArgs(name: "O-Kalpha", kind: ffGauss)]
  ##   buildFitFunc(cEpicFunc, cEpicF)

  let
    thVarId = ident(name.strVal & "F_Mangle")
    funcId = ident(name.strVal & "Func")
  var ffSeq = nnkBracket.newTree()

  for s in stmts:
    expectKind(s, nnkCall)
    echo s.treeRepr
    let fkind = s[0]
    let ffName = s[1][0].strVal
    ffSeq.add quote do:
      FitFuncArgs(name: `ffName`, kind: `fKind`)

  result = quote do:
    const `thVarId` = `ffSeq`
    buildFitFunc(`funcId`, `thVarId`)


declareFitFunc(cuNi15):
  ffExpGauss: "Cu-esc"
  ffExpGauss: "Cu-Kalpha"
declareFitFunc(mnCr12):
  ffExpGauss: "Mn-esc"
  ffExpGauss: "Mn-Kalpha"
declareFitFunc(tiTi9):
  ffGauss: "Ti-esc-alpha"
  ffGauss: "Ti-esc-beta"
  ffExpGauss: "Ti-Kalpha"
  ffGauss: "Ti-Kbeta"
declareFitFunc(agAg6):
  ffexpGauss: "Ag-Lalpha"
  ffGauss: "Ag-Lbeta"
declareFitFunc(alAl4):
  ffexpGauss: "Al-Kalpha"
declareFitFunc(cuEpic2):
  ffGauss: "Cu-Lalpha"
declareFitFunc(cuEpic0_9):
  ffGauss: "O-Kalpha"
  ffGauss: "C-Kalpha"
  ffGauss: "Fe-Lalphabeta"
  ffGauss: "Ni-Lalphabeta"
declareFitFunc(cEpic0_6):
  ffGauss: "C-Kalpha"
  ffGauss: "O-Kalpha"
    #gmu = 17.0 ##maybe declare some fixed parameter here later


# old code for reference and understanding. TODOs still relevant
# TODO: add additional constants for other pairs
#const cEpicF = @[FitFuncArgs(name: "C-Kalpha", kind: ffGauss),
#                 FitFuncArgs(name: "O-Kalpha", kind: ffGauss)]

# TODO: call buildFitProc for the other pairs
# buildFitProc(cEpicFuncAlt, cEpicF)
# buildFitFunc(cEpicFunc, cEpicF)

# TODO: call the template above with some custom functions mainly consisting
# e.g. of linears and parabola functions and see if they evaluate correctly!
# just generate some dots and plot them and compare with manual defintion

proc serialize(parts: seq[FitFuncArgs]): seq[float] =
  for p in parts:
    case p.kind
    of ffGauss:
      result.add @[p.gN, p.gmu, p.gs]
    of ffExpGauss:
      result.add @[p.ea, p.eb, p.eN, p.emu, p.es]
    else: discard


macro genTfToFitFunc(pname: untyped): untyped =
  let tfkind = getType(TargetFilterKind)
  # first generate the string combinations
  #echo tfKind.treeRepr
  var funcNames: seq[string]
  for x in tfKind:
    if x.kind != nnkEmpty:
      let xStr = ($(x.getImpl))
        .toLowerAscii
        .replace("-", "")
        .replace(".", "")
        .replace("kv", "") & "Func"
      funcNames.add xStr
  # given the names, write a proc that returns the function
  let
    # arg1 = ident"target"
    # argt1 = ident"TargetKind"
    # arg2 = ident"filter"
    # argt2 = ident"FilterKind"
    ##now with target and filter combined
    arg = ident"tfKind"
    argType = ident"TargetFilterKind"
    cdf = ident"CdlFitFunc"
    tfNameNode = ident"n"
    resIdent = ident"result"
  var caseStmt = nnkCaseStmt.newTree(tfNameNode)
  for n in funcNames:
    let retId = ident(n)
    let retval = quote do:
      `resIdent` = `retId`
    caseStmt.add nnkOfBranch.newTree(newLit n, retval)
  result = quote do:
    proc `pname`(`arg`: `argType`): `cdf` =
      let `tfNameNode` = ($`arg`)
        .toLowerAscii
        .replace("-", "")
        .replace(".", "")
        .replace("kv", "") & "Func"
      `caseStmt`
  echo result.repr

# generate the =getCdlFitFunc= used to get the correct fit function
# based on a `TargetKind` and `FilterKind`
genTfToFitFunc(getCdlFitFunc)


proc histoCdl(data: seq[SomeInteger], binSize: float = 3.0): (seq[float], seq[float]) =

  let low = 0.0#-0.5
  var high = max(data).float# + 0.5
  let nbins = (ceil((high - low) / binSize)).int
  # using correct nBins, determine actual high
  high = low + binSize * nbins.float
  let bin_edges = linspace(low, high, nbins)# + 1)
  let hist = data.histogram(bins = nbins + 1, range = (low, high))

  echo bin_edges.len
  echo hist.len
  echo bin_edges

  result[0] = hist.mapIt(it.float)
  result[1] = bin_edges# [0 .. ^1]

proc fitCdlImpl(hist, binedges: seq[float], tfKind: TargetFilterKind): seq[float] =
  let lines = getLines(hist, binedges, tfKind)
  let params = lines.serialize

  let err = hist.mapIt(1.0 / sqrt(it))

  let (pRes, res) = fit(getCdlFitFunc(tfKind),
                        params,
                        binedges,
                        hist,
                        err)

  echoResult(pRes, res=res)
  result = pRes


const cuni15FuncCharge = cuNi15Func
const mnCr12FuncCharge = mnCr12Func
const tiTi9FuncCharge = tiTi9Func
const agAg6FuncCharge = agAg6Func
const alAl4FuncCharge = alAl4Func
const cuEpic2FuncCharge = cuEpic2Func
const cuEpic0_9FuncCharge = cuEpic0_9Func
const cEpic0_6FuncCharge = cEpic0_6Func

# TODO: impl rest of functions + Charge functions

# proc getCdlFits(): Table[string, CdlFitFunc] =
  # result = {"C-EPIC-0.6kV" : cEpicFunc}.toTable
          #{ "C-EPIC-0.6kV" : cEpicFunc,
             #"Cu-EPIC-0.9kV" : cuEpicLowFunc,
             #"Cu-EPIC-2kV" : cuEpicHighFunc,
             #"Al-Al-4kV" : alAlFunc,
             #"Ag-Ag-6kV" : agAgFunc,
             #"Ti-Ti-9kV" : tiTiFunc,
             #"Mn-Cr-12kV" : mnCrFunc,
             #"Cu-Ni-15kV" : cuNiFunc}.toTable

proc toCutStr(run: CdlRun): string =
  let hv = block:
    if run.hv > 1.0 and run.hv < 10.0:
      &"{run.hv:1}"
    elif run.hv > 10.0:
      &"{run.hv:2}"
    else:
      &"{run.hv:1.1f}"
  result = &"{run.target}-{run.filter}-{hv}kV"

proc readRuns(fname: string): seq[CdlRun] =
  var s = newFileStream(fname, fmRead)
  var parser: CsvParser
  if not s.isNil:
    parser.open(s, fname, separator = '|')
    parser.readHeaderRow()
    discard parser.readRow()
    while parser.readRow:
      let row = parser.row
      let run = CdlRun(number: row[1].strip.parseInt,
                       runType: parseEnum[RunTypeKind](row[2].strip, rtNone),
                       hasFadc: row[3].strip.parseBool,
                       target: parseEnum[TargetKind](row[4].strip, tEmpty),
                       filter: parseEnum[FilterKind](row[5].strip, fEmpty),
                       hv: if row[6].strip.len > 0: row[6].strip.parseFloat else: 0.0)
      result.add run
      #echo run.toCutStr

proc totfkind(run: CdlRun): TargetFilterKind =
  result = parseEnum[TargetFilterKind](&"{toCutStr(run)}")

proc main =

  let args = docopt(doc)
  echo args

  let h5file = $args["<h5file>"]

  const filename = "../../resources/cdl_runs_2019.org"
  const cutparams = "../../resources/params.org"

  let runs = readRuns(filename)

  var h5f = H5file(h5file, "rw")
  defer: discard h5f.close()
  let cutTab = getXraySpectrumCutVals()
  for r in runs:
    if r.number != 347:
     continue
    sleep 500
    case r.runType
    of rtXrayFinger:
      let grp = h5f[(recoDataChipBase(r.number) & "3").grp_str]
      let cut = cutTab[r.toCutStr]
      #echo cut
      let passIdx = cutOnProperties(h5f,
                                    grp,
                                    cut.cutTo,
                                    ("rmsTransverse", cut.minRms, cut.maxRms),
                                    ("length", 0.0, cut.maxLength),
                                    ("hits", cut.minPix, Inf),
                                    ("eccentricity", 0.0, cut.maxEccentricity))
      let nevents = passIdx.len
      proc writeDset(dsetWrite, dsetRead: string) =
        var
          dset = h5f.create_dataset(grp.name / dsetWrite, (nevents, 1),
                                    int64)
        if dsetWrite == "CdlSpectrumIndices":
          dset[dset.all] = passIdx
        else:
          let read = h5f[grp.name / dsetRead, int64]
          dset[dset.all] = passIdx.mapIt(read[it])
        dset.attrs["Target"] = $r.target
        dset.attrs["Filter"] = $r.filter
        dset.attrs["HV"] = $r.hv
      writeDset("CdlSpectrumIndices", "")
      writeDset("CdlSpectrum", "hits")
      writeDset("CdlSpectrumEvents", "eventNumber")

      # TODO: here we can now:
      # - create the CdlFit object for this CdlRun object making use of the
      #   `getCdlFitFunc` procedure
      # - call a `fitSpectrum` function, which performs the steps outlined in
      #   the Org file in the "Calling functions at runtime" chapter
      ##done with histoCdl and fitCdlImpl

      let hitsRawData = h5f[grp.name / "hits", int64]
      let hitsCdl = h5f[grp.name / "CdlSpectrum", int64]
      let (histdata, bins) = histoCdl(hitsCdl, binSize = 1.0)
      let tfk = r.totfkind
      let pRes = fitCdlImpl(histdata, bins, tfk)
      let fitres = bins.mapIt(getCdlFitFunc(tfk)(pRes, it))
      let cdlplot = scatterPlot(bins, fitres).mode(PlotMode.Lines)


      let hitsRaw = histPlot(hitsRawData.mapIt(it.float64))
        .binSize(1.0)
        .binRange(0.0, 400.0)
      let hitsCut = histPlot(hitsCDL.mapIt(it.float64))
        .binSize(1.0)
        .binRange(0.0, 400.0)
      hitsRaw.layout.barMode = BarMode.Overlay
      let plt = hitsRaw.addTrace(hitsCut.traces[0])
        .addTrace(cdlplot.traces[0])
      plt.layout.title = &"run number: {r.number} target: {r.toCutStr}"
      plt.layout.showlegend = true
      plt.traces[1].opacity = 0.5
      plt.traces[0].name = "raw data"
      plt.traces[1].name = "data with cuts"
      plt.traces[2].name = "fit curve"
      plt.layout.yaxis.title = "Occurence"
      plt.layout.xaxis.title = "Number of pixels"
      plt.show()

    else:
      discard

proc testdeclarefitmacro()=
  const
    mu1 = 2.0
    mu2 = 10.0
    d1  = 5.0
    d2  = 1.0
    N1  = 20.0
    N2  = 5.0

  let xa = linspace(0.0,20.0,1000)
  let p  = @[N1, mu1, d1, N2, mu2, d2]
  let y  = xa.mapIt(cEpic0_6Func(p, it))
  scatterPlot(xa,y).show()



when isMainModule:
  main()
