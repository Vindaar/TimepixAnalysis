import parsecsv, os, streams, strutils, strformat, nimhdf5, tables, sequtils, macros
import seqmath, algorithm
import plotly, mpfit, nlopt, nimpy
import ingrid / [ingrid_types, tos_helpers, calibration]
import docopt
import helpers / utils
import chroma


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

const Color1 = color(1.0, 0.0, 102.0 / 256.0)
const Color2 = color(0.0, 153.0 / 256.0, 204 / 256.0)
const black = color(0.0, 0.0, 0.0)

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

const fixed = NaN

func getLines(hist, binning: seq[float], tfKind: TargetFilterKind): seq[FitFuncArgs] =
  ## this is a runtime generator for the correct fitting function prototype,
  ## i.e. it returns a seq of parts, which need to be combined to the complete
  ## function at runtime
  let muIdx = argmax(hist)
  case tfKind
  of tfCuNi15:
    result.add FitFuncArgs(name: "Cu-esc",
                           kind: ffExpGauss,
                           ea: hist[muIdx] * 1e-10,#fixed, ##1e-10 on laptop
                           eb: hist[muIdx] * 1e-12,#fixed,
                           eN: hist[muIdx] / 4.0,#10.0, ##4 on laptop
                           emu: 170.0,#binning[muIdx] ,#fixed, ##classic binning[muIdx] on laptop
                           es: 20.0)#hist[muIdx] / 30.0)#40.00)# ##classic hist[muIdx] on laptop
    result.add FitFuncArgs(name: "Cu-Kalpha",
                           kind: ffExpGauss,
                           ea: hist[muIdx] * 1e-10,#fixed, ##1e-10 on laptop
                           eb: hist[muIdx] * 1e-12,#fixed,
                           eN: hist[muIdx],#80.0
                           emu: binning[muIdx],#fixed,
                           es: hist[muIdx] / 20.0)#18.0)
  of tfMnCr12:
    result.add FitFuncArgs(name: "Mn-esc",
                           kind: ffExpGauss,
                           ea: -hist[muIdx] * 1e-10,
                           eb: -hist[muIdx] * 1e-12,
                           eN: hist[muIdx] / 9.0,#50.0,
                           emu: binning[muIdx] / 2.0,#160.0,
                           es: hist[muIdx] / 30.0)#16.0)
    result.add FitFuncArgs(name: "Mn-Kalpha",
                           kind: ffExpGauss,
                           ea: -hist[muIdx] * 1e-10,
                           eb: -hist[muIdx] * 1e-12,
                           eN: hist[muIdx],#350.0,
                           emu: binning[muIdx],#200.0,
                           es: hist[muIdx] / 15.0)#13.0)
  of tfTiTi9:
    result.add FitFuncArgs(name: "Ti-esc-alpha",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx],
                           gs: hist[muIdx] / 15.0) ##30 on laptop
    result.add FitFuncArgs(name: "Ti-esc-beta",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx],
                           gs: hist[muIdx] / 15.0) ##30 on laptop
    result.add FitFuncArgs(name: "Ti-Kalpha",
                           kind: ffExpGauss,
                           ea: hist[muIdx] * 1e-12,
                           eb: -hist[muIdx] * 1e-12,
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
                          ea: hist[muIdx] * 1e-10,
                          eb: -hist[muIdx] * 1e-12,
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
                          ea: hist[muIdx]* 1e-10,
                          eb: -hist[muIdx] * 1e-12,
                          eN: hist[muIdx],
                          emu: binning[muIdx],
                          es: hist[muIdx] )
  of tfCuEpic2:
   result.add FitFuncArgs(name: "Cu-Lalpha",
                          kind: ffGauss,
                          gmu: binning[muIdx],
                          gN: hist[muIdx],
                          gs: hist[muIdx] )

   result.add FitFuncArgs(name: "Cu-Lbeta",
                          kind: ffExpGauss,
                          ea: hist[muIdx]* 1e-10,
                          eb: -hist[muIdx] * 1e-12,
                          eN: hist[muIdx],
                          emu: binning[muIdx],
                          es: hist[muIdx] / 30.0)
  of tfCuEpic0_9:
   result.add FitFuncArgs(name: "O-Kalpha",
                          kind: ffGauss,
                          gmu: binning[muIdx],
                          gN: hist[muIdx],
                          gs: hist[muIdx] )
   #result.add FitFuncArgs(name: "C-Kalpha",
   #                       kind: ffGauss,
   #                       gmu: binning[muIdx],
   #                       gN: hist[muIdx],
   #                       gs: hist[muIdx] )
   result.add FitFuncArgs(name: "O-Lbeta",
                          kind: ffExpGauss,
                          ea: hist[muIdx]* 1e-10,
                          eb: -hist[muIdx] * 1e-12,
                          eN: hist[muIdx],
                          emu: binning[muIdx],
                          es: hist[muIdx] / 20.0)
   #result.add FitFuncArgs(name: "Fe-Lalphabeta",
   #                       kind: ffGauss,
   #                       gmu: binning[muIdx],
   #                       gN: hist[muIdx],
   #                       gs: hist[muIdx] )
   #result.add FitFuncArgs(name: "Ni-Lalphabeta",
   #                       kind: ffGauss,
   #                       gmu: binning[muIdx],
   #                       gN: hist[muIdx],
   #                       gs: hist[muIdx] )
  of tfCEpic0_6:
    result.add FitFuncArgs(name: "C-Kalpha",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx] / 2.0,
                           gs: hist[muIdx] / 2.0)
    result.add FitFuncArgs(name: "O-Kalpha",
                           kind: ffGauss,
                           gN: hist[muIdx] ,
                           gmu: binning[muIdx] , #
                           gs: hist[muIdx] )#fixed) #


func getLinesCharge(hist, binning: seq[float], tfKind: TargetFilterKind): seq[FitFuncArgs] =
  ## this is a runtime generator for the correct fitting function prototype,
  ## i.e. it returns a seq of parts, which need to be combined to the complete
  ## function at runtime
  let muIdx = argmax(hist)
  case tfKind
  of tfCuNi15:
    result.add FitFuncArgs(name: "Cu-esc",
                           kind: ffGauss,
                           gN: hist[muIdx] / 5.0,
                           gmu: binning[muIdx] / 2.0,
                           gs: hist[muIdx] * 1e3)
    result.add FitFuncArgs(name: "Cu-Kalpha",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx],
                           gs: hist[muIdx] * 1e3)
  of tfMnCr12:
    result.add FitFuncArgs(name: "Mn-esc",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx],
                           gs: hist[muIdx] / 15.0)
    result.add FitFuncArgs(name: "Mn-Kalpha",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx],
                           gs: hist[muIdx] / 15.0)
    result.add FitFuncArgs(name: "p0",
                           kind: ffConst,
                           c: hist[muIdx])
    result.add FitFuncArgs(name: "p1",
                           kind: ffPol1,
                           cp: hist[muIdx])
    result.add FitFuncArgs(name: "p2",
                           kind: ffPol2,
                           cpp: hist[muIdx])
  of tfTiTi9:
    result.add FitFuncArgs(name: "Ti-esc-alpha",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx],
                           gs: hist[muIdx] / 15.0)
    result.add FitFuncArgs(name: "Ti-esc-beta",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx],
                           gs: hist[muIdx] / 15.0)
    result.add FitFuncArgs(name: "Ti-Kalpha",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx],
                           gs: hist[muIdx] / 15.0)
    result.add FitFuncArgs(name: "Ti-Kbeta",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx],
                           gs: hist[muIdx] / 30.0)
  of tfAgAg6:
    result.add FitFuncArgs(name: "Ag-Lalpha",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx],
                           gs: hist[muIdx] / 15.0)
    result.add FitFuncArgs(name: "Ag-Lbeta",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx],
                           gs: hist[muIdx] / 15.0)
    result.add FitFuncArgs(name: "p0",
                           kind: ffConst,
                           c: hist[muIdx])
    result.add FitFuncArgs(name: "p1",
                           kind: ffPol1,
                           cp: hist[muIdx])
    result.add FitFuncArgs(name: "p2",
                           kind: ffPol2,
                           cpp: hist[muIdx])
  of tfAlAl4:
    result.add FitFuncArgs(name: "Ag-Kaplha",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx],
                           gs: hist[muIdx] / 15.0)
    result.add FitFuncArgs(name: "p0",
                           kind: ffConst,
                           c: hist[muIdx])
    result.add FitFuncArgs(name: "p1",
                           kind: ffPol1,
                           cp: hist[muIdx])
    result.add FitFuncArgs(name: "p2",
                           kind: ffPol2,
                           cpp: hist[muIdx])
  of tfCuEpic2:
   result.add FitFuncArgs(name: "Cu-Lalpha",
                          kind: ffGauss,
                          gmu: binning[muIdx],
                          gN: hist[muIdx],
                          gs: hist[muIdx] )
  of tfCuEpic0_9:
   result.add FitFuncArgs(name: "O-Kalpha",
                          kind: ffGauss,
                          gmu: binning[muIdx],
                          gN: hist[muIdx],
                          gs: hist[muIdx] )
  of tfCEpic0_6:
    result.add FitFuncArgs(name: "C-Kalpha",
                           kind: ffGauss,
                           gN: hist[muIdx],
                           gmu: binning[muIdx] / 2.0,
                           gs: hist[muIdx] / 2.0)
    result.add FitFuncArgs(name: "O-Kalpha",
                           kind: ffGauss,
                           gN: hist[muIdx] ,
                           gmu: binning[muIdx],
                           gs: hist[muIdx] )

proc genFitFuncImpl(resultNode, idx, xNode, pFitNode: NimNode, paramsNode: seq[NimNode]): NimNode =
  ## the compilet time procedure that creates the implementation lines for
  ## the <target><Filter>Funcs that we create, which calls the correct functions,
  ## e.g.
  ##   result += p_ar[i] * gauss(x, p[i + 1], p[i + 2])
  ##   inc i, 3
  expectKind(pFitNode, nnkObjConstr)
  let fkind = parseEnum[FitFuncKind](pFitNode[2][1].strVal)
  for x in paramsNode:
    echo x.repr
  case fKind
  of ffConst:
    let p0 = paramsNode[0]
    result = quote do:
      `resultNode` += `p0`
  of ffPol1:
    let p0 = paramsNode[0]
    result = quote do:
      `resultNode` += `p0` * `xNode`
  of ffPol2:
    let p0 = paramsNode[0]
    result = quote do:
      `resultNode` += `p0` * `xNode` * `xNode`
  of ffGauss:
    let
      p0 = paramsNode[0]
      p1 = paramsNode[1]
      p2 = paramsNode[2]
    result = quote do:
      `resultNode` += `p0` * gauss(`xNode`,
                                   `p1`,
                                   `p2`)

  of ffExpGauss:
    let
      p0 = paramsNode[0]
      p1 = paramsNode[1]
      p2 = paramsNode[2]
      p3 = paramsNode[3]
      p4 = paramsNode[4]
    result = quote do:
      `resultNode` += expGauss(@[`p0`, `p1`, `p2`, `p3`, `p4`], `xNode`)


proc idOf(x: int): NimNode =
  let param = ident"p_ar"
  result = quote do:
    `param`[`x`]

proc drop(p: NimNode, frm: int): NimNode =
  result = copy(p)
  result.del(0, frm)
  echo result.treeRepr

proc incAndAdd(s: var seq[NimNode], idx: var int, cTab: var Table[string, NimNode],
               key: string) =
  if key in cTab:
    s.add cTab[key]
  else:
    cTab[key] = idOf(idx)
    s.add cTab[key]
    inc idx

proc parseParamsNodes(p, paramsNode: NimNode, idx: var int): seq[NimNode] =
  ##
  echo "XXX ", p.treeRepr
  let fkind = parseEnum[FitFuncKind](p[2][1].strVal)
  var cTab = initTable[string, NimNode]()
  var toReplace: NimNode
  if p.len > 3:
    toReplace = p.drop(3)
    for x in toReplace:
      let id = x[0].strVal
      cTab[id] = x[1]
  #var idx = 0
  case fKind
  of ffConst:
    incAndAdd(result, idx, cTab, "c")
  of ffPol1:
    incAndAdd(result, idx, cTab, "cp")
  of ffPol2:
    incAndAdd(result, idx, cTab, "cpp")
    inc idx
  of ffGauss:
    incAndAdd(result, idx, cTab, "gN")
    incAndAdd(result, idx, cTab, "gmu")
    incAndAdd(result, idx, cTab, "gs")
  of ffExpGauss:
    incAndAdd(result, idx, cTab, "ea")
    incAndAdd(result, idx, cTab, "eb")
    incAndAdd(result, idx, cTab, "eN")
    incAndAdd(result, idx, cTab, "emu")
    incAndAdd(result, idx, cTab, "es")
  echo result.repr

proc buildFitFunc(name, parts: NimNode): NimNode =
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

  var i = 0
  for p in parts:
    # add the lines for the function calls
    let parsedParams = parseParamsNodes(p, paramsNode, i)
    procBody.add genFitFuncImpl(resultNode, idx, xNode, p, parsedParams)

  # now define the result variable as a new proc
  result = newProc(name = name,
                   params = [retType, retParNode, retXNode],
                   body = procBody,
                   procType = nnkFuncDef)

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
    case s[1].len
    of 1:
      let ffName = s[1][0].strVal
      ffSeq.add quote do:
        FitFuncArgs(name: `ffName`, kind: `fKind`)
    else:
      var ffName = ""
      var ffArg = nnkObjConstr.newTree(ident"FitFuncArgs")
      for x in s[1]:
        expectKind(x, nnkAsgn)
        if x[0].basename.ident == toNimIdent"name":
          ffArg.add nnkExprColonExpr.newTree(ident"name", x[1])
        else:
          ffArg.add nnkExprColonExpr.newTree(x[0], x[1])
      ffArg.insert(2, nnkExprColonExpr.newTree(ident"kind", fKind))

      ffSeq.add ffArg

  result = buildFitFunc(funcId, ffSeq)
  echo result.repr

declareFitFunc(cuNi15):
  ffExpGauss: "Cu-esc"
    #name = "Cu-esc"
    #ea = -1.644
    #eb = 0.023
    #eN = 28.47
    #emu = 187.0
    #es = 13.5
  ffExpGauss: "Cu-Kalpha"
    #name = "Cu-Kalpha"
    #ea = -8.305
    #eb = 0.08444
    #eN = 195.2
    #emu = 286.7
    #es = 16.23
declareFitFunc(cuNi15Charge):
  ffGauss: "Cu-esc"
  ffGauss: "Cu-Kalpha"
declareFitFunc(mnCr12):
  ffExpGauss: "Mn-esc"
    #name = "Mn-esc"
    #ea = -0.06356
    #eb = 0.0467
    #emu = 102.1
    #es = 8.44
  ffExpGauss: "Mn-Kalpha"
    #name = "Mn-Kalpha"
    #ea = -1.834
    #eb = 0.04455
    #emu = 200.3
    #es = 12.7
declareFitFunc(mnCr12Charge):
  ffGauss: "Mn-esc"
  ffGauss: "Mn-Kalpha"
  ffConst: "p0"
  ffPol1: "p1"
  ffPol2: "p2"
declareFitFunc(tiTi9):
  ffGauss: "Ti-esc-alpha"
  ffGauss: "Ti-esc-beta"
    #name = "Ti-esc-beta"
    #gmu = ??
    #gs = ??
  ffExpGauss: "Ti-Kalpha"
  ffGauss: "Ti-Kbeta"
    #name = "Ti-Kbeta"
    #gmu = ??
    #gs = ??
declareFitFunc(tiTi9Charge):
  ffGauss: "Ti-esc-aplha"
  ffGauss: "Ti-esc-beta"
  ffGauss: "Ti-Kalpha"
  ffGauss: "Ti-Kbeta"
declareFitFunc(agAg6):
  ffexpGauss: "Ag-Lalpha"
  ffGauss: "Ag-Lbeta"
    #name = "Ag-Lbeta"
    #gmu = ??
    #gs = ??
declareFitFunc(agAg6Charge):
  ffGauss: "Ag-esc"
  ffGauss: "Ag-Kalpha"
  ffConst: "p0"
  ffPol1: "p1"
  ffPol2: "p2"
declareFitFunc(alAl4):
  ffexpGauss: "Al-Kalpha"
declareFitFunc(alAl4Charge):
  ffGauss: "Al-Kalpha"
  ffConst: "p0"
  ffPol1: "p1"
  ffPol2: "p2"
declareFitFunc(cuEpic2):
  ffGauss: "Cu-Lalpha"
  ffexpGauss: "Cu-Lbeta"
declareFitFunc(cuEpic2Charge):
  ffGauss: "Cu-Lalpha"
declareFitFunc(cuEpic0_9):
  ffGauss: "O-Kalpha"
  #ffGauss: "C-Kalpha"
  ffGauss: "O-Lbeta"
    #name = "C-Kalpha"
    #gmu = ??
    #gs = ??
  #ffGauss: "Fe-Lalphabeta"
    #name = "Fe-Lalphabeta"
    #gmu = ??
    #gs = ??
  #ffGauss: "Ni-Lalphabeta"
    #name = "Ni-Lalphabeta"
    #gmu = ??
    #gs = ??
declareFitFunc(cuEpic0_9Charge):
  ffGauss: "O-Kalpha"
declareFitFunc(cEpic0_6):
  ffGauss: "C-Kalpha"
  ffGauss: "O-Kalpha"
    #name = "O-Kalpha"
    #gmu = 15.0
    #gs = -12.0
declareFitFunc(cEpic0_6Charge):
  ffGauss: "C-Kalpha"
  ffGauss: "O-Kalpha"

func filterNaN(s: openArray[float]): seq[float] =
  result = newSeqOfCap[float](s.len)
  for x in s:
    if classify(x) != fcNaN:
      result.add x

proc serialize(parts: seq[FitFuncArgs]): seq[float] =
  for p in parts:
    case p.kind
    of ffConst:
      result.add filterNaN([p.c])
    of ffPol1:
      result.add filterNaN([p.cp])
    of ffPol2:
      result.add filterNaN([p.cpp])
    of ffGauss:
      result.add filterNaN([p.gN, p.gmu, p.gs])
    of ffExpGauss:
      result.add filterNaN([p.ea, p.eb, p.eN, p.emu, p.es])
    #else:
      # TODO: Finish other cases!
      #discard

macro genTfToFitFunc(pname: untyped): untyped =
  let tfkind = getType(TargetFilterKind)
  # first generate the string combinations
  var funcNames: seq[string]
  #var funcNamesC: seq[string]
  for x in tfKind:
    if x.kind != nnkEmpty:
      let xStr = ($(x.getImpl))
        .toLowerAscii
        .replace("-", "")
        .replace(".", "")
        .replace("kv", "") #& "Func"
      funcNames.add xStr# & "Func"
      #funcNamesC.add xStr & "ChargeFunc"
  # given the names, write a proc that returns the function
  let
    ##now with target and filter combined
    arg = ident"tfKind"
    argType = ident"TargetFilterKind"
    cdf = ident"CdlFitFunc"
    tfNameNode = ident"n"
    resIdent = ident"result"
  var caseStmt = nnkCaseStmt.newTree(tfNameNode)
  var caseStmtC = nnkCaseStmt.newTree(tfNameNode)
  for n in funcNames:
    let retId = ident($n & "Func")
    let retval = quote do:
      `resIdent` = `retId`
    caseStmt.add nnkOfBranch.newTree(newLit $retId, retval)
    let retIdC = ident($n & "ChargeFunc")
    let retvalC = quote do:
      `resIdent` = `retIdC`
    caseStmtC.add nnkOfBranch.newTree(newLit $retIdC, retvalC)
  let hitsname = pname
  let chargename = ident($pname & "Charge")
  result = quote do:
    #var hitsname = $pname
    #var chargename = $pname & "charge"
    proc `hitsname`(`arg`: `argType`): `cdf` =
      let `tfNameNode` = ($`arg`)
        .toLowerAscii
        .replace("-", "")
        .replace(".", "")
        .replace("kv", "") & "Func"
      `caseStmt`
    proc `chargename`(`arg`: `argType`): `cdf` =
      let `tfNameNode` = ($`arg`)
        .toLowerAscii
        .replace("-", "")
        .replace(".", "")
        .replace("kv", "") & "ChargeFunc"
      `caseStmtC`

  echo result.repr

# generate the =getCdlFitFunc= used to get the correct fit function
# based on a `TargetKind` and `FilterKind`
genTfToFitFunc(getCdlFitFunc)
#genTfToFitFunc(getCdlFitFuncCharge)


proc histoCdl(data: seq[SomeNumber], binSize: float = 3.0): (seq[float], seq[float]) =
  let low = -0.5
  var high = max(data).float + 0.5
  let nbins = (ceil((high - low) / binSize)).int
  # using correct nBins, determine actual high
  high = low + binSize * nbins.float
  let bin_edges = linspace(low, high, nbins + 1)
  let hist = data.histogram(bins = nbins, range = (low, high))

  result[0] = hist.mapIt(it.float)
  result[1] = bin_edges[1 .. ^1]

proc fitCdlImpl(hist, binedges: seq[float], tfKind: TargetFilterKind):
               (seq[float], seq[float], seq[float]) =
  let lines = getLines(hist, binedges, tfKind)
  let params = lines.serialize
  #echo params
  let testcum = cumsum(hist)
  let testsum = sum(hist)
  let testquo = testcum.mapIt(it/testsum)
  let lowIdx = testquo.lowerBound(0.0005)
  let highIdx = testquo.lowerBound(0.98)
  let passbin = binedges[lowIdx .. highIdx]
  let passhist = hist[lowIdx .. highIdx]

  let passIdx = toSeq(0 .. passhist.high).filterIt(passhist[it] > 0)
  let fitBins = passIdx.mapIt(passbin[it])
  let fitHist = passIdx.mapIt(passhist[it])
  let err = fitHist.mapIt(1.0)# / sqrt(it))

  let (pRes, res) = fit(getCdlFitFunc(tfKind),
                        params,
                        fitBins,
                        fitHist,
                        err)

  echoResult(pRes, res=res)
  result = (pRes, fitBins, fitHist)

proc fitCdlChargeImpl(hist, binedges: seq[float], tfKind: TargetFilterKind):
               (seq[float], seq[float], seq[float]) =
  let linesCharge = getLinesCharge(hist, binedges, tfKind)
  echo linesCharge, "??\n"
  let params = linesCharge.serialize
  echo params
  let testcum = cumsum(hist)
  let testsum = sum(hist)
  let testquo = testcum.mapIt(it/testsum)
  let lowIdx = testquo.lowerBound(0.0005)
  let highIdx = testquo.lowerBound(0.98)
  let passbin = binedges[lowIdx .. highIdx]
  let passhist = hist[lowIdx .. highIdx]

  let passIdx = toSeq(0 .. passhist.high).filterIt(passhist[it] > 0)
  let fitBins = passIdx.mapIt(passbin[it])
  let fitHist = passIdx.mapIt(passhist[it])
  let err = fitHist.mapIt(1.0)# / sqrt(it))

  let (pRes, res) = fit(getCdlFitFuncCharge(tfKind),
                        params,
                        fitBins,
                        fitHist,
                        err)

  echoResult(pRes, res=res)
  result = (pRes, fitBins, fitHist)


#const cuni15FuncCharge = cuNi15Func
#const mnCr12FuncCharge = mnCr12Func
#const tiTi9FuncCharge = tiTi9Func
#const agAg6FuncCharge = agAg6Func
#const alAl4FuncCharge = alAl4Func
#const cuEpic2FuncCharge = cuEpic2Func
#const cuEpic0_9FuncCharge = cuEpic0_9Func
#const cEpic0_6FuncCharge = cEpic0_6Func

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
  #echo args
  let h5file = $args["<h5file>"]

  const filename = "../../resources/cdl_runs_2019.org"
  const cutparams = "../../resources/cutparams.org"
  ##actually cutparams isn't necessary since cuts are choosen in tos helpers

  let runs = readRuns(filename)
  var h5f = H5file(h5file, "rw")
  defer: discard h5f.close()
  let cutTab = getXraySpectrumCutVals()


  proc recoAndWrite (h5file: string) =
    for r in runs:
      #if r.number != 315:
        #continue
      sleep 500
      case r.runType
      of rtXrayFinger:
        let grp = h5f[(recoDataChipBase(r.number) & "3").grp_str]
        #echo grp
        let cut = cutTab[r.toCutStr]
        let passIdx = cutOnProperties(h5f,
                                     grp,
                                     cut.cutTo,
                                     ("rmsTransverse", cut.minRms, cut.maxRms),
                                     ("length", 0.0, cut.maxLength),
                                     ("hits", cut.minPix, Inf),
                                     ("eccentricity", 0.0, cut.maxEccentricity))
        let nevents = passIdx.len

        proc writeDset(dsetWrite, dsetRead: string, datatype: typedesc) =
          var
            dset = h5f.create_dataset(grp.name / dsetWrite, (nevents, 1),
                                      datatype)
          if dsetWrite == "CdlSpectrumIndices":
            dset[dset.all] = passIdx
          else:
            let read = h5f[grp.name / dsetRead, datatype]
            dset[dset.all] = passIdx.mapIt(read[it])
          dset.attrs["Target"] = $r.target
          dset.attrs["Filter"] = $r.filter
          dset.attrs["HV"] = $r.hv
        writeDset("CdlSpectrumIndices", "", int64)
        writeDset("CdlSpectrum", "hits", int64)
        writeDset("CdlSpectrumEvents", "eventNumber", int64)
        writeDset("CdlSpectrumCharge", "totalCharge", float64)

        let runnum = h5f[(recoBase() & $r.number).grp_str]
        #echo runnum
        runnum.attrs["tfKind"] = $r.toCutStr
      else:
        discard

  recoAndWrite(h5file)

  proc fitAndPlot(h5file: string, tfKind: TargetFilterKind) =
    var h5file = h5file
    let targetFilter = tfkind
    var rawseq: seq[int64]
    var cutseq: seq[int64]
    var rawChargeseq: seq[float64]
    var cutChargeseq: seq[float64]



    for r in runs:
      #if r.number != 347:
        #continue
      sleep 500
      case r.runType
      of rtXrayFinger:
        let grp = h5f[(recoDataChipBase(r.number) & "3").grp_str]
        let tfk = r.totfkind
        #echo tfk
        if tfk == targetFilter:
            let hitsRawDataseq = h5f[grp.name / "hits", int64]
            let hitsCdlseq = h5f[grp.name / "CdlSpectrum", int64]
            let chargeRawDataseq = h5f[grp.name / "totalCharge", float64]
            let chargeCdlseq = h5f[grp.name / "CdlSpectrumCharge", float64]
            #echo chargeCdlseq

            rawseq.add(hitsRawDataseq)
            cutseq.add(hitsCdlseq)
            rawChargeseq.add(chargeRawDataseq)
            cutChargeseq.add(chargeCdlseq)
        #echo rawseq
      else:
         discard


    #let hitsRawData = h5f[grp.name / "hits", int64]
    #let hitsCdl = h5f[grp.name / "CdlSpectrum", int64]
    let (histdata, bins) = histoCdl(cutseq, binSize = 1.0)
    let (pRes, fitBins, fitHist) = fitCdlImpl(histdata, bins, targetFilter)
    let (histdataC, binsC) = histoCdl(cutChargeseq, binSize = 10000.0)
    let (pResC, fitBinsC, fitHistC) = fitCdlChargeImpl(histdataC, binsC, targetFilter)

    let fitres = fitBins.mapIt(getCdlFitFunc(targetFilter)(pRes, it))
    let cdlplot = scatterPlot(fitBins, fitres).mode(PlotMode.Lines)
    let fitresC = fitBins.mapIt(getCdlFitFuncCharge(targetFilter)(pResC, it))
    let cdlplotC = scatterPlot(fitBinsC, fitresC).mode(PlotMode.Lines)


    let ff = getCdlFitFunc(targetFilter)
    let ffC = getCdlFitFuncCharge(targetFilter)

    #let scipy = pyImport("scipy.optimize")
     #proc convertff(fitObject: PyObject, p0, p1, p2, p3, p4,p5, p6, p7, p8, p9: float): seq[float] =
           #for x in fitObject:
          #let xNim = x.to(float)
          #result.add ff(@[p0, p1, p2, p3, p4, p5, p6, p7, p8, p9], xNim)

    fitForNlopt(convertNlopt, ff)
    var opt = newNloptOpt("LN_BOBYQA", pRes.len)
    var fitObj = FitObject(x: fitBins, y: fitHist, yErr: fitHist.mapIt(sqrt(it)))
    var vstruct = newVarStruct(convertNlopt, fitObj)
    opt.setFunction(vstruct)

    opt.xtol_rel = 1e-10
    opt.ftol_rel = 1e-10
    opt.maxtime  = 5.0
    let (paramsN, minN) = opt.optimize(pRes)
    #echo opt.status
    #echo paramsN
    nlopt_destroy(opt.optimizer)


    fitForNlopt(convertNloptC, ffC)
    var optC = newNloptOpt("LN_BOBYQA", pResC.len)
    var fitObjC = FitObject(x: fitBinsC, y: fitHistC, yErr: fitHistC.mapIt(sqrt(it)))
    var vstructC = newVarStruct(convertNloptC, fitObjC)
    optC.setFunction(vstructC)

    optC.xtol_rel = 1e-10
    optC.ftol_rel = 1e-10
    optC.maxtime  = 5.0
    let lChargeStart = getLinesCharge(histDataC, binsC, targetFilter).serialize
    let (paramsNC, minNC) = optC.optimize(lChargeStart)
    #echo opt.status
    #echo paramsN
    nlopt_destroy(optC.optimizer)

    #let res = scipy.curve_fit(convertff, fitBins, fitHist, p0 = pRes,
    #                          sigma = fitHist.mapIt(sqrt(it)),
    #                          `method` = "trf", diff_step = 1e-5)
    #let popt = res[0].mapIt(it.to(float))
    #echo popt
    #let pcov = res[1]
    #let fitResPy = fitBins.mapIt(getCdlFitFunc(tfk)(popt, it))
    let fitResNlopt = fitBins.mapIt(getCdlFitFunc(targetFilter)(paramsN, it))
    #let cdlplotPy = scatterPlot(fitBins, fitResPy).mode(PlotMode.Lines)
    let cdlplotNlopt = scatterPlot(fitBins, fitResNlopt).mode(PlotMode.Lines)

    let fitResNloptC = fitBins.mapIt(getCdlFitFuncCharge(targetFilter)(paramsNC, it))
    let cdlplotNloptC = scatterPlot(fitBinsC, fitResNloptC).mode(PlotMode.Lines)

    ##plot of hits
    let hitsRaw = histPlot(rawseq.mapIt(it.float64))
      .binSize(1.0)
      .binRange(0.0, 400.0)
    let hitsCut = histPlot(cutseq.mapIt(it.float64))
      .binSize(1.0)
      .binRange(0.0, 400.0)
    hitsRaw.layout.barMode = BarMode.Overlay
    let plt = hitsRaw.addTrace(hitsCut.traces[0])
      #.addTrace(cdlplot.traces[0])
      #.addTrace(cdlPlotPy.traces[0])
      .addTrace(cdlPlotNlopt.traces[0])
    plt.layout.title = &"target: {targetFilter}"
    plt.layout.showlegend = true
    plt.traces[1].opacity = 0.5
    plt.traces[0].name = "raw data"
    plt.traces[0].marker = Marker[float](color: @[Color1])
    plt.traces[1].name = "data with cuts"
    plt.traces[1].marker = Marker[float](color: @[Color2])
    plt.traces[1].opacity = 1.0
    #plt.traces[1].name = "fit curve"
    #plt.traces[3].name = "fit curve py"
    plt.traces[2].name = "fit curve nlopt"
    plt.traces[2].marker = Marker[float](color: @[black])
    plt.layout.yaxis.title = "Occurence"
    plt.layout.xaxis.title = "Number of pixels"
    plt.show(&"{targetFilter}-2019.svg")

    ##plot of charge
    let ChargeRaw = histPlot(rawChargeseq.mapIt(it.float64))
      .binSize(10000.0)
      .binRange(0.0, 3500000.0)
    let ChargeCut = histPlot(cutChargeseq.mapIt(it.float64))
      .binSize(10000.0)
      .binRange(0.0, 3500000.0)
    ChargeRaw.layout.barMode = BarMode.Overlay
    let pltC = ChargeRaw.addTrace(ChargeCut.traces[0])
      .addTrace(cdlplotC.traces[0])
      .addTrace(cdlPlotNloptC.traces[0])
    pltC.layout.title = &"target: {targetFilter}"
    pltC.layout.showlegend = true
    pltC.traces[1].opacity = 0.5
    pltC.traces[0].name = "raw data"
    pltC.traces[0].marker = Marker[float](color: @[Color1])
    pltC.traces[1].name = "data with cuts"
    pltC.traces[1].marker = Marker[float](color: @[Color2])
    pltC.traces[1].opacity = 1.0
    pltC.traces[2].name = "fit curve"
    pltC.traces[2].marker = Marker[float](color: @[black])
    pltC.layout.yaxis.title = "Occurence"
    pltC.layout.xaxis.title = "Charge"
    pltC.show(&"{targetFilter}Charge-2019.svg")



  for tfkind in TargetFilterKind:
    fitAndPlot(h5file, tfkind)


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
