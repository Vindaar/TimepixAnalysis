import parsecsv, os, streams, strutils, strformat, nimhdf5, tables, sequtils, macros
import seqmath, algorithm, times, strscans, typeinfo
import plotly, mpfit, nlopt, nimpy
import ingrid / [ingrid_types, tos_helpers, calibration]
import docopt
import helpers / utils
import chroma
import nimsvg, ospaths, random, xmlparser, xmltree

const docStr = """
Usage:
  cdl_spectrum_creation <h5file> [options]
  cdl_spectrum_creation -h | --help
  cdl_spectrum_creation --version
  cdl_spectrum_creation <h5file> --dumpAccurate [options]
  cdl_spectrum_creation <h5file> --cutcdl [options]
  cdl_spectrum_creation <h5file> --genRefFile --year=YEAR [options]
  cdl_spectrum_creation <h5file> --genCdlFile --year=YEAR [options]

Options:
  -h, --help      Show this help
  --version       Show the version number
  --cutcdl        Creates CDL data in h5
  --dumpAccurate  If set will dump the fit parameters to a
                  `fitParameters_<timestamp>.txt` file with
                  higher accuracy (4 decimal places instead of 2).
  --genRefFile    Generates the X-ray reference data file. Basically
                  TargetFilterKinds filtered by charge cut on peaks.
  --genCdlFile    Generate the combined CDL calibration file.
                  Mainly input file regrouped by TargetFilterKind
                  instead of run numbers.
  --outfile=NAME  Name of the output file. Optional.
  --year=YEAR     Year to add to output filenames.
"""
const doc = withDocopt(docStr)

##some constants depending on the run
const filename = "../../resources/cdl_runs_2019.org"
#const cutparams = "../../resources/cutparams.org"
#actually cutparams isn't necessary since cuts are choosen in tos helpers
const outdate = &"2019"
const chipnumber = "3"

## some different color definitions
const Color1 = color(1.0, 0.0, 102.0 / 256.0)
const Color2 = color(0.0, 153.0 / 256.0, 204 / 256.0)
const black = color(0.0, 0.0, 0.0)
const ColorTGelb = parseHex("FFC107")
const ColorTBlau = parseHex("0288D1")
const ColorTDBlau = parseHex("0288D1")
const ColorTBGrau = parseHex("BDBDBD")
const ColorTGrau = parseHex("757575")
const ColorTHGrau = color(0.92, 0.92, 0.92)


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

  TargetFilterKind* = enum
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

  DataKind = enum
    Dhits = "hits"
    Dcharge = "charge"

  YearKind = enum
    yr2014 = "2014"
    yr2018 = "2018"

const fixed = NaN

func getLines(hist, binning: seq[float], tfKind: TargetFilterKind): seq[FitFuncArgs] =
  ## this is a runtime generator for the correct fitting function prototype,
  ## i.e. it returns a seq of parts, which need to be combined to the complete
  ## function at runtime
  let muIdx = argmax(hist)
  case tfKind
  of tfCuNi15:
    result.add FitFuncArgs(name: "Cu-Kalpha",
                          kind: ffExpGauss,
                          ea: hist[muIdx] * 1e-10,
                          eb: hist[muIdx] * 1e-12,
                          eN: hist[muIdx] / 30.0,
                          emu: 250.0, #binning[muIdx], #200.0
                          es: hist[muIdx] / 30.0) #13.0
    result.add FitFuncArgs(name: "Cu-esc",
                          kind: ffExpGauss,
                          ea: hist[muIdx] *  1e-10,
                          eb: hist[muIdx] * 1e-12,
                          eN: hist[muIdx] / 10.0, #100.0
                          emu: 150.0, #binning[muIdx],# / 2.0,  #100.0
                          es: hist[muIdx]) # / 2.0) #/ 15.0) #16.0
  of tfMnCr12:
    result.add FitFuncArgs(name: "Mn-Kalpha",
                          kind: ffExpGauss,
                          ea: hist[muIdx] * 1e-10,
                          eb: hist[muIdx] * 1e-12,
                          eN: hist[muIdx],# / 10.0, #30.0
                          emu: 197.0, #binning[muIdx], #200.0
                          es: hist[muIdx] / 15.0) # 13.0
    result.add FitFuncArgs(name: "Mn-esc",
                          kind: ffExpGauss,
                          ea: hist[muIdx]  * 1e-10,
                          eb: hist[muIdx]  * 1e-12,
                          eN: hist[muIdx] / 10.0, #120
                          emu: 95.0, #binning[muIdx], #100.0
                          es: hist[muIdx] / 15.0) #16.0
  of tfTiTi9:
    result.add FitFuncArgs(name: "Ti-Kalpha",
                          kind: ffExpGauss,
                          ea: hist[muIdx] * 1e-10,
                          eb: hist[muIdx] * 1e-12,
                          eN: hist[muIdx],# / 10.0,
                          emu: binning[muIdx],
                          es: hist[muIdx] / 30.0)
    result.add FitFuncArgs(name: "Ti-esc-alpha",
                          kind: ffGauss,
                          gN: hist[muIdx] / 20.0,
                          gmu: binning[muIdx] / 2.0,
                          gs: hist[muIdx] / 15.0)
    result.add FitFuncArgs(name: "Ti-esc-beta",
                          kind: ffGauss,
                          gN: hist[muIdx],# / 20.0,
                          gmu: fixed, #binning[muIdx],
                          gs: fixed) #hist[muIdx] / 15.0)
    result.add FitFuncArgs(name: "Ti-Kbeta",
                          kind: ffGauss,
                          gN: hist[muIdx],# / 10.0,
                          gmu: fixed, #binning[muIdx],
                          gs: fixed) #hist[muIdx] / 15.0)
  of tfAgAg6:
    result.add FitFuncArgs(name: "Ag-Lalpha",
                          kind: ffExpGauss,
                          ea: hist[muIdx] * 1e-10,
                          eb: hist[muIdx] * 1e-12,
                          eN: hist[muIdx], # / 10.0,
                          emu: binning[muIdx],
                          es: hist[muIdx] / 30.0) ##30.0
    result.add FitFuncArgs(name: "Ag-Lbeta",
                          kind: ffGauss,
                          gN: fixed, #binning[muIdx] / 3.0,
                          gmu: fixed, #hist[muIdx],
                          gs: fixed) #hist[muIdx] / 10.0)
  of tfAlAl4:
    result.add FitFuncArgs(name: "Al-Kalpha",
                          kind: ffExpGauss,
                          ea: hist[muIdx] * 1e-10,
                          eb: hist[muIdx] * 1e-12,
                          eN: hist[muIdx],# / 10.0,
                          emu: binning[muIdx],
                          es: hist[muIdx] / 80.0) ##30.0
  of tfCuEpic2:
    result.add FitFuncArgs(name: "Cu-Lalpha",
                          kind: ffGauss,
                          gN: hist[muIdx],
                          gmu: binning[muIdx],# / 4.0,
                          gs: hist[muIdx])# / 10.0)
    result.add FitFuncArgs(name: "Cu-Lbeta",
                          kind: ffGauss,
                          gN: hist[muIdx] / 10.0,
                          gmu: binning[muIdx],
                          gs: hist[muIdx]  / 40.0)
  of tfCuEpic0_9:
    result.add FitFuncArgs(name: "O-Kalpha",
                          kind: ffGauss,
                          gN: hist[muIdx],# / 1.5,
                          gmu: binning[muIdx],# * 2.0,
                          gs: hist[muIdx] / 100.0)
    #result.add FitFuncArgs(name: "C-Kalpha",
    #                      kind: ffGauss,
    #                      gN: hist[muIdx] / 100.0,
    #                      gmu: fixed, #binning[muIdx], # / 2.0,
    #                      gs: fixed) #hist[muIdx] / 30.0)
    #result.add FitFuncArgs(name: "Fe-Lalphabeta",
    #                      kind: ffGauss,
    #                      gN: hist[muIdx] / 10.0,
    #                      gmu: fixed, #binning[muIdx],
    #                      gs: fixed) #hist[muIdx] )
    #result.add FitFuncArgs(name: "Ni-Lalphabeta",
    #                      kind: ffGauss,
    #                      gN: hist[muIdx] / 10.0,
    #                      gmu: fixed, #binning[muIdx],
    #                      gs: fixed) #hist[muIdx] )
  of tfCEpic0_6:
    result.add FitFuncArgs(name: "C-Kalpha",
                          kind: ffGauss,
                          gN: hist[muIdx],# / 100.0,# / 4.0,
                          gmu: binning[muIdx], # * 2.0,
                          gs: hist[muIdx])# *  2.0)
    result.add FitFuncArgs(name: "O-Kalpha",
                          kind: ffGauss,
                          gN: hist[muIdx] / 100.0,# / 10.0,
                          gmu: fixed, #binning[muIdx],
                          gs: fixed) #hist[muIdx])

func getbounds(tfKind:TargetFilterKind): seq[tuple[l, u:float]] =
  case tfKind
  of tfCuNi15:
    result = @[(l: -Inf, u:Inf),
               (l: -Inf, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: -Inf, u:Inf),
               (l: -Inf, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf)]
  of tfMnCr12:
    result = @[(l: -Inf, u:Inf),
               (l: -Inf, u:Inf),
               (l: 80.0, u:Inf),
               (l: 190.0, u:Inf),
               (l: 1.5, u:20.0),
               (l: -Inf, u:Inf),
               (l: -Inf, u:Inf),
               (l: 5.0, u:Inf),
               (l: 80.0, u:110.0),
               (l: 1.0, u:15.0)]
  of tfTiTi9:
    result = @[(l: -Inf, u:Inf),
               (l: -Inf, u:Inf),
               (l: 260.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:100.0),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf)]
  of tfAgAg6:
    result = @[(l: -Inf, u:Inf),
               (l: -Inf, u:Inf),
               (l: 310.0, u:380.0),
               (l: 50.0, u:116.0),
               (l: 10.0, u:13.5)]
  of tfAlAl4:
    result = @[(l: -Inf, u:Inf),
               (l: -Inf, u:Inf),
               (l: 480.0, u:535.0),
               (l: 10.0, u:57.5),
               (l: 5.0, u:7.5)]
  of tfCuEpic2:
    result = @[(l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf)]
  of tfCuEpic0_9:
    result = @[(l: 1.0, u:Inf),
               (l: 1.0, u:23.0),
               (l: 1.0, u:7.0)]
               #(l: 1.0, u:60.0)]
               #(l: 5.0, u:60.0),
               #(l: 5.0, u:110.0)]
  of tfCEpic0_6:
    result = @[(l: 600.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:100.0)]


func getLinesCharge(hist, binning: seq[float], tfKind: TargetFilterKind): seq[FitFuncArgs] =
  ## this is a runtime generator for the correct fitting function prototype,
  ## i.e. it returns a seq of parts, which need to be combined to the complete
  ## function at runtime
  let muIdx = argmax(hist)
  case tfKind
  of tfCuNi15:
    result.add FitFuncArgs(name: "Cu-Kalpha",
                          kind: ffGauss,
                          gN: hist[muIdx],
                          gmu: 1700.0e3, #binning[muIdx],
                          gs: hist[muIdx] * 1.5e3)
    result.add FitFuncArgs(name: "Cu-esc",
                          kind: ffGauss,
                          gN: hist[muIdx] / 10.0,
                          gmu: 1070.0e3, #binning[muIdx] / 2.0,
                          gs: hist[muIdx] * 1e3)
  of tfMnCr12:
    result.add FitFuncArgs(name: "Mn-Kalpha",
                          kind: ffGauss,
                          gN: hist[muIdx] / 4.0,
                          gmu: binning[muIdx],
                          gs: hist[muIdx] * 2.0e3)
    result.add FitFuncArgs(name: "Mn-esc",
                          kind: ffGauss,
                          gN: hist[muIdx] / 10.0,
                          gmu: binning[muIdx] / 2.0,
                          gs: hist[muIdx] * 1.0e3)
    #result.add FitFuncArgs(name: "p0",
    #                      kind: ffConst,
    #                      c: -hist[muIdx])
    #result.add FitFuncArgs(name: "p1",
    #                      kind: ffPol1,
    #                      cp: hist[muIdx])
    #result.add FitFuncArgs(name: "p2",
    #                      kind: ffPol2,
    #                      cpp: -hist[muIdx])
  of tfTiTi9:
    result.add FitFuncArgs(name: "Ti-Kalpha",
                          kind: ffGauss,
                          gN: hist[muIdx],# / 10.0,
                          gmu: binning[muIdx] , #1140.0e3,
                          gs: hist[muIdx] * 4e2)
    result.add FitFuncArgs(name: "Ti-esc-alpha",
                          kind: ffGauss,
                          gN: hist[muIdx] / 20.0,
                          gmu: binning[muIdx] / 3.0, #400.0e3,
                          gs: hist[muIdx] * 1e2)
    result.add FitFuncArgs(name: "Ti-esc-beta",
                          kind: ffGauss,
                          gN: hist[muIdx] / 30.0,
                          gmu: fixed, #binning[muIdx] ,
                          gs: fixed) #hist[muIdx] * 1e3)
    result.add FitFuncArgs(name: "Ti-Kbeta",
                          kind: ffGauss,
                          gN: hist[muIdx] / 30.0,
                          gmu: fixed, #binning[muIdx], #* 1e5,
                          gs: fixed) #hist[muIdx] * 1e3)
  of tfAgAg6:
    result.add FitFuncArgs(name: "Ag-Lalpha",
                          kind: ffGauss,
                          gN: hist[muIdx],# / 2.0,
                          gmu: binning[muIdx],# * 1e3,
                          gs: hist[muIdx] * 1e2)# * 10)
    result.add FitFuncArgs(name: "Ag-Lbeta",
                          kind: ffGauss,
                          gN: fixed, #hist[muIdx] / 10.0,
                          gmu: fixed, #binning[muIdx],
                          gs: fixed) #hist[muIdx] * 1.5e3)
    #result.add FitFuncArgs(name: "p0",
    #                      kind: ffConst,
    #                      c: hist[muIdx] * 1e-3)
    #result.add FitFuncArgs(name: "p1",
    #                      kind: ffPol1,
    #                      cp: hist[muIdx] * 1e-3)
    #result.add FitFuncArgs(name: "p2",
    #                      kind: ffPol2,
    #                      cpp: hist[muIdx] * 1e-3)
  of tfAlAl4:
    result.add FitFuncArgs(name: "Ag-Kalpha",
                          kind: ffGauss,
                          gN: hist[muIdx],# / 2.0,
                          gmu: binning[muIdx],# 350.0e3
                          gs: hist[muIdx] * 1e2)# * 10)
    #result.add FitFuncArgs(name: "p0",
    #                      kind: ffConst,
    #                      c: -hist[muIdx])
    #result.add FitFuncArgs(name: "p1",
    #                      kind: ffPol1,
    #                      cp: hist[muIdx])
    #result.add FitFuncArgs(name: "p2",
    #                      kind: ffPol2,
    #                      cpp: -hist[muIdx])
  of tfCuEpic2:
    result.add FitFuncArgs(name: "Cu-Lalpha",
                          kind: ffGauss,
                          gN: hist[muIdx],# / 20.0,
                          gmu: binning[muIdx], #* 0.5e3,
                          gs: hist[muIdx] * 2e2)
    result.add FitFuncArgs(name: "Cu-Lbeta",
                          kind: ffGauss,
                          gN: hist[muIdx],# / 50.0,
                          gmu: binning[muIdx],
                          gs: hist[muIdx] * 0.5e2)
  of tfCuEpic0_9:
    result.add FitFuncArgs(name: "O-Kalpha",
                          kind: ffGauss,
                          gN: hist[muIdx],# / 4.0,
                          gmu: binning[muIdx],
                          gs: hist[muIdx] * 1e2)
    result.add FitFuncArgs(name: "C-Kalpha",
                          kind: ffGauss,
                          gN: hist[muIdx],# / 4.0,
                          gmu: binning[muIdx],
                          gs: hist[muIdx] * 1e3)
  of tfCEpic0_6:
    result.add FitFuncArgs(name: "C-Kalpha",
                          kind: ffGauss,
                          gN: hist[muIdx],# / 4.0,
                          gmu: binning[muIdx],
                          gs: hist[muIdx]  * 1e2)
    result.add FitFuncArgs(name: "O-Kalpha",
                          kind: ffGauss,
                          gN: hist[muIdx],# / 10.0,
                          gmu: fixed, #binning[muIdx], #* 1e3,
                          gs: fixed) #hist[muIdx]  * 1e3)

func getboundsCharge(tfKind:TargetFilterKind): seq[tuple[l, u:float]] =
  case tfKind
  of tfCuNi15:
    result = @[(l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf)]
  of tfMnCr12:
    result = @[(l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf)]
  of tfTiTi9:
    result = @[(l: 260.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:100.0),
               (l: 1.0, u:450.0e3),
               (l: 1.0, u:Inf),
               (l: 5.0, u:Inf),
               (l: 5.0, u:Inf)]
  of tfAgAg6:
    result = @[(l: 390.0, u:450.0),
               (l: 1.0, u:Inf),
               (l: 1.0, u:110.0e3)]
  of tfAlAl4:
    result = @[(l: 630.0, u:Inf),
               (l: 1.0, u:360.0e3),
               (l: 1.0, u:70.0e3)]
  of tfCuEpic2:
    result = @[(l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf)]
  of tfCuEpic0_9:
    result = @[(l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf)]
  of tfCEpic0_6:
    result = @[(l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf)]

proc genFitFuncImpl(idx: var int, xNode,
                    pFitNode: NimNode, paramsNode: seq[NimNode]): NimNode =
  ## the compile time procedure that creates the implementation lines for
  ## the <target><Filter>Funcs that we create, which calls the correct functions,
  ## e.g.
  ##   result += p_ar[i] * gauss(x, p[i + 1], p[i + 2])
  ##   inc i, 3
  let resultNode = ident"result"
  expectKind(pFitNode, nnkObjConstr)
  let fkind = parseEnum[FitFuncKind](pFitNode[2][1].strVal)
  case fKind
  of ffConst:
    let p0 = paramsNode[idx]
    result = quote do:
      `resultNode` += `p0`
    inc idx
  of ffPol1:
    let p0 = paramsNode[idx]
    result = quote do:
      `resultNode` += `p0` * `xNode`
    inc idx
  of ffPol2:
    let p0 = paramsNode[idx]
    result = quote do:
      `resultNode` += `p0` * `xNode` * `xNode`
    inc idx
  of ffGauss:
    let
      p0 = paramsNode[idx]
      p1 = paramsNode[idx + 1]
      p2 = paramsNode[idx + 2]
    result = quote do:
      `resultNode` += `p0` * gauss(`xNode`, `p1`, `p2`)
    inc idx, 3
  of ffExpGauss:
    let
      p0 = paramsNode[idx]
      p1 = paramsNode[idx + 1]
      p2 = paramsNode[idx + 2]
      p3 = paramsNode[idx + 3]
      p4 = paramsNode[idx + 4]
    result = quote do:
      `resultNode` += expGauss(@[`p0`, `p1`, `p2`, `p3`, `p4`], `xNode`)
    inc idx, 5

proc idOf(x: int): NimNode =
  let param = ident"p_ar"
  result = quote do:
    `param`[`x`]

proc drop(p: NimNode, frm: int): NimNode =
  result = copy(p)
  result.del(0, frm)
  echo result.treeRepr


proc resolveCall(tab: OrderedTable[string, NimNode],
                 n: NimNode): NimNode =
  expectKind(n, nnkCall)
  const paramIdent = "p_ar"
  if n[1].strVal == paramIdent:
    # call is an nnkOpenSymChoice from the `p_ar[x]`
    result = n
  else:
    # call is a reference to a different parameter
    let key = n[0].strVal & "_" & n[1].strVal
    result = tab[key]

proc resolveNodes(tab: OrderedTable[string, NimNode]): seq[NimNode] =
  ## resolve the potentially still ambiguous nodes (e.g. a reference to a
  ## parameter of a different line) and convert the `OrderedTable` into a
  ## serialized sequence of the nodes
  for k, v in tab:
    case v.kind
    of nnkCall:
      result.add resolveCall(tab, v)
    else:
      if v.len == 0:
        result.add v
      else:
        # have to walk the tree to replace references
        var node = copyNimTree(v)
        for i in 0 ..< node.len:
          case node[i].kind
          of nnkCall:
            # replace the call with the resolved node
            node[i] = resolveCall(tab, node[i])
          else:
            # else leave node unchanged
            discard
        result.add node

proc incAndFill(tab: var OrderedTable[string, NimNode],
                idx: var int, keyBase, lineName: string) =
  let key = keyBase & "_" & lineName
  if key notin tab:
    tab[key] = idOf(idx)
    inc idx
  else:
    # if key in the table, remove it and re-add so that the position
    # is the correct one in our OrderedTable!
    let tmp = tab[key]
    tab.del(key)
    tab[key] = tmp

proc fillParamsTable(tab: var OrderedTable[string, NimNode],
                     p: NimNode, idx: var int) =
  ## creates a table of NimNodes, which stores the NimNode that will show up
  ## in the created fitting procedure under a unique string, which consists
  ## of `FitParameter_LineName`. This table can then later be used to
  ## reference arbitrary parameters of the fitting functions in a different
  ## or same line
  # get the correct part of the parameters
  let fkind = parseEnum[FitFuncKind](p[2][1].strVal)
  #var cTab = initOrderedTable[string, NimNode]()
  var toReplace: NimNode
  let lineName = p[1][1].strVal
  if p.len > 3:
    # if p.len has more than 3 children, one ore more parameters have
    # been fixed to some constant or relative value. Fill the `tab` with
    # the NimNodes corresponding to those values
    toReplace = p.drop(3)
    for x in toReplace:
      case x[1].kind
      of nnkIdent, nnkIntLit .. nnkFloatLit, nnkInfix, nnkCall:
        let id = x[0].strVal & "_" & lineName
        tab[id] = x[1]
      else: error("Unsupported kind " & $x[1].kind)

  # Now walk over the table with the parameters of the current line
  # and either enter the correct `p_ar[idx]` field as the value in
  # the table or put the value to be replaced at the correct position
  # in the table
  case fKind
  of ffConst:
    incAndFill(tab, idx, "c", lineName)
  of ffPol1:
    incAndFill(tab, idx, "cp", lineName)
  of ffPol2:
    incAndFill(tab, idx, "cpp", lineName)
    inc idx
  of ffGauss:
    incAndFill(tab, idx, "gN", lineName)
    incAndFill(tab, idx, "gmu", lineName)
    incAndFill(tab, idx, "gs", lineName)
  of ffExpGauss:
    incAndFill(tab, idx, "ea", lineName)
    incAndFill(tab, idx, "eb", lineName)
    incAndFill(tab, idx, "eN", lineName)
    incAndFill(tab, idx, "emu", lineName)
    incAndFill(tab, idx, "es", lineName)

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

  var paramsTab = initOrderedTable[string, NimNode]()
  for p in parts:
    # create the table which maps the parameter of each line to the
    # correct NimNode. May either be
    # - `p_ar[idx]`
    # - some literal int / float
    # - some calculation involving a reference to some other parameter
    paramsTab.fillParamsTable(p, i)
  # now that we have the whole tableresolve all still remaining references
  # to other parameters
  let params = resolveNodes(paramsTab)
  i = 0
  for p in parts:
    procBody.add genFitFuncImpl(i, xNode, p, params)

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
  ffExpGauss: "Cu-Kalpha"
  ffExpGauss: "Cu-esc"
declareFitFunc(cuNi15Charge):
  ffGauss: "Cu-Kalpha"
  ffGauss: "Cu-esc"
declareFitFunc(mnCr12):
  ffExpGauss: "Mn-Kalpha"
  ffExpGauss: "Mn-esc"
declareFitFunc(mnCr12Charge):
  ffGauss: "Mn-Kalpha"
  ffGauss: "Mn-esc"
declareFitFunc(tiTi9):
  ffExpGauss: "Ti-Kalpha"
  ffGauss: "Ti-esc-alpha"
    #name = "Ti-esc-alpha"
    #gmu = emu("Ti-Kalpha") * (1.537/4.511)
  ffGauss: #"Ti-esc-beta"
    name = "Ti-esc-beta"
    gmu = emu("Ti-Kalpha") * (1.959/4.511)
    gs = gs("Ti-esc-alpha")
  ffGauss: #"Ti-Kbeta"
    name = "Ti-Kbeta"
    gmu = emu("Ti-Kalpha") * (4.932/4.511)
    gs = es("Ti-Kalpha")
declareFitFunc(tiTi9Charge):
  ffGauss: "Ti-Kalpha"
  ffGauss: "Ti-esc-alpha"
    #name = "Ti-esc-alpha"
    #gmu = emu("Ti-Kalpha") * (1.537/4.511)
  ffGauss: #"Ti-esc-beta"
    name = "Ti-esc-beta"
    gmu = gmu("Ti-Kalpha") * (1.959/4.511)
    gs = gs("Ti-esc-alpha")
  ffGauss: #"Ti-Kbeta"
    name = "Ti-Kbeta"
    gmu = gmu("Ti-Kalpha") * (4.932/4.511)
    gs = gs("Ti-Kalpha")
declareFitFunc(agAg6):
  ffExpGauss: "Ag-Lalpha"
  ffGauss: #"Ag-Lbeta"
    name = "Ag-Lbeta"
    gN = eN("Ag-Lalpha") * 0.1
    gmu = emu("Ag-Lalpha") * (3.151/2.984)
    gs = es("Ag-Lalpha")
declareFitFunc(agAg6Charge):
  ffGauss: "Ag-Lalpha"
  ffGauss: #"Ag-Lbeta"
    name = "Ag-Lbeta"
    gN = gN("Ag-Lalpha") * 0.1
    gmu = gmu("Ag-Lalpha") * (3.151/2.984)
    gs = gs("Ag-Lalpha")
  #ffConst: "p0"
  #ffPol1: "p1"
  #ffPol2: "p2"
declareFitFunc(alAl4):
  ffExpGauss: "Al-Kalpha"
declareFitFunc(alAl4Charge):
  ffGauss: "Al-Kalpha"
  #ffConst: "p0"
  #ffPol1: "p1"
  #ffPol2: "p2"
declareFitFunc(cuEpic2):
  ffGauss: "Cu-Lalpha"
  ffGauss: "Cu-Lbeta"
declareFitFunc(cuEpic2Charge):
  ffGauss: "Cu-Lalpha"
  ffGauss: "Cu-Lbeta"
declareFitFunc(cuEpic0_9):
  ffGauss: "O-Kalpha"
  #ffGauss: #"C-Kalpha"
  #  name = "C-Kalpha"
  #  gmu = gmu("O-Kalpha") * (0.277/0.525)
  #  gs = gs("O-Kalpha")
  #ffGauss: #"Fe-Lalphabeta"
  #  name = "Fe-Lalphabeta"
  #  gmu = gmu("O-Kalpha") * (0.71/0.525)
  #  gs = gs("O-Kalpha")
  #ffGauss: #"Ni-Lalphabeta"
  #  name = "Ni-Lalphabeta"
  #  gmu = gmu("O-Kalpha") * (0.86/0.525)
  #  gs = gs("O-Kalpha")
declareFitFunc(cuEpic0_9Charge):
  ffGauss: "O-Kalpha"
  ffGauss: "C-Kalpha"
declareFitFunc(cEpic0_6):
  ffGauss: "C-Kalpha"
  ffGauss: #"O-Kalpha"
    name = "O-Kalpha"
    gmu = gmu("C-Kalpha") * (0.525/0.277)
    gs = gs("C-Kalpha")
declareFitFunc(cEpic0_6Charge):
  ffGauss: "C-Kalpha"
  ffGauss: #"O-Kalpha"
    name = "O-Kalpha"
    gmu = gmu("C-Kalpha") * (0.525/0.277)
    gs = gs("C-Kalpha")

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


macro genTfToFitFunc(pname: untyped): untyped =
  let tfkind = getType(TargetFilterKind)
  #first generate the string combinations
  var funcNames: seq[string]
  for x in tfKind:
    if x.kind != nnkEmpty:
      let xStr = ($(x.getImpl))
        .toLowerAscii
        .replace("-", "")
        .replace(".", "")
        .replace("kv", "")
      funcNames.add xStr
  #given the names, write a proc that returns the function
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


proc histoCdl(data: seq[SomeNumber], binSize: float = 3.0, dKind: DataKind): (seq[float], seq[float]) =
  let low = -0.5 * binSize
  var high = max(data).float + (0.5 * binSize)
  let nbins = (ceil((high - low) / binSize)).round.int
  # using correct nBins, determine actual high
  high = low + binSize * nbins.float
  let (hist, bin_edges) = data.histogram(bins = nbins, range = (low, high))
  #let bin_edges = linspace(low, high, nbins + 1)
  #let hist = data.histogram(bins = nbins, range = (low, high))

  result[0] = hist.mapIt(it.float)
  case dKind
  of Dhits:
    result[1] = bin_edges[1 .. ^1]
  of Dcharge:
    result[1] = bin_edges[0 .. ^2]
  #echo "Bin edges len ", bin_edges.len
  #echo "Result len ", result[1].len


proc fitCdlImpl(hist, binedges: seq[float], tfKind: TargetFilterKind, dKind: DataKind):
               (seq[float], seq[float], seq[float], seq[float], seq[float]) =
  ##generates the fit paramter
  var lines: seq[FitFuncArgs]
  var fitfunc: CdlFitFunc
  case dKind
  of Dhits:
    lines = getLines(hist, binedges, tfKind)
    fitfunc = getCdlFitFunc(tfKind)
  of Dcharge:
    lines = getLinesCharge(hist, binedges, tfKind)
    fitfunc = getCdlFitFuncCharge(tfKind)

  let params = lines.serialize
  let cumu = cumsum(hist)
  let sum = sum(hist)
  let quotient = cumu.mapIt(it/sum)
  let lowIdx = quotient.lowerBound(0.0005)
  let highIdx = quotient.lowerBound(0.98)
  let passbin = binedges[lowIdx .. highIdx]
  let passhist = hist[lowIdx .. highIdx]

  let passIdx = toSeq(0 .. passhist.high).filterIt(passhist[it] > 0)
  let fitBins = passIdx.mapIt(passbin[it])
  let fitHist = passIdx.mapIt(passhist[it])
  let err = fitHist.mapIt(1.0)# / sqrt(it))

  let (pRes, res) = fit(fitfunc,
                        params,
                        fitBins,
                        fitHist,
                        err)

  echoResult(pRes, res=res)
  echo "error futcdlimpl ", res.error
  result = (params, pRes, fitBins, fitHist, res.error)

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

proc totfkind(run: CdlRun): TargetFilterKind =
  result = parseEnum[TargetFilterKind](&"{toCutStr(run)}")

iterator tfRuns(h5f: var H5FileObj, tfKind: TargetFilterKind): H5Group =
  ## Yields the center chip group of all runs from `filename`,
  ## which match `tfKind`
  let runs = readRuns(filename)
  for r in runs:
    case r.runType
    of rtXrayFinger:
      let tfk = r.totfkind
      if tfk == tfKind:
        let runGrp = h5f[recoRunGrpStr(r.number)]
        let centerChip = runGrp.attrs["centerChip", int]
        let chpGrp = h5f[(runGrp.name / "chip_" & $centerChip).grp_str]
        yield chpGrp
    else:
      # nothing to yield if not an "XrayFinger" (read CDL) run
      discard

proc cutAndWrite(h5file: string) =
  let runs = readRuns(filename)
  var h5f = H5file(h5file, "rw")
  defer: discard h5f.close()
  let cutTab = getXraySpectrumCutVals()
  for r in runs:
    #if r.number != 315:
      #continue
    sleep 500
    case r.runType
    of rtXrayFinger:
      let grp = h5f[(recoDataChipBase(r.number) & chipnumber).grp_str]
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
      runnum.attrs["tfKind"] = $r.toCutStr
    else:
      discard

proc energycurve(energyResHits: seq[float], energyResCharge: seq[float], errorHits: seq[float], errorCharge: seq[float]) =
  #var energyval = @[15.0, 12.0, 9.0, 6.0, 4.0, 2.0, 0.9, 0.6]
  var energyval = @[8.048, 5.899, 4.511, 2.984, 1.487, 0.930, 0.525, 0.277]
  let energyplotHits = scatterPlot(energyval, energyResHits)
  let energyplotCharge = scatterPlot(energyval, energyResCharge)
  let errorbarHits = newErrorbar(errorHits)
  let errorbarCharge = newErrorbar(errorCharge)
  let energyplot = energyplotHits.addTrace(energyplotCharge.traces[0])
    .legendLocation(x = 0.8, y = 0.9)
    #.legendBgColor(ColorTHGrau)
    .backgroundColor(ColorTHGrau)
    .gridColor(color())
  energyplot.layout.showlegend = true
  energyplot.layout.title = "Energy resolution plot"
  energyplot.layout.xaxis.title = "Energy in keV"
  energyplot.layout.yaxis.title = "Energy resolution"
  energyplot.traces[0].name = "number of hits"
  energyplot.traces[0].marker = Marker[float](color: @[ColorTDBlau])#, size: @[8.0])
  energyplot.traces[0].ys_err = errorbarHits
  energyplot.traces[1].name = "total charge"
  energyplot.traces[1].marker = Marker[float](color: @[ColorTGelb])#, size: @[8.0])
  energyplot.traces[1].ys_err = errorbarCharge
  energyplot.show(&"energyresoplot-{outdate}.svg")

proc peakfit(peakpos: seq[float], name: string, error: seq[float]) =
  var energyval = @[8.048, 5.899, 4.511, 2.984, 1.487, 0.930, 0.525, 0.277]
  let peakfit = scatterPlot(energyval, peakpos)
    .legendLocation(x = 0.1, y = 0.9)
    #.legendBgColor(ColorTHGrau)
    .backgroundColor(ColorTHGrau)
    .gridColor(color())
  let errorbar = newErrorbar(error)
  peakfit.layout.showlegend = true
  peakfit.layout.title = name
  peakfit.traces[0].name = &"{name}"
  peakfit.traces[0].marker = Marker[float](color: @[ColorTDBlau])
  peakfit.traces[0].ys_err = errorbar
  peakfit.layout.xaxis.title = "Energy in keV"
  peakfit.layout.yaxis.title = &"Peakposition for {name}"
  peakfit.show(&"{name}.svg")

proc dumpFitParameters(outfile, svgFname: string,
                       params: seq[float], errors: seq[float],
                       tfKind: TargetFilterKind, dKind: DataKind) =
  ## dumps the fit paramters and their names, plus the filename of the corresponding SVG
  ## to a txt file
  var outf = open(outfile, fmAppend)
  outf.write(&"svg: {svgFname}\n")
  outf.write(&"tfKind: {tfKind}\n")
  # now get the correct names for the fit parameters via a call to getLines
  var fitLines: seq[FitFuncArgs]
  case dKind
  of DHits:
    fitLines = getLines(@[0'f64], @[0'f64], # dummy values
                        tfKind)
  of DCharge:
    fitLines = getLinesCharge(@[0'f64], @[0'f64], # dummy values
                              tfKind)
  let maxLineWidth = fitLines.mapIt(it.name.len).max

  # determine max length of params and errors fields as a string
  func toString(f: float): string =
    if f >= 1e5:
      if "accurate" in outfile:
        result = &"{f:.2e}"
      else:
        result = &"{f:.1e}"
    else:
      if "accurate" in outfile:
        result = &"{f:.4f}"
      else:
        result = &"{f:.2f}"
  func maxWidth(s: seq[float]): int =
    let strs = s.mapIt(it.toString)
    result = max(strs.mapIt(it.len))
  func paddedStrs(s: seq[float]): seq[string] =
    let maxW = s.maxWidth
    result = s.mapIt(it.toString)
    for i in 0 ..< result.len:
      result[i] = result[i].align(maxW)

  let pStrs = params.paddedStrs
  let eStrs = errors.paddedStrs
  # iterate the lines and then unroll the FitFuncArgs object
  var i = 0
  var lineName: string
  for line in fitLines:
    for field, val in fieldPairs(line):
      when type(val) is string:
        lineName = val
        outf.write(&"{lineName}:\n")
      elif type(val) is float:
        let fstr = $field
        if classify(val) != fcNaN:
          outf.write(&"{fstr:<3} = {pStrs[i]} \\pm {eStrs[i]}\n")
          inc i
  outf.close()

proc calcfitcurve(minbin: float, maxbin: float,
                  cdlFitFunc: CdlFitFunc,
                  fitparams: seq[float]): (seq[float], seq[float]) =
  let
    minvalue = minbin
    maxvalue = maxbin
    range = linspace(minvalue.float, maxvalue.float, 1500)
    yvals = range.mapIt(cdlFitFunc(fitparams, it))
  result = (range, yvals)

proc fitAndPlot[T: SomeNumber](h5file, fitParamsFname: string,
                               tfKind: TargetFilterKind, dKind: DataKind):
               (seq[float], seq[float], seq[float], seq[float]) =
  var h5f = H5file(h5file, "rw")
  defer: discard h5f.close()

  let runs = readRuns(filename)
  var ploterror: seq[float]
  var rawseq: seq[T]
  var cutseq: seq[T]
  var binsizeplot: float
  var binrangeplot: float
  var xtitle: string
  var outname: string
  var fitfunc: CdlFitFunc
  var lines: seq[FitFuncArgs]
  var dummy = @[1.0, 2.0]
  case dKind
  of Dhits:
    fitfunc = getCdlFitFunc(tfKind)
    binsizeplot = 1.0
    binrangeplot = 400.0
    xtitle = "Number of pixles"
    outname = &"{tfKind}"
    lines = getLines(dummy, dummy, tfKind)
  of Dcharge:
    fitfunc = getCdlFitFuncCharge(tfKind)
    binsizeplot = 10000.0
    binrangeplot = 3500000.0
    xtitle = "Charge"
    outname = &"{tfKind}Charge"
    lines = getLinesCharge(dummy, dummy, tfKind)

  for grp in tfRuns(h5f, tfKind):
    case dKind
    of Dhits:
      let RawDataSeq = h5f[grp.name / "hits", T]
      let Cdlseq = h5f[grp.name / "CdlSpectrum", T]
      rawseq.add(RawDataseq)
      cutseq.add(Cdlseq)
    of Dcharge:
      let RawDataSeq = h5f[grp.name / "totalCharge", T]
      let Cdlseq = h5f[grp.name / "CdlSpectrumCharge", T]
      rawseq.add(RawDataseq)
      cutseq.add(Cdlseq)

  proc calcfit(dataseq: seq[SomeNumber],
               cdlFitFunc: CdlFitFunc,
               binSize: float,
               tfKind: TargetFilterKind,
               dKind: DataKind): (seq[float], seq[float], float, float) =
    let (histdata, bins) = histoCdl(dataseq, binSize, dKind)
    let (pStart, pRes, fitBins, fitHist, errorres) = fitCdlImpl(histdata, bins, tfKind, dKind)
    #fitForNlopt(convertNlopt, cdlFitFunc)
    #fitForNloptLnLikelihood(convertNlopt, cdlFitFunc)
    fitForNloptLnLikelihoodGrad(convertNlopt, cdlFitFunc)
    var bounds: seq[tuple[l, u:float]]
    case dKind
    of Dhits:
      bounds = getbounds(tfKind)
    of Dcharge:
      bounds = getboundsCharge(tfKind)
    doassert pStart.len == bounds.len
    #echo "P start len ", pStart.len
    #echo "Bounds len ", bounds.len
    #var opt = newNloptOpt(LD_TNEWTON_PRECOND, pStart.len, bounds)
    var fitObj = FitObject(x: fitBins, y: fitHist) #, yErr: fitHist.mapIt(sqrt(it)))
    var vstruct = newVarStruct(convertNlopt, fitObj)

    var opt = newNloptOpt[type(fitObj)](LD_MMA, pStart.len, bounds)
    opt.setFunction(vstruct)
    opt.xtol_rel = 1e-10
    opt.ftol_rel = 1e-10
    #opt.xtol_abs = 1e-14
    #opt.ftol_abs = 1e-14
    opt.maxtime  = 20.0
    #opt.maxEval  = 2000
    #opt.initialStep = 1
    let (paramsN, minN) = opt.optimize(pStart)
    nlopt_destroy(opt.optimizer)

    let minbin = fitBins.min
    let maxbin = fitBins.max
    result = (pRes, paramsN, minbin, maxbin)
    #echo "fit start params ", pStart
    #echo "fit pRes ", pRes
    #echo "fit paramsN ", paramsN

  let (histdata, bins) = histoCdl(cutseq, binsizeplot, dKind)
  let (pStart, pRes, fitBins, fitHist, errorsres) = fitCdlImpl(histdata, bins, tfKind, dKind)
  let fitresults = calcfit(cutseq, fitfunc, binsizeplot, tfKind, dKind)
  echo "fitresults nlopt", fitresults[1]
  ##get the interesting fit params

  ##some sketchy work around for errors
  let cumu = cumsum(histdata)
  let sum = sum(histdata)
  let quotient = cumu.mapIt(it/sum)
  let lowIdx = quotient.lowerBound(0.0005)
  let highIdx = quotient.lowerBound(0.98)
  let passbin = bins[lowIdx .. highIdx]
  let passhist = histdata[lowIdx .. highIdx]
  let passIdx = toSeq(0 .. passhist.high).filterIt(passhist[it] > 0)
  let fitBinserr = passIdx.mapIt(passbin[it])
  let fitHisterr = passIdx.mapIt(passhist[it])
  let err = fitHist.mapIt(1.0)# / sqrt(it))

  let (FitError, mpfitError) = fit(fitfunc, fitresults[1], fitBinserr, fitHisterr, err)
  echo "fit errors: ", mpfitError.error
  ploterror = mpfitError.error


  var fitmu: float
  var fitsig: float
  var fitmuerr: float
  var fitsigerr: float
  var energyres: float
  var museq: seq[float]
  var sigseq: seq[float]
  var energyseq: seq[float]
  var energyreserr: float
  var energyseqerr: seq[float]
  #echo "lines ", lines[0].kind
  case lines[0].kind
  of ffGauss:
    fitmu = fitresults[1][1]
    fitsig = fitresults[1][2]
    fitmuerr = ploterror[1]
    fitsigerr = ploterror[2]
    museq.add(fitmu)
  of ffExpGauss:
    fitmu = fitresults[1][3]
    fitsig = fitresults[1][4]
    fitmuerr = ploterror[3]
    fitsigerr = ploterror[4]
    museq.add(fitmu)
  else:
    discard
  #echo "fitmu ", fitmu
  #echo "fitsig ", fitsig
  #echo "fitmuseq ", museq
  energyres = fitsig / fitmu
  energyreserr = sqrt(pow((fitsigerr / fitmu), 2) + pow((fitsig * fitmuerr / pow(fitmu,2)),2) )
  energyseqerr.add(energyreserr)
  energyseq.add(energyres)
  echo "energyres ", energyres
  result[0].add(fitmu)
  result[1].add(energyres)
  result[2].add(fitmuerr)
  result[3].add(energyreserr)

  let mpfitres = calcfitcurve(fitresults[2], fitresults[3], fitfunc, fitresults[0])
  let nloptres = calcfitcurve(fitresults[2], fitresults[3], fitfunc, fitresults[1])
  let startval = calcfitcurve(fitresults[2], fitresults[3], fitfunc, pStart)
  let cdlPlot = scatterPlot(mpfitres[0], mpfitres[1]).mode(PlotMode.Lines)
  let cdlPlotNlopt = scatterPlot(nloptres[0], nloptres[1]).mode(PlotMode.Lines)
  let startPlot = scatterPlot(startval[0], startval[1]).mode(PlotMode.Lines)

  ##test of annotations
  #echo "testforparams ", testmu
  #echo test8
  #let teststring = test8.string
  #let testanno = Annotation(x: 350,
  #                          y: 100,
  #                          text: &"test: " & test8)

  ##plot of hits and charge
  let hitsRaw = histPlot(rawseq.mapIt(it.float64))
    .binSize(binsizeplot)
    .binRange(0.0, binrangeplot)
  let hitsCut = histPlot(cutseq.mapIt(it.float64))
    .binSize(binsizeplot)
    .binRange(0.0, binrangeplot)
  hitsRaw.layout.barMode = BarMode.Overlay
  let plt = hitsRaw.addTrace(hitsCut.traces[0])
    #.addTrace(cdlPlot.traces[0])
    .addTrace(cdlPlotNlopt.traces[0])
    #.addTrace(startPlot.traces[0])
    .legendLocation(x = 0.8, y = 0.9)
    #.legendBgColor(ColorTHGrau)
    .backgroundColor(ColorTHGrau)
    .gridColor(color())
  plt.layout.title = &"target: {tfKind}"
  plt.layout.showlegend = true
  #plt.legendBgColor(ColorTB)
  plt.traces[0].opacity = 1.0
  plt.traces[0].name = "raw data"
  plt.traces[0].marker = Marker[float](color: @[ColorTGrau])
  plt.traces[1].name = "data with cuts"
  plt.traces[1].marker = Marker[float](color: @[ColorTDBlau])
  plt.traces[1].opacity = 1.0
  plt.traces[2].name = "fit curve nlopt"
  plt.traces[2].marker = Marker[float](color: @[ColorTGelb])
  #plt.traces[3].name = "fit start"
  #plt.traces[3].marker = Marker[float](color: @[black])
  plt.layout.yaxis.title = "Occurence"
  plt.layout.xaxis.title = xtitle
  #plt.layout.annotations.add [testanno]

  let fname = &"{outname}-{outdate}.svg"
  plt.show(fname)

  # now dump the fit results, SVG filename and correct parameter names to a file
  dumpFitParameters(fitParamsFname, fname, fitresults[1], ploterror, tfKind, dKind)

proc cdlToXrayTransform(h5fout: var H5FileObj,
                        passedData: seq[float],
                        tfKindStr: string,
                        outname: string,
                        year: YearKind) =
  ## performs whole conversion from from calibration-cdl data to Xray cut
  ## if the charge cut has been applied (`passIdx` indicates passing indices).
  # means we have a 1:1 map, so get the required binning
  var
    numBins: int
    minVal: float
    maxVal: float
  case year
  of yr2014:
    (numBins, minVal, maxVal) = cdlToXrayBinning2014(outname)
  of yr2018:
    (numBins, minVal, maxVal) = cdlToXrayBinning2018(outname) # not implemented yet
  # given passing data, calculate histogram and write to file
  let (hist, bins) = histogram(passedData,
                               numBins,
                               range = (minVal, maxVal),
                               upperRangeBinRight = false)
  # combine the hist bins data to a seq2D
  let histBins = @[bins[0 .. ^2], hist.mapIt(it.float)].transpose
  # create dataset
  let dsetToWrite = h5fout.create_dataset((tfKindStr / outname),
                                          histBins.shape,
                                          float)
  dsetToWrite[dsetToWrite.all] = histBins

proc readAndFilter(h5f: var H5FileObj,
                   dsetName: string,
                   passIdx: seq[int]): seq[float] =
  ## reads the dataset given by `dsetName` as a `float` seq and filters
  ## by `passIdx`.
  let data = h5f.readAs(dsetName, float)
  # apply passIdx
  result = passIdx.mapIt(data[it])

proc generateCdlCalibrationFile(h5file: string, year: YearKind,
                                outfile = "calibration-cdl") =
  ## generates the CDL calibration data file from a HDF5 file containing
  ## all CDL runs. Supports either 2014 CDL data or 2019 CDL data.
  # walk all runs corresponding to a single `TargetFilterKind` and
  # combine the datasets into the output files
  var h5f = H5file(h5file, "r")
  let runs = readRuns(filename)
  for tfKind in TargetFilterKind:
    for grp in tfRuns(h5f, tfKind):
      # will not iterate the datasets and write to the outfile
      # via hyperslabs, potentially appending to the existing
      # dataset
      discard
  discard h5f.close()

proc generateXrayReferenceFile(h5file: string, year: YearKind,
                               outfile = "XrayReferenceFile") =
  ## generates the X-ray reference data file
  # this is achieved by taking the raw TargetFilterKind runs, combining them
  # into a the full CDL calibration file (if the input is the raw run based
  # file, else it skips the first step)
  # then we apply the charge cuts and bin the data by N bins
  # the result is written to the new file as (N, 2) datasets
  var h5f = H5file(h5file, "r")
  if "reconstruction" in h5f:
    # is a raw file, first create the CDL calibration file
    discard h5f.close()
    let cdlOut = "auto_calibration-cdl_" & $year
    generateCdlCalibrationFile(h5file, year, cdlOut)
    h5f = H5file(cdlOut, "r")

  # now walk all groups in root of h5f, read the datasets required for
  # charge cuts, write all passing indices back to file as binned
  # datasets
  var
    date: string
    tfKindStr: string
  const xrayRefTab = getXrayRefTable()
  let xrayRefCuts = getEnergyBinMinMaxVals() #XraySpectrumCutVals()

  var h5fout = H5file(outfile & $year & ".h5", "rw")

  for group in h5f:
    var mgrp = group
    echo group.name
    if scanf(group.name, "/calibration-cdl-$+-$+kV", date, tfKindStr):
      case date
      of "apr2014":
        doAssert year == yr2014
        # now perform cuts on datasets
        doAssert tfKindStr & "kV" in xrayRefCuts
        let cut = xrayRefCuts[tfKindStr & "kV"]
        let passIdx = cutOnProperties(h5f,
                                      group,
                                      crSilver, # try cutting to silver
                                      ("RmsTransverse", cut.minRms, cut.maxRms),
                                      ("Length", 0.0, cut.maxLength),
                                      ("NumberOfPixels", cut.minPix, Inf),
                                      ("TotalCharge", cut.minCharge, cut.maxCharge))
        echo "Number of passing indices ", passIdx.len
        #let pix = h5f.readAs(group.name / "NumberOfPixels", float)
        #let pixReduced = passIdx.mapIt(pix[it])
        #let dfRaw = seqsToDf({"pix" : pix})
        #let dfReduced = seqsToDf({"pix" : pixReduced})
        #ggplot(dfRaw, aes("pix")) + geom_histogram(bins = 300) + ggsave("test_" & $tfKindStr & ".pdf")
        #ggplot(dfReduced, aes("pix")) + geom_histogram(bins = 300) + ggsave("test_cut_" & $tfKindStr & ".pdf")
        # given passIdx, now read each dataset iteratively and apply cuts
        for dset in mgrp:
          case dset.dtypeAnyKind
          of akSequence:
            # variable length data, x, y, charge, will be dropped in conversion
            discard
          else:
            # get the name of the dataset in the output
            let outname = cdlToXray2014(dset.name.extractFilename)
            if outname.len > 0:
              # float32 data
              let passedData = readAndFilter(h5f, dset.name, passIdx)
              cdlToXrayTransform(h5fout, passedData, tfKindStr & "kV", outname, year)
        # 2014 dataset does not contain `lengthdivbyrmsy` in the calibration file
        # create that manually
        let length = readAndFilter(h5f, mgrp.name / "Length", passIdx)
        let rmsTrans = readAndFilter(h5f, mgrp.name / "RmsTransverse", passIdx)
        let lengthDivRmsTrans = toSeq(0 ..< length.len).mapIt(
          length[it] / rmsTrans[it]
        )
        # also write this dataset
        cdlToXrayTransform(h5fout,
                           lengthDivRmsTrans,
                           tfKindStr & "kV",
                           "lengthdivbyrmsy",
                           year)
      of "feb2019":
        doAssert year == yr2018
        # now perform cuts on datasets
        doAssert tfKindStr & "kV" in xrayRefCuts
        let cut = xrayRefCuts[tfKindStr & "kV"]
        let passIdx = cutOnProperties(h5f,
                                      group,
                                      ("rmsTransverse", cut.minRms, cut.maxRms),
                                      ("length", 0.0, cut.maxLength),
                                      ("hits", cut.minPix, Inf),
                                      ("totalCharge", cut.minCharge, cut.maxCharge))
        echo "Number of passing indices ", passIdx.len
        # given passIdx, now read each dataset iteratively and apply cuts
        for dset in mgrp:
          case dset.dtypeAnyKind
          of akSequence:
            # variable length data, x, y, charge, will be dropped in conversion
            discard
          else:
            # output name as as input for 2018
            let outname = dset.name.extractFilename
            # float64 data
            let passedData = readAndFilter(h5f, dset.name, passIdx)
            cdlToXrayTransform(h5fout, passedData, tfKindStr & "kV", outname, year)
      else:
        raise newException(ValueError, "Invalid year string for calibration: " &
          $date)
    else:
      raise newException(ValueError, "Could not match calibration group in " &
        $h5file & ", group name was " & $group.name)
  discard h5fout.close()
  discard h5f.close()

proc main =
  #echo "OK"
  let args = docopt(doc)
  echo "ARGS", args

  let h5file = $args["<h5file>"]
  let year = parseEnum[YearKind]($args["--year"])
  let reco_order = args["--cutcdl"].toBool
  let genRefFile = args["--genRefFile"].toBool
  let genCdlFile = args["--genCdlFile"].toBool
  let dumpAccurate = if args["--dumpAccurate"].toBool:
                       true
                     else:
                       false

  if genRefFile:
    generateXrayReferenceFile(h5file, year)
  if genCdlFIle:
    generateCdlCalibrationFile(h5file, year)
  if reco_order:
    cutAndWrite(h5file)

  if not genRefFile and not genCdlFile:
    # only perform CDL fits if neither CDL calibration file nor
    # reference file created
    var fitParamsFname = ""
    if dumpAccurate:
      fitParamsFname = "fitparams_accurate_" & $(epochTime().round.int) & ".txt"
    else:
      fitParamsFname = "fitparams_" & $(epochTime().round.int) & ".txt"


    var peakposHits: seq[float]
    var peakHitsErr: seq[float]
    var energyResHits: seq[float]
    var energyHitsErr: seq[float]
    var peakposCharge: seq[float]
    var peakChargeErr: seq[float]
    var energyResCharge: seq[float]
    var energyChargeErr: seq[float]
    #let a = fitAndPlot[int64](h5file, tfCuEpic0_9, Dhits)
    #let b = fitAndPlot[float64](h5file, tfCuEpic0_9, Dcharge)
    for tfkind in TargetFilterKind:
      let energyHits = fitAndPlot[int64](h5file, fitParamsFname, tfkind, Dhits)
      peakposHits.add(energyHits[0])
      energyResHits.add(energyHits[1])
      peakHitsErr.add(energyHits[2])
      energyHitsErr.add(energyHits[3])
      #echo "energyres ", energyResHit
      let energyCharge = fitAndPlot[float64](h5file, fitParamsFname, tfkind, Dcharge)
      peakposCharge.add(energyCharge[0])
      energyResCharge.add(energyCharge[1])
      peakChargeErr.add(energyCharge[2])
      energyChargeErr.add(energyCharge[3])
    energycurve(energyResHits, energyResCharge, energyHitsErr, energyChargeErr)
    peakfit(peakposHits, "Hits", peakHitsErr)
    peakfit(peakposCharge, "Charge", peakChargeErr)

when isMainModule:
  main()
