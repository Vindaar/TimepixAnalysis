import parsecsv, os, streams, strutils, strformat, nimhdf5, tables, sequtils
import plotly
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
    tfCuNi = "Cu-Ni"
    tfMnCr = "Mn-Cr"
    tfTiTi = "Ti-Ti"
    tfAgAg = "Ag-Ag"
    tfAlAl = "Al-Al"
    tfCuEpic = "Cu-Epic"
    tfCEpic =  "C-Epic"

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

  cutsB = object
    target: TargetKind
    filter: FilterKind
    hv: float
    ck: int
    length: float
    rms_min: float
    rms_max: float
    eccentricity: float

  cutsC = object
    target: TargetKind
    filter: FilterKind
    hv: float
    ck: int
    charge_min: float
    charge_max: float
    length: float
    rms_min: float
    rms_max: float

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
    of ffExpGauss:
      ea: float
      eb: float
      eN: float
      emu: float
      es: float
    of ffGauss:
      gN: float
      gmu: float
      gs: float

func getLines(hist, binning: seq[float], tfKind: TargetFilterKind): seq[FitFuncArgs] =
  let muIdx = argmax(hist)
  case tfKind
  of tfCuNi: #need to add for rest combinatons
     discard
  of tfMnCr:
     discard
  of tfTiTi:
     discard
  of tfAgAg:
     discard
     ## result.add FitFuncArgs(name: "Ag-Lalpha",
     ##                        kind: ffExpGauss)
     ## result.add FitFuncArgs(name: "Ag-Lbeta",
     ##                        kind: ffGauss,
     ##                        gmu: binning[muIdx],
     ##                        gN: hist[muIdx],
     ##                        gs: hist[muIdx] / 10.0)
  of tfAlAl:
    discard
  of tfCuEpic:
     discard
  of tfCEpic:
    #var ffConst = 0.0
    result.add FitFuncArgs(name: "C-Kalpha",
                           kind: ffGauss,
                           gmu: binning[muIdx],
                           gN: hist[muIdx],
                           gs: hist[muIdx] / 10.0)
    result.add FitFuncArgs(name: "O-Kalpha",
                           kind: ffGauss,
                           gmu: binning[muIdx],
                           gN: hist[muIdx],
                           gs: hist[muIdx] / 10.0)
  #else:
    #discard


proc cEpicFunc(p_ar: seq[float], x: float): float =
  var p: seq[float]
  result = p[0] * gauss(x, p[1], p[2]) + p[3] * gauss(x, p[4], p[5])

template buildFitProc(name: untyped, parts: seq[FitFuncArgs]): untyped =
  var params: seq[float]
  var funcBody: float
  for p in parts:
    case p.kind
    of ffGauss:
      funcBody += p.gN * gauss(x, p.gmu, p.gs)
    of ffExpGauss:
      funcBody += expGauss(@[p.ea, p.eb, p.eN, p.emu, p.es], x)
  proc `name`(p_ar: seq[float], x: float): float =
    result = funcBody
  `name`


const cEpicFuncCharge = cEpicFunc
# TODO: impl rest of functions + Charge functions

proc getCdlFits(): Table[string, CdlFit] =
  result = {"C-EPIC-0.6kV" : cEpicFunc}.toTable
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
      echo run.toCutStr

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
    if r.number != 342:
      continue
    case r.runType
    of rtXrayFinger:
      let grp = h5f[(recoDataChipBase(r.number) & "3").grp_str]
      let cut = cutTab[r.toCutStr]
      echo cut
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

      let hitsRaw = histPlot(h5f[grp.name / "hits", int64])
        .binSize(2.0)
        .binRange(0.0, 400.0)
      let hitsCut = histPlot(h5f[grp.name / "CdlSpectrum", int64])
        .binSize(2.0)
        .binRange(0.0, 400.0)
      hitsRaw.layout.barMode = BarMode.Overlay
      let plt = hitsRaw.addTrace(hitsCut.traces[0])
      plt.traces[1].opacity = 0.5
      plt.show()

    else:
      discard


when isMainModule:
  main()
