import parsecsv, os, streams, strutils, strformat, tables, sequtils, macros, fenv
import seqmath, algorithm, times, strscans, typeinfo
import mpfit, nlopt, nimhdf5
import ingrid / [ingrid_types, tos_helpers]
import ingrid / calibration
import ingrid / calibration / [fit_functions, calib_fitting]
import cdlFitting / cdlFitMacro
import docopt
import helpers / utils
import chroma

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
import projectDefs
const filename = TpxDir / "resources/cdl_runs_2014.org"
#const cutparams = "../../resources/cutparams.org"
#actually cutparams isn't necessary since cuts are choosen in tos helpers
const outdate = &"2019"
const chipnumber = "3"
# this global can be changed to create the XrayReferenceFile.h5 using the
# 2014 naming scheme (default TPA naming)
const XrayRefGenerateNamingScheme = fkTpa
# this global can be changed to create the calibration-cdl.h5 using the
# 2014 naming scheme (default TPA naming)
const CdlGenerateNamingScheme = fkTpa


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

  DataKind = enum
    Dhits = "hits"
    Dcharge = "charge"

  YearKind = enum
    yr2014 = "2014"
    yr2018 = "2018"

func getLines(hist, binning: seq[float], tfKind: TargetFilterKind): seq[FitFuncArgs] =
  ## this is a runtime generator for the correct fitting function prototype,
  ## i.e. it returns a seq of parts, which need to be combined to the complete
  ## function at runtime
  let (mu_main, sigma_main, n_main, _, _, _) = getLines(hist, binning)
  case tfKind
  of tfCuNi15:
    result.add FitFuncArgs(name: "Cu-Kalpha",
                          kind: ffExpGauss,
                          ea: n_main * 1e-10,
                          eb: n_main * 1e-12,
                          eN: n_main / 1.5,
                          emu: mu_main, #200.0
                          es: sigma_main) #13.0
    result.add FitFuncArgs(name: "Cu-esc",
                          kind: ffExpGauss,
                          ea: n_main *  1e-10,
                          eb: n_main * 1e-12,
                          eN: n_main / 10.0, #100.0
                          emu: mu_main / 1.5,  #100.0
                          es: sigma_main) # / 2.0) #/ 15.0) #16.0
  of tfMnCr12:
    result.add FitFuncArgs(name: "Mn-Kalpha",
                          kind: ffExpGauss,
                          ea: n_main * 1e-10,
                          eb: n_main * 1e-12,
                          eN: n_main,# / 10.0, #30.0
                          emu: mu_main, #200.0
                          es: sigma_main) # 13.0
    result.add FitFuncArgs(name: "Mn-esc",
                          kind: ffExpGauss,
                          ea: n_main  * 1e-10,
                          eb: n_main  * 1e-12,
                          eN: n_main / 10.0, #120
                          emu: mu_main / 2.0, #100.0
                          es: sigma_main) #16.0
  of tfTiTi9:
    result.add FitFuncArgs(name: "Ti-Kalpha",
                          kind: ffExpGauss,
                          ea: n_main * 1e-10,
                          eb: n_main * 1e-12,
                          eN: n_main / 2.0,
                          emu: mu_main,
                          es: sigma_main)
    result.add FitFuncArgs(name: "Ti-esc-alpha",
                          kind: ffGauss,
                          gN: n_main / 20.0,
                          gmu: mu_main / 2.0,
                          gs: sigma_main)
    result.add FitFuncArgs(name: "Ti-esc-beta",
                          kind: ffGauss,
                          gN: n_main / 20.0,# / 20.0,
                          gmu: fixed, #mu_main,
                          gs: fixed) #n_main / 15.0)
    result.add FitFuncArgs(name: "Ti-Kbeta",
                          kind: ffGauss,
                          gN: n_main / 2.0,# / 10.0,
                          gmu: fixed, #mu_main,
                          gs: fixed) #n_main / 15.0)
  of tfAgAg6:
    result.add FitFuncArgs(name: "Ag-Lalpha",
                          kind: ffExpGauss,
                          ea: n_main * 1e-10,
                          eb: n_main * 1e-12,
                          eN: n_main, # / 10.0,
                          emu: mu_main,
                          es: sigma_main) ##30.0
    result.add FitFuncArgs(name: "Ag-Lbeta",
                          kind: ffGauss,
                          gN: fixed, #mu_main / 3.0,
                          gmu: fixed, #n_main,
                          gs: fixed) #n_main / 10.0)
  of tfAlAl4:
    result.add FitFuncArgs(name: "Al-Kalpha",
                          kind: ffExpGauss,
                          ea: n_main * 1e-10,
                          eb: n_main * 1e-12,
                          eN: n_main,# / 10.0,
                          emu: mu_main,
                          es: sigma_main) ##30.0
  of tfCuEpic2:
    result.add FitFuncArgs(name: "Cu-Lalpha",
                          kind: ffGauss,
                          gN: n_main / 1.5,
                          gmu: mu_main,# / 4.0,
                          gs: sigma_main)# / 10.0)
    result.add FitFuncArgs(name: "Cu-Lbeta",
                          kind: ffGauss,
                          gN: n_main / 10.0,
                          gmu: mu_main,
                          gs: sigma_main)
  of tfCuEpic0_9:
    result.add FitFuncArgs(name: "O-Kalpha",
                          kind: ffGauss,
                          gN: n_main,# / 1.5,
                          gmu: mu_main,# * 2.0,
                          gs: sigma_main)
    #result.add FitFuncArgs(name: "C-Kalpha",
    #                      kind: ffGauss,
    #                      gN: n_main / 100.0,
    #                      gmu: fixed, #mu_main, # / 2.0,
    #                      gs: fixed) #n_main / 30.0)
    #result.add FitFuncArgs(name: "Fe-Lalphabeta",
    #                      kind: ffGauss,
    #                      gN: n_main / 10.0,
    #                      gmu: fixed, #mu_main,
    #                      gs: fixed) #n_main )
    #result.add FitFuncArgs(name: "Ni-Lalphabeta",
    #                      kind: ffGauss,
    #                      gN: n_main / 10.0,
    #                      gmu: fixed, #mu_main,
    #                      gs: fixed) #n_main )
  of tfCEpic0_6:
    result.add FitFuncArgs(name: "C-Kalpha",
                          kind: ffGauss,
                          gN: n_main,# / 100.0,# / 4.0,
                          gmu: mu_main, # * 2.0,
                          gs: sigma_main)# *  2.0)
    result.add FitFuncArgs(name: "O-Kalpha",
                          kind: ffGauss,
                          gN: n_main / 100.0,# / 10.0,
                          gmu: fixed, #mu_main,
                          gs: fixed) #n_main)

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
  let (mu_main, sigma_main, n_main, _, _, _) = getLines(hist, binning)
  case tfKind
  of tfCuNi15:
    result.add FitFuncArgs(name: "Cu-Kalpha",
                          kind: ffGauss,
                          gN: n_main,
                          gmu: mu_main,
                          gs: sigma_main)
    result.add FitFuncArgs(name: "Cu-esc",
                          kind: ffGauss,
                          gN: n_main / 5.0,
                          gmu: mu_main / 1.5,
                          gs: sigma_main)
  of tfMnCr12:
    result.add FitFuncArgs(name: "Mn-Kalpha",
                          kind: ffGauss,
                          gN: n_main,
                          gmu: mu_main,
                          gs: sigma_main)
    result.add FitFuncArgs(name: "Mn-esc",
                          kind: ffGauss,
                          gN: n_main / 10.0,
                          gmu: mu_main / 2.0,
                          gs: sigma_main)
    #result.add FitFuncArgs(name: "p0",
    #                      kind: ffConst,
    #                      c: -n_main)
    #result.add FitFuncArgs(name: "p1",
    #                      kind: ffPol1,
    #                      cp: n_main)
    #result.add FitFuncArgs(name: "p2",
    #                      kind: ffPol2,
    #                      cpp: -n_main)
  of tfTiTi9:
    result.add FitFuncArgs(name: "Ti-Kalpha",
                          kind: ffGauss,
                          gN: n_main,# / 10.0,
                          gmu: mu_main , #1140.0e3,
                          gs: sigma_main)
    result.add FitFuncArgs(name: "Ti-esc-alpha",
                          kind: ffGauss,
                          gN: n_main / 20.0,
                          gmu: mu_main / 3.0, #400.0e3,
                          gs: sigma_main)
    result.add FitFuncArgs(name: "Ti-esc-beta",
                          kind: ffGauss,
                          gN: n_main / 30.0,
                          gmu: fixed, #mu_main ,
                          gs: fixed) #n_main * 1e3)
    result.add FitFuncArgs(name: "Ti-Kbeta",
                          kind: ffGauss,
                          gN: n_main / 30.0,
                          gmu: fixed, #mu_main, #* 1e5,
                          gs: fixed) #n_main * 1e3)
  of tfAgAg6:
    result.add FitFuncArgs(name: "Ag-Lalpha",
                          kind: ffGauss,
                          gN: n_main,# / 2.0,
                          gmu: mu_main,# * 1e3,
                          gs: sigma_main)# * 10)
    result.add FitFuncArgs(name: "Ag-Lbeta",
                          kind: ffGauss,
                          gN: fixed, #n_main / 10.0,
                          gmu: fixed, #mu_main,
                          gs: fixed) #n_main * 1.5e3)
    #result.add FitFuncArgs(name: "p0",
    #                      kind: ffConst,
    #                      c: n_main * 1e-3)
    #result.add FitFuncArgs(name: "p1",
    #                      kind: ffPol1,
    #                      cp: n_main * 1e-3)
    #result.add FitFuncArgs(name: "p2",
    #                      kind: ffPol2,
    #                      cpp: n_main * 1e-3)
  of tfAlAl4:
    result.add FitFuncArgs(name: "Ag-Kalpha",
                          kind: ffGauss,
                          gN: n_main,# / 2.0,
                          gmu: mu_main,# 350.0e3
                          gs: sigma_main)# * 10)
    #result.add FitFuncArgs(name: "p0",
    #                      kind: ffConst,
    #                      c: -n_main)
    #result.add FitFuncArgs(name: "p1",
    #                      kind: ffPol1,
    #                      cp: n_main)
    #result.add FitFuncArgs(name: "p2",
    #                      kind: ffPol2,
    #                      cpp: -n_main)
  of tfCuEpic2:
    result.add FitFuncArgs(name: "Cu-Lalpha",
                          kind: ffGauss,
                          gN: n_main / 2.0,# / 20.0,
                          gmu: mu_main, #* 0.5e3,
                          gs: sigma_main)
    result.add FitFuncArgs(name: "Cu-Lbeta",
                          kind: ffGauss,
                          gN: n_main / 2.0,# / 50.0,
                          gmu: mu_main,
                          gs: sigma_main)
  of tfCuEpic0_9:
    result.add FitFuncArgs(name: "O-Kalpha",
                          kind: ffGauss,
                          gN: n_main / 2.0,# / 4.0,
                          gmu: mu_main,
                          gs: sigma_main)
    result.add FitFuncArgs(name: "C-Kalpha",
                          kind: ffGauss,
                          gN: n_main / 2.0,# / 4.0,
                          gmu: mu_main,
                          gs: sigma_main)
  of tfCEpic0_6:
    result.add FitFuncArgs(name: "C-Kalpha",
                          kind: ffGauss,
                          gN: n_main,# / 4.0,
                          gmu: mu_main,
                          gs: sigma_main)
    result.add FitFuncArgs(name: "O-Kalpha",
                          kind: ffGauss,
                          gN: n_main,# / 10.0,
                          gmu: fixed, #mu_main, #* 1e3,
                          gs: fixed) #sigma_main  * 1e3)

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

## Declaration of all fit functions
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

proc energyResolution(energyResHits, energyResCharge,
                      errorHits, errorCharge: seq[float]) =
  ## Creates the plot of the energy resolution for both the pixel and charge
  ## spectra
  const energies = @[8.048, 5.899, 4.511, 2.984, 1.487, 0.930, 0.525, 0.277]
  let dfPix = seqsToDf({ "Energy" : energies, "Resolution" : energyResHits,
                         "ResErr" : errorHits })
  let dfCh = seqsToDf({ "Energy" : energies, "Resolution" : energyResCharge,
                        "ResErr" : errorCharge })
  let df = bind_rows([("Pixels", dfPix), ("Charge", dfCh)], id = "Type")
  ggplot(df, aes("Energy", "Resolution", color = "Type")) +
    geom_point() +
    geom_errorbar(aes(yMin = f{`Resolution` - `ResErr`},
                      yMax = f{`Resolution` + `ResErr`})) +
    xlab("Energy / keV") +
    ylab("Energy resolution / %") +
    ggtitle("Energy resolution plot") +
    ggsave(&"energyresoplot-{outdate}.pdf")

proc peakFit(peakPos: seq[float], name: string, error: seq[float]) =
  const energies = @[8.048, 5.899, 4.511, 2.984, 1.487, 0.930, 0.525, 0.277]
  let df = seqsToDf({ "Energy" : energies, "PeakPos" : peakPos,
                      "PeakErr" : error })
  ggplot(df, aes("Energy", "PeakPos")) +
    geom_point() +
    geom_errorbar(aes(yMin = f{`PeakPos` - `PeakErr`},
                      yMax = f{`PeakPos` + `PeakErr`})) +
    xlab("Energy / keV") +
    ylab(&"Peak position for {name}") +
    ggtitle("Peak position of all targets") +
    ggsave(&"{name}.pdf")

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

proc calcFit[T: SomeNumber](dataseq: seq[T],
                            cdlFitFunc: CdlFitFunc,
                            binSize: float,
                            tfKind: TargetFilterKind,
                            dKind: DataKind): (seq[float], seq[float], float, float) =
  let (histdata, bins) = histoCdl(dataseq, binSize, dKind)
  let (pStart, pRes, fitBins, fitHist, errorres) = fitCdlImpl(histdata, bins, tfKind, dKind)
  var bounds: seq[tuple[l, u:float]]
  case dKind
  of Dhits:
    bounds = getbounds(tfKind)
  of Dcharge:
    bounds = getboundsCharge(tfKind)
  doassert pStart.len == bounds.len
  #var opt = newNloptOpt(LD_TNEWTON_PRECOND, pStart.len, bounds)
  fitForNlopt(convertNlopt, cdlFitFunc,
              nfKind = nfMle,
              toExport = false)
  # have to mixin the `convertNlopt` name, since generic symbol resolution
  # happens before macro call (so `convertNlopt` isn't known yet)
  mixin convertNlopt
  var fitObj = FitObject(x: fitBins, y: fitHist) #, yErr: fitHist.mapIt(sqrt(it)))
  var vstruct = newVarStruct(convertNlopt, fitObj)

  #var opt = newNloptOpt[type(fitObj)](LD_MMA, pStart.len, bounds)
  var opt = newNloptOpt[type(fitObj)](LN_COBYLA, pStart.len, bounds)
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

proc fitAndPlot[T: SomeNumber](h5f: var H5FileObj, fitParamsFname: string,
                               tfKind: TargetFilterKind, dKind: DataKind):
               (seq[float], seq[float], seq[float], seq[float]) =
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

  let (histdata, bins) = histoCdl(cutseq, binsizeplot, dKind)
  let (pStart, pRes, fitBins, fitHist, errorsres) = fitCdlImpl(histdata, bins, tfKind, dKind)
  let fitresults = calcfit(cutseq, fitfunc, binsizeplot, tfKind, dKind)
  echo "fitresults nlopt", fitresults[1]
  ##get the interesting fit params

  ## TODO: FIX ME!
  let
    cumu = cumsum(histdata)
    sum = sum(histdata)
    quotient = cumu.mapIt(it/sum)
    lowIdx = quotient.lowerBound(0.0005)
    highIdx = quotient.lowerBound(0.98)
    passbin = bins[lowIdx .. highIdx]
    passhist = histdata[lowIdx .. highIdx]
    passIdx = toSeq(0 .. passhist.high).filterIt(passhist[it] > 0)
    fitBinserr = passIdx.mapIt(passbin[it])
    fitHisterr = passIdx.mapIt(passhist[it])
    err = fitHist.mapIt(1.0)# / sqrt(it))

  let (FitError, mpfitError) = fit(fitfunc, fitresults[1], fitBinserr, fitHisterr, err)
  echo "fit errors: ", mpfitError.error
  ploterror = mpfitError.error

  ## TODO: FIX ME!
  var
    fitmu: float
    fitsig: float
    fitmuerr: float
    fitsigerr: float
    energyres: float
    museq: seq[float]
    sigseq: seq[float]
    energyseq: seq[float]
    energyreserr: float
    energyseqerr: seq[float]

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
  energyres = fitsig / fitmu
  energyreserr = sqrt(pow((fitsigerr / fitmu), 2) + pow((fitsig * fitmuerr / pow(fitmu,2)),2) )
  energyseqerr.add(energyreserr)
  energyseq.add(energyres)

  result[0].add(fitmu)
  result[1].add(energyres)
  result[2].add(fitmuerr)
  result[3].add(energyreserr)

  let mpfitres = calcfitcurve(fitresults[2], fitresults[3], fitfunc, fitresults[0])
  let nloptres = calcfitcurve(fitresults[2], fitresults[3], fitfunc, fitresults[1])
  let startval = calcfitcurve(fitresults[2], fitresults[3], fitfunc, pStart)
  let df = bind_rows([("MPFIT", toTab({ "Energy" : mpfitres[0],
                                       "Counts" : mpfitres[1] })),
                      ("NLopt", toTab({ "Energy" : nloptres[0],
                                       "Counts" : nloptres[1] })),
                      ("Start", toTab({ "Energy" : startval[0],
                                       "Counts" : startval[1] }))],
                     id = "Type")

  ##plot of hits and charge
  # modify bin range if necessary
  if fitmu < binrangeplot / 3.0:
    binrangeplot = binrangeplot / 1.5
  elif fitmu > binrangeplot / 2.0:
    binrangeplot = binrangeplot * 1.5

  let dfRaw = bind_rows([("RawData", toTab({ "Counts" : rawSeq })),
                         ("CutData", toTab({ "Counts" : cutSeq }))],
                         id = "Cut")
  let fname = &"{outname}-{outdate}.pdf"
  ggplot(df, aes("Energy", "Counts")) +
    geom_histogram(data = dfRaw.filter(f{float: `Counts` < binRangePlot}),
                   aes = aes("Counts", fill = "Cut"),
                   position = "identity",
                   alpha = some(0.5),
                   binWidth = binSizePlot) +
    geom_line(aes(color = "Type")) +
    xlab(xtitle) +
    #xlim(0.0, binrangeplot) +
    ylab("Counts") +
    ggtitle(&"target: {tfKind}") +
    ggsave(fname)

  # now dump the fit results, SVG filename and correct parameter names to a file
  ## TODO: add fit parameters as annotation to plot!!
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
    (numBins, minVal, maxVal) = cdlToXrayBinning2018(outname)
  # given passing data, calculate histogram and write to file
  if numBins > 0: # means its a dataset that will show up in the resulting file
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
  echo "INFO: Apply filtering to dset: ", dsetName
  let data = h5f.readAs(dsetName, float)
  # apply passIdx
  result = passIdx.mapIt(data[it])

template datasetCreation(h5f: untyped, name, dlen, `type`: untyped): untyped =
  ## inserts the correct data set creation parameters
  h5f.create_dataset(name,
                     dlen,
                     dtype = `type`,
                     chunksize = @[5000, 1], # some hardcoded chunk size
                     maxshape = @[int.high, 1])

template createAndWrite(h5read, h5write, `type`, dset, outname: untyped): untyped =
  let vlenType = special_type(`type`)
  let outDset = h5write.datasetCreation(outname, dset.shape, vlenType)
  outDset[outDset.all] = h5read[dset.name, vlenType, `type`]

template getAndAdd(h5read, h5write, dset, `type`, outDset: untyped): untyped =
  let vlenType = special_type(`type`)
  outDset.add h5read[dset.name, vlenType, `type`]

proc generateCdlCalibrationFile(h5file: string, year: YearKind,
                                outfile = "calibration-cdl") =
  ## generates the CDL calibration data file from a HDF5 file containing
  ## all CDL runs. Supports either 2014 CDL data or 2019 CDL data.
  # walk all runs corresponding to a single `TargetFilterKind` and
  # combine the datasets into the output files
  var h5f = H5file(h5file, "r")
  var h5fout = H5file(&"{outfile}-{year}.h5", "rw")
  let runs = readRuns(filename)
  for tfKind in TargetFilterKind:
    for grp in tfRuns(h5f, tfKind):
      # will not iterate the datasets and write to the outfile
      # via hyperslabs, potentially appending to the existing
      # dataset
      var mgrp = grp
      for dset in mgrp:
        if dset.shape.len == 1 or dset.shape[1] > 1:
          # skip all datasets in the input, which are 2 dimensional. That includes
          # {polya, polyaFit, FeSpectrum, CDL...}
          # dset.shape == 1 is for `FeSpectrumCharge`, since it's the only dataset
          # atm that is actually 1D (see TPA issue #39).
          continue
        # TODO: add appropriate conversion of TPA naming to Marlin naming for
        # 2014 raw input or allow `genRefFile` option to deal with TPA naming
        # even for 2014 input!
        let outname = cdlGroupName($tfKind, $year, dset.name)
        if outname notin h5fout:
          # create dataset first
          case dset.dtypeAnyKind
          of dkSequence:
            # keep their base type and create special type
            case dset.dtypeBaseKind
            of dkUint8:
              createAndWrite(h5f, h5fout, uint8, dset, outname)
            of dkUint16:
              createAndWrite(h5f, h5fout, uint16, dset, outname)
            of dkFloat64:
              # charge
              createAndWrite(h5f, h5fout, float64, dset, outname)
            else:
              raise newException(Exception, "??? " & $dset)
          else:
            let outDset = h5fout.datasetCreation(outname, dset.shape, float)
            outDset[outDset.all] = h5f.readAs(dset.name, float)
        else:
          # append to dataset
          var outDset = h5fout[outname.dset_str]
          case outDset.dtypeAnyKind
          of dkSequence:
            # keep their base type and create special type
            case outDset.dtypeBaseKind
            of dkUint8:
              getAndAdd(h5f, h5fout, dset, uint8, outDset)
            of dkUint16:
              getAndAdd(h5f, h5fout, dset, uint16, outDset)
            of dkFloat64:
              getAndAdd(h5f, h5fout, dset, float64, outDset)
            else:
              raise newException(Exception, "??? " & $dset)
          else:
            outDset.add h5f.readAs(dset.name, float)

  h5fout.attrs["FrameworkKind"] = $CdlGenerateNamingScheme
  discard h5f.close()
  discard h5fout.close()

proc generateXrayReferenceFile(h5file: string, year: YearKind,
                               outfile = "XrayReferenceFile") =
  ## generates the X-ray reference data file
  # this is achieved by taking the raw TargetFilterKind runs, combining them
  # into a the full CDL calibration file (if the input is the raw run based
  # file, else it skips the first step)
  # then we apply the charge cuts and bin the data by N bins
  # the result is written to the new file as (N, 2) datasets
  var h5f = H5file(h5file, "r")
  # read framework kind from `h5f`
  var frameworkKind = fkMarlin
  if "FrameworkKind" in h5f.attrs:
    frameworkKind = parseEnum[FrameworkKind](h5f.attrs["FrameworkKind", string])

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
  var xrayRefCuts: Table[string, Cuts]
  case year
  of yr2014:
    xrayRefCuts = getEnergyBinMinMaxVals2014() #XraySpectrumCutVals()
  of yr2018:
    xrayRefCuts = getEnergyBinMinMaxVals2018()

  var h5fout = H5file(outfile & $year & ".h5", "rw")
  for group in h5f:
    var mgrp = group
    echo group.name
    if scanf(group.name, "/calibration-cdl-$+-$+kV", date, tfKindStr):
      doAssert year == yr2014 or year == yr2018
      # now perform cuts on datasets
      doAssert tfKindStr & "kV" in xrayRefCuts
      let cut = xrayRefCuts[tfKindStr & "kV"]
      let passIdx = cutOnProperties(
        h5f,
        group,
        crSilver, # try cutting to silver
        (toDset(igRmsTransverse, frameworkKind), cut.minRms, cut.maxRms),
        (toDset(igLength, frameworkKind), 0.0, cut.maxLength),
        (toDset(igHits, frameworkKind), cut.minPix, Inf),
        (toDset(igTotalCharge, frameworkKind), cut.minCharge, cut.maxCharge))
      echo "Number of passing indices ", passIdx.len
      # given passIdx, now read each dataset iteratively and apply cuts
      for dset in mgrp:
        case dset.dtypeAnyKind
        of dkSequence:
          # variable length data, x, y, charge, will be dropped in conversion
          discard
        else:
          # get the name of the dataset in the output
          let dsetName = dset.name.extractFilename
          var outname = ""
          case frameworkKind
          of fkMarlin:
            outname = cdlToXray2014(dsetName)
          of fkTpa:
            # if 2014 data reconstructed with TPA, we'll end up here
            outname = dsetName # use dataset name as is
          doAssert outname.len > 0, "Not possible anymore"
          # parse outname and check if it is in `XrayReferenceDsets`
          if outname.toInGridDset in XrayReferenceDsets:
            let passedData = readAndFilter(h5f, dset.name, passIdx)
            cdlToXrayTransform(h5fout, passedData, tfKindStr & "kV", outname, year)
      if frameworkKind == fkMarlin:
        # data generated by Marlin does not contain `lengthdivbyrmsy` in the calibration file
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
    else:
      raise newException(ValueError, "Could not match calibration group in " &
        $h5file & ", group name was " & $group.name)
  # finally add the type of reference file to the root group of the output,
  # which indicates which naming scheme is used
  h5fout.attrs["FrameworkKind"] = $XrayRefGenerateNamingScheme
  discard h5fout.close()
  discard h5f.close()

proc main =
  let args = docopt(doc)
  echo "ARGS", args

  let h5file = $args["<h5file>"]
  let reco_order = args["--cutcdl"].toBool
  let genRefFile = args["--genRefFile"].toBool
  let genCdlFile = args["--genCdlFile"].toBool
  let genOutfile = if $args["--outfile"] != "nil":
                     $args["--outfile"]
                   else:
                     ""
  let dumpAccurate = if args["--dumpAccurate"].toBool:
                       true
                     else:
                       false

  template callGenFunc(fn: untyped): untyped =
    let year = parseEnum[YearKind]($args["--year"])
    if genOutFile.len > 0:
      fn(h5file, year, genOutFile)
    else:
      fn(h5file, year)
  if genRefFile:
    callGenFunc(generateXrayReferenceFile)
  if genCdlFIle:
    callGenFunc(generateCdlCalibrationFile)
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
    var h5f = H5file(h5file, "rw")
    for tfkind in TargetFilterKind:
      let energyHits = fitAndPlot[int64](h5f, fitParamsFname, tfkind, Dhits)
      peakposHits.add(energyHits[0])
      energyResHits.add(energyHits[1])
      peakHitsErr.add(energyHits[2])
      energyHitsErr.add(energyHits[3])
      #echo "energyres ", energyResHit
      let energyCharge = fitAndPlot[float64](h5f, fitParamsFname, tfkind, Dcharge)
      peakposCharge.add(energyCharge[0])
      energyResCharge.add(energyCharge[1])
      peakChargeErr.add(energyCharge[2])
      energyChargeErr.add(energyCharge[3])
    discard h5f.close()
    energyResolution(energyResHits, energyResCharge, energyHitsErr, energyChargeErr)
    peakFit(peakposHits, "Hits", peakHitsErr)
    peakFit(peakposCharge, "Charge", peakChargeErr)

when isMainModule:
  main()
