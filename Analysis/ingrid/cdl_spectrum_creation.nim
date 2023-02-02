#[ WARNING:
 Compiling this file requires a compiler patch! Without it you'll be greeted by a "generic instantiation
 too nested" error.

diff --git a/compiler/seminst.nim b/compiler/seminst.nim
index bd5eb1ec3..170061aaf 100644
--- a/compiler/seminst.nim
+++ b/compiler/seminst.nim
@@ -321,7 +321,7 @@ proc generateInstance(c: PContext, fn: PSym, pt: TIdTable,
   # no need to instantiate generic templates/macros:
   internalAssert c.config, fn.kind notin {skMacro, skTemplate}
   # generates an instantiated proc
-  if c.instCounter > 50:
+  if c.instCounter > 1000:
     globalError(c.config, info, "generic instantiation too nested")
   inc(c.instCounter)
   # careful! we copy the whole AST including the possibly nil body!
]#

import std / [parsecsv, os, streams, strutils, strformat, tables, sequtils,
              macros, fenv, algorithm, times, strscans, typeinfo]
import ingrid / [ingrid_types, tos_helpers]
import ingrid / calibration
import ingrid / calibration / [fit_functions, calib_fitting]
import cdlFitting / cdlFitMacro
import helpers / utils
import cligen / macUt
import pkg / [mpfit, nlopt, nimhdf5, parsetoml, seqmath, measuremancer]

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

when defined(linux):
  const commitHash = staticExec("git rev-parse --short HEAD")
else:
  const commitHash = ""
# get date using `CompileDate` magic
const compileDate = CompileDate & " at " & CompileTime
const versionStr = "Version: $# built on: $#" % [commitHash, compileDate]


##some constants depending on the run
import projectDefs
const filename = TpxDir / "resources/cdl_runs_2019.org"
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

  FitResult = object
    pStart: seq[float]
    pRes: seq[float]
    pErr: seq[float]
    χ²dof: float

  FitData = object
    bins: seq[float]
    hist: seq[float]
    errs: seq[float]

  MainPeak = object
    tfKind: TargetFilterKind
    runNumber: int # If any, if `fitByRun` is false, this remains 0
    dKind: DataKind
    fit_μ: Measurement[float]
    fit_σ: Measurement[float]
    energyRes: Measurement[float]

proc removeSuffix(s, suff: string): string =
  result = s
  result.removeSuffix(suff)

# helper converters to get targets, filters and kV values from a `TargetFilterKind`
proc toTarget(tfKind: TargetFilterKind): TargetKind = ($tfKind).split("-")[0].parseEnum[:TargetKind]()
proc toFilter(tfKind: TargetFilterKind): FilterKind = ($tfKind).split("-")[1].parseEnum[:FilterKind]()
proc toHV(tfKind: TargetFilterKind, withSuffix = true): string =
  if withSuffix: ($tfKind).split("-")[2]
  else: ($tfKind).split("-")[2].removeSuffix("kV")

proc readFitByRun(config: string): bool =
  ## Reads the `fitByRun` field from the configuration file, deciding whether to
  ## perform our fits by run or by target/filter kind only.
  let configPath = if config.len == 0:
                     const sourceDir = currentSourcePath().parentDir
                     sourceDir / "config.toml"
                   else:
                     config
  let config = parseToml.parseFile(configPath)
  result = config["CDL"]["fitByRun"].getBool

func getLines(hist, binning: seq[float], tfKind: TargetFilterKind): seq[FitFuncArgs] =
  ## this is a runtime generator for the correct fitting function prototype,
  ## i.e. it returns a seq of parts, which need to be combined to the complete
  ## function at runtime
  let (mu_main, sigma_main, n_main, _, _, _) = getLines(hist, binning)
  # Note: It is important that the *first* line is the *main peak* we fit!
  case tfKind
  of tfCuNi15:
    result.add FitFuncArgs(name: "Cu-Kalpha",
                          kind: ffExpGauss,
                          ea: n_main * 1e-10,
                          eb: n_main * 1e-12,
                          eN: n_main / 1.5,
                          emu: mu_main, #200.0
                          es: sigma_main) #13.0
    result.add FitFuncArgs(name: "Cu-Kalpha-esc",
                          kind: ffExpGauss,
                          ea: n_main *  1e-10,
                          eb: n_main * 1e-12,
                          eN: n_main / 10.0, #100.0
                          emu: mu_main / 1.5,  #100.0
                          es: sigma_main) # / 2.0) #/ 15.0) #16.0
  of tfMnCr12:
    result.add FitFuncArgs(name: "Mn-Kalpha",
                          kind: ffExpGauss,
                          ea: -5.0,
                          eb: 0.049,#1e-3,
                          eN: n_main,# / 10.0, #30.0
                          emu: mu_main, #200.0
                          es: 15.0) #sigma_main / 2.0) # 13.0
    result.add FitFuncArgs(name: "Mn-Kalpha-esc",
                          kind: ffGauss,
                          #ea: 5.0,
                          #eb: 1e-3,
                          gN: n_main / 10.0, #120
                          gmu: mu_main / 2.0, #100.0
                          gs: 10.0) # sigma_main / 2.0) #16.0
  of tfTiTi9:
    result.add FitFuncArgs(name: "Ti-Kalpha",
                          kind: ffExpGauss,
                          ea: n_main * 1e-10,
                          eb: n_main * 1e-12,
                          eN: n_main / 2.0,
                          emu: mu_main,
                          es: sigma_main)
    result.add FitFuncArgs(name: "Ti-Kalpha-esc",
                          kind: ffGauss,
                          gN: n_main / 20.0,
                          gmu: mu_main / 3.0,
                          gs: sigma_main)
    result.add FitFuncArgs(name: "Ti-Kbeta",
                          kind: ffGauss,
                          gN: n_main / 2.0,# / 10.0,
                          gmu: fixed, #mu_main,
                          gs: fixed) #n_main / 15.0)
    result.add FitFuncArgs(name: "Ti-Kbeta-esc",
                          kind: ffGauss,
                          gN: n_main / 20.0,# / 20.0,
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
                          gN: fixed, # n_main / 10.0,# / 50.0,
                          gmu: fixed,
                          gs: fixed)
    result.add FitFuncArgs(name: "O-Kalpha",
                          kind: ffGauss,
                          gN: fixed,
                          gmu: fixed,
                          gs: fixed)
    result.add FitFuncArgs(name: "unknown",
                           kind: ffGauss,
                           gN: n_main / 3.0,
                           gmu: mu_main * 1.75,
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
    result.add FitFuncArgs(name: "unknown",
                           kind: ffGauss,
                           gN: n_main / 3.0,
                           gmu: mu_main * 1.75,
                           gs: sigma_main)
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

func getBounds(tfKind:TargetFilterKind): seq[tuple[l, u:float]] =
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
               (l: 175.0, u:Inf),
               (l: 1.5, u:20.0),
               #(l: 0.0, u:Inf), # force a, b for escape to > 0, as otherwise it "takes over" and breaks the fit
               #(l: 0.0, u:Inf),
               (l: 5.0, u:50),
               (l: 50.0, u:120.0),
               (l: 5.0, u:15.0)]
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
               (l: 100.0, u:400.0),
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
               (l: 1.0, u:20.0),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:20.0)]
  of tfCuEpic0_9:
    result = @[(l: 1.0, u:Inf),
               (l: 1.0, u:23.0),
               (l: 1.0, u:7.0),
               (l: 1.0, u:Inf),
               (l: 5.0, u:60.0),
               (l: 5.0, u:110.0)]
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
  # Note: It is important that the *first* line is the *main peak* we fit!
  case tfKind
  of tfCuNi15:
    result.add FitFuncArgs(name: "Cu-Kalpha",
                          kind: ffGauss,
                          gN: n_main,
                          gmu: mu_main,
                          gs: sigma_main)
    result.add FitFuncArgs(name: "Cu-Kalpha-esc",
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
    result.add FitFuncArgs(name: "Mn-Kalpha-esc",
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
    result.add FitFuncArgs(name: "Ti-Kalpha-esc",
                          kind: ffGauss,
                          gN: n_main / 8.0,
                          gmu: mu_main / 3.0, #400.0e3,
                          gs: sigma_main / 3.0)
    result.add FitFuncArgs(name: "Ti-Kbeta",
                          kind: ffGauss,
                          gN: n_main / 30.0,
                          gmu: fixed, #mu_main, #* 1e5,
                          gs: fixed) #n_main * 1e3)
    result.add FitFuncArgs(name: "Ti-Kbeta-esc",
                          kind: ffGauss,
                          gN: n_main / 30.0,
                          gmu: fixed, #mu_main ,
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
                          gN: n_main,# / 20.0,
                          gmu: mu_main, #* 0.5e3,
                          gs: sigma_main / 4.0)
    result.add FitFuncArgs(name: "Cu-Lbeta",
                          kind: ffGauss,
                          gN: fixed, # n_main / 10.0,# / 50.0,
                          gmu: fixed,
                          gs: fixed)
                          #gmu: mu_main / 3.0,
                          #gs: sigma_main / 4.0)
    result.add FitFuncArgs(name: "O-Kalpha",
                          kind: ffGauss,
                          gN: fixed,
                          gmu: fixed,
                          gs: fixed)
    result.add FitFuncArgs(name: "unknown",
                          kind: ffGauss,
                          gN: n_main / 2.0,# / 4.0,
                          gmu: mu_main * 2.0,
                          gs: sigma_main)
  of tfCuEpic0_9:
    result.add FitFuncArgs(name: "O-Kalpha",
                          kind: ffGauss,
                          gN: n_main,# / 4.0,
                          gmu: mu_main,
                          gs: sigma_main)
    result.add FitFuncArgs(name: "C-Kalpha",
                          kind: ffGauss,
                          gN: fixed, #n_main / 10.0,# / 4.0,
                          gmu: fixed,
                          gs: fixed)
    result.add FitFuncArgs(name: "unknown",
                          kind: ffGauss,
                          gN: n_main / 2.0,# / 4.0,
                          gmu: mu_main * 2.0,
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

func getBoundsCharge(tfKind:TargetFilterKind): seq[tuple[l, u:float]] =
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
    result = @[(l: 500.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 1.0, u:Inf),
               (l: 50.0, u:100.0),
               (l: 1.0, u:3e5),
               (l: 1.0, u:Inf),
               (l: 5.0, u:50.0),
               (l: 5.0, u:50.0)]
  of tfAgAg6:
    result = @[(l: 390.0, u:Inf),
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
               #(l: 1.0, u:Inf)]
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
# Note: It is important that the *first* line is the *main peak* we fit!
declareFitFunc(cuNi15):
  ffExpGauss: "Cu-Kalpha"
  ffExpGauss: "Cu-Kalpha-esc"
declareFitFunc(cuNi15Charge):
  ffGauss: "Cu-Kalpha"
  ffGauss: "Cu-Kalpha-esc"
declareFitFunc(mnCr12):
  ffExpGauss: "Mn-Kalpha"
  ffGauss: "Mn-Kalpha-esc"
declareFitFunc(mnCr12Charge):
  ffGauss: "Mn-Kalpha"
  ffGauss: "Mn-Kalpha-esc"
declareFitFunc(tiTi9):
  ffExpGauss: "Ti-Kalpha"
  ffGauss: "Ti-Kalpha-esc"
    #name = "Ti-esc-alpha"
    #gmu = emu("Ti-Kalpha") * (1.537/4.511)
  ffGauss: #"Ti-Kbeta"
    name = "Ti-Kbeta"
    gmu = emu("Ti-Kalpha") * (4.932/4.511)
    gs = es("Ti-Kalpha")
  ffGauss: #"Ti-esc-beta"
    name = "Ti-Kbeta-esc"
    gmu = emu("Ti-Kalpha") * (1.959/4.511)
    gs = gs("Ti-Kalpha-esc")
declareFitFunc(tiTi9Charge):
  ffGauss: "Ti-Kalpha"
  ffGauss: "Ti-Kalpha-esc"
    #name = "Ti-esc-alpha"
    #gmu = emu("Ti-Kalpha") * (1.537/4.511)
  ffGauss: #"Ti-Kbeta"
    name = "Ti-Kbeta"
    gmu = gmu("Ti-Kalpha") * (4.932/4.511)
    gs = gs("Ti-Kalpha")
  ffGauss: #"Ti-esc-beta"
    name = "Ti-Kbeta-esc"
    gmu = gmu("Ti-Kalpha") * (1.959/4.511)
    gs = gs("Ti-Kalpha-esc")
declareFitFunc(agAg6):
  ffExpGauss: "Ag-Lalpha"
  ffGauss: #"Ag-Lbeta"
    name = "Ag-Lbeta"
    gN = eN("Ag-Lalpha") * 0.56
    gmu = emu("Ag-Lalpha") * (3.151/2.984)
    gs = es("Ag-Lalpha")
declareFitFunc(agAg6Charge):
  ffGauss: "Ag-Lalpha"
  ffGauss: #"Ag-Lbeta"
    name = "Ag-Lbeta"
    gN = gN("Ag-Lalpha") * 0.56
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
  #ffGauss: "Cu-Lbeta"
  ffGauss:
    name = "Cu-Lbeta"
    gN = gN("Cu-Lalpha") * (0.65 / 1.11)
    gmu = gmu("Cu-Lalpha") * (0.9498 / 0.9297) # energies in keV of the lines
    gs = gs("Cu-Lalpha")
  ffGauss:
    name = "O-Kalpha"
    gN = gN("Cu-Lalpha") / 3.5 ## this has no physical motivation!
    gmu = gmu("Cu-Lalpha") * (0.5249 / 0.9297)
    gs = gs("Cu-Lalpha") / 2.0
  ffGauss: "unknown"
declareFitFunc(cuEpic2Charge):
  ffGauss: "Cu-Lalpha" ## XXX: replace two lines by mix of one, adjust energy to relative importance
  ffGauss:
    name = "Cu-Lbeta"
    gN = gN("Cu-Lalpha") * (0.65 / 1.11) #/ 5.0
    gmu = gmu("Cu-Lalpha") * (0.9498 / 0.9297) # energies in keV of the lines
    gs = gs("Cu-Lalpha")
  ffGauss:
    name = "O-Kalpha"
    gN = gN("Cu-Lalpha") / 3.5 ## this has no physical motivation!
    gmu = gmu("Cu-Lalpha") * (0.5249 / 0.9297)
    gs = gs("Cu-Lalpha") / 2.0
  ffGauss: "unknown"
  #ffGauss:
  #  name = "O-Kalpha"
  #  gN = gN("Cu-Lalpha") / 4.0
  #  gmu = gmu("Cu-Lalpha") * 0.5249 / 0.9297
  #  gs = gs("Cu-Lalpha") / 2.0
declareFitFunc(cuEpic0_9):
  ffGauss: "O-Kalpha"
  ffGauss: "unknown"
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
  ffGauss:
    name = "C-Kalpha"
    gN = gN("O-Kalpha") / 10.0 ## this has no physical motivation!
    gmu = gmu("O-Kalpha") * (277.0 / 524.9)
    gs = gs("O-Kalpha")
  ffGauss: "unknown"
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
  #echo result.repr

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

proc getFitData(hist, bins: seq[float]): FitData =
  ## Returns the `FitData` (i.e. the real bins, counts and errors we fit to)
  ## based on the cut CDL data. This mainly implies removing any outliers still left.
  ##
  ## Note that this returns the data such that the bins correspond to bin *centers*!
  ## That's because that is what we actually fit to.
  ## XXX: can't we simplify this?
  let cumu = cumsum(hist)
  let sum = sum(hist)
  let quotient = cumu.mapIt(it/sum)
  let lowIdx = quotient.lowerBound(0.001)
  let highIdx = quotient.lowerBound(0.98)
  result = FitData(bins: newSeqOfCap[float](hist.len),
                   hist: newSeqOfCap[float](hist.len),
                   errs: newSeqOfCap[float](hist.len))
  doAssert bins.len > 0
  # we know all bin widths are the same, so this is fine:
  let binWidth = bins[1] - bins[0]
  for idx in lowIdx .. highIdx:
    if hist[idx] > 0.0:
      result.bins.add bins[idx] + (binWidth / 2.0) # add half bin width, as bins are left edges
      result.hist.add hist[idx]
      result.errs.add sqrt(hist[idx])

proc getFuncAndStartParams(fitData: FitData,
                           tfKind: TargetFilterKind, dKind: DataKind): (seq[float], CdlFitFunc) =
  ## Returns the start parameters for the desired fit and data as well as the `CdlFitFunc`
  ## (the function we actually fit)
  case dKind
  of Dhits:
    let lines = getLines(fitData.hist, fitData.bins, tfKind)
    result = (lines.serialize, getCdlFitFunc(tfKind))
  of Dcharge:
    let lines = getLinesCharge(fitData.hist, fitData.bins, tfKind)
    result = (lines.serialize, getCdlFitFuncCharge(tfKind))

proc fitCdlMpfit(fitData: FitData, tfKind: TargetFilterKind, dKind: DataKind): FitResult =
  ## Performs the fit using `mpfit`.
  let (params, fitFunc) = getFuncAndStartParams(fitData, tfKind, dKind)
  var bounds: seq[tuple[l, u:float]]
  case dKind
  of Dhits:
    bounds = getBounds(tfKind)
  of Dcharge:
    bounds = getBoundsCharge(tfKind)

  let (pRes, res) = fit(fitFunc,
                        params,
                        fitData.bins,
                        fitData.hist,
                        fitData.errs)

  echoResult(pRes, res=res)
  echo "error fitcdlimpl ", res.error
  result = FitResult(pStart: params,
                     pRes: pRes,
                     pErr: res.error,
                     χ²dof: res.reducedChiSq)

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

proc toTfKind(run: CdlRun): TargetFilterKind =
  result = parseEnum[TargetFilterKind](&"{toCutStr(run)}")

iterator tfRuns(h5f: H5File, tfKind: TargetFilterKind): (int, H5Group) =
  ## Yields the center chip group of all runs from `filename`,
  ## which match `tfKind`
  let runs = readRuns(filename)
  for r in runs:
    case r.runType
    of rtXrayFinger:
      let tfk = r.toTfKind
      if tfk == tfKind:
        let runGrp = h5f[recoRunGrpStr(r.number)]
        let centerChip = runGrp.attrs["centerChip", int]
        let chpGrp = h5f[(runGrp.name / "chip_" & $centerChip).grp_str]
        yield (r.number, chpGrp)
    else:
      # nothing to yield if not an "XrayFinger" (read CDL) run
      discard

proc getCdlCutIdxs(h5f: H5File, runNumber: int, tfKind: TargetFilterKind): seq[int] =
  let cutTab = getXrayCleaningCuts()
  let grp = h5f[(recoDataChipBase(runNumber) & chipnumber).grp_str]
  let cut = cutTab[$tfKind]
  result = cutOnProperties(h5f,
                           grp,
                           cut.cutTo,
                           ("rmsTransverse", cut.minRms, cut.maxRms),
                           ("length", 0.0, cut.maxLength),
                           ("hits", cut.minPix, Inf),
                           ("eccentricity", 0.0, cut.maxEccentricity))

proc readCutCDL[T](h5f: H5File, runNumber: int, dset: string,
                   tfKind: TargetFilterKind, _: typedesc[T]): seq[T] =
  let passIdx = h5f.getCdlCutIdxs(runNumber, tfKind)
  let data = h5f.readAs(recoDataChipBase(runNumber) & $chipnumber / dset, T)
  result = passIdx.mapIt(data[it])

proc readCutCDL[T](h5f: H5File, runNumber: int, dset: string,
                   tfKind: TargetFilterKind, passIdx: seq[int], _: typedesc[T]): seq[T] =
  let data = h5f.readAs(recoDataChipBase(runNumber) & $chipnumber / dset, T)
  result = passIdx.mapIt(data[it])

proc cutAndWrite(h5file: string) =
  ## XXX: this is imo not really very useful. Instead we should make sure to always cut the data
  ## on the fly. The cuts are cheap after all!
  let runs = readRuns(filename)
  var h5f = H5open(h5file, "rw")
  defer: discard h5f.close()
  for r in runs:
    #if r.number != 315:
      #continue
    case r.runType
    of rtXrayFinger:
      let grp = h5f[(recoDataChipBase(r.number) & chipnumber).grp_str]
      let passIdx = h5f.getCdlCutIdxs(r.number, r.toTfKind())
      let nevents = passIdx.len

      proc writeDset(r: CDLRun, dsetWrite, dsetRead: string, datatype: typedesc) =
        var
          dset = h5f.create_dataset(grp.name / dsetWrite, nevents,
                                    datatype)
        if dsetWrite == "CdlSpectrumIndices":
          dset[dset.all] = passIdx
        else:
          let read = h5f[grp.name / dsetRead, datatype]
          dset[dset.all] = passIdx.mapIt(read[it])
        dset.attrs["Target"] = $r.target
        dset.attrs["Filter"] = $r.filter
        dset.attrs["HV"] = $r.hv
      writeDset(r, "CdlSpectrumIndices", "", int64)
      writeDset(r, "CdlSpectrum", "hits", int64)
      writeDset(r, "CdlSpectrumEvents", "eventNumber", int64)
      writeDset(r, "CdlSpectrumCharge", "totalCharge", float64)

      let runnum = h5f[(recoBase() & $r.number).grp_str]
      runnum.attrs["tfKind"] = $r.toCutStr
    else:
      discard

proc energyResolution(peaksHits, peaksCharge: seq[MainPeak],
                      pathPrefix: string) =
  ## Creates the plot of the energy resolution for both the pixel and charge
  ## spectra
  ## XXX: implement splits by `run` if `fitByRun`
  const peakEnergies = getXrayFluorescenceLines().reversed
  # hits & charge peaks in same order, so energies then same
  let energies = peaksHits.mapIt(peakEnergies[ord(it.tfKind)])
  let eH = peaksHits.mapIt(it.energyRes.value)
  let eHErr = peaksHits.mapIt(it.energyRes.error)
  let dfPix = toDf({ "Energy" : energies, "Resolution" : eH,
                     "ResErr" : eHErr })
  let eC = peaksCharge.mapIt(it.energyRes.value)
  let eCErr = peaksCharge.mapIt(it.energyRes.error)
  let dfCh = toDf({ "Energy" : energies, "Resolution" : eC,
                    "ResErr" : eCErr })
  let df = bind_rows([("Pixels", dfPix), ("Charge", dfCh)], id = "Type")
  ggplot(df, aes("Energy", "Resolution", color = "Type")) +
    geom_point() +
    geom_errorbar(aes(yMin = f{`Resolution` - `ResErr`},
                      yMax = f{`Resolution` + `ResErr`})) +
    xlab("Energy [keV]") +
    ylab("Energy resolution [%]") +
    ggtitle("Energy resolution depending on energy") +
    ggsave(&"{pathPrefix}/energyresoplot-{outdate}.pdf", width = 800, height = 480)

proc peakFit(mainPeaks: seq[MainPeak], name: string, pathPrefix: string) =
  ## XXX: implement splits by `run` if `fitByRun`
  const peakEnergies = getXrayFluorescenceLines().reversed
  let peaks = mainPeaks.mapIt(it.fit_μ.value)
  let pErr = mainPeaks.mapIt(it.fit_μ.error)
  let energies = mainPeaks.mapIt(peakEnergies[ord(it.tfKind)])
  let df = toDf({ "Energy" : energies, "PeakPos" : peaks, "PeakErr" : pErr })
  ggplot(df, aes("Energy", "PeakPos")) +
    geom_point() +
    geom_errorbar(aes(yMin = f{`PeakPos` - `PeakErr`},
                      yMax = f{`PeakPos` + `PeakErr`})) +
    xlab("Energy [keV]") +
    ylab(&"Peak position for {name}") +
    ggtitle("Peak position of all targets") +
    ggsave(&"{pathPrefix}/{name}.pdf", width = 800, height = 480)

proc serializeFitParameters(fitResult: FitResult,
                            tfKind: TargetFilterKind, dKind: DataKind,
                            dumpAccurate: bool): string =
  result.add &"tfKind: {tfKind}\n"
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
      if dumpAccurate:
        result = &"{f:.2e}"
      else:
        result = &"{f:.1e}"
    else:
      if dumpAccurate:
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

  let pStrs = fitResult.pRes.paddedStrs
  let eStrs = fitResult.pErr.paddedStrs
  # iterate the lines and then unroll the FitFuncArgs object
  var i = 0
  var lineName: string
  for line in fitLines:
    var allNaN = true
    for field, val in fieldPairs(line):
      when type(val) is string:
        lineName = val
        result.add &"{lineName}:"
      elif type(val) is float:
        let fstr = $field
        if classify(val) != fcNaN:
          result.add &"\n  {fstr:<3} = {pStrs[i]} ± {eStrs[i]}"
          inc i
          allNaN = false
    if allNaN:
      result.add " fixed\n"
    else:
      result.add "\n"
  result.add &"χ²/dof = {fitResult.χ²dof:.2f}"

proc dumpFitParameters(outfile, svgFname: string,
                       fitResult: FitResult,
                       tfKind: TargetFilterKind, dKind: DataKind) =
  ## dumps the fit paramters and their names, plus the filename of the corresponding SVG
  ## to a txt file
  var outf = open(outfile, fmAppend)
  outf.write &"svg: {svgFname}\n"
  outf.write(serializeFitParameters(fitResult, tfKind, dKind, "accurate" in outfile))
  outf.close()

proc calcFitCurve(minbin: float, maxbin: float,
                  cdlFitFunc: CdlFitFunc,
                  params: seq[float]): DataFrame =
  ## Computes the line of the fit parameters for the given function and returns it
  ## as a DataFrame.
  let
    minvalue = minbin
    maxvalue = maxbin
    range = linspace(minvalue.float, maxvalue.float, 1500)
    yvals = range.mapIt(cdlFitFunc(params, it))
  result = toDf({"Energy" : range, "Counts" : yvals})

proc calcFitCurve(fitData: FitData,
                  fitRes: FitResult,
                  cdlFitFunc: CdlFitFunc): DataFrame =
  ## Wrapper around above.
  result = calcFitCurve(fitData.bins.min, fitData.bins.max, cdlFitFunc, fitRes.pRes)

proc chiSq(fitData: FitData, pRes: seq[float], fitFunc: CdlFitFunc): float =
  let dof = fitData.hist.len - pRes.len
  for i in 0 ..< fitData.hist.len:
    result += ( (fitData.hist[i] - fitFunc(pRes, fitData.bins[i])) / fitData.errs[i] )^2
  result /= dof.float

proc fitCdlNlopt(fitData: FitData,
                 tfKind: TargetFilterKind,
                 dKind: DataKind): FitResult =
  ## Performs the fit using `Nlopt`
  let (params, fitFunc) = getFuncAndStartParams(fitData, tfKind, dKind)
  var bounds: seq[tuple[l, u:float]]
  case dKind
  of Dhits:
    bounds = getBounds(tfKind)
  of Dcharge:
    bounds = getBoundsCharge(tfKind)

  # now make sure all start parameters are within bounds, else clamp
  var pStart = params
  doAssert pStart.len == bounds.len
  for i in 0 ..< pStart.len:
    pStart[i] = clamp(pStart[i], bounds[i].l, bounds[i].u)
  #var opt = newNloptOpt(LD_TNEWTON_PRECOND, pStart.len, bounds)
  fitForNlopt(convertNlopt, fitFunc,
              nfKind = nfChiSqGrad,
              toExport = false)
  # have to mixin the `convertNlopt` name, since generic symbol resolution
  # happens before macro call (so `convertNlopt` isn't known yet)
  mixin convertNlopt
  var fitObj = FitObject(x: fitData.bins, y: fitData.hist, yErr: fitData.errs)
  var vstruct = newVarStruct(convertNlopt, fitObj)

  var opt = newNloptOpt[type(fitObj)](LD_MMA, pStart.len, bounds)
  #var opt = newNloptOpt[type(fitObj)](LN_COBYLA, pStart.len, bounds)
  #var opt = newNloptOpt[type(fitObj)](LN_BOBYQA, pStart.len, bounds)
  opt.setFunction(vstruct)
  opt.xtol_rel = 1e-10
  opt.ftol_rel = 1e-10
  opt.xtol_abs = 1e-14
  opt.ftol_abs = 1e-14
  opt.maxtime  = 5.0
  opt.maxEval  = 20000
  #opt.initialStep = 1.0
  let (paramsN, minN) = opt.optimize(pStart)
  if opt.status == NLOPT_INVALID_ARGS:
    echo "NLOPT Fit failed. Likely cause: satrt parameter outside bounds:"
    for i in 0 ..< pStart.len:
      echo "pStart[", i, "] =", pStart[i], " with bound: ", bounds[i]
    quit()
  elif opt.status < NLOPT_SUCCESS:
    echo "NLOPT fit failed with code ", opt.status
    quit()
  echo "=============== NLOPT parameters ==============="
  for i in 0 ..< pStart.len:
    echo "p[", i, "] = ", paramsN[i]
  let χ²dof = chiSq(fitData, paramsN, fitFunc)
  echo "χ²dof = ", χ²dof
  echo "===============       end        ==============="

  nlopt_destroy(opt.optimizer)
  result = FitResult(pStart: pStart,
                     pRes: paramsN,
                     pErr: newSeq[float](), ## Do not have classical errors. Could compute manually though.
                     χ²dof: χ²dof)

proc getMuSigma(lines: seq[FitFuncArgs],
                fitResult: FitResult
               ): (Measurement[float], Measurement[float]) =
  ##get the interesting fit params
  case lines[0].kind
  of ffGauss: # in a gaussian, parameters are `[0]: N, [1]: μ, [2]: σ`
    result = (fitResult.pRes[1] ± fitResult.pErr[1],
              fitResult.pRes[2] ± fitResult.pErr[2])
  of ffExpGauss: # in an exp gaussian, parameters are `[0]: a, [1]: b, [2]: N, [3]: μ, [4]: σ`
    result = (fitResult.pRes[3] ± fitResult.pErr[3],
              fitResult.pRes[4] ± fitResult.pErr[4])
  else:
    discard

proc getAmplitude(lines: seq[FitFuncArgs],
                  fitResult: FitResult
                 ): Measurement[float] =
  case lines[0].kind
  of ffGauss: # in a gaussian, parameters are `[0]: N, [1]: μ, [2]: σ`
    result = fitResult.pRes[0] ± fitResult.pErr[0]
  of ffExpGauss: # in an exp gaussian, parameters are `[0]: a, [1]: b, [2]: N, [3]: μ, [4]: σ`
    result = fitResult.pRes[2] ± fitResult.pErr[2]
  else:
    discard

proc getFitFuncLinesAndInfo(tfKind: TargetFilterKind, dKind: DataKind): (CdlFitFunc, seq[FitFuncArgs], float, string, string) =
  var binSize: float
  var xtitle: string
  var outname: string
  var fitFunc: CdlFitFunc
  var lines: seq[FitFuncArgs]
  let dummy = @[1.0, 2.0]
  case dKind
  of Dhits:
    fitFunc = getCdlFitFunc(tfKind)
    binSize = 1.0
    xtitle = "Number of pixels"
    outname = &"{tfKind}"
    lines = getLines(dummy, dummy, tfKind)
  of Dcharge:
    fitFunc = getCdlFitFuncCharge(tfKind)
    binSize = if ord(tfKind) >= ord(tfAlAl4): 5000.0 else: 10000.0
    xtitle = "Charge [e⁻]"
    outname = &"{tfKind}Charge"
    lines = getLinesCharge(dummy, dummy, tfKind)
  result = (fitFunc, lines, binSize, xtitle, outname)

proc fitAndPlotImpl(h5f: H5File, dfU: DataFrame, runNumber: int, fitParamsFname: string,
                    tfKind: TargetFilterKind, dKind: DataKind,
                    showStartParams, hideNloptFit: bool,
                    plotPath: string): MainPeak =
  # get the correct fit function, lines & plot settings
  let (fitFunc, lines, binSize, xtitle, outname) = getFitFuncLinesAndInfo(tfKind, dKind)
  # histogram of the CDL data to have data to fit to
  let (histdata, bins) = histoCdl(dfU.filter(f{idx("Cut?") == "Cut"})["Counts", float].toSeq1D, binSize, dKind)
  let fitData = getFitData(histdata, bins)
  # perform fits
  let fitResultMpfit = fitCdlMpfit(fitData, tfKind, dKind)
  let fitResultNlopt = fitCdlNlopt(fitData, tfKind, dKind)
  # and compute smooth fit curves
  let mpfitres = calcFitCurve(fitData, fitResultMpfit, fitfunc)
  let nloptres = calcFitCurve(fitData, fitResultNlopt, fitfunc)
  let startval = calcFitCurve(fitData.bins.min, fitData.bins.max, fitfunc, fitResultMpfit.pStart)
  var df = bind_rows([("MPFIT", mpfitres),
                      ("NLopt", nloptres),
                      ("Start", startval)],
                     id = "Type")
  if not showStartParams:
    df = df.filter(f{string: `Type` != "Start"})
  if hideNloptFit:
    df = df.filter(f{string: `Type` != "NLopt"})

  ## Extract μ and σ parameters incl errors. Need mpfit for errors at the moment
  let (fit_μ, fit_σ) = getMuSigma(lines, fitResultMpfit)
  # calculate energy resolution. Error propagated automatically.
  echo "FIT σ ", fit_σ, " and fit μ ", fit_μ, " / ", fit_σ / fit_μ
  # compute energy resolution & assign result
  let energyRes = fit_σ / fit_μ
  result = MainPeak(tfKind: tfKind, dKind: dKind, runNumber: runNumber,
                    fit_μ: fit_μ, fit_σ: fit_σ, energyRes: energyRes)

  # define range of plots using mean of main peak
  let binRangePlot = fit_μ.value * 3.0

  let transparent = color(0.0, 0.0, 0.0, 0.0)
  let fname = &"{plotPath}/{outname}-{outdate}_run_{runNumber}.pdf"
  # gather unbinned data for plot
  let cLow  = (fit_μ - (3 * fit_σ))
  let cHigh = (fit_μ + (3 * fit_σ))
  let cN = getAmplitude(lines, fitResultMpfit)
  let cLV = cLow.value; let cHV = cHigh.value; let cNV = cN.value
  let cLE = cLow.error; let cHE = cHigh.error
  var plt = ggplot(df, aes("Energy")) +
    geom_histogram(data = dfU.filter(f{float: `Counts` < binRangePlot}),
                   aes = aes("Counts", fill = "Cut?"),
                   position = "identity",
                   alpha = some(0.5),
                   hdKind = hdOutline,
                   binWidth = binSize) +
    # add error bar bands
    geom_tile(aes = aes(x = cLV - cLE/2.0, y = 0.0, width = cLE, height = cNV),
              alpha = 0.1, fillColor = "grey", color = transparent) +
    geom_tile(aes = aes(x = cHV - cHE/2.0, y = 0.0, width = cHE, height = cNV),
              alpha = 0.1, fillColor = "grey", color = transparent)
  if not showStartParams and hideNloptFit: # if only MPFIT, just draw black
    plt = plt + geom_line(aes(y = "Counts"))
  else:
    plt = plt + geom_line(aes(y = "Counts", color = "Type"))
  # continue rest of plot
  plt +
    # add lines of the lower / upper cut
    geom_linerange(aes = aes(x = cLV, y = cNV / 2.0, yMin = 0.0, yMax = cNV)) +
    geom_linerange(aes = aes(x = cHV, y = cNV / 2.0, yMin = 0.0, yMax = cNV)) +
    annotate(serializeFitParameters(fitResultMpfit, tfKind, dKind, true),
             0.55, 0.6,
             font = font(12.0, family = "monospace")) +
    xlab(xtitle) +
    xlim(0.0, binRangePlot) +
    ylab("Counts") +
    ggtitle(&"target: {tfKind}, run: {runNumber}") +
    ggsave(fname, width = 800, height = 480)

  # now dump the fit results, SVG filename and correct parameter names to a file
  dumpFitParameters(fitParamsFname, fname, fitResultMpfit, tfKind, dKind)

proc calcEnergyFromFits(df: DataFrame, mainPeak: MainPeak): DataFrame =
  ## Given the fit result of this data type & target/filter combination compute the energy
  ## of each cluster by using the mean position of the main peak and its known energy
  result = df
  result["Type"] = $(mainPeak.dKind)
  result["Target"] = $(mainPeak.tfKind)
  let invTab = getInverseXrayRefTable()
  let energies = getXrayFluorescenceLines()
  let lineEnergy = energies[invTab[$(mainPeak.tfKind)]]
  result = result.mutate(f{float: "Energy" ~ `Counts` / mainPeak.fit_μ.value * lineEnergy})

proc readCdlRunTfKind(h5f: H5File, grp: H5Group, run: int, tfKind: TargetFilterKind,
                      dKind: DataKind): DataFrame =
  case dKind
  of Dhits:
    let RawDataSeq = h5f[grp.name / "hits", int64]
    let Cdlseq = h5f.readCutCDL(run, "hits", tfKind, int64)
    result = bind_rows([("Raw", toDf({"Counts": RawDataSeq})),
                        ("Cut", toDf({"Counts": Cdlseq}))],
                       id = "Cut?")
    result["runNumber"] = run
  of Dcharge:
    let RawDataSeq = h5f[grp.name / "totalCharge", float]
    let Cdlseq = h5f.readCutCDL(run, "totalCharge", tfKind, float)
    result = bind_rows([("Raw", toDf({"Counts": RawDataSeq})),
                        ("Cut", toDf({"Counts": Cdlseq}))],
                       id = "Cut?")
    result["runNumber"] = run

iterator getCdlData(h5f: H5File, tfKind: TargetFilterKind, dKind: DataKind, fitByRun: bool): (int, DataFrame) =
  ## Yields the run number (if `fitByRun`) and the corresponding data frame. If `fitByRun` is
  ## `false` just yields `0` for run number.
  if fitByRun:
    var df = newDataFrame()
    for (run, grp) in tfRuns(h5f, tfKind):
      yield (run, h5f.readCdlRunTfKind(grp, run, tfKind, dKind))
  else:
    var df = newDataFrame()
    for (run, grp) in tfRuns(h5f, tfKind):
      df.add h5f.readCdlRunTfKind(grp, run, tfKind, dKind)
    yield (0, df)

proc fitAndPlot(h5f: H5File, fitParamsFname: string,
                tfKind: TargetFilterKind, dKind: DataKind,
                showStartParams, hideNloptFit, fitByRun: bool,
                plotPath: string):
                  (DataFrame, seq[MainPeak]) =
  # DataFrame for the unbinned data
  var dfU = newDataFrame()
  var mainPeaks = newSeq[MainPeak]()
  for (runNumber, dfLoc) in getCdlData(h5f, tfKind, dKind, fitByRun):
    let peak = h5f.fitAndPlotImpl(dfLoc, runNumber, fitParamsFname,
                                  tfKind, dKind, showStartParams, hideNloptFit,
                                  plotPath)
    # calibrate energy using `fit_μ` for unbinned data `dfLoc`
    dfU.add calcEnergyFromFits(dfLoc, peak)
    mainPeaks.add peak

  # first plot of only cut data by
  let (_, _, binSize, xtitle, outname) = getFitFuncLinesAndInfo(tfKind, dKind)
  let fnameByRun = &"{plotPath}/{outname}-{outdate}_by_run.pdf"
  let dfUCut = dfU.filter(f{idx("Cut?") == "Cut"})
  ggplot(dfUCut, aes("Counts", fill = factor("runNumber"))) +
    geom_histogram(binWidth = binSize,
                   hdKind = hdOutline,
                   position = "identity",
                   alpha = 0.5) +
    ggtitle("Cleaned data of & " & $dKind & " " & $tfKind & " split by run") +
    xlab(xtitle) +
    ggsave(fnameByRun, width = 800, height = 480)

  result = (dfU, mainPeaks)

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
                     chunksize = @[5000], # some hardcoded chunk size
                     maxshape = @[int.high])

template createAndWrite(h5read, h5write, `type`, dset, outname: untyped): untyped =
  let vlenType = special_type(`type`)
  let outDset = h5write.datasetCreation(outname, dset.shape, vlenType)
  outDset[outDset.all] = h5read[dset.name, vlenType, `type`]

template getAndAdd(h5read, h5write, dset, `type`, outDset: untyped): untyped =
  let vlenType = special_type(`type`)
  outDset.add h5read[dset.name, vlenType, `type`]

proc getFitParamsFname(dumpAccurate: bool): string =
  if dumpAccurate:
    result = "fitparams_accurate_" & $(epochTime().round.int) & ".txt"
  else:
    result = "fitparams_" & $(epochTime().round.int) & ".txt"

proc writeChargeCutBounds(h5f: H5File, peak: MainPeak,
                          tfKind: TargetFilterKind, year: YearKind, fitByRun: bool, runNumber: int) =
  let grpName = cdlGroupName($tfKind, $year, "", fitByRun, runNumber)
  let grp = h5f.getOrCreateGroup(grpName) # create or get
  ## Define the cut ranges around the main peak of ±3σ
  grp.attrs["ChargeLow"]  = (peak.fit_μ - (3 * peak.fit_σ)).value
  grp.attrs["ChargeHigh"] = (peak.fit_μ + (3 * peak.fit_σ)).value

proc generateCdlCalibrationFile(h5file: string, year: YearKind, fitByRun: bool,
                                outfile = "calibration-cdl") =
  ## generates the CDL calibration data file from a HDF5 file containing
  ## all CDL runs. Supports either 2014 CDL data or 2019 CDL data.

  # walk all runs corresponding to a single `TargetFilterKind` and
  # combine the datasets into the output files
  var h5f = H5open(h5file, "r")
  var h5fout = H5open(&"{outfile}-{year}.h5", "rw")
  let runs = readRuns(filename)
  let plotPath = h5f.attrs[PlotDirPrefixAttr, string]
  let fitParamsFname = getFitParamsFname(true)
  for tfKind in TargetFilterKind:
    var df = newDataFrame()
    var peak: MainPeak
    for (run, grp) in tfRuns(h5f, tfKind):
      # first get data for the CDL fit
      if fitByRun:
        df = h5f.readCdlRunTfKind(grp, run, tfKind, Dcharge)
        # fit directly, in case of `fitByRun = false` we fit _after_ this loop!
        peak = h5f.fitAndPlotImpl(df, run, fitParamsFname,
                                  tfKind, Dcharge, showStartParams = false, hideNloptFit = true,
                                  plotPath)
        # calibrate energy using `fit_μ` for unbinned data `df`
        ## XXX: write energy to output
        df = calcEnergyFromFits(df, peak)
        # write cuts for the charge based on fit
        h5fout.writeChargeCutBounds(peak, tfKind, year, fitByRun, run)
      else:
        df.add h5f.readCdlRunTfKind(grp, run, tfKind, Dcharge)

      # will not iterate the datasets and write to the outfile
      # via hyperslabs, potentially appending to the existing
      # dataset
      var mgrp = grp
      for dset in mgrp:
        if dset.shape.len > 1:
          # skip all datasets in the input, which are 2 dimensional. That includes
          # {polya, polyaFit, FeSpectrum, CDL...}
          # dset.shape == 1 is for `FeSpectrumCharge`, since it's the only dataset
          # atm that is actually 1D (see TPA issue #39).
          continue
        # TODO: add appropriate conversion of TPA naming to Marlin naming for
        # 2014 raw input or allow `genRefFile` option to deal with TPA naming
        # even for 2014 input!
        let outname = cdlGroupName($tfKind, $year, dset.name, fitByRun, run)
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
          of dkInt .. dkUint64: # includes all float
            let outDset = h5fout.datasetCreation(outname, dset.shape, float)
            outDset[outDset.all] = h5f.readAs(dset.name, float)
          else:
            echo "Skipping dataset: ", dset.name
            continue
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
          of dkInt .. dkUint64:
            outDset.add h5f.readAs(dset.name, float)
          else:
            continue
    if not fitByRun: ## If not fitting by run, then use the built DF to fit everything now
      peak = h5f.fitAndPlotImpl(df, runNumber = 0, fitParamsFname,
                                tfKind, Dcharge, showStartParams = false, hideNloptFit = true,
                                plotPath)
      # calibrate energy using `fit_μ` for unbinned data `dfLoc`
      ## XXX: write energy to output
      df = calcEnergyFromFits(df, peak)
      # write cuts for the charge based on fit
      h5fout.writeChargeCutBounds(peak, tfKind, year, fitByRun, 0)


  h5fout.attrs["FrameworkKind"] = $CdlGenerateNamingScheme
  h5fout.attrs["fitByRun"] = $fitByRun
  discard h5f.close()
  discard h5fout.close()

proc generateXrayReferenceFile(h5file: string, year: YearKind, fitByRun: bool,
                               outfile = "XrayReferenceFile") =
  ## generates the X-ray reference data file
  # this is achieved by taking the raw TargetFilterKind runs, combining them
  # into a the full CDL calibration file (if the input is the raw run based
  # file, else it skips the first step)
  # then we apply the charge cuts and bin the data by N bins
  # the result is written to the new file as (N, 2) datasets
  var h5f = H5open(h5file, "r")
  # read framework kind from `h5f`
  var frameworkKind = fkMarlin
  if "FrameworkKind" in h5f.attrs:
    frameworkKind = parseEnum[FrameworkKind](h5f.attrs["FrameworkKind", string])

  if "reconstruction" in h5f:
    # is a raw file, first create the CDL calibration file
    discard h5f.close()
    let cdlOut = "auto_calibration-cdl_" & $year
    generateCdlCalibrationFile(h5file, year, fitByRun, cdlOut)
    h5f = H5open(cdlOut, "r")

  # now walk all groups in root of h5f, read the datasets required for
  # charge cuts, write all passing indices back to file as binned
  # datasets
  var
    date: string
    tfKindStr: string
  const xrayRefTab = getXrayRefTable()
  var xrayRefCuts: OrderedTable[string, Cuts]
  case year
  of yr2014:
    xrayRefCuts = getEnergyBinMinMaxVals2014() #XraySpectrumCutVals()
  of yr2018:
    xrayRefCuts = getEnergyBinMinMaxVals2018()

  var h5fout = H5open(outfile & $year & ".h5", "rw")
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

proc readGasGain(h5f: H5File, tfKind: TargetFilterKind): (DataFrame, DataFrame) =
  var gainDf = newDataFrame()
  var tempDf = newDataFrame()
  for (run, grp) in tfRuns(h5f, tfKind):
    var gains = newSeq[float]()
    for gainSlice in iterGainSlicesFromDset(h5f, grp):
      gains.add gainSlice.G
    gainDf.add toDf({ "Gain" : gains,
                      "Target" : $tfKind,
                      "runNumber" : run })
    # check if this run has temperature data, if so read it too
    if grp.name.parentDir / "temperatures" in grp:
      let temps = h5f[grp.name.parentDir / "temperatures", TemperatureLogEntry]
      tempDf.add(toDf({ "Timestamp" : temps.mapIt(it.timestamp),
                        "IMB" : temps.mapIt(it.imb),
                        "runNumber" : run,
                        "Septem" : temps.mapIt(it.septem) })
        .mutate(f{"Time since run start" ~ `Timestamp` - min(col("Timestamp"))}))
  result = (gainDf, tempDf)

proc energyHistograms(df: DataFrame, plotPath: string) =
  let dfC = df.filter(f{`Type` == "charge" and idx("Cut?") == "Cut" and `Energy` < 10.0})
  let fname = "calibrated_cdl_energy"

  ggplot(dfC, aes("Energy", fill = "Target")) +
    geom_histogram(bins = 100, hdKind = hdOutline, density = true, position = "identity", alpha = 0.6, color = "black", lineWidth = 1.0) +
    ggtitle("Normalized histograms of energy calibrated CDL data by target / filter") +
    xlab("Energy [keV]") +
    ggsave(&"{plotPath}/{fname}_histos.pdf", width = 800, height = 480)

  ggplot(dfC, aes("Energy", fill = "Target")) +
    geom_density(alpha = 0.6, color = "black", size = 1.0, normalize = true) +
    ggtitle("Normalized KDE of energy calibrated CDL data by target / filter") +
    xlab("Energy [keV]") +
    ggsave(&"{plotPath}/{fname}_kde.pdf", width = 800, height = 480)

proc plotRunGains(gainDf, tempDf: DataFrame, plotPath: string) =
  ggplot(gainDf, aes("runNumber", "Gain", color = factor("Target"))) +
    geom_point() +
    ggsave(&"{plotPath}/gas_gain_by_run_and_tfkind.pdf", width = 800, height = 480)

  let minGain = gainDf["Gain", float].min
  let maxGain = gainDf["Gain", float].max
  let tempMean = tempDf.group_by("runNumber").summarize(f{float: "Temp" << mean(`Septem`)})
    .mutate(f{"TempNorm" ~ `Temp` / max(col("Temp")) * maxGain}) # (maxGain - minGain) *
  let maxTemp = tempMean["Temp", float].max

  ggplot(gainDf, aes("runNumber", "Gain", color = factor("Target"))) +
    geom_point() +
    geom_point(data = tempMean, aes = aes("runNumber", "TempNorm"), marker = mkCross) +
    scale_y_continuous(secAxis = secAxis(name = "Temp [°C]", trans = f{maxTemp / maxGain})) +
    margin(right = 7) +
    legendPosition(0.8, 0.0) +
    ggsave(&"{plotPath}/gas_gain_by_run_and_tfkind_with_cdl_temp.pdf", width = 800, height = 480)

  ggplot(tempDf, aes("Timestamp", "Septem", color = factor("runNumber"))) +
    geom_point(size = 0.5) +
    ggsave(&"{plotPath}/septem_temperature_cdl.pdf")

  let title = "Temperatures of CDL runs which contain temperature data"
  ggplot(tempDf, aes("Time since run start", "Septem", color = factor("runNumber"))) +
    facet_wrap("runNumber", scales = "free") +
    facetMargin(1.0) +
    geom_point(size = 0.5) +
    xlab(rotate = -45, alignTo = "right") +
    ggtitle(title) +
    ggsave(&"{plotPath}/septem_temperature_facet_cdl_time_since_start.pdf", width = 1200, height = 1000)

  proc toDate(v: float): string =
    result = v.int.fromUnix.format("dd-MM HH:mm")
  ggplot(tempDf, aes("Timestamp", "Septem", color = factor("runNumber"))) +
    facet_wrap("runNumber", scales = "free") +
    facetMargin(1.0) +
    geom_point(size = 0.5) +
    xlab(rotate = -45, alignTo = "right") +
    scale_x_continuous(labels = toDate) +
    ggtitle(title) +
    ggsave(&"{plotPath}/septem_temperature_facet_cdl.pdf", width = 1200, height = 1000)

  ggplot(tempDf.gather(["Septem", "IMB"], "Type", "Temp"), aes("Timestamp", "Temp", color = factor("Type"))) +
    facet_wrap("runNumber", scales = "free") +
    facetMargin(1.0) +
    geom_point(size = 0.5) +
    xlab(rotate = -45, alignTo = "right") +
    scale_x_continuous(labels = toDate) +
    ggtitle(title) +
    ggsave(&"{plotPath}/septem_imb_temperature_facet_cdl.pdf", width = 1200, height = 1000)

proc plotIngridProperties(h5f: H5File, tfKind: TargetFilterKind, plotPath: string) =
  var dfAll = newDataFrame()
  var dsets = newSeq[string]()
  for dset in InGridDsetKind:
    if dset in {igInvalid, igEnergyFromCharge, igEnergyFromPixel, igLikelihood, igNumClusters,
                igFractionInHalfRadius, igRadiusDivRmsTrans, igRadius, igBalance, igLengthDivRadius, igEventNumber}:
      continue
    var df = newDataFrame()
    let dsetStr = dset.toDset()
    for (run, grp) in tfRuns(h5f, tfKind):
      df.add toDf({ dsetStr : h5f.readCutCDL(run, dsetStr, tfKind, float), "run" : run })
    echo df
    let maxVal = df[dsetStr, float].max
    dfAll[dsetStr] = df[dsetStr, float].map_inline(x / maxVal)
    dfAll["run"] = df["run", int]
    dsets.add dsetStr
    ggplot(df, aes(dsetStr, fill = factor("run"))) +
      geom_histogram(bins = 100, hdKind = hdOutline, alpha = 0.5, position = "identity", density = true) +
      ggtitle("Dataset " & $dsetStr & " for " & $tfKind) +
      ggsave(&"{plotPath}/{$tfKind}_{dsetStr}_histogram_by_run.pdf")
    when false:
      ggplot(df, aes(dsetStr, fill = factor("run"))) +
        geom_density(normalize = true) +
        ggtitle("Dataset " & $dsetStr & " for " & $tfKind) +
        ggsave(&"{plotPath}/{$tfKind}_{dsetStr}_kde_by_run.pdf")
  dfAll = dfAll.gather(dsets, "Dset", "Value")
  ggplot(dfAll, aes("Value", fill = factor("run"))) +
    ggridges("Dset", overlap = 4.0) +
    geom_density(normalize = true, alpha = 0.6, color = "black", size = 0.5) +
    ggtitle("Ridgeline plot of " & $tfKind & " for all InGrid properties in arbitrary units: x/max(x)") +
    margin(left = 4) +
    ggsave(&"{plotPath}/{$tfKind}_ridgeline_kde_by_run.pdf", width = 900, height = 600)

proc plotsAndEnergyResolution(input: string,
                              dumpAccurate, showStartParams, hideNloptFit, fitByRun: bool) =
  let fitParamsFname = getFitParamsFname(dumpAccurate)
  var peaksHits: seq[MainPeak]
  var peaksCharge: seq[MainPeak]
  var h5f = H5open(input, "rw")
  var df = newDataFrame()
  let plotPath = h5f.attrs[PlotDirPrefixAttr, string]
  var gainDf = newDataFrame()
  var tempDf = newDataFrame()
  for tfkind in TargetFilterKind:
    echo "\n\n------------------------------ Working on : ", tfKind, " ----------------------------------------"
    let (dfH, pHits) = fitAndPlot(h5f, fitParamsFname, tfkind, Dhits,
                                  showStartParams, hideNloptFit, fitByRun, plotPath)
    df.add dfH
    peaksHits.add pHits
    echo "and now charge:"
    let (dfC, pCharge) = fitAndPlot(h5f, fitParamsFname, tfkind, Dcharge,
                                    showStartParams, hideNloptFit, fitByRun, plotPath)
    peaksCharge.add pCharge
    df.add dfC
    #if tfKind == tfCuEpic2:
    #  quit()

    let (gDf, tDf) = readGasGain(h5f, tfKind)
    gainDf.add gDf
    tempDf.add tDf

    # plot the ingrid properties by run
    h5f.plotIngridProperties(tfKind, plotPath)


  # plot the gas gain of all data
  plotRunGains(gainDf, tempDf, plotPath)

  discard h5f.close()
  energyResolution(peaksHits, peaksCharge, plotPath)
  peakFit(peaksHits, "Hits", plotPath)
  peakFit(peaksCharge, "Charge", plotPath)

  # finally make a histogram of all data containing the energies of the (raw/cut) clusters
  # 'calibrated' by the charge of the known main peak energy.
  energyHistograms(df, plotPath)

import std / [parseutils, macrocache]
const FuncsTeX = CacheTable"TeXFunctions"
proc getStr(n: NimNode): string = n.strVal
proc upcaseWord(s: string): string =
  result = s
  result[0] = s[0].toUpperAscii()

proc fnNameToTfKind(fn: string): TargetFilterKind =
  # not the prettiest parsing, but gets the job done without too much
  # specifics of the details...
  var tmp = fn.replace("_", ".").removeSuffix("Func").removeSuffix("Charge")
  var target = ""
  var filter = ""
  var idx = tmp.parseUntil(target, until = {'A' .. 'Z'}, start = 0) # skip first char!
  target = target.upcaseWord()
  idx += tmp.parseUntil(filter, until = {'0' .. '9'}, start = idx)
  let kV = tmp[idx .. ^1] & "kV"
  let tfKindStr = @[target, filter, kV].join("-")
  result = parseEnum[TargetFilterKind](tfKindStr)

proc getFuncs(charge = true): string =
  var arr: array[TargetFilterKind, string]
  for k, v in FuncsTeX:
    try:
      let tfKind = k.fnNameToTfKind()
      if charge and "Charge" in k:
        arr[tfKind] = v.getStr()
      elif not charge and "Charge" notin k:
        arr[tfKind] = v.getStr()
    except ValueError:
      echo "Cannot parse " & $k & " into a TargetFilterKind. Skipping it."
      # expect this for any other defined func not for CDL, e.g. feSpectrum(Charge)
  result = "| Target | Filter | HV [kV] | Fit functions |\n"
  for k, x in arr:
    result.add "| " & [$k.toTarget(), $k.toFilter(), $k.toHV(false), x].join("|") & "|\n"

proc raiseNoInput(input: string) =
  if input.len == 0:
    raise newException(OSError, "In order to process data we need an input, but none given.")

proc main(input = "",
          cutcdl = false,
          dumpAccurate = false,
          genCdlFile = false,
          genRefFile = false,
          showStartParams = false,
          hideNloptFit = false,
          outfile = "",
          year: string = "2018",
          config = "",
          printFunctions = false
         ) =
  ##
  docCommentAdd(versionStr)
  let year = parseEnum[YearKind](year)
  let fitByRun = readFitByRun(config)

  template callGenFunc(fn: untyped): untyped =
    raiseNoInput(input)
    if outfile.len > 0:
      fn(input, year, fitByRun, outfile)
    else:
      fn(input, year, fitByRun)
  if genRefFile:
    callGenFunc(generateXrayReferenceFile)
  if genCdlFile:
    callGenFunc(generateCdlCalibrationFile)
  if cutcdl:
    raiseNoInput(input)
    cutAndWrite(input)

  if not genRefFile and not genCdlFile and not printFunctions:
    # only perform CDL fits if neither CDL calibration file nor
    # reference file created
    raiseNoInput(input)
    plotsAndEnergyResolution(input, dumpAccurate, showStartParams, hideNloptFit, fitByRun)

  if printFunctions:
    echo "Charge functions:"
    const chargeFns = getFuncs()
    const pixelFns = getFuncs(charge = false)
    echo chargeFns
    echo "Pixel functions:"
    echo pixelFns


when isMainModule:
  import cligen
  dispatch main
