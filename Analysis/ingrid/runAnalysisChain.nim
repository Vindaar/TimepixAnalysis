import shell
import std / [strutils, strformat, os, options, algorithm, sequtils, math]
from sequtils import toSeq
import helpers / utils
import logging

const docStr = """
Usage:
  runAnalysisChain <dataPath> (--2014 | --2017 | --2018) [--noBack --noCalib --noRaw --noReco] [options]

Options:
  --2014       Run 2014 analysis chain
  --2017       Run 2017 analysis chain
  --2018       Run 2018_2 analysis chain
  --noBack     Do not perform analysis chain on background runs
  --noCalib    Do not perform analysis chain on calibration runs
  --noRaw      Do not perform raw data manipulation
  --noReco     Do not perform reconstruction
  -h, --help   Show this help
  --version    Show the version number
"""

from ingrid / ingrid_types import RunTypeKind
from reconstruction import RecoFlags, RecoConfig, initRecoConfig

const subPaths = ["2014_15", "2017", "2018_2"]
const recoOptions = @[{rfNone},
                      {rfOnlyFadc},
                      {rfOnlyCharge},
                      {rfOnlyGasGain},
                      {rfOnlyGainFit},
                      {rfOnlyEnergyElectrons}]
const relIDPath = "../../InGridDatabase/src/resources"

type
  AnaFlags = enum
    afBack, afCalib, afRaw, afReco

  InputKind = enum
    ikFileH5, # refers to a single H5 input file
    ikRunFolder, # refers to a folder of a run (possibly containing multiple runs)
    ikDataFolder # refers to a folder containing directories required for 2017/2018
                 # CAST data analysis (sub dirs for 2017, 2018_2 with each calibration & background)

  DataYear = enum
    dy2014 = 2014
    dy2017 = 2017
    dy2018 = 2018

  # let's test those juicy default values
  Config = object ## Configuration for the analysis
    input: string # the input file / directory
    path: string  # the path of the input. If input is a file `path` is the parent dir. Else it's `input`
    inputKind: InputKind
    anaFlags: set[AnaFlags] = {}
    runNumber: Option[int]
    recoCfg: RecoConfig
    runType: RunTypeKind = rtNone
    outpath = "" # the output path in which the generated H5 files will be placed
    outfile = "" # ability to override output file setitngs. Using this file exactly,
                 # placed into `outpath` if given
    # prefixes for the different resulting files
    tpx3Prefix = "tpx3" # the prefix of the output file after `readTpx3`
    # the following
    rawStep = "raw" # the prefix of the output file after `raw_data_manipulation`
    recoStep = "reco" # the prefix of the output file after `reconstruction`
    # the type of data
    calibType = "calibration"
    backType = "data"
    # whether multiple or single run
    singleRun = "run"
    mulitpleRun = "runs"
    # template to use for files with multiple runs
    nameTmpl = "$type$runs$year_$step.h5"

# set up the logger
var L = newConsoleLogger()
if not dirExists("logs"):
  createDir("logs")
var fL = newFileLogger("logs/runAnalysisChain.log", fmtStr = verboseFmtStr)
when isMainModule:
  addHandler(L)
  addHandler(fL)

proc toStr(rt: RunTypeKind): string =
  case rt
  of rtNone: "none"
  of rtBackground: "data"
  of rtCalibration: "calibration"
  of rtXrayFinger: "xray"

proc toRunStr(cfg: Config): string =
  result = if cfg.runNumber.isSome: "Run_" & $cfg.runNumber.get
           else: "Runs"

## XXX: fix me: use run type or use `afCalib/afBack`? Differentiate CAST & else
proc toOutfile(cfg: Config, runType: RunTypeKind, flag: AnaFlags, year: string): string =
  doAssert flag in {afRaw, afReco}
  let step = if flag == afRaw: "raw"
             else: "reco"
  result = cfg.nameTmpl % ["type", runType.toStr.capitalizeAscii,
                           "runs", cfg.toRunStr,
                           "year", year,
                           "step", step.capitalizeAscii]

proc initConfig(input: string,
                outpath: string,
                inputKind: InputKind,
                recoCfg: RecoConfig,
                runType: RunTypeKind): Config =
  let path = if inputKind == ikFileH5: input.parentDir()
             else: input
  let outpath = if outpath.len == 0: path
                else: outpath
  result = Config(input: input, path: path, inputKind: inputKind,
                  outpath: outpath,
                  recoCfg: recoCfg, runType: runType,
                  runNumber: recoCfg.runNumber)

proc toDataPath(year: DataYear): string =
  case year
  of dy2014: result = "2014_15"
  of dy2017: result = "2017"
  of dy2018: result = "2018_2"

proc rawDataManipulation(path, runType, outName: string): bool =
  let runTypeStr = &"--runType={runType}"
  let outStr = &"--out={outName}"
  info "Running raw_data_manipulation on " & $path & $runTypeStr & $outStr
  let res = shellVerbose:
    ./raw_data_manipulation ($path) ($runTypeStr) ($outStr)
  #info "Last commands output: " & $res[0]
  info "Last commands exit code: " & $res[1]
  result = res[1] == 0

proc toStr(opt: set[RecoFlags]): string = opt.toSeq.join(" ")

proc reconstruction(rawName: string, recoName: string, options: set[RecoFlags],
                    cfg: Config): bool =
  let option = options.toStr()
  let runNumber = if cfg.recoCfg.runNumber.isSome: "--runNumber " & $cfg.recoCfg.runNumber.get
                  else: ""
  info "Running reconstruction on " & $rawName & $option
  var res: (string, int)
  if option.len > 0: # has any options, implying not raw->reco step
    res = shellVerbose:
      ./reconstruction ($recoName) ($runNumber) ($option)
  else:
    res = shellVerbose:
      ./reconstruction ($rawName) "--out" ($recoName) ($runNumber)
  #info "Last commands output: " & $res[0]
  info "Last commands exit code: " & $res[1]
  result = res[1] == 0

func likelihood() = discard

template tc(cmd: untyped): untyped {.dirty.} =
  ## convenience template to wrap a command in toContinue checks
  if toContinue:
    toContinue = cmd

proc runCastChain(dYear: DataYear, cfg: Config): bool =
  ## performs the whole chain of the given dataset
  let dataPath = path / (yearKind.toDataPath())
  let outpath = if cfg.outpath.len > 0: cfg.outpath
                else: dataPath
  var toContinue = true

  ## XXX: Make use of file name template!
  # raw data for calibration
  let calibOut = outPath / cfg.toOutfile(rtCalibration, afRaw, $dYear)
  if afCalib in cfg.anaFlags and afRaw in cfg.anaFlags:
    tc(rawDataManipulation(dataPath / "CalibrationRuns", "calib",
                           calibOut))
  # raw data for background
  let backOut = outPath / cfg.toOutfile(rtBackground, afRaw, $dYear)
  if afBack in cfg.anaFlags and afRaw in cfg.anaFlags:
    tc(rawDataManipulation(dataPath / "DataRuns", "back",
                           backOut))

  # if any reco flags given, use only those instead of our predefined set
  let recoOptions = if cfg.recoCfg.flags.card > 0: cfg.recoCfg.flags.toSeq.sorted.mapIt({it})
                    else: recoOptions
  let recoCalibOut = cfg.path / cfg.toOutfile(rtCalibration, afReco, $dYear)
  let recoBackOut = cfg.path / cfg.toOutfile(rtBackground, afReco, $dYear)
  for opt in recoOptions:
    if afCalib in cfg.anaFlags and afReco in cfg.anaFlags:
      tc(reconstruction(calibOut,
                        recoCalibOut,
                        opt, cfg))
    if rfOnlyGainFit notin opt and afBack in cfg.anaFlags and afReco in cfg.anaFlags:
      tc(reconstruction(backOut,
                        recoBackOut,
                        opt, cfg))
  result = toContinue

proc runChain(cfg: Config) =
  ## Performs analysis on the given input. Either a directory (possibly with subs)
  var toContinue = true

  ##
  if afRaw in cfg.anaFlags:
    tc(rawDataManipulation(
      cfg.path, $cfg.runType,
      cfg.outpath / cfg.toOutfile(cfg.runType, afRaw, "")) # no explicit year
    )

  # if any reco flags given, use only those instead of our predefined set
  let recoOptions = if cfg.recoCfg.flags.card > 0: cfg.recoCfg.flags.toSeq.sorted.mapIt({it})
                    else: recoOptions.filterIt(rfOnlyGainFit notin it) # gain fit does not make sense in single input case

  if afReco in cfg.anaFlags:
    tc(reconstruction(cfg.path / cfg.toOutfile(cfg.runType, afRaw, ""), # no explicit year
                      cfg.path / cfg.toOutfile(cfg.runType, afReco, ""),
                      {},
                      cfg))

proc main(
  path: string,
  outpath = "",
  years: seq[int] = @[],
  inputKind: InputKind = ikDataFolder,
  back = false,
  calib = false,
  raw = false,
  reco = false,
  tpx3 = false,
  runType = rtNone,
  runNumber = -1,
  create_fe_spec = false,
  only_fadc = false,
  only_fe_spec = false,
  only_charge = false,
  only_gas_gain = false,
  only_gain_fit = false,
  only_energy_from_e = false,
  only_energy = NaN,
  config = ""
     ) =
  var toContinue = true

  ## Parse full analysis chain related parameters
  let anaFlags = block:
    var flags: set[AnaFlags]
    if back:
      flags.incl afBack
    if calib:
      flags.incl afCalib
    if raw:
      flags.incl afRaw
    if reco:
      flags.incl afReco
    flags

  ## Parse `reconstruction` related parameters
  let recoCfg = block:
    var
      recoFlags: set[RecoFlags]
      calibFactor = NaN
      runNumberArg: Option[int]
    if runNumber < 0:
      recoFlags.incl rfReadAllRuns
    else:
      runNumberArg = some(runNumber)
    if classify(only_energy) != fcNaN:
      recoFlags.incl rfOnlyEnergy
      calibFactor = only_energy
    if create_fe_spec:
      recoFlags.incl rfCreateFe
    if onlyCharge:
      recoFlags.incl rfOnlyCharge
    if onlyFadc:
      recoFlags.incl rfOnlyFadc
    if onlyFeSpec:
      recoFlags.incl rfOnlyFeSpec
    if onlyGasGain:
      recoFlags.incl rfOnlyGasGain
    if onlyGainFit:
      recoFlags.incl rfOnlyGainFit
    if onlyEnergyFromE:
      recoFlags.incl rfOnlyEnergyElectrons
    initRecoConfig(recoFlags, runNumberArg, calibFactor)

  let inputKind = if path.endsWith(".h5"): ikFileH5 else: inputKind
  var cfg = initConfig(path, outpath, inputKind, recoCfg, runType)


  ## XXX: if the input is a H5 file open it and determine what type of
  ## file it is (raw, reco) and thus adjust what is the "input" in the
  ## following procedures. E.g. if it's `raw` then we need to use that
  ## as input to `reconstruction`.

  ## XXX: now differentiate between cases of:
  ## 1. full reconstruction of all runs from raw data
  ## 2. continuation of a "full run" setup based on a H5 file
  ## 3. reconstruction from a _single_ run folder
  ## 4. continuation of a single run file

  case cfg.inputKind
  of ikDataFolder:
    for year in years:
      doAssert year.uint16 in {2014, 2017, 2018}, "Years supported are 2014, 2017, 2018."
      let yearKind = DataYear(year)
      tc(runCastChain(yearKind, cfg))
  else:
    doAssert cfg.runType != rtNone, "For individual runs a `runType` is required."
    runChain(cfg)

when isMainModule:
  import cligen
  # we could do `multiDispatch`, but the majority of arguments apply to all
  dispatch(main, help = {
    "path"           : "The path containing the data to be processed.",
    "inputKind"      : """Type of input given. If not given assumes full CAST data, but overwritten
to ikFileH5 if explicit H5 file given. For a single run handing this is required.""",
    "years"          : "If input is full CAST data (ikDataFolder) the datasets to reconstruct.",
    "outpath"        : "The output path in which all files will be placed. If none given input path is used.",
    "back"           : "(CAST) If set perform analysis of background data. Only relevant for full CAST data.",
    "calib"          : "(CAST) If set perform analysis of calibration data. Only relevant for full CAST data.",
    "raw"            : "If set performs parsing of data and storing in HDF5.",
    "reco"           : "If set performs geometric clustering and calculation of geometric properties.",
    "tpx3"           : "If set indicates input is a Tpx3 based data file.",
#    "plot"           : "If set generates all plots as described desired via `plotDataSuffix`.",
#    "outpath"        : "The path in which all generated H5 files will be stored.",
    "config"         : "The path to the config file for this tool.",
    "runType"        : "The run type of the input files {rtCalibration : ⁵⁵Fe, rtBackground : background data, rtXrayFinger : X-ray finger run}.",
    "runNumber"          : "Only work on this run",
    "create_fe_spec"     : """(reco) Toggle to create Fe calibration spectrum based on cuts
Takes precedence over --calib_energy if set!""",
    "only_fadc"          : """(reco) If this flag is set, the reconstructed FADC data is used to calculate
FADC values such as rise and fall times among others, which are written
to the H5 file.""",
    "only_fe_spec"       : """(reco) Toggle to /only/ create the Fe spectrum for this run and perform the
fit of it. Will try to perform a charge calibration, if possible.""",
    "only_charge"        : """(reco) Toggle to /only/ calculate the charge for each TOT value based on
the TOT calibration. The `ingridDatabase.h5` needs to be present.""",
    "only_gas_gain"      : """(reco) Toggle to /only/ calculate the gas gain for the runs in the
input file based on the polya fits to time slices defined by `gasGainInterval`. `ingridDatabase.h5`
needs to be present.""",
    "only_gain_fit"      : """(reco) Toggle to /only/ calculate the fit mapping the energy calibration
factors of the 55Fe runs to the gas gain values for each time slice. Required to calculate the
energy in any run using `only_energy_from_e`.""",
    "only_energy_from_e" : """(reco) Toggle to /only/ calculate the energy for each cluster based on
the Fe charge spectrum vs gas gain calibration""",
    "only_energy"        : """(reco) Toggle to /only/ perform energy calibration using the given factor.
Takes precedence over --create_fe_spec if set.
If no runNumber is given, performs energy calibration on all runs
in the HDF5 file.""",

  })
