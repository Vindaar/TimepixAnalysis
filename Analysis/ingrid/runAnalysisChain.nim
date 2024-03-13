import shell
import std / [strutils, strformat, options, algorithm, sequtils, math, times]
import std / os except FileInfo, getFileInfo
from sequtils import toSeq
import helpers / utils
import logging

import projectDefs

import ingrid / [ingrid_types]
import ingrid / private / [hdf5_utils, pure]
from reconstruction import RecoConfig, initRecoConfig

const subPaths = ["2014_15", "2017", "2018_2"]
const recoOptions = @[{rfOnlyFadc},
                      {rfOnlyCharge},
                      {rfOnlyFeSpec}, # for FADC 55Fe spectra!
                      {rfOnlyGasGain},
                      {rfOnlyGainFit},
                      {rfOnlyEnergyElectrons}]
const relIDPath = TpxDir / "InGridDatabase/src/resources"

const TrackingLogs = TpxDir / "resources/LogFiles/tracking-logs"

type
  AnaFlags = enum
    afBack, afCalib, afCDL, afRaw, afReco, afLogL, afTracking

  InputKind = enum
    ikFileH5, # refers to a single H5 input file
    ikRunFolder, # refers to a folder of a run (possibly containing multiple runs)
    ikDataFolder # refers to a folder containing directories required for 2017/2018
                 # CAST data analysis (sub dirs for 2017, 2018_2 with each calibration & background)

  StepKind = enum
    skRaw  # run folder / data folder (but data folder is _not_ equivalent. can contain finished H5 files!)
    skReco # ikFileH5 as input

  DataYear = enum
    dy2014 = 2014
    dy2017 = 2017
    dy2018 = 2018
    dy2019 = 2019 # for CDL data

  # let's test those juicy default values
  Config = object ## Configuration for the analysis
    input: string # the input file / directory
    path: string  # the path of the input. If input is a file `path` is the parent dir. Else it's `input`
    trackingLogs: string = TrackingLogs # path to the tracking files
    case inputKind: InputKind
    of ikFileH5: fileInfo: FileInfo
    else: discard
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
    nameTmpl = r"$type$runs$year_$step.h5"
    # default names for data, calib and CDL file
    dataRuns = "DataRuns"
    calibrationRuns = "CalibrationRuns"
    calibCdlFilename = "calibration-cdl-2018.h5"

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
  result = if cfg.runNumber.isSome: "Run_" & $cfg.runNumber.get & "_"
           else: "Runs"

## XXX: fix me: use run type or use `afCalib/afBack`? Differentiate CAST & else
proc toOutfile(cfg: Config, runType: RunTypeKind, flag: AnaFlags, year: string): string =
  doAssert flag in {afRaw, afReco}
  let step = if flag == afRaw: "raw"
             else: "reco"
  let typ = if flag == afCDL: "CDL_"
            else: runType.toStr.capitalizeAscii
  result = cfg.nameTmpl % ["type", typ,
                           "runs", cfg.toRunStr,
                           "year", year & "_", ## XXX: the `_` is dropped for some reason
                           "step", step.capitalizeAscii]

proc initConfig(input: string,
                outpath: string,
                inputKind: InputKind,
                recoCfg: RecoConfig,
                runType: RunTypeKind,
                anaFlags: set[AnaFlags],
                trackingLogs: string): Config =
  let path = if inputKind == ikFileH5: input.parentDir()
             else: input
  let outpath = if outpath.len == 0: path
                else: outpath
  result = Config(input: input, path: path,
                  inputKind: inputKind,
                  outpath: outpath,
                  recoCfg: recoCfg, runType: runType,
                  runNumber: recoCfg.runNumber,
                  anaFlags: anaFlags)
  if trackingLogs.len > 0:
    result.trackingLogs = trackingLogs

  case inputKind
  of ikFileH5:
    result.fileInfo = getFileInfo(input)
    # if only single run in file, set that as run number to get correct naming
    if result.runNumber.isNone and result.fileInfo.runs.len == 1:
      result.runNumber = some(result.fileInfo.runs[0])
    result.runType = result.fileInfo.runType
  of ikRunFolder:
    let (is_rf, runNumber, _, _) = isTosRunFolder(input)
    doAssert is_rf, "Input is not a run folder! Input: " & $input
    result.runNumber = some(runNumber)
  else:
    discard

proc toDataPath(year: DataYear): string =
  case year
  of dy2014: result = "2014_15"
  of dy2017: result = "2017"
  of dy2018: result = "2018_2"
  of dy2019: result = "CDL_2019"

proc rawDataManipulation(path, outName: string, runType: RunTypeKind): bool =
  let runTypeStr = &"--runType={runType}"
  let outStr = &"--out={outName}"
  info "Running raw_data_manipulation on " & $path & $runTypeStr & $outStr
  let res = shellVerbose:
    raw_data_manipulation -p ($path) ($runTypeStr) ($outStr)
  #info "Last commands output: " & $res[0]
  info "Last commands exit code: " & $res[1]
  result = res[1] == 0

proc toStr(opt: set[RecoFlags], overwrite: bool): string =
  result = opt.toSeq.join(" ")
  if overwrite:
    result.add " --overwrite"

proc reconstruction(input: string, options: set[RecoFlags],
                    cfg: Config,
                    output: string = ""): bool =
  let option = options.toStr(cfg.recoCfg.overwrite)
  let runNumber = if cfg.recoCfg.runNumber.isSome: "--runNumber " & $cfg.recoCfg.runNumber.get
                  else: ""
  info "Running reconstruction on " & $input & $option
  var res: (string, int)
  if option.len > 0: # has any options, implying not raw->reco step
    res = shellVerbose:
      reconstruction -i ($input) ($runNumber) ($option)
  elif output.len > 0:
    res = shellVerbose:
      reconstruction -i ($input) "--out" ($output) ($runNumber)
  #info "Last commands output: " & $res[0]
  info "Last commands exit code: " & $res[1]
  result = res[1] == 0

proc cdl_spectrum_creation(input: string): bool =
  info "Running cdl_spectrum_creation on " & $input & " to produce `calibration-cdl-2018.h5` file."
  var res: (string, int)
  res = shellVerbose:
    cdl_spectrum_creation ($input) "--genCdlFile --year=2018"
  info "Last commands exit code: " & $res[1]
  result = res[1] == 0

proc likelihood(input: string,
                cdlFile: string,
                cfg: Config,
                output: string = ""): bool =
  info "Running likelihood on " & $input & " to compute `likelihood` dataset"
  ## XXX: make adjustable
  let cdlData = &"--cdlYear 2018 --cdlFile {cdlFile}"
  let logLOpt = "--computeLogL"
  var res: (string, int)
  res = shellVerbose:
    likelihood -f ($input) ($cdlData) ($logLOpt)
  #info "Last commands output: " & $res[0]
  info "Last commands exit code: " & $res[1]
  result = res[1] == 0

proc trackingInfo(input: string,
                  cfg: Config,
                  output: string = ""): bool =
  info "Running likelihood on " & $input & " to compute `likelihood` dataset"
  var res: (string, int)
  # Note: to be generic we include the whole data taking campaign of CAST (and then some)
  # to make sure we don't lose tracking info on either end
  let logs = cfg.trackingLogs
  res = shellVerbose:
    cast_log_reader tracking "-p" ($logs) --startTime "2017/01/01" --endTime "2018/12/31" --h5out ($input)
  #info "Last commands output: " & $res[0]
  info "Last commands exit code: " & $res[1]
  result = res[1] == 0

template tc(cmd: untyped): untyped {.dirty.} =
  ## convenience template to wrap a command in toContinue checks
  if toContinue:
    toContinue = cmd

type
  Files = object
    rawPath: string # path to all runs of this kind
    raw: string     # path to the `_raw.h5` file
    reco: string    # path to the `_reco.h5` file

  DataFiles = object
    background: Files
    calibration: Files
    cdl: string # path to the `calibration-cdl-2018.h5` file

proc fillFiles(cfg: Config, path, outpath: string, runType: RunTypeKind, dYear: DataYear): Files =
  let rawPath = case runType
                of rtBackground: path / cfg.dataRuns
                of rtXrayFinger: path # CDL runs are inside `CDL_2019` path
                else: path / cfg.calibrationRuns
  result = Files(rawPath: path,
                 raw: outpath / cfg.toOutfile(runType, afRaw, $(ord(dYear))),
                 reco: outpath / cfg.toOutfile(runType, afReco, $(ord(dYear))))

proc fillDataFiles(cfg: Config, dYear: DataYear): DataFiles =
  let dataPath = cfg.path / (dYear.toDataPath())
  let outpath = if cfg.outpath.len > 0: cfg.outpath
                else: dataPath
  result = DataFiles(background: cfg.fillFiles(dataPath, outpath, rtBackground, dYear),
                     calibration: cfg.fillFiles(dataPath, outpath, rtCalibration, dYear),
                     cdl: cfg.path / (dy2019.toDataPath()) / cfg.calibCdlFilename)

proc fillCdlFiles(cfg: Config): Files =
  let dataPath = cfg.path / (dy2019.toDataPath)
  result = cfg.fillFiles(dataPath, dataPath, rtXrayFinger, dy2019)

proc runCastChain(cfg: Config, data: DataFiles): bool =
  ## performs the whole chain of the given dataset
  var toContinue = true

  # raw data for calibration
  if afCalib in cfg.anaFlags and afRaw in cfg.anaFlags:
    tc(rawDataManipulation(data.calibration.rawPath, data.calibration.raw, rtCalibration))
  # raw data for background
  if afBack in cfg.anaFlags and afRaw in cfg.anaFlags:
    tc(rawDataManipulation(data.background.rawPath, data.background.raw, rtBackground))

  # if any reco flags given, use only those instead of our predefined set
  let recoOptions = if cfg.recoCfg.flags.card > 0: cfg.recoCfg.flags.toSeq.sorted.mapIt({it})
                    else: recoOptions
  if afCalib in cfg.anaFlags and afReco in cfg.anaFlags:
    tc(reconstruction(data.calibration.raw,
                      {}, cfg,
                      output = data.calibration.reco
    ))
  if afBack in cfg.anaFlags and afReco in cfg.anaFlags:
    tc(reconstruction(data.background.raw,
                      {}, cfg,
                      output = data.background.reco
    ))
  # else apply reco options
  # if `raw` in flags but reco not we don't want to run options! Means input was a folder or
  # raw data file, but we didn't reconstruct it.
  ## XXX: check using tpa file kind
  if afReco in cfg.anaFlags:
    for opt in recoOptions:
      if afCalib in cfg.anaFlags:
        tc(reconstruction(data.calibration.reco, opt, cfg))
      if afBack in cfg.anaFlags:
        tc(reconstruction(data.background.reco, opt, cfg))

  ## XXX: We should rerun `--create_fe_spec` for the FADC spectrum fit. Because we only have
  ## `minVal` after the FADC reco has run nowadays. This means the default Fe fitting does not
  ## perform the fit anymore.
  if afReco in cfg.anaFlags:
    if {rfOnlyFadc} in recoOptions and afCalib in cfg.anaFlags:
      tc(reconstruction(data.calibration.reco, {rfOnlyFeSpec}, cfg))

  ## XXX: Technically for this step we also have to perform the generation of the CDL file
  ## in the first place, which is not part of this!

  ## As a final step compute the likelihood dataset in the file and add the tracking information
  if afLogL in cfg.anaFlags:
    tc(likelihood(data.calibration.reco, data.cdl, cfg))
    tc(likelihood(data.background.reco, data.cdl, cfg))

  ## And now the tracking information
  if afTracking in cfg.anaFlags:
    tc(trackingInfo(data.background.reco, cfg))

  result = toContinue

proc runCdlChain(cfg: Config, data: Files): bool =
  ## Performs the analysis of the CAST CDL 2019 data and produces the `calibration-cdl-2018.h5` file.
  var toContinue = true

  if afRaw in cfg.anaFlags:
    tc(rawDataManipulation(data.rawPath, data.raw, rtXrayFinger))
  if afReco in cfg.anaFlags:
    tc(reconstruction(data.raw,
                      {}, cfg,
                      output = data.reco
    ))
    const recoOptions = @[{rfOnlyFadc},
                          {rfOnlyCharge},
                          {rfOnlyGasGain}]
    for opt in recoOptions:
      tc(reconstruction(data.reco, opt, cfg))

  # now produce the cdl file
  tc(cdl_spectrum_creation(data.reco))
  result = toContinue

proc runChain(cfg: Config) =
  ## Performs analysis on the given input. Either a directory (possibly with subs)
  var toContinue = true

  if cfg.inputKind == ikRunFolder and afRaw notin cfg.anaFlags:
    raise newException(Exception, "Input does not do anything. Input is a run folder, but `--raw` " &
      "flag not provided. Cannot continue after.")

  var input = cfg.input

  ## dispatch to a generic procedure? This same thing here with different args? or whatever
  if cfg.inputKind == ikRunFolder: # skip if H5 input
    let output = cfg.outpath / cfg.toOutfile(cfg.runType, afRaw, "") # no explicit year
    if afRaw in cfg.anaFlags:
      tc(rawDataManipulation(
        input, output, cfg.runType)
      )
    input = output

  if cfg.inputKind in {ikRunFolder, ikFileH5}:
    # if any reco flags given, use only those instead of our predefined set
    let recoOptions = if cfg.recoCfg.flags.card > 0: cfg.recoCfg.flags.toSeq.sorted.mapIt({it})
                      else: recoOptions.filterIt(rfOnlyGainFit notin it) # gain fit does not make sense in single input case

    let output = cfg.outpath / cfg.toOutfile(cfg.runType, afReco, "")
    if afReco in cfg.anaFlags and
      (cfg.inputKind == ikRunFolder or
       (cfg.inputKind == ikFileH5 and cfg.fileInfo.tpaFileKind == tpkRawData)):
      tc(reconstruction(input, # no explicit year
                        {},
                        cfg,
                        output = output))
      input = output
    elif afReco in cfg.anaFlags and
       cfg.inputKind == ikFileH5 and cfg.fileInfo.tpaFileKind != tpkRawData:
      echo "WARNING: Gave `--reco` flag, but input: ", cfg.input, " is ", cfg.fileInfo.tpaFileKind
      echo "\tReconstruction step will be ignored."
    # else apply reco options
    # if `raw` in flags but reco not we don't want to run options! Means input was a folder or
    # raw data file, but we didn't reconstruct it.
    ## XXX: check using tpa file kind
    if afReco in cfg.anaFlags or (afReco notin cfg.anaFlags and afRaw notin cfg.anaFlags):
      for opt in recoOptions:
        tc(reconstruction(input, opt, cfg)) # no output

proc main(
  input: string,
  outpath = "",
  years: seq[int] = @[],
  inputKind: InputKind = ikDataFolder,
  trackingLogs = "", # if CAST chain, path to the tracking log files
  back = false,
  calib = false,
  cdl = false,
  raw = false,
  reco = false,
  logL = false,
  tracking = false,
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
  overwrite = false,
  config = "" # XXX: not yet implemented!
     ) =
  var toContinue = true

  ## Parse full analysis chain related parameters
  let anaFlags = block:
    var flags: set[AnaFlags]
    if back:     flags.incl afBack
    if calib:    flags.incl afCalib
    if cdl:      flags.incl afCDL
    if raw:      flags.incl afRaw
    if reco:     flags.incl afReco
    if logL:     flags.incl afLogL
    if tracking: flags.incl afTracking
    flags

  ## Parse `reconstruction` related parameters
  let recoCfg = block:
    var
      recoFlags: set[RecoFlags]
      calibFactor = none(float)
      runNumberArg: Option[int]
    #if runNumber < 0:
    #  recoFlags.incl rfReadAllRuns
    if runNumber >= 0: # else:
      runNumberArg = some(runNumber)
    if classify(only_energy) != fcNaN:
      recoFlags.incl rfOnlyEnergy
      calibFactor = some(only_energy)
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
    initRecoConfig(recoFlags, runNumberArg, calibFactor, overwrite)

  let (is_rf, _, _, _) = isTosRunFolder(input)
  let inputKind = if is_rf: ikRunFolder
                  elif input.endsWith(".h5"): ikFileH5
                  else: inputKind
  var cfg = initConfig(input, outpath, inputKind, recoCfg, runType, anaFlags, trackingLogs)


  ## XXX: if the input is a H5 file open it and determine what type of
  ## file it is (raw, reco) and thus adjust what is the "input" in the
  ## following procedures. E.g. if it's `raw` then we need to use that
  ## as input to `reconstruction`.

  ## XXX: now differentiate between cases of:
  ## 1. full reconstruction of all runs from raw data
  ## 2. continuation of a "full run" setup based on a H5 file
  ## 3. reconstruction from a _single_ run folder
  ## 4. continuation of a single run file
  let t0 = epochTime()
  case cfg.inputKind
  of ikDataFolder:
    if afCDL in cfg.anaFlags: # perform CDL parsing / reconstruction
      let cdlFiles = fillCdlFiles(cfg)
      tc(cfg.runCdlChain(cdlFiles))

    for year in years:
      doAssert year.uint16 in [2014'u16, 2017, 2018], "Years supported are 2014, 2017, 2018."
      let dYear = DataYear(year)

      let dataFiles = fillDataFiles(cfg, dYear)
      if not existsFile(dataFiles.cdl):
        raise newException(ValueError, "The CDL file " & dataFiles.cdl & " does not seem to exist yet. Did you forget " &
          "the `--cdl` option?")
      tc(cfg.runCastChain(dataFiles))
  else:
    doAssert cfg.runType != rtNone, "For individual runs a `runType` is required."
    runChain(cfg)
  echo "The entire analysis chain took: ", (epochTime() - t0) / 60.0, " min"

when isMainModule:
  import cligen
  # we could do `multiDispatch`, but the majority of arguments apply to all
  dispatch(main, help = {
    "input"          : "The path containing the data to be processed.",
    "inputKind"      : """Type of input given. If not given assumes full CAST data, but overwritten
to ikFileH5 if explicit H5 file given. For a single run handing this is required.""",
    "years"          : "If input is full CAST data (ikDataFolder) the datasets to reconstruct.",
    "outpath"        : "The output path in which all files will be placed. If none given input path is used.",
    "trackingLogs"   : "(CAST) Can be used to overwrite default path for tracking log files.",
    "back"           : "(CAST) If set perform analysis of background data. Only relevant for full CAST data.",
    "calib"          : "(CAST) If set perform analysis of calibration data. Only relevant for full CAST data.",
    "cdl"            : "(CAST) If ste perform analysis of the CDL runs. Only relevant for full CAST data.",
    "raw"            : "If set performs parsing of data and storing in HDF5.",
    "reco"           : "If set performs geometric clustering and calculation of geometric properties.",
    "logL"           : "If set computes the `likelihood` dataset for the CAST data chain",
    "tracking"       : "If set adds the tracking information to the CAST data files",
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
    "overwrite"          : """(reco) If set will overwrite the results of the given calibration step, e.g.
to redo the energy or charge calibration."""

  })
