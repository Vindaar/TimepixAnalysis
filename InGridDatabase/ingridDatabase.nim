import strutils, strformat, sequtils, re, os, times, sugar
import ingrid/tos_helpers
import helpers/utils
import arraymancer, parsetoml, seqmath, nimhdf5, zero_functional

import ingridDatabase/databaseDefinitions
import ingridDatabase/databaseRead
import ingridDatabase/databaseWrite
import ingridDatabase/databaseUtils

# cannot import `Chip` due to clash with this modules `Chip`
import ingrid/ingrid_types except Chip
import ingrid / calibration / calib_fitting

when defined(linux):
  const commitHash = staticExec("git rev-parse --short HEAD")
else:
  const commitHash = ""

# get date using `CompileDate` magic
const currentDate = CompileDate & " at " & CompileTime

const vsn = "$# built on: $#" % [commitHash, currentDate]
const docTmpl = """
Version: $#

The InGrid database management tool.
  This tool and library aims to provide two things. Running it as a main
  module, it is a tool to add / manage the InGrid database HDF5 file, which
  stores information like FSR, thresholds, TOT calibration etc. for each chip
  with the corresponding name and location (e.g. Septem H).
  As a library it provides convenience functions to access the database to
  read and write to it, e.g.
  getTot(chipName)
  to get the ToT values of that chip.
"""
const doc = docTmpl % [vsn]
import macros
macro insertDocComment(s: static string): untyped =
  result = newNimNode(nnkCommentStmt)
  result.strVal = s

proc calibrateType(tcKind: TypeCalibrate, runPeriod: string) =
  ## calibrates the given type and stores the fit results as attributes
  ## of each chip with valid data for calibration
  # open H5 file with write access
  var h5f = H5open(dbPath, "rw")
  case tcKind
  of TotCalibrate:
    # perform TOT calibration for all chips, which contain a TOT
    # calibration
    # - get groups in root of file == chips in file
    # - get TOT object of each chip
    # - perform Fit as in plotCalibration
    # - take `mpfit` result and save fit parameters as attrs
    for runPeriodGrp in h5f.items(start_path = "/", depth = 1):
      if runPeriod.len > 0 and not runPeriodGrp.name.startsWith(runPeriod.formatName()): continue
      for chipGrp in h5f.items(start_path = runPeriodGrp.name, depth = 1):
        # group name is the chip name
        let chip = chipGrp.name
        let runPeriod = runPeriodGrp.name
        let tot = h5f.getTotCalib(chip, runPeriod)
        let totCalib = fitToTCalib(tot, 0.0)
        # given fit result write attributes:
        h5f.writeTotCalibAttrs(chipGrp, totCalib)

  of SCurveCalibrate:
    discard

  let err = h5f.close()
  if err < 0:
    echo "Could not properly close database! err = ", err


proc parseCalibrateType(typeStr: string): TypeCalibrate =
  ## checks the given `typeStr` whether it's a valid type to
  ## calibrate for. If so returns that type, else raises an
  ## exception
  case typeStr.toLower
  of $TotCalibrate:
    result = TotCalibrate
  of $SCurveCalibrate:
    result = ScurveCalibrate
  else:
    raise newException(KeyError, "Given type is not a valid type to calibrate " &
      """for. Valid types are ["TOT", "SCurve"] (case insensitive)""")

proc removePrefix(s, pr: string): string =
  result = s
  result.removePrefix(pr)

proc parseChipInfo(filename: string): Chip =
  ## parses the `chipInfo` file and returns a chip object from that
  let info = readFile(filename).splitLines
  doAssert info[0].startsWith(ChipNameLineReg), "The chipInfo file needs to contain the " &
    "`chipName: ` as the first line"
  doAssert info[1].startsWith(RunPeriodLineReg), "The chipInfo file needs to contain the " &
    "`runPeriod: ` as the second line"
  let chipStr = info[0].removePrefix(ChipNameLineReg).strip
  # fill the `name` and `run` field
  result.name = parseChipName(chipStr)
  result.run = info[1].removePrefix(RunPeriodLineReg).strip
  ## TODO: possibly use the `version` field instead of using the `info` table
  ## as we currently do for things like TimepixVersion
  #result.version = parseEnum[TimepixVersion](info[2].removePrefix(TimepixVersionReg).strip)
  # parse the optional content of further notes on the chip
  # e.g. board, chip number on that board etc.
  result.info = initTable[string, string]()
  for line in info[2 .. info.high]:
    if line.len > 0:
      if ':' in line:
        let
          lineSplit = line.split(':')
          n = lineSplit[0].strip
          c = lineSplit[1].strip
        result.info[n] = c
      else:
        # support for lines without a key
        result.info[line.strip] = ""

proc parseFsr(filename: string): FSR =
  ## parses the contents of a (potential) given FSR file
  ## uses regex `FsrContentReg` to parse the file
  echo "Filename is ", filename
  if existsFile(filename) == true:
    let fLines = readFile(filename).splitLines
    var dacVal: array[2, string]
    result = initTable[string, int]()
    for line in fLines:
      if match(line.strip, FsrContentReg, dacVal):
        result[dacVal[0]] = dacVal[1].parseInt

proc parseScurveFolder(folder: string): SCurveSeq =
  ## parses all SCurve files (matching `SCurveReg`) in the folder
  ## and returns an `SCurveSeq` from it
  var scurveMatch: array[1, string]
  echo folder
  for f in walkFiles(joinPath(folder, SCurvePattern)):
    if match(f, SCurveReg, scurveMatch):
      result.files.add f
      let voltage = scurveMatch[0].parseInt
      let scurve = readScurveVoltageFile(f)
      result.curves.add scurve

proc parseTotFile(filename: string): Tot =
  ## checks for the existence of a TOT calibration file, reads it
  ## and returns a `Tot` object
  if existsFile(filename) == true:
    let (chip, tot) = readToTFile(filename,
                                  startRead = StartTotRead,
                                  totPrefix = TotPrefix)
    result.pulses = tot.pulses
    result.mean = tot.mean
    result.std = tot.std

proc parseThresholdFile(filename: string): Threshold =
  ## parses a Threshold(Means) file and returns it
  echo filename
  let flat = readFile(filename).splitLines -->
    map(it.splitWhitespace) -->
    flatten() -->
    map(it.parseInt)
  result = flat.toTensor().reshape([256, 256])

proc addChip(folder: string) =
  ## handles adding the data for a chip from a folder
  ## it parses all available data in the folder and adds it to the correct group
  ## in the database file
  ## the given folder needs to conform to a few things to be parsed properly:
  ## Folder structure:
  ## - chipInfo.txt
  ##   Contains:
  ##     - chipName: <name as given in TOS>
  ##     - arbitrary lines, which are individually added as string attributes for
  ##       this chips group. If lines have a colon, the string before that is the
  ##       name of the attribute
  ## Optional: if the following exists, it will be added
  ## - fsr?.txt
  ##   The FSR file of that chip
  ## - threshold?.txt
  ##   The threshold equalization of that chip
  ## - thresholdMeans?.txt
  ##   The mean values of threshold equalization of that chip
  ## - SCurve/
  ##   A folder containing the SCurves for that chip
  ## - TOTCalib?.txt
  ##   the TOT calibrration of that chip
  doAssert dirExists(folder), "Please hand an existing folder!"
  doAssert existsFile joinPath(folder, ChipInfo), "The given folder needs to contain a `chipInfo.txt`!"

  # read chip info file
  let chip = parseChipInfo joinPath(folder, ChipInfo)
  echo "Chip is ", chip

  # after chip check for fsr
  var fsrFileName = ""
  for f in walkFiles(joinPath(folder, FsrPattern)):
    echo f
    fsrFileName = f
  let fsr = parseFsr fsrFileName
  if fsr.len > 0:
    echo fsr

  # check for SCurves
  var scurves: SCurveSeq
  if dirExists(joinPath(folder, SCurveFolder)) == true:
    scurves = parseScurveFolder(joinPath(folder, SCurveFolder))

  # check for TOT
  var tot: Tot
  for f in walkFiles(joinPath(folder, TotPattern)):
    tot = parseTotFile(f)
  for f in walkFiles(joinPath(folder, TotPatternTpx3)):
    tot = parseTotFile(f)


  # check for threshold / threshold means
  var
    threshold: Threshold
    thresholdMeans: ThresholdMeans
  for f in walkFiles(joinPath(folder, ThresholdPattern)):
    threshold = parseThresholdFile(f)
  #for f in walkFiles(joinPath(folder, ThresholdMeansPattern)):
    #thresholdMeans = parseThresholdFile(f)
  #  assert false, "Threshold means writing not yet implemented!"

  # given all data, add it to the H5 file
  addChipToH5(chip, fsr, scurves, tot, threshold, thresholdMeans)

proc addRunPeriod(name: string) =
  ## Adds a new (or updates) a run period to the Ingrid Database.
  ## A run period is a set of runs with timestamps of the first and last
  ## (TODO: and if desired for all runs?) which are used to decide which
  ## calibration data is read/written for a given chip.
  ## The input file is simply an ASCII file using TOML syntax (i.e. a TOML file).
  ##
  ## The run period file needs to have the following structure:
  ## .. code-block::
  ##   # optional title
  ##   title = "Run period 2 of CAST, 2017/18"
  ##   # list of the run periods defined in the file
  ##   runPeriods = ["Run2"]
  ##
  ##   [Run2] # name must match the periods listed above in `runPeriods`
  ##   start = 2017-10-30 # from when it will be valid, dates as full days will be fully inclusive
  ##   stop = 2018-04-11 # up to when it will be valid
  ##   # run numbers (used to look up which run period a run belongs to)
  ##   # either as a sequence of run numbers
  ##   validRuns = [ 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89
  ##     90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,101,102,103
  ##     104,105,106,107,108,109,110,111,112,113,114,115,116,117
  ##     118,119,120,121,122,123,124,125,126,127,128,145,146,147
  ##     148,149,150,151,152,153,154,155,156,157,158,159,160,161
  ##     162,163,164,165,166,167,168,169,170,171,172,173,174,175
  ##     176,177,178,179,180,181,182,183,184,185,186,187,188
  ##   ]
  ##   # or as simply a range given as start and stop values
  ##   firstRun = 76
  ##   lastRun = 188
  ##   additionalFields = "will be added to the database as written here"
  ##
  ## The name of the file is unimportant, but it is recommended to
  ## call it runPeriod.toml
  ## As mentioned in the comment on the example, the file requires a set of run
  ## numbers to be handed. These will be used to perform a lookup on the correct run
  ## period in combination with the chip name. Alternatively a lookup can be done using
  ## chip name + a timestamp.
  let run = parseToml.parseFile(name)
  # check if the data is valid
  if "runPeriods" notin run:
    raise newException(ValueError, "Given TOML file does not contain a run period " &
      "variable about, which run periods are described in the file. Add a `runPeriods` " &
      "variable at top level, which has an array of strings as values. The strings " &
      "must correspond to tables defined in the TOML file, e.g.\n" &
      """
runPeriods = ["Run123"]
[Run123]
...
""")

  var runPeriods = newSeq[RunPeriod]()
  for k in run["runPeriods"].getElems:
    let tab = run[k.getStr]
    let startInFile = "start" in tab
    let stopInFile = "stop" in tab
    let runsAvailable = ("validRuns" in tab and tab["validRuns"].getElems.len > 0) or
      ("firstRun" in tab and "lastRun" in tab)
    if not (startInFile and stopInFile and runsAvailable):
      raise newException(ValueError, "Given TOML file does not contain `start`, " &
        "`stop`, and either `validRuns` or `firstRun` and `lastRun` in run period " &
        k.getStr)
    var runP = RunPeriod(name: k.getStr)
    for key, v in pairs(tab.tableVal):
      case key
      of RunPeriodStart:
        runP.start = v.parseTomlTime()
      of RunPeriodStop:
        runP.stop = v.parseTomlTime()
      of RunPeriodRunsAvailable:
        runP.validRuns = v.mapIt(it.checkAndGetInt)
      of RunPeriodFirstRun:
        runP.firstRun = v.checkAndGetInt
      of RunPeriodLastRun:
        runP.lastRun = v.checkAndGetInt
      else:
        runP.additionalInfo[key] = v
    if runP.validRuns.len == 0:
      runP.validRuns = toSeq(runP.firstRun .. runP.lastRun)
    elif runP.firstRun == runP.lastRun:
      runP.firstRun = min(runP.validRuns)
      runP.lastRun = max(runP.validRuns)
    else:
      doAssert(runP.firstRun == min(runP.validRuns) and
               runP.lastRun == max(runP.validRuns),
               "firstRun/lastRun and availableRuns do not agree on the covered " &
               "runs!\n" & &"firstRun = {runP.firstRun}, lastRun = {runP.lastRun}, " &
               &"availableRuns = {min(runP.validRuns)} .. {max(runP.validRuns)}")
    echo runP
    runPeriods.add runP
  # all checks of the file passed, add to database
  addRunPeriod(runPeriods)

proc main(addPeriod = "", # add run periods from given `runPeriod.toml` file
          add = "", # add the chip described by data in given directory
          calibrate = "", # calibrate all chips (unless restricted by runPeriod)
          runPeriod = "", # only calibrate chips in this period
          modify = false,
          delete = "", # delete chip with given name from DB
         ) =
  insertDocComment(doc)
  if modify:
    doAssert false, "Modify support is not implmented yet."
  elif delete.len > 0:
    doAssert false, "Delete support is not implmented yet."
  elif calibrate.len > 0:
    calibrateType calibrate.parseCalibrateType, runPeriod
  elif add.len > 0:
    # call proc to add data from folder to database
    addChip(add)
  elif addPeriod.len > 0:
    addRunPeriod(addPeriod)
  else:
    echo "Noting to do."

when isMainModule:
  import cligen
  clCfg.version = vsn
  dispatch(main, doc = doc, short = {"version" : 'V'}, help = {
    "addPeriod" : """Adds the run periods stored in the given `runPeriod.toml`
file to the database to which chips can be assigned""",
    "add" : "Add the chip contained in the given folder to the database",
    "calibrate" : """Perform given calibration {TOT, SCurve} and save
fit parameters as attributes. Calibrates all chips with suitable data.""",
    "runPeriod" : """Restrict calibration to all chips of the given run period.
Must match a run period in the database.""",
     "modify" : """Start a CLI interface to allow viewing the files
content and modify attributes / delete elements.""",
    "delete" : "Delete the chip with the given name from the database.",
    "version" : "Show the version number."
    })
