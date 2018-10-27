import shell
import strutils, sequtils, strscans, times
from ingrid/tos_helpers import isTosRunFolder, getRunInfo, formatAsOrgDate, parseRunType
from helpers/utils import getListOfFiles, getDaysHoursMinutes
from ingrid/ingrid_types import RunTypeKind, RunInfo
import ospaths, os
import docopt
import zero_functional

when defined(linux):
  const commitHash = staticExec("git rev-parse --short HEAD")
else:
  const commitHash = ""
# get date using `CompileDate` magic
const currentDate = CompileDate & " at " & CompileTime

const docTmpl = """
Version: $# built on: $#
InGrid run backup tool

Usage:
  backupRun <folder> --runType=<type> [options]

Options:
  <folder>            Name of the run folder we should backup
  --runType=<type>    Select run type (Calib | Back | Xray)
                      The following are parsed case insensetive:
                      Calib = {"calib", "calibration", "c"}
                      Back = {"back", "background", "b"}
                      Xray = {"xray", "xrayfinger", "x"}
  -h --help           Show this help
  --version           Show version.

"""
const doc = docTmpl % [commitHash, currentDate]

const
  # define constants needed for backup
  tosCAST = "/home/ingrid/TOS/data/runs"
  dataCAST = "/home/ingrid/Data/AutoBackup/"
  krb5 = "kinit basti90"
  baf = "basti90@baf.physik.uni-bonn.de:~/AutoBackup"
  
  runFile = expandTilde("~/org/auto_run_list.org")
  # a set that stores all backed up run numbers
  runSetFile = expandTilde("~/org/auto_run_set.org")

type
  hostNameKind = enum
    hnCast = "InGrid-DAQ"
    hnTpc18 = "tpc18"
    hnVoid = "void"
    hnTpc00 = "tpc00"
    hnBaf = "baf03.physik.uni-bonn.de" # we won't always be on baf03!
    hnHome = "basti-MS-7885"

func findLastRun(runFolder: string, runSet: set[uint16]): set[uint16] =
  ## takes a look at the directory in which the run folders should be stores
  ## `tosCAST` and returns the run with the highest run number
  var
    dummy: string
    runNumber: string
    runs: set[uint16]
  
  for pcKind, path in walkDir(runFolder):
    case pcKind
    of pcDir:
      if scanf(path, "$*Run_$*_", dummy, runNumber):
        runs.incl runNumber.parseInt.uint16
      debugecho "cc ", path
    else:
      discard

  # now cross check with already backed up runs, if there is
  # any new run
  var
    runMax: uint16
    backedUpMax: uint16
  if runs.card > 0:
    runMax = max(toSeq(items(runs)))
  if runSet.card > 0:
    backedUpMax = max(toSeq(items(runSet)))
  if backedUpMax < runMax:
    # need to extract ``all`` run numbers bigger than runMax
    let runsTodo = (toSeq(items(runs))) --> filter(it > backedUpMax.uint16) --> to(seq[uint16])
    for r in runs:
      if r > backedUpMax:
        result.incl r

func getRunPath(runNumber: int, runFolder: string): string =
  ## returns the path to the run with `runNumber` in `runFolder`
  var
    dummy: string
    readNum: string
  for pcKind, path in walkDir(runFolder):
    case pcKind
    of pcDir:
      if scanf(path, "$*Run_$*_", dummy, readNum):
        if readNum.parseInt == runNumber:
          result = path
    else: discard

func determineType(runInfo: RunInfo): RunTypeKind =
  ## depending on run length and ratio of FADC to normal InGrid events,
  ## will deterime run type
  const
    evRatio = 0.5
    calibShort = 1800 # min length of calibration run 1/2 h
    calibLong = 12 * 3600 # max length of calibration run 12 h
    backShort = 7200 # min length of background run 2 h

  let ratio = runInfo.nFadcEvents.float / runInfo.nEvents.float
  let runLen = runInfo.timeInfo.t_length
  if runLen >= initDuration(seconds = calibShort) and
     runLen <= initDuration(seconds = calibLong) and
     ratio > evRatio:
    result = rtCalibration
  elif runLen > initDuration(seconds = backShort) and
    ratio < evRatio:
    # even in super noisy runs a ratio of 0.5 (half as many FADC events)
    # will not occur
    result = rtBackground
  else:
    result = rtNone

proc genRunSet(): set[uint16] =
  ## generates the backed up run set from the run list file
  for line in lines(runSetFile):
    if line.len > 0:
      result.incl line.splitWhitespace()[0].parseInt.uint16
  echo result

proc recordRun(runInfo: RunInfo) =
  ## records the given run folder in the run_list.org file.
  ## - reads the start and end time of the run
  ## - reads the duration of the run
  ## - marks the run as backed up
  ## - outputs the table row for the Org table and writes them to the file
  let
    parsed_first = formatAsOrgDate(runInfo.timeInfo.t_start)
    parsed_last  = formatAsOrgDate(runInfo.timeInfo.t_end)
    lenStr = getDaysHoursMinutes(runInfo.timeInfo.t_length)
    ratio = runInfo.nFadcEvents.float / runInfo.nEvents.float

  var f = open(runFile, fmAppend)
  let entry = &"| {runInfo.runNumber} | {runInfo.runType} | {runInfo.rfKind} " &
      &"| <{parsed_first}> | <{parsed_last}> | {lenStr} " &
      &"| {runInfo.nEvents} | {runInfo.nFadcEvents} | y |\n"
  echo entry
  f.write(entry)
  f.close()

  # now append to run file
  f = open(runSetFile, fmAppend)
  f.write(&"{runInfo.runNumber}\t{runInfo.runType}\n")
  f.close()

proc castPC(runFolder: string, runType: RunTypeKind) =
  ## commands to run on the InGrid CAST PC
  shell:
    # first get Kerberos ticket
    # `$krb5`
    # start bash first...
    # bash
    one:
      cd `$tosCAST`
      tar -czf `$runFolder`_`$runType`.tar.gz `$runFolder`
    mv `$tosCAST`/`$runFolder`.tar.gz `$dataCAST`
    scp `$dataCAST`/`$runFolder`.tar.gz `$baf`

  # then depending on whether this is a calibration or background run
  # put it into the appropriate folder
  var dir = ""
  case runType
  of rtBackground:
    dir = "DataRuns"
  of rtCalibration:
    dir = "CalibrationRuns"
  else:
    echo "Unsupported run type: ", runType
    return
  shell:
    one:
      cd `$dataCAST`
      mv `$runFolder` "2018_2"/`$dir`

proc bafNode(runFolder: string, runType: RunTypeKind) =
  ## commands to run on the BAF to backup the run(s) created by
  ## `castPC` to copy them over to `tpc18` and `tpc00`
  # on the baf we first need to 
  #shell:
  discard

proc main =
  let
    args = docopt(doc)
    runFolder = $args["<folder>"]
    runTypeStr = $args["--runType"]

  var runType: RunTypeKind
  if runTypeStr != "nil":
    runType = parseRunType(runTypeStr)

  var hostname = ""
  shellAssign:
    hostname = hostname

  let hnKind = parseEnum[hostNameKind](hostname, hnBaf)

  let runSet = genRunSet()
  let runs = findLastRun(runFolder, runSet)
  for run in runs:
    var runInfo = getRunInfo(getRunPath(run.int, runFolder))
    runInfo.runType = runInfo.determineType
    
    case hnKind
    of hnCast:
      castPC(runFolder, runType)
    of hnTpc18:
      recordRun(runInfo)
    of hnTpc00:
      discard
    of hnVoid:
      discard
    of hnBaf:
      bafNode(runFolder, runType)
    of hnHome:
      recordRun(runInfo)

when isMainModule:
  main()
