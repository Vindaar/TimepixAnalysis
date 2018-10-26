import shell
import strutils, sequtils
from ingrid/tos_helpers import isTosRunFolder, getRunTimeInfo, formatAsOrgDate, parseRunType
from helpers/utils import getListOfFiles, getDaysHoursMinutes
from ingrid/ingrid_types import RunTypeKind
import ospaths
import docopt

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
  dataCAST = "/home/ingrid/Data"
  krb5 = "kinit basti90"
  baf = "basti90@baf.physik.uni-bonn.de:~/"

type
  hostNameKind = enum
    hnCast = "InGrid-DAQ"
    hnTpc18 = "tpc18"
    hnVoid = "void"
    hnTpc00 = "tpc00"
    hnBaf = "baf03.physik.uni-bonn.de" # we won't always be on baf03!

proc recordRun(runFolder: string, runType: RunTypeKind) =
  ## records the given run folder in the run_list.org file.
  ## - reads the start and end time of the run
  ## - reads the duration of the run
  ## - marks the run as backed up
  ## - outputs the table row for the Org table and writes them to the file

  const runFile = expandTilde("~/org/auto_run_list.org")

  let regex = r"^/([\w-_]+/)*data\d{6}\.txt$"
  let (is_run_folder, runNumber, rfKind, contains_run_folder) = isTosRunFolder(run_folder)
  let files = getListOfFiles(run_folder, regex)
  let rt_info = getRunTimeInfo(files)

  let
    parsed_first = formatAsOrgDate(rt_info.t_start)
    parsed_last  = formatAsOrgDate(rt_info.t_end)
    lenStr = getDaysHoursMinutes(rt_info.t_length)

  var f = open(runFile, fmAppend)
  f.write(&"| {runNumber} | {runType} | <{parsed_first}> | <{parsed_last}> | {lenStr} | y |     |\n")
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
      tar -czf `$runFolder`.tar.gz `$runFolder`
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
  case hnKind
  of hnCast:
    castPC(runFolder, runType)
  of hnTpc18:
    recordRun(runFolder, runType)
  of hnTpc00:
    discard
  of hnVoid:
    discard
  of hnBaf:
    bafNode(runFolder, runType)

when isMainModule:
  main()
