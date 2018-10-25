import shell
import strutils, sequtils
from ingrid/utils/pure import parseRunType
from ingrid/ingrid_types import RunTypeKind
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

proc baf(runFolder: string, runType: RunTypeKind) =
  ## commands to run on the BAF to backup the run(s) created by
  ## `castPC` to copy them over to `tpc18` and `tpc00`
  shell:


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
    discard
  of hnTpc00:
    discard
  of hnVoid:
    discard
  of hnBaf:
    baf(runFolder, runType)

when isMainModule:
  main()
