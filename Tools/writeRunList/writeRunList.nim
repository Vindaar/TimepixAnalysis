import times, sequtils, algorithm, strutils, os, strformat
import ingrid / [tos_helpers, ingrid_types]
from helpers/utils import getDaysHoursMinutes
import docopt
import nimhdf5

when defined(linux):
  const commitHash = staticExec("git rev-parse --short HEAD")
else:
  const commitHash = ""

# get date using `CompileDate` magic
const currentDate = CompileDate & " at " & CompileTime

const docTmpl = """
Version: $# built on: $#
A tool to write the run list including # trackings

Usage:
  writeRunList <H5file> <H5file2> [options]

Options:
  -h, --help             Show this help
  --version              Show the version number

The first file should be the `DataRuns.h5` and the second the `CalibrationRuns.h5`.
"""
const doc = docTmpl % [commitHash, currentDate]

const HeaderElements = ["Run #",
                        "Type",
                        "DataType",
                        "Start",
                        "End",
                        "Length",
                        "# trackings",
                        "# frames",
                        "# FADC",
                        "Backup?",
                        "Notes"]
const HeaderLine = HeaderElements.foldl(a & " | " & b, "") & " |\n"

type
  RunList = seq[ExtendedRunInfo]

func determineType(runInfo: RunInfo): RunTypeKind =
  ## depending on run length and ratio of FADC to normal InGrid events,
  ## will deterime run type
  const
    evRatio = 0.5
    calibShort = 1800 # min length of calibration run 1/2 h
    calibLong = 24 * 3600 # max length of calibration run 12 h
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

proc createRunList(fname: string, runType: RunTypeKind): RunList =
  var h5f = H5open(fname, "r")
  let fInfo = getFileInfo(h5f)
  for r in fInfo.runs:
    result.add getExtendedRunInfo(h5f, r, runType)

proc `<`(a, b: ExtendedRunInfo): bool =
  result = a.timeInfo.t_start < b.timeInfo.t_start

proc sortRuns(runs: varargs[RunList]): RunList =
  var allRuns: RunList
  for r in runs:
    allRuns.add r
  result = allRuns.sorted

proc write(fname: string, runs: RunList) =
  template createOrAppend(f, runFile: untyped): untyped =
    if existsFile(runFile):
      f = open(runFile, fmAppend)
    else:
      f = open(runFile, fmWrite)
      f.write(HeaderLine)
  var f: File
  createOrAppend(f, fname)

  for r in runs:
    let
      start = r.timeInfo.t_start.formatAsOrgDate
      stop = r.timeInfo.t_end.formatAsOrgDate
      lenStr = getDaysHoursMinutes(r.timeInfo.t_length)
    let entry = &"| {r.runNumber} | {r.runType} | {r.rfKind} " &
        &"| <{start}> | <{stop}> | {lenStr} | {r.trackings.len} " &
        &"| {r.nEvents} | {r.nFadcEvents} | y |\n"
    f.write(entry)
  f.close()

proc writeSummary(runs: RunList, runType: RunTypeKind) =
  ## writes a summary of the total time stored in the current type
  echo "Type: ", runType
  var trackingDuration = initDuration()
  var nonTrackingDuration = initDuration()
  for r in runs:
    trackingDuration = trackingDuration + r.trackingDuration
    nonTrackingDuration = nonTrackingDuration + r.nonTrackingDuration
  echo "\t trackingDuration: ", trackingDuration
  echo "\t nonTrackingDuration: ", nonTrackingDuration

proc main =
  let args = docopt(doc)
  let backFname = $args["<H5file>"]
  let calibFname = $args["<H5file2>"]

  let pltBack = createRunList(backFname, rtBackground)
  let pltCalib = createRunList(calibFname, rtCalibration)

  let sortedRuns = sortRuns(pltBack, pltCalib)
  echo sortedRuns

  writeSummary(pltBack, rtBackground)
  writeSummary(pltCalib, rtCalibration)

  write("runList.org", sortedRuns)

when isMainModule:
  main()
