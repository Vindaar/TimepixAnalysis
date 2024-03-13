import times, sequtils, algorithm, strutils, os, strformat
import ingrid / [tos_helpers, ingrid_types]
from helpers/utils import getDaysHoursMinutes
import nimhdf5
import cligen / macUt # for `docCommentAdd`

when defined(linux):
  const commitHash = staticExec("git rev-parse --short HEAD")
else:
  const commitHash = ""
# get date using `CompileDate` magic
const compileDate = CompileDate & " at " & CompileTime
const versionStr = "Version: $# built on: $#" % [commitHash, compileDate]

const docTmpl = """
A tool to write the run list including # trackings
"""
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
  withH5(fname, "r"):
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
  var activeTrackingTime = initDuration()
  var activeNonTrackingTime = initDuration()
  var totalTime = initDuration()
  var activeTime = initDuration()
  var activeFromIndices = initDuration()
  for r in runs:
    trackingDuration += r.trackingDuration
    nonTrackingDuration += r.nonTrackingDuration
    activeTrackingTime += r.activeTrackingTime
    activeNonTrackingTime += r.activeNonTrackingTime
    echo "Run: ", r.runNumber, " activeTrackingTime = ", r.activeTrackingTime.inSeconds(), " compare tracking: ", r.trackingDuration.inSeconds()
    for i, tr in r.trackings:
      echo "Tracking: ", i, " duration from eventDuration = ", tr.durations.sum
      activeFromIndices += initDuration(microseconds = (tr.durations.sum * 1e6).round.int)
    totalTime += r.totalTime
    activeTime += r.activeTime

  proc asHours(d: Duration): float =
    ## Returns the time as exact float hours
    result = d.inMicroSeconds().float / 3600e6

  echo "\t total duration: ", totalTime
  echo "\t   In hours: ", totalTime.asHours()
  echo "\t   active duration: ", activeTime.asHours()
  echo "\t trackingDuration: ", trackingDuration
  echo "\t   In hours: ", trackingDuration.asHours()
  echo "\t   active tracking duration: ", activeTrackingTime.asHours()
  echo "\t   active tracking duration from event durations: ", activeFromIndices.asHours()
  echo "\t nonTrackingDuration: ", nonTrackingDuration
  echo "\t   In hours: ", nonTrackingDuration.asHours()
  echo "\t   active background duration: ", activeNonTrackingTime.asHours()
  let df = toDf({ "Solar tracking [h]" : trackingDuration.asHours(),
                  "Background [h]" : nonTrackingDuration.asHours(),
                  "Active tracking [h]" : activeTrackingTime.asHours(),
                  "Active tracking (eventDuration) [h]" : activeFromIndices.asHours(),
                  "Active background [h]" : activeNonTrackingTime.asHours(),
                  "Total time [h]" : totalTime.asHours(),
                  "Active time [h]" : activeTime.asHours() })
  echo df.toOrgTable(precision = 6)

proc main(back, calib: string,
          runList = "runList.org") =
  ## Given a path to a `DataRuns.h5` (`back`) and `CalibrationRuns.h5` file (`calib`) outputs a
  ## run list and the total tracking and non tracking durations.
  docCommentAdd(versionStr)

  let pltBack = createRunList(back, rtBackground)
  let pltCalib = createRunList(calib, rtCalibration)

  let sortedRuns = sortRuns(pltBack, pltCalib)
  #echo sortedRuns

  writeSummary(pltBack, rtBackground)
  writeSummary(pltCalib, rtCalibration)

  write(runList, sortedRuns)

when isMainModule:
  import cligen
  dispatch main
