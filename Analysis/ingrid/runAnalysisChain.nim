import shell, docopt
import strutils, strformat, os
import helpers / utils
import logging

const docStr = """
Usage:
  runAnalysisChain <dataPath> (--2014 | --2017 | --2018) [options]

Options:
  --2014       Run 2014 analysis chain
  --2017       Run 2017 analysis chain
  --2018       Run 2018_2 analysis chain
  -h, --help   Show this help
  --version    Show the version number
"""
const doc = withDocopt(docStr)

const subPaths = ["2014_15", "2017", "2018_2"]
const recoOptions = ["", "--only_fadc", "--only_charge", "--only_gas_gain",
                     "--only_gain_fit", "--only_energy_from_e"]
const relIDPath = "../../InGridDatabase/src/resources"

type
  DataYear = enum
    dy2014 = "2014"
    dy2017 = "2017"
    dy2018 = "2018"

# set up the logger
var L = newConsoleLogger()
if not dirExists("logs"):
  createDir("logs")
var fL = newFileLogger("logs/runAnalysisChain.log", fmtStr = verboseFmtStr)
when isMainModule:
  addHandler(L)
  addHandler(fL)

proc rawDataManipulation(path, runType, outName: string): bool =
  let runTypeStr = &"--runType={runType}"
  let outStr = &"--out={outName}"
  info "Running raw_data_manipulation on " & $path & $runTypeStr & $outStr
  let res = shellVerbose:
    ./raw_data_manipulation `$path` `$runTypeStr` `$outStr`
  info "Last commands output: " & $res[0]
  info "Last commands exit code: " & $res[1]
  result = res[1] == 0

proc reconstruction(path: string, option = ""): bool =
  info "Running reconstruction on " & $path & $option
  let res = shellVerbose:
    ./reconstruction `$path` `$option`
  info "Last commands output: " & $res[0]
  info "Last commands exit code: " & $res[1]
  result = res[1] == 0

func likelihood() = discard

template tc(cmd: untyped): untyped {.dirty.} =
  ## convenience template to wrap a command in toContinue checks
  if toContinue:
    toContinue = cmd

proc runChain(path: string, dYear: DataYear): bool =
  ## performs the whole chain of the given dataset
  var toContinue = true
  case dYear
  of dy2017, dy2018:
    # copy the correct ingrid Database file
    # TODo: if we only copy the backup databases, we always overwrite potentially
    # new fit parameters that were added to the database during the reconstruction!
    # Need to copy back after running the reconstruction in some form.
    copyFile(relIDPath / &"ingridDatabase{$dYear}.h5", relIdPath / "ingridDatabase.h5")
  else: discard

  # raw data for calibration
  tc(rawDataManipulation(path / "CalibrationRuns", "calib", path / &"CalibrationRuns{$dYear}.h5"))
  # raw data for background
  #tc(rawDataManipulation(path / "DataRuns", "back", path / &"DataRuns{$dYear}.h5"))

  # reconstruction for both
  #for opt in recoOptions:
  #  tc(reconstruction(path / &"CalibrationRuns{$dYear}.h5", opt))
  #  if opt != "--only_gain_fit":
  #    tc(reconstruction(path / &"DataRuns{$dYear}.h5", opt))
  #result = toContinue

proc main =
  let args = docopt(doc)
  let dataPath = $args["<dataPath>"]
  var toContinue = true

  if args["--2014"].toBool:
    tc(runChain(dataPath / "2014_15", dy2014))
  if toContinue and args["--2017"].toBool:
    tc(runChain(dataPath / "2017", dy2017))
  if toContinue and args["--2018"].toBool:
    tc(runChain(dataPath / "2018_2", dy2018))

when isMainModule:
  main()
