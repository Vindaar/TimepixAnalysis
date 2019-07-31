import shell, docopt
import strutils, strformat, os
import helpers / utils

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

proc rawDataManipulation(path, runType, outName: string) =
  let runTypeStr = &"--runType={runType}"
  let outStr = &"--out={outName}"
  shell:
    ./raw_data_manipulation `$path` `$runTypeStr` `$outStr`

proc reconstruction(path: string,
                    option = "") =
  shell:
    ./reconstruction `$path` `$option`

func likelihood() = discard

proc runChain(path: string, dYear: DataYear) =
  ## performs the whole chain of the given dataset
  case dYear
  of dy2017, dy2018:
    # copy the correct ingrid Database file
    copyFile(relIDPath / &"ingridDatabase{$dYear}.h5", relIdPath / "ingridDatabase.h5")
  else: discard

  # raw data for calibration
  rawDataManipulation(path / "CalibrationRuns", "calib", path / &"CalibrationRuns{$dYear}.h5")
  # raw data for background
  rawDataManipulation(path / "DataRuns", "back", path / &"DataRuns{$dYear}.h5")

  # reconstruction for both
  for opt in recoOptions:
    reconstruction(path / &"CalibrationRuns{$dYear}.h5", opt)
    if opt != "--only_gain_fit":
      reconstruction(path / &"DataRuns{$dYear}.h5", opt)

proc main =
  let args = docopt(doc)
  let dataPath = $args["<dataPath>"]

  if args["--2014"].toBool:
    runChain(dataPath / "2014_15", dy2014)
  if args["--2017"].toBool:
    runChain(dataPath / "2017", dy2017)
  if args["--2018"].toBool:
    runChain(dataPath / "2018_2", dy2018)

when isMainModule:
  main()
