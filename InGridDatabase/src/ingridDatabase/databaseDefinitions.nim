import strutils, ospaths, re, tables, macros

# helper proc to remove the ``src`` which is part of `nimble path`s output
# this is a bug, fix it.
proc removeSuffix(s: string, rm: string): string {.compileTime.} =
  result = s
  result.removeSuffix(rm)

# TODO: this is way too complicated. For some reason I didn't realize yesterday
# that `dirExists` DOES work at compile time :S Fix this

const dbDir = "resources"
const path1 = staticExec("nimble path ingridDatabase") / dbDir
const path2 = staticExec("nimble path ingridDatabase").removeSuffix("src") / dbDir
var tmpPath {.compileTime.} = ""

# check whether path exists to check whether we need the `src` or not
when defined(windows):
  const doesExist = "True"
  const doesNotExist = "False"
  # Need to be run in PowerShell as far as I'm aware
  const exists1 = gorge("Test-Path " & path1)
  const exists2 = gorge("Test-Path " & path2)
  static:
    when exists1 == doesExist:
      tmpPath = path1
    elif exists2 == doesExist:
      tmpPath = path2
    else:
      # else write a warning and put path to local folder
      {.fatal: "Could not find valid path to ingridDatabase.h5 file! Did you forget" &
        "to install the `ingridDatabase` nim module?".}
else:
  const doesExist = 0
  const doesNotExist = 1
  const existsTup1 = gorgeEx("stat " & path1)
  const existsTup2 = gorgeEx("stat " & path2)
  static:
    when existsTup1[1] == doesExist:
      tmpPath = path1
    elif existsTup2[1] == doesExist:
      # else check without src
      tmpPath = path2
    else:
      # else write a warning and put path to local folder
      {.fatal: "Could not find valid path to ingridDatabase.h5 file! Did you forget" &
        "to install the `ingridDatabase` nim module?".}

# if we haven't quit we found the path
const ingridPath* = tmpPath
static:
  hint("Found path " & ingridPath)
const dbPath* = joinPath(ingridPath, "ingridDatabase.h5")

const
  ChipInfo* = "chipInfo.txt"
  FsrPrefix* = "fsr"
  FsrPattern* = "fsr*.txt"
  TotPrefix* = "TOTCalib"
  TotPattern* = TotPrefix & r"*\.txt"
  ThresholdPrefix* = "threshold"
  ThresholdPattern* = r"threshold[0-9]\.txt"
  ThresholdMeansPrefix* = ThresholdPrefix & "Means"
  ThresholdMeansPattern* = ThresholdPrefix & r"Means*\.txt"
  SCurvePrefix* = "voltage_"
  SCurveRegPrefix* = r".*" & SCurvePrefix
  SCurvePattern* = r"voltage_*\.txt"
  SCurveFolder* = "SCurve/"
  StartTotRead* = 20.0
let
  ChipNameLineReg* = re(r"chipName:")
  ChipNameReg* = re(r".*([A-Z])\s*([0-9]+)\s*W\s*([0-9]{2}).*")
  FsrReg* = re(FsrPrefix & r"([0-9])\.txt")
  FsrContentReg* = re(r"(\w+)\s([0-9]+)")
  TotReg* = re(TotPrefix & r"([0-9])\.txt")
  SCurveReg* = re(SCurveRegPrefix & r"([0-9]+)\.txt")

type
  ChipName* = object
    col*: char
    row*: int
    wafer*: int

  Chip* = object
    name*: ChipName
    info*: Table[string, string]

  TotType* = object
    pulses*: int
    mean*: float
    std*: float

  TypeCalibrate* {.pure.} = enum
    TotCalibrate = "tot"
    SCurveCalibrate = "scurve"

  # type to store results of fitting with mpfit
  FitResult* = object
    x*: seq[float]
    y*: seq[float]
    pRes*: seq[float]
    pErr*: seq[float]
    redChiSq*: float

proc `$`*(chip: ChipName): string =
  result = $chip.col & $chip.row & " W" & $chip.wafer
