import strutils, ospaths, re, tables, macros, os, times
import parsetoml

# helper proc to remove the ``src`` which is part of `nimble path`s output
# this is a bug, fix it.
proc removeSuffix(s: string, rm: string): string {.compileTime.} =
  result = s
  result.removeSuffix(rm)

# TODO: this is way too complicated. For some reason I didn't realize yesterday
# that `dirExists` DOES work at compile time :S Fix this

const dbDir = "resources"
const path1 = staticExec("nimble path ingridDatabase").strip / dbDir
const path2 = staticExec("nimble path ingridDatabase").strip.removeSuffix("src") / dbDir
var tmpPath {.compileTime.} = ""
# check whether path exists to check whether we need the `src` or not
static:
  # needs to be in static, otherwise tmpPath won't be set
  when dirExists(path1):
    tmpPath = path1
  elif dirExists(path2):
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
  ChargeCalibGasGain* = "chargeCalibGasGain"

  RunPeriodStart* = "start"
  RunPeriodStartTimestamp* = "startTimestamp"
  RunPeriodStop* = "stop"
  RunPeriodStopTimestamp* = "stopTimestamp"
  RunPeriodFirstRun* = "firstRun"
  RunPeriodLastRun* = "lastRun"
  RunPeriodRunsAvailable* = "validRuns"
  RunPeriodRunDset* = "runs"
  RunPeriodChipsAttr* = "chipsInRunPeriod"
  RunPeriodAttr* = "runPeriod" # run period of chip, just the parent group essentially

  # defines the "center" chips of different detectors, which are natively supported
  # by the ingrid database. However, this is only for reference.
  centerChip2014* = "D03W63"
  centerChip2017* = "H10W69"

let
  ChipNameLineReg* = "chipName:"
  RunPeriodLineReg* = "runPeriod:"
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
    run*: string # name of run period (has to exist in file!)
    info*: Table[string, string]

  TotType* = object
    pulses*: int
    mean*: float
    std*: float

  TypeCalibrate* {.pure.} = enum
    TotCalibrate = "tot"
    SCurveCalibrate = "scurve"

  RunPeriod* = object
    name*: string
    start*: DateTime
    stop*: DateTime
    validRuns*: seq[int]
    firstRun*: int
    lastRun*: int
    additionalInfo*: Table[string, TomlValueRef]

proc `$`*(chip: ChipName): string =
  result = $chip.col & $chip.row & " W" & $chip.wafer

#####################################################
#### Procs specifically realted to hardware #########
#####################################################

proc getSeptemHChip*(chipNumber: int): string =
  ## returns the name of a given SeptemH chip
  const names = ["E6 W69",
                 "K6 W69",
                 "H9 W69",
                 "H10 W69",
                 "G10 W69",
                 "D9 W69",
                 "L8 W69"]
  result = names[chipNumber]
