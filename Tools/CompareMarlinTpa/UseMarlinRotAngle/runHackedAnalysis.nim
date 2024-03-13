import shell, strutils, os

proc checkRes[T](res: T, msg = "Compilation failed") =
  when T is tuple:
    if res[1] != 0:
      quit(msg)
  else:
    if res != 0:
      quit(msg)

proc copyConfigs(recoPath: string): bool =
  if fileExists($recoPath / "config.nims"):
    let res = shellVerbose:
      cp ($recoPath)/config.nims ($recoPath)/config.nims_backup
    checkRes res, "backing up config.nims file failed!"
    result = true
  let res = shellVerbose:
    cp runHackedAnalysis.nims ($recoPath)/config.nims
  checkRes res, "copying confing.nims failed"

proc doUndoBackup(recoPath: string, undo: bool) =
  if undo:
    let res = shellVerbose:
      cp ($recoPath)/config.nims_backup ($recoPath)/config.nims
    checkRes res, "reverting config.nims from backup failed"
  else:
    let res = shellVerbose:
      rm ($recoPath)/config.nims
    checkRes res, "removing config.nims failed"

#const commonArgs = "-f -d:danger -d:release --threads:on -d:activateHijack"
const commonArgs = "-f -d:danger -d:release --threads:on -d:activateHijack -d:hijackRotationAngle"
proc compileBackground(recoPath: string): int =
  let res = shellVerbose:
    nim c ($commonArgs) "-d:hijackBackground" ($recoPath)/reconstruction.nim
  result = res[1]

proc compileCalibration(recoPath: string): int =
  let res = shellVerbose:
    nim c ($commonArgs) "-d:hijackCalibration" ($recoPath)/reconstruction.nim
  result = res[1]

proc getFile(prefix: string, raw = false): string =
  if raw:
    result = &"/mnt/1TB/CAST/{prefix}_2014_Raw.h5"
  else:
    result = &"/mnt/1TB/CAST/{prefix}_2014_MarlinHijackedCluster_RotAng.h5"

proc runRawToReco(prefix: string) =
  let infile = getFile(prefix, raw = true)
  let outfile = getFile(prefix, raw = false)
  let res = shellVerbose:
    "../../../Analysis/ingrid/reconstruction" ($infile) "--out" ($outfile)
  checkRes res

proc runRecoFlags() =
  const recoOptions = ["--only_charge", "--only_gas_gain",
                       # "--only_gain_fit",
                       "--only_energy_from_e"]
  const filePrefixes = ["DataRuns"] # ["CalibrationRuns", "DataRuns"]
  for f in filePrefixes:
    let infile = getFile(f, raw = false)
    for opt in recoOptions:
      let res = shellVerbose:
        "../../../Analysis/ingrid/reconstruction" ($infile) ($opt)
      checkRes res, "running reconstruction on " & $infile & " with opt " & $opt & " failed!"

proc runLogReader() =
  let file = "/mnt/1TB/CAST/DataRuns_2014_MarlinHijackedCluster_RotAng.h5"
  let logFolder = "../../../resources/LogFiles/tracking-logs"
  template run(args: string): untyped =
    block:
      let res = shellVerbose:
        "../../../LogReader/cast_log_reader" ($logFolder) "--h5out" ($file) ($args)
      checkRes res, "running log reader failed with args " & $args
  # first potentially delete log files
  run("--delete")
  run("")

proc runLikelihood() =
  let lPath = "../../../Analysis/ingrid/likelihood"
  let file = "/mnt/1TB/CAST/DataRuns_2014_MarlinHijackedCluster_RotAng.h5"
  let outname = "/mnt/1TB/CAST/LHood_2014_MarlinHijackedCluster_RotAng.h5"
  let res = shellVerbose:
    ($lPath) ($file) "--h5out" ($outname)

proc runFull() =
  runRawToReco("CalibrationRuns")
  runRawToReco("DataRuns")

when isMainModule:
  let isBackground = true
  let recoPath = "../../../Analysis/ingrid/"
  let toUndo = copyConfigs(recoPath)
  var res: int
  if isBackground:
    res = compileBackground(recoPath)
  else:
    res = compileCalibration(recoPath)
  checkRes res

  let prefix = block:
    if isBackground:
      "DataRuns"
    else:
      "CalibrationRuns"
  runRawToReco(prefix)
  runRecoFlags()
  #runLogReader()
  runLikelihood()

  var undoBackup = false
  if undoBackup:
    doUndoBackup(recoPath, toUndo)
