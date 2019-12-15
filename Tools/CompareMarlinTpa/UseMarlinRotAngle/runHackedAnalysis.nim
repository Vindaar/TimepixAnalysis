import shell, strutils, os

const isBackground = false
when isBackground:
  let prefix = "DataRuns"
else:
  let prefix = "CalibrationRuns"


let recoPath = "../../../Analysis/ingrid/"
var undoBackup = false
if fileExists($recoPath / "config.nims"):
  shell:
    cp ($recoPath)/config.nims ($recoPath)/config.nims_backup
  undoBackup = true
shell:
  cp runHackedAnalysis.nims ($recoPath)/config.nims
when isBackground:
  let res = shellVerbose:
    nim c "-f -d:danger -d:release --threads:on -d:activateHijack -d:hijackBackground" ($recoPath)/reconstruction.nim
else:
  let res = shellVerbose:
    nim c "-f -d:danger -d:release --threads:on -d:activateHijack -d:hijackCalibration" ($recoPath)/reconstruction.nim
if res[1] != 0:
  quit("Compilation failed")
if undoBackup:
  shell:
    cp ($recoPath)/config.nims_backup ($recoPath)/config.nims
else:
  shell:
    rm ($recoPath)/config.nims
let outfile = &"/mnt/1TB/CAST/{prefix}_2014_MarlinHijackedCluster.h5"
let infile = &"/mnt/1TB/CAST/{prefix}_2014_Raw.h5"
echo "INFILE ", infile
echo "OUTFILE ", outfile
shell:
  #"../../../Analysis/ingrid/reconstruction ../../../Tests/run_245_2014.h5" "--out" testfile.h5
  "../../../Analysis/ingrid/reconstruction" ($infile) "--out" ($outfile)
  #"../../../Analysis/ingrid/reconstruction /mnt/1TB/CAST/CalibrationRuns_2014_Raw.h5" "--out" ($outfile) "--runNumber 300"
  #"../../../Analysis/ingrid/reconstruction" ($outfile) "--runNumber 300 --only_fe_spec"
