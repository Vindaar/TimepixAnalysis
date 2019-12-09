import shell, strutils, os

let recoPath = "../../../Analysis/ingrid/"
var undoBackup = false
if fileExists($recoPath / "config.nims"):
  shell:
    cp ($recoPath)/config.nims ($recoPath)/config.nims_backup
  undoBackup = true
shell:
  cp runHackedAnalysis.nims ($recoPath)/config.nims
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
let outfile = "/mnt/1TB/CAST/CalibrationRuns_2014_MarlinHijackedCluster.h5"
shell:
  #"../../../Analysis/ingrid/reconstruction ../../../Tests/run_245_2014.h5" "--out" testfile.h5
  "../../../Analysis/ingrid/reconstruction /mnt/1TB/CAST/CalibrationRuns_2014_Raw.h5" "--out" ($outfile) "--runNumber 256"
