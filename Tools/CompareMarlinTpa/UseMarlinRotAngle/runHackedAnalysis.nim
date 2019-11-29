import shell, strutils, os

let recoPath = "../../../Analysis/ingrid/"
var undoBackup = false
if fileExists($recoPath / "config.nims"):
  shell:
    cp ($recoPath)/config.nims ($recoPath)/config.nims_backup
  undoBackup = true
shell:
  cp runHackedAnalysis.nims ($recoPath)/config.nims
shell:
  nim c "-f --threads:on -d:danger -d:activateHijack" ($recoPath)/reconstruction.nim
if undoBackup:
  shell:
    cp ($recoPath)/config.nims_backup ($recoPath)/config.nims
else:
  shell:
    rm ($recoPath)/config.nims
shell:
  "../../../Analysis/ingrid/reconstruction ../../../Tests/run_245_2014.h5" "--out" testfile.h5
