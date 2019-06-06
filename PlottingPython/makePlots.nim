import shell
import strutils, os, strformat, sequtils

const kinds = ["all", "bronze", "silver", "gold"]

template onlyFadc(files: seq[string]): untyped {.dirty.} =
  let fadc = files.filterIt("fadc" in it and "scinti" notin it)[0]
  shell:
    "./Plotting/PyS_createBackgroundRate.py" `$none` "--files" `$fadc` "--ck_binning --fancy --title \"\""
  shell:
    "./Plotting/PyS_createBackgroundRate.py" `$none` "--files" `$fadc` "--ck_binning --fancy --title \"\" --log"     
template onlyScinti(files: seq[string]): untyped {.dirty.} =
  let scinti = files.filterIt("fadc" notin it and "scinti" in it)[0]
  shell:
    "./Plotting/PyS_createBackgroundRate.py" `$none` "--files" `$scinti` "--ck_binning --fancy --title \"\""
  shell:
    "./Plotting/PyS_createBackgroundRate.py" `$none` "--files" `$scinti` "--ck_binning --fancy --title \"\" --log"      
template onlyBoth(files: seq[string]): untyped {.dirty.} = 
  let both = files.filterIt("fadc" in it and "scinti" in it)[0]
  shell:
    "./Plotting/PyS_createBackgroundRate.py" `$none` "--files" `$both` "--ck_binning --fancy --title \"\""
  shell:
    "./Plotting/PyS_createBackgroundRate.py" `$none` "--files" `$both` "--ck_binning --fancy --title \"\" --log"    
template all(files: seq[string]): untyped {.dirty.} =
  shell:
    "./Plotting/PyS_createBackgroundRate.py" `$none` "--files" `$fadc` `$scinti` `$both` "--ck_binning --fancy --title \"\""
  shell:
    "./Plotting/PyS_createBackgroundRate.py" `$none` "--files" `$fadc` `$scinti` `$both` "--ck_binning --fancy --title \"\" --log"  

for k in kinds:
  echo "-----"
  var files: seq[string]
  for f in walkFiles(&"./*{k}*.h5"):
    files.add f

  let none = files.filterIt("fadc" notin it and "scinti" notin it)[0]
  onlyFadc(files)
  onlyScinti(files)
  onlyBoth(files)
  all(files)
  
