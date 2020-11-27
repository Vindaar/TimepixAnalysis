import os, strutils, sequtils

for file in walkDirRec(".", yieldFilter = {pcFile}):
  let fname = file.extractFilename
  if fname.startsWith("chipInfo") and not fname.endsWith("txt~"):
    let data = file.readFile().splitLines()
    let runPeriod = file.parentDir.parentDir
    let toAdd = "runPeriod: " & runPeriod
    let toWrite = concat(@[data[0], toAdd], data[1 .. ^1]).join("\n")
    echo "Updating ", toAdd, " with ", runPeriod
    writeFile(file, toWrite)
