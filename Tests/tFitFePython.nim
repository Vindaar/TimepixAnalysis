## test to check whether calling python from Nim works in
## calibration.nim to fit the Fe spectrum

import os
import nimhdf5
import strutils
import ingrid / [calibration, tos_helpers]

when isMainModule:
  if paramCount() < 1:
    quit("Please hand a h5 file!")

  let h5file = paramStr(1)
  var h5f = H5open(h5file, "rw")
  h5f.visitFile()

  for runStr, grp in runs(h5f):
    let runNumber = runStr.parseInt
    h5f.fitToFeSpectrum(runNumber, 3)
