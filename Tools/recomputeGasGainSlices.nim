import shell, strutils, os

#[
Note, this can crash in case two processes running in parallel try to open the ingrid
database at the same time. Instead of failing on access to the db, we should retry
accessing it multiple times with a short wait in between.
]#

const path = "/mnt/1TB/CAST/"
const files = ["2017/CalibrationRuns2017_Reco.h5",
               "2017/DataRuns2017_Reco.h5",
               "2018_2/CalibrationRuns2018_Reco.h5",
               "2018_2/DataRuns2018_Reco.h5"]

proc runFile(f: string) =
  let (res, err, code) = shellVerboseErr:
    #reconstruction ($f) "--only_gas_gain"
    reconstruction ($f) "--only_energy_from_e"
    #reconstruction ($f) "--only_gain_fit"
  if code != 0:
    raise newException(Exception, "Error running the given command for " & $f)

proc main =
  createDir("out")
  var thr: array[files.len, Thread[string]]
  for i in 0 ..< files.len:
    createThread(thr[i], runFile, path / files[i])
  joinThreads(thr)

when isMainModule:
  main()
