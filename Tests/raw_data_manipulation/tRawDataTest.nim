#import ingrid / ingrid_types
#import ingrid / tos_helpers
#import ingridDatabase / databaseDefinitions
import sequtils, strutils, os, algorithm, strformat
import unittest
#import json, sugar
#from ingrid / reconstruction import recoEvent
import nimhdf5
import shell

from ggplotnim import almostEqual
import seqmath

const pwd = currentSourcePath().parentDir
const dataPwd = pwd / "../../resources/TPAresources/raw_data_manipulation/"

## WARNING: this test assumes that the `raw_data_manipulation` file was compiled
## on the current commit!
## We could add a `--commit` option to the main programs, which simply outputs the
## commit it was compiled with. We could then compare that with the current commit
## of the repo and fail if it doesn't match.

proc checkRun_wo_FADC(number: int, name: string): bool =
  ## helper proc, which checks whether our reduced run 240 HDF5 file actually
  ## looks the way we expect it to, if created with `--nofadc` flag.
  var h5f = H5file(name, "r")
  defer: discard h5f.close()
  # first check that "reconstruction" is no longer part of this file
  result = "reconstruction" notin h5f
  # "runs" is in file
  result = "runs" in h5f
  template getDset(name, number: untyped): untyped =
    h5f[("runs/run_" & $number / name).dset_str]
  template checkDims(name, `type`: untyped): untyped =
    let dset = getDset(name, number)
    result = dset.shape == @[1001, 1]
    result = dset.dtype == $`type`
  checkDims("eventDuration", float64)
  checkDims("eventNumber", int64)
  checkDims("fadcReadout", int64)
  checkDims("fadcTriggerClock", int64)
  checkDims("szint1ClockInt", int64)
  checkDims("szint2ClockInt", int64)
  checkDims("szint1ClockInt", int64)
  checkDims("timestamp", int64)
  checkDims("useHvFadc", int64)
  let meanDuration = sum(h5f["runs/run_" & $number / "eventDuration", float]) / 1001.0
  echo meanDuration
  if number == 240:
    result = almostEqual(meanDuration, 2.121233475324665)
  else:
    result = almostEqual(meanDuration, 0.03090720322177823)
  var outf = open("tot.txt", fmAppend)
  const occSum240 = @[559605,
                      564600,
                      683400,
                      753375,
                      731115,
                      506625,
                      470745]
  const occSum241 = @[42790  ,
                      37884  ,
                      79640  ,
                      5604577,
                      112464 ,
                      32890  ,
                      45936  ]
  template checkChips(chip, numTot: int): untyped =
    let hits = getDset("chip_" & $chip / "Hits", number)
    result = hits.shape == @[1001, 1]
    let occ = getDset("chip_" & $chip / "Occupancy", number)
    result = occ.shape == @[256, 256]
    if number == 240:
      result = occSum240[chip] == occ[int64].sum
      echo "??? is ", result, " from ", occSum240[chip], " and ", occ[int64].sum
    else:
      echo "ok ", chip
      result = occSum241[chip] == occ[int64].sum
      echo "??? is ", result
    let tot = getDset("chip_" & $chip / "ToT", number)
    outf.write("Tot " & ($chip / $number) & " : " & $tot[uint16].sum & "\n")

    #result = hits.shape == @[1001, 1]
    #let raw_x = getDset("chip_" & $chip / "raw_x", number)
    #result = hits.shape == @[1001, 1]
    #let raw_y = getDset("chip_" & $chip / "raw_y", number)
    #result = hits.shape == @[1001, 1]
    #let raw_ch = getDset("chip_" & $chip / "raw_ch", number)
    #result = hits.shape == @[1001, 1]
  for i in 0 .. 6:
    checkChips(i, 0)
  outf.close()


suite "raw data manipulation":
  ## these tests check whether the raw data manipulation produces HDF5
  ## files as we expect them given certain command line arguments
  # first run raw data manipulation
  const runs = [(run: "Run_240_181021-14-54", outName: "run_240.h5",
                 runType: "rtCalibration", num: 240),
                (run: "Run_241_181022-16-16", outName: "run_241.h5",
                 runType: "rtBackground", num: 241)]
  # first raw data of both
  test "Without fadc: --nofadc":
    for r in runs:
      let res = shellVerbose:
        "../../Analysis/ingrid/raw_data_manipulation" ($(dataPwd/r.run)) "--nofadc" "--out" ($r.outName) "--runType" ($r.runType)
      check fileExists(r.outName)
      check checkRun_wo_FADC(r.num, r.outName)
