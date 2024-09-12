import sequtils, strutils, os, algorithm, strformat
import unittest
import nimhdf5
import shell
import seqmath

import helpers/testUtils

const pwd = currentSourcePath().parentDir
const dataPwd = pwd / "../../resources/TPAresources/raw_data_manipulation/"

## WARNING: this test assumes that the `raw_data_manipulation` file was compiled
## on the current commit!
## We could add a `--commit` option to the main programs, which simply outputs the
## commit it was compiled with. We could then compare that with the current commit
## of the repo and fail if it doesn't match.

## In addition this test requires the `raw_data_manipulation` tool to be in `PATH`

#template checkReturn(arg: untyped): untyped {.dirty.} =
#  result = arg
#  if not result:
#    return

template checkRun(number: int, name: string, withFadc = false): untyped =
  ## helper proc, which checks whether our reduced run 240 HDF5 file actually
  ## looks the way we expect it to, if created with `--nofadc` flag.
  var h5f = H5open(name, "r")
  defer: discard h5f.close()
  # first check that "reconstruction" is no longer part of this file
  check "reconstruction" notin h5f
  # "runs" is in file
  check "runs" in h5f
  template getDset(n, num: untyped): untyped =
    h5f[("runs/run_" & $num / n).dset_str]
  template checkDims(n, `type`: untyped): untyped =
    let dset = getDset(n, number)
    check dset.shape == @[1001, 1]
    check dset.dtype == $`type`
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
    check almostEq(meanDuration, 2.121233475324665)
  else:
    check almostEq(meanDuration, 0.03090720322177823)
  const occSum240 = @[37307,
                      37640,
                      45560,
                      50225,
                      48741,
                      33775,
                      31383]
  const occSum241 = @[3890  ,
                      3444  ,
                      7240  ,
                      509507,
                      10224 ,
                      2990  ,
                      4176  ]
  const totSum240 = @[53602,
                      45025,
                      28457,
                      55978,
                      25170,
                      13452,
                      23118]
  const totSum241 = @[37094,
                      36949,
                      52573,
                      37036,
                      43633,
                      45386,
                      917  ]
  const averageHits240 = @[37.26973026973027,
                           37.6023976023976 ,
                           45.51448551448551,
                           50.17482517482517,
                           48.69230769230769,
                           33.741258741258741,
                           31.35164835164835]
  const averageHits241 = @[3.886113886113886,
                           3.44055944055944 ,
                           7.232767232767233,
                           508.998001998002 ,
                           10.21378621378621,
                           2.987012987012987,
                           4.171828171828172]
  const attrs = ["runNumber", "runTime", "runTimeFrames", "numChips", "shutterTime",
                 "runMode", "fastClock", "externalTrigger", "pathName", "dateTime", "shutterMode"]
  template checkChips(chip, numTot: int): untyped =
    let hits = getDset("chip_" & $chip / "Hits", number)
    check hits.shape == @[1001]
    let occ = getDset("chip_" & $chip / "Occupancy", number)
    check occ.shape == @[256, 256]
    let tot = getDset("chip_" & $chip / "ToT", number)
    check hits.shape == @[1001]
    let raw_x = getDset("chip_" & $chip / "raw_x", number)
    let raw_y = getDset("chip_" & $chip / "raw_y", number)
    let raw_ch = getDset("chip_" & $chip / "raw_ch", number)
    if number == 240:
      check occSum240[chip] == occ[int64].sum
      check totSum240[chip].uint16 == tot[uint16].sum
      check almostEq(averageHits240[chip], raw_x[special_type(uint8), uint8].mapIt(it.len).sum.float / 1001.0, 1e-8)
      check almostEq(averageHits240[chip], raw_y[special_type(uint8), uint8].mapIt(it.len).sum.float / 1001.0, 1e-8)
      check almostEq(averageHits240[chip], raw_ch[special_type(uint16), uint16].mapIt(it.len).sum.float / 1001.0, 1e-8)
    else:
      check occSum241[chip] == occ[int64].sum
      check totSum241[chip].uint16 == tot[uint16].sum
      check almostEq(averageHits241[chip], raw_x[special_type(uint8), uint8].mapIt(it.len).sum.float / 1001.0, 1e-8)
      check almostEq(averageHits241[chip], raw_y[special_type(uint8), uint8].mapIt(it.len).sum.float / 1001.0, 1e-8)
      check almostEq(averageHits241[chip], raw_ch[special_type(uint16), uint16].mapIt(it.len).sum.float / 1001.0, 1e-8)
    let grp = h5f[("runs/run_" & $number).grp_str]
    let rawGrp = h5f["runs".grp_str]
    # compare the attributes
    let rawAttrsJson = rawGrp.attrsToJson(withType = true)
    let grpAttrsJson = grp.attrsToJson(withType = true)
    ## Run the following line to update the JSON files (e.g. after a new attribute
    ## was added)
    #writeFile(pwd / "run_" & $number & "_attrs_rawgrp.json", rawAttrsJson.pretty)
    #writeFile(pwd / "run_" & $number & "_attrs_rungrp.json", grpAttrsJson.pretty)
    check compareJson(
      rawAttrsJson,
      parseFile(pwd / "run_" & $number & "_attrs_rawgrp.json"),
      exceptKeys = @["raw_data_manipulation_version", "raw_data_manipulation_compiled_on"]
    )
    check compareJson(
      grpAttrsJson,
      parseFile(pwd / "run_" & $number & "_attrs_rungrp.json"),
      exceptKeys = @["raw_data_manipulation_version", "raw_data_manipulation_compiled_on"]
    )

  for i in 0 .. 6:
    checkChips(i, 0)
  if withFadc:
    # TODO: create some check for the raw FADC data!
    let fadcPath = "runs/run_" & $number / "fadc"
    check fadcPath in h5f
    var fadcgrp = h5f[fadcPath.grp_str]
    check fadcPath / "eventNumber" in h5f
    check fadcPath / "raw_fadc" in h5f
    check fadcPath / "trigger_record" in h5f
    let expFadcAttrs = %* {
      "posttrig": 80,
      "n_channels": 0,
      "pedestal_run": 0,
      "sampling_mode": 0,
      "frequency": 2,
      "channel_mask": 15,
      "pretrig": 15000
    }
    let jattrs = fadcgrp.attrsToJson
    check compareJson(expFadcAttrs, jattrs)

suite "raw data manipulation":
  ## these tests check whether the raw data manipulation produces HDF5
  ## files as we expect them given certain command line arguments
  # first run raw data manipulation
  const runs = [(run: "Run_240_181021-14-54", outName: "run_240.h5",
                 runType: "rtBackground", num: 240),
                (run: "Run_241_181022-16-16", outName: "run_241.h5",
                 runType: "rtCalibration", num: 241)]
  # first raw data of both
  test "Without fadc: --nofadc":
    for r in runs:
      ## First remove existing file, if any
      if fileExists(r.outName):
        removeFile(r.outName)
      let res = shellVerbose:
        raw_data_manipulation -p ($(dataPwd/r.run)) "--nofadc" "--out" ($r.outName) "--runType" ($r.runType)
      check fileExists(r.outName)
      checkRun(r.num, r.outName)
      removeFile(r.outName)

  test "With fadc":
    for r in runs:
      ## First remove existing file, if any
      if fileExists(r.outName):
        removeFile(r.outName)
      let res = shellVerbose:
        raw_data_manipulation -p ($(dataPwd/r.run)) "--out" ($r.outName) "--runType" ($r.runType)
      check fileExists(r.outName)
      checkRun(r.num, r.outName, withFadc = true)
      # test does not delete file, is input for ../reconstruction/tReconstruction.nim
      #shell:
      #  rm ($r.outName)
