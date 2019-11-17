import sequtils, strutils, os, algorithm, strformat
import unittest
import nimhdf5
import shell

from ggplotnim import almostEqual
import seqmath

import json

const pwd = currentSourcePath().parentDir
# `dataInPath` contains the H5 files created by the tRawDataManipulation.nim test,
# which serve as the input for this test
const dataInPath = pwd / "raw_data_manipulation/"

suite "reconstruction":
  const runs = [(inName: "run_240.h5", outName: "reco_240.h5",
                 runType: "rtCalibration", num: 240),
                (inName: "run_241.h5", outName: "reco_241.h5",
                 runType: "rtBackground", num: 241)]
  test "Default args":
    for r in runs:
      let res = shellVerbose:
        "../../Analysis/ingrid/reconstruction" ($(dataInPath/r.inName)) "--out" ($r.outName) "--runType" ($r.runType)
      check fileExists(r.outName)
