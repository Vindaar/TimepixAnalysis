import sequtils, strutils, os, algorithm, strformat, sets
import unittest
import nimhdf5
import shell

from ggplotnim import almostEqual
import seqmath

import json

const pwd = currentSourcePath().parentDir
# `dataInPath` contains the H5 files created by the tRawDataManipulation.nim test,
# which serve as the input for this test
const dataInPath = pwd / "../raw_data_manipulation/"

const DatasetSet = ["skewnessTransverse",
                    "length",
                    "hits",
                    "centerY",
                    "fractionInTransverseRms",
                    "centerX",
                    "eventNumber",
                    "width",
                    "ToT",
                    "rmsTransverse",
                    "skewnessLongitudinal",
                    "x",
                    "kurtosisLongitudinal",
                    "eccentricity",
                    "rmsLongitudinal",
                    "y",
                    "kurtosisTransverse",
                    "lengthDivRmsTrans",
                    "rotationAngle",
                    "sumTot"].toHashSet
const ChipGroupsSet = (toSeq(0 .. 6).mapIt("chip_" & $it)).toHashSet
const RunGroupsSet = toHashSet(["run_240", "run_241"])
const FadcDatasetSet = toHashSet(["fadc_data",
                                  "minvals",
                                  "noisy",
                                  "eventNumber"])

proc checkContent(h5f: H5FileOBj, runNumber: int, withFadc = false): bool =
  template check(cond: untyped): untyped =
    if cond:
      result = true
    else:
      echo "Failed: ", astToStr(cond), " was ", cond
      return false

  check "/reconstruction" in h5f
  let r = "run_" & $runNumber
  let runAttrs = attrsToJson(h5f[("reconstruction" / r).grp_str], withType = true)
  # screw it, we compare by string. For some reason it appears `==` for Json doesn't
  # properly handle comparisons?!
  #writeFile(&"run_{runNumber}.json", runAttrs.pretty)
  check runAttrs.pretty == parseJson(readFile(&"run_{runNumber}.json")).pretty
  check "/reconstruction" / r in h5f
  for ch in ChipGroupsSet:
    check "/reconstruction" / r / ch in h5f
    for dset in DatasetSet:
      check  "/reconstruction" / r / ch / dset in h5f
    if withFadc:
      for dset in FadcDatasetSet:
        check  "/reconstruction" / r / "fadc" / dset in h5f

suite "reconstruction":
  const runs = [(inName: "run_240.h5", outName: "reco_240.h5",
                 runType: "rtCalibration", num: 240),
                (inName: "run_241.h5", outName: "reco_241.h5",
                 runType: "rtBackground", num: 241)]
  test "Default args":
    for r in runs:
      var res = shellVerbose:
        "../../Analysis/ingrid/reconstruction" ($(dataInPath/r.inName)) "--out" ($r.outName)
      check res[1] == 0
      check fileExists(r.outName)

      # get all groups and datasets in the files
      var h5f = H5file(r.outName, "r")
      check checkContent(h5f, r.num, withFadc = true)
      check h5f.close() >= 0
      removeFile(r.inName)

      # now run different command line options
      res = shellVerbose:
        "../../Analysis/ingrid/reconstruction" ($(dataInPath/r.inName)) "--out" ($r.outName)

      removeFile(r.outName)
