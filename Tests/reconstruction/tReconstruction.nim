import sequtils, strutils, os, algorithm, strformat, sets, os, typeinfo
import std / sha1
import unittest
import nimhdf5
import shell
import seqmath

import helpers/testUtils

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
  check compareJObjects(runAttrs, parseFile(&"run_{runNumber}.json"))
  check "/reconstruction" / r in h5f
  for ch in ChipGroupsSet:
    check "/reconstruction" / r / ch in h5f
    for dset in DatasetSet:
      check  "/reconstruction" / r / ch / dset in h5f
    if withFadc:
      for dset in FadcDatasetSet:
        check  "/reconstruction" / r / "fadc" / dset in h5f
  if runNumber == 241:
    # check Fe Spectrum datasets and attributes
    # for the dataset we're going a different route. We're going to read all data
    # into a seq[T] and then (if float) round all values to nearest int? 2 digits?
    # not sure and convert to string. Then we simply calc the sha1 hash as a
    # checksum for the data and compare that.
    let FeDsets = ["", "Charge", "ChargePlot", "ChargePlotFit",
                   "Events", "Indices", "Plot", "PlotFit"].mapIt("FeSpectrum" & it)
    let runGrp = h5f[(&"/reconstruction/run_{runNumber}").grp_str]
    let centerChip = runGrp.attrs["centerChip", int]
    check centerChip == 3
    for chip in 0 ..< 6:
      if chip == centerChip:
        # check for FeSpec related dsets and content
        check FeDsets.filterIt(
          runGrp.name / "chip_" & $chip / it in h5f
        ).len == FeDsets.len
        # check attributes
        var feAttrs = newJObject()
        for dsetName in FeDsets:
          let dset = h5f[(runGrp.name / "chip_" & $chip / dsetName).dset_str]
          feAttrs[dsetName] = dset.attrsToJson(withType = true)
        #gwriteFile(&"fe_spectrum_attributes_run_{runNumber}.json", feAttrs.pretty)
        check compareJObjects(
          feAttrs,
          parseFile(&"fe_spectrum_attributes_run_{runNumber}.json")
        )
      else:
        check FeDsets.filterIt(
          runGrp.name / "chip_" & $chip / it notin h5f
        ).len == FeDsets.len

  # TODO: Write total charge test


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
      withH5(r.outName, "r"):
        check checkContent(h5f, r.num, withFadc = true)
      removeFile(r.outName)
      # now run different command line options
      # first of all just check whether all options actually work the way they should


      # giving run number should be exactly the same as above
      #res = shellVerbose:
      #  "../../Analysis/ingrid/reconstruction" ($(dataInPath/r.inName)) "--out" ($r.outName) "--runNumber" ($r.num)
      #withH5(r.outName, "r"):
      #  check checkContent(h5f, r.num, withFadc = true)
      #removeFile(r.outName)

      #res = shellVerbose:
      #  "../../Analysis/ingrid/reconstruction" ($(dataInPath/r.inName)) "--out" ($r.outName) "--only_charge"
      #withH5(r.outName, "r"):
      #  check checkContent(h5f, r.num, withFadc = true)

      #removeFile(r.inName)
