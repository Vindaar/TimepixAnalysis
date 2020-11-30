import sequtils, strutils, os, algorithm, strformat, sets, os
import std / sha1
import unittest
import nimhdf5
import shell
import seqmath

import helpers / testUtils

import ggplotnim
from ginger import ggColorHue

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

proc customFloatRepr(s: var string, f: float) =
  s.add &"{f:.1f}"

proc checkContent(h5f: H5FileObj, runNumber: int, withFadc = false): bool =
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
        var feDsetHashes = newJObject()
        for dsetName in FeDsets:
          let dset = h5f[(runGrp.name / "chip_" & $chip / dsetName).dset_str]
          feAttrs[dsetName] = dset.attrsToJson(withType = true)
          # read dset, hash content
          case dset.dtypeAnyKind
          of dkInt64:
            let data = dset[int64].mapIt(&"{it}")
            feDsetHashes[dsetName] = % $secureHash($(% data))
          of dkFloat64:
            let data = dset[float64].mapIt(&"{it:.1f}")
            feDsetHashes[dsetName] = % $secureHash($(% data))
          else:
            doAssert false, "what " & $dset
        #writeFile(&"hashes_fe_spetrum_run_{runNumber}.json", feDsetHashes.pretty)
        check compareJObjects(
          feDsetHashes,
          parseFile(&"hashes_fe_spetrum_run_{runNumber}.json")
        )
        #writeFile(&"fe_spectrum_attributes_run_{runNumber}.json", feAttrs.pretty)
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
                 runType: "rtBackground", num: 240),
                (inName: "run_241.h5", outName: "reco_241.h5",
                 runType: "rtCalibration", num: 241)]
  test "Default args":
    for r in runs:
      check fileExists(dataInPath/r.inName)
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
      res = shellVerbose:
        "../../Analysis/ingrid/reconstruction" ($(dataInPath/r.inName)) "--out" ($r.outName) "--runNumber" ($r.num)
      withH5(r.outName, "r"):
        check checkContent(h5f, r.num, withFadc = true)

      res = shellVerbose:
        "../../Analysis/ingrid/reconstruction" ($r.outName) "--only_charge"
      withH5(r.outName, "r"):
        check checkContent(h5f, r.num, withFadc = true)
      ## TODO: write test to check the fit of the Fe spectrum!
      ## TODO: write tests for remaining options
      removeFile(r.inName)

  test "Gas gain":
    let r = runs[1]
    var res = shellVerbose:
      "../../Analysis/ingrid/raw_data_manipulation ../../resources/TPAresources/gas_gain/Run_241_181022-16-16 --out raw_241_full.h5 --runType rtCalibration"
      "../../Analysis/ingrid/reconstruction raw_241_full.h5 --out reco_241_full.h5"
      "../../Analysis/ingrid/reconstruction reco_241_full.h5 --only_charge"
      "../../Analysis/ingrid/reconstruction reco_241_full.h5 --only_gas_gain"

    var h5f = H5open("reco_241_full.h5", "r")
    let shapePolya = h5f[("reconstruction/run_" & $241 / "chip_3/polya").dset_str].shape
    let polya = h5f["reconstruction/run_" & $r.num / "chip_3/polya", float64].reshape2D(shapePolya)
    let polyaF = h5f["reconstruction/run_" & $r.num / "chip_3/polyaFit", float64].reshape2D(shapePolya)

    let (x, p) = polya.split(SplitSeq.Seq2Col)
    let (xFit, pFit) = polyaF.split(SplitSeq.Seq2Col)
    let dfR = seqsToDf({ "x" : x,
                        "polya" : p })
    let dfFit = seqsToDf({ "x" : xFit,
                           "polya" : pFit })
    let dfAlt = bind_rows([("Polya", dfR), ("Fit", dfFit)],
                          id = "From")
    ggplot(dfAlt, aes("x", "polya")) +
      geom_histogram(data = dfAlt.filter(fn {c"From" == "Polya"}),
                     stat = "identity",
                     color = some(ggColorHue(2)[1])) +
      geom_line(data = dfAlt.filter(fn {c"From" == "Fit"}),
                color = some(ggColorHue(2)[0])) +
      ggsave("gasgain.pdf")
    discard h5f.close()
   # now check the gas gain
   # now remove the files
   # removeFile(r.outName)

  test "Gas gain 2014 data compare":
    var res = shellVerbose:
      "../../Analysis/ingrid/raw_data_manipulation ../../resources/TPAresources/gas_gain/525-Run151111_21-31-47 --out raw_525_full.h5 --runType back" #rtBackground"
      "../../Analysis/ingrid/reconstruction raw_525_full.h5 --out reco_525_full.h5"
      "../../Analysis/ingrid/reconstruction reco_525_full.h5 --only_charge"
      "../../Analysis/ingrid/reconstruction reco_525_full.h5 --only_gas_gain"

    var h5f = H5open("reco_525_full.h5", "r")
    let shapePolya = h5f[("reconstruction/run_" & $525 / "chip_0/polya").dset_str].shape
    let polya = h5f["reconstruction/run_" & $525 / "chip_0/polya", float64].reshape2D(shapePolya)
    let polyaF = h5f["reconstruction/run_" & $525 / "chip_0/polyaFit", float64].reshape2D(shapePolya)

    let (x, p) = polya.split(SplitSeq.Seq2Col)
    let (xFit, pFit) = polyaF.split(SplitSeq.Seq2Col)
    let dfR = seqsToDf({ "x" : x,
                        "polya" : p })
    let dfFit = seqsToDf({ "x" : xFit,
                           "polya" : pFit })
    let dfAlt = bind_rows([("Polya", dfR), ("Fit", dfFit)],
                          id = "From")
      # filter to max 2e4 electrons
      .filter(fn {c"x" <= 2.0e4})
    let dsetFit = h5f[("reconstruction/run_" & $525 / "chip_0/polyaFit").dset_str]
    let attrs = dsetFit.attrsToJson
    echo attrs.pretty
    let G = attrs["G"].getFloat
    let G_fit = attrs["G_fit"].getFloat
    let G_fitmean = attrs["G_fitmean"].getFloat


    ggplot(dfAlt, aes("x", "polya")) +
      geom_histogram(data = dfAlt.filter(fn {c"From" == "Polya"}),
                     stat = "identity",
                     color = some(ggColorHue(2)[1])) +
      geom_line(data = dfAlt.filter(fn {c"From" == "Fit"}),
                color = some(ggColorHue(2)[0])) +
      ggtitle(&"Polya fit of run 525; G = {G:.1f}, G_fit = {G_fit:.1f}, " &
        &"G_fitMean = {G_fitmean:.1f}") +
      ggsave("gasgain_2014.pdf")
    discard h5f.close()
    # now check the gas gain
    # now remove the files
    # removeFile(r.outName)

    ## TODO: make this a proper test!
