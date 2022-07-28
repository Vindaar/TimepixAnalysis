import sequtils, strutils, os, algorithm, strformat, sets, os
import std / sha1
import unittest
import nimhdf5
import shell
import seqmath

import helpers / testUtils
from ingrid / calibration / fit_functions import polyaImpl
from ingrid / calibration import calcMeanGainFit

import ggplotnim
from ginger import ggColorHue

const pwd = currentSourcePath().parentDir
# `dataInPath` contains the H5 files created by the tRawDataManipulation.nim test,
# which serve as the input for this test
const dataInPath = pwd / "../raw_data_manipulation/"
const resourcePath = pwd / "../../resources/TPAresources/"

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

proc toDsetSuffix(idx, interval: int): string =
  result = &"_{idx}_{interval}"

proc plotGasGain(h5f: H5File, grp: H5Group, run: int) =
  let polyaGrp = h5f[(grp.name / "polyaDsets").grp_str]
  let numIntervals = polyaGrp.attrs["numGasGainSlices", int]
  let interval = polyaGrp.attrs["gasGainInterval", int]

  for idx in 0 ..< numIntervals:
    let polyaDset = h5f[(polyaGrp.name / "polya" & $toDsetSuffix(idx, interval)).dset_str]
    let shapePolya = polyaDset.shape
    let polya = polyaDset[float64].reshape2D(shapePolya)

    # Read fit parameters from polya
    let N = polyaDset.attrs["N", float]
    let G_fit = polyaDset.attrs["G_fit", float]
    let theta = polyaDset.attrs["theta", float]

    let (x, p) = polya.split(SplitSeq.Seq2Col)
    let xFit = linspace(x.min, x.max, 1000)
    let pFit = xFit.mapIt(polyaImpl(@[N, G_fit, theta], it))
    let dfR = toDf({ "x" : x,
                        "polya" : p })
    let dfFit = toDf({ "x" : xFit,
                           "polya" : pFit })
    let dfAlt = bind_rows([("Polya", dfR), ("Fit", dfFit)],
                          id = "From")
      # filter to max 2e4 electrons
      .filter(fn {c"x" <= 2.0e4})
    ## Compute the mean of the polya data & data described by fit
    let G = histMean(x, p)
    let G_fitmean = calcMeanGainFit(xFit, pFit)
    ggplot(dfAlt, aes("x", "polya")) +
      geom_histogram(data = dfAlt.filter(fn {c"From" == "Polya"}),
                     stat = "identity",
                     fillColor = ggColorHue(2)[1]) +
      geom_line(data = dfAlt.filter(fn {c"From" == "Fit"}),
                color = ggColorHue(2)[0]) +
      ggtitle(&"Polya fit of run {run}; G = {G:.1f}, G_fit = {G_fit:.1f}, " &
        &"G_fitMean = {G_fitmean:.1f}") +
      ggsave(&"gasgain_run_{run}_slice_{idx}.pdf")

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
  check compareJson(runAttrs, parseFile(pwd / &"run_{runNumber}.json"))
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
        ## Run the following two commands to regenaret the JSON containing the SHA
        ## hashes, if the reconstruction changes
        #writeFile(&"{pwd}/hashes_fe_spetrum_run_{runNumber}.json", feDsetHashes.pretty)
        #writeFile(&"{pwd}/fe_spectrum_attributes_run_{runNumber}.json", feAttrs.pretty)
        check compareJson(
          feDsetHashes,
          parseFile(pwd / &"hashes_fe_spetrum_run_{runNumber}.json")
        )
        check compareJson(
          feAttrs,
          parseFile(pwd / &"fe_spectrum_attributes_run_{runNumber}.json")
        )
      else:
        check FeDsets.filterIt(
          runGrp.name / "chip_" & $chip / it notin h5f
        ).len == FeDsets.len

      ## XXX: generate plots for the Fe55 fit & store it. THat way if we do something
      ## that changes the plot we can visually inspect whether it's a regression

  # TODO: Write total charge test

suite "reconstruction":
  const runs = [(inName: "run_240.h5", outName: pwd / "reco_240.h5",
                 runType: "rtBackground", num: 240),
                (inName: "run_241.h5", outName: pwd / "reco_241.h5",
                 runType: "rtCalibration", num: 241)]
  test "Default args":
    for r in runs:
      check fileExists(dataInPath/r.inName)
      # remove existing file
      shell:
        rm ($r.outName)
      var res = shellVerbose:
        reconstruction ($(dataInPath/r.inName)) "--out" ($r.outName)
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
        reconstruction ($(dataInPath/r.inName)) "--out" ($r.outName) "--runNumber" ($r.num)
      withH5(r.outName, "r"):
        check checkContent(h5f, r.num, withFadc = true)

      res = shellVerbose:
        reconstruction ($r.outName) "--only_charge"
      withH5(r.outName, "r"):
        check checkContent(h5f, r.num, withFadc = true)
      ## TODO: write test to check the fit of the Fe spectrum!
      ## TODO: write tests for remaining options
      removeFile(r.inName)

  test "Gas gain":
    # delete existing raw / reco files
    shell:
      rm ($pwd"/raw_241_full.h5")
      rm ($pwd"/reco_241_full.h5")
    let r = runs[1]
    var res = shellVerbose:
      raw_data_manipulation ($resourcePath"gas_gain/Run_241_181022-16-16") "--out" ($pwd"/raw_241_full.h5") "--runType rtCalibration"
      reconstruction ($pwd"/raw_241_full.h5") "--out" ($pwd"/reco_241_full.h5")
      reconstruction ($pwd"/reco_241_full.h5") "--only_charge"
      reconstruction ($pwd"/reco_241_full.h5") "--only_gas_gain"

    var h5f = H5open(pwd / "reco_241_full.h5", "r")
    let grp = h5f[("reconstruction/run_" & $241 / "chip_3").grp_str]
    plotGasGain(h5f, grp, r.num)

    discard h5f.close()
   # now check the gas gain
   # now remove the files
   # removeFile(r.outName)

  test "Gas gain 2014 data compare":
    ## Remove the raw files if they exist
    shell:
      rm ($pwd"/raw_525_full.h5")
      rm ($pwd"/reco_525_full.h5")
    var res = shellVerbose:
      raw_data_manipulation ($resourcePath"/gas_gain/525-Run151111_21-31-47") "--out" ($pwd"/raw_525_full.h5") "--runType back" #rtBackground"
      reconstruction ($pwd"/raw_525_full.h5") "--out" ($pwd"/reco_525_full.h5")
      reconstruction ($pwd"/reco_525_full.h5") "--only_charge"
      reconstruction ($pwd"/reco_525_full.h5") "--only_gas_gain"

    var h5f = H5open($pwd / "reco_525_full.h5", "r")

    let grp = h5f[("reconstruction/run_" & $525 / "chip_0").grp_str]
    plotGasGain(h5f, grp, 525)

    discard h5f.close()
    # now check the gas gain
    # now remove the files
    # removeFile(r.outName)

    ## TODO: make this a proper test!
