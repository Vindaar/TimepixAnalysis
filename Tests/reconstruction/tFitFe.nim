import ggplotnim, sequtils, strutils, strformat, ginger, json
import nimhdf5
import os
import nimpy
import unittest

import ingrid / private / python_utils
import ingrid / calibration / fit_functions
import ingrid / calibration / [calib_fitting, calib_plotting]

suite "Comparison of Fe spectrum fitting in Nim and Python":
  const runNumber = 241
  const chipNumber = 3
  var h5f = H5file("reco_241_full.h5", "r")
  # only nim
  let groupName = "/reconstruction/run_" & $runNumber & "/chip_" & $chipNumber & "/FeSpectrum"
  let groupNameCh = "/reconstruction/run_" & $runNumber & "/chip_" & $chipNumber & "/FeSpectrumCharge"
  var feDset = h5f[groupName.dsetStr]
  let feData = feDset[int64]
  var feDsetCh = h5f[groupNameCh.dsetStr]
  let feDataCh = feDsetCh[float64]

  #test "Fitting Nim":
  #  # call python function with data
  #  let resNim = fitFeSpectrum(feData)
  #  let data = resNim[0].mapIt(feSpectrumFunc(resNim[2], it))
  #  echo data
  #  echo type(data)
  #  let dfNim = seqsToDf({ "hist" : resNim[1],
  #                         "bins" : resNim[0],
  #                         "fit" : data})
  #  #echo dfNim.pretty(-1)
  #  ggplot(dfNim, aes("bins", "hist")) +
  #    geom_histogram(stat = "identity") +
  #    geom_line(aes(y = "fit"), color = some(ggColorHue(4)[3])) +
  #    ggsave(&"fe_spec_run_{runNumber}_chip_{chipNumber}_nim.pdf")
  #
  #  check true
  #
  #test "Fitting Python":
  #  # only Python
  #  let pyFitFe = pyImport("ingrid.fit_fe_spectrum")
  #  # call python function with data
  #  let res = pyFitFe.fitAndPlotFeSpectrum([feData], "", ".", runNumber,
  #                                         true)
  #  let hist = res[0].hist.toNimSeq(float)
  #  let bins = res[0].binning.toNimSeq(float)
  #
  #  let xFit = res[0].x_pl.toNimSeq(float)
  #  let yFit = res[0].y_pl.toNimSeq(float)
  #
  #  let df = seqsToDf({ "hist" : hist,
  #                      "bins" : bins,
  #                      "xFit" : xFit,
  #                      "yFit" : yFit })
  #  let params = res[0].popt.toNimSeq(float)
  #  let calcedY = pyFitFe.feSpectrumFunc(xFit, res[0].popt).toNimSeq(float)
  #  check calcedY == yFit
  #  ggplot(df, aes("bins", "hist")) +
  #    geom_histogram(stat = "identity") +
  #    geom_line(aes(y = "yFit"), color = some(ggColorHue(4)[1])) +
  #    ggsave(&"fe_spec_run_{runNumber}_chip_{chipNumber}_python.pdf")
  #
  #  check true

  test "Fitting both":
    # call python function with data
    let resNim = fitFeSpectrum(feData)

    # only Python
    let pyFitFe = pyImport("ingrid.fit_fe_spectrum")
    # call python function with data
    let res = pyFitFe.fitAndPlotFeSpectrum([feData], "", ".", runNumber,
                                           true)
    let hist = res[0].hist.toNimSeq(float)
    let bins = res[0].binning.toNimSeq(float)

    let xFit = res[0].x_pl.toNimSeq(float)
    let yFit = res[0].y_pl.toNimSeq(float)

    let dfPython = seqsToDf({ "hist" : hist,
                              "bins" : bins,
                              "xFit" : xFit,
                              "yFit" : yFit })
    let dfNim = seqsToDf({ "hist" : resNim.hist,
                           "bins" : resNim.binning,
                           "xFit" : resNim.xFit,
                           "yFit" : resNim.yFit})
    let df = bind_rows(("Nim", dfNim), ("Python", dfPython), id = "Lang")
    ggplot(df.filter(fn {`Lang` == "Nim"}), aes("bins", "hist")) +
      geom_histogram(stat = "identity", position = "identity") +
      geom_line(data = df, aes = aes(x = "xFit", y = "yFit", color = "Lang")) +
      ggsave(&"fe_spec_run_{runNumber}_chip_{chipNumber}_both.pdf")


    # call energy calib
    let resEnergy = fitEnergyCalib(resNim)
    plotFeEnergyCalib(resEnergy, 241, isPixel = true)

    # roughly compare energy caibration fits
    # fit parameters for fe spec are "very" different
    for i in 0 ..< resNim.pRes.len:
      check abs(resNim.pRes[i] - res[0].popt[i].to(float)) < 5.0
    # location of peaks is pretty good though, so energy calib yields satisfying
    # results
    for i in 0 ..< resEnergy.pRes.len:
      check abs(resEnergy.pRes[i] - res[1].popt[i].to(float)) < 0.1

  test "Compare binned data from Python and Nim":
    # compare the json files, produced from Run 241 (inserted:
    # Nim:
    # var jnode = newJObject()
    # jnode["hist"] = % data_tofit
    # jnode["bins"] = % bins_tofit
    # writeFile("/tmp/run_241_fe_spec_histo.json", pretty(jnode))
    # Python:
    # import json
    # f = open("/tmp/run_241_fe_spec_histo_python.json", "w")
    # json.dump({"hist" : data_tofit.tolist(),
    #            "bins" : bipns_tofit.tolist()}, f)
    # f.close()
    let jNim = parseJson(readFile("run_241_fe_spec_histo.json"))
    let jPython = parseJson(readFile("run_241_fe_spec_histo_python.json"))
    var idx = 0
    check jNim["hist"].getElems.mapIt(it.getFloat.int) == jPython["hist"].getElems.mapIt(it.getInt)
    check jNim["bins"] == jPython["bins"]

  #test "Compare energy calibration":
  #  let res = fitEnergyCalib(

  test "Fitting both charge":
    # call python function with data
    let resNim = fitFeSpectrumCharge(feDataCh)

    # only Python
    let pyFitFe = pyImport("ingrid.fit_fe_spectrum")
    # call python function with data
    let res = pyFitFe.fitAndPlotFeSpectrumCharge([feDataCh], "", ".", runNumber,
                                                 true)
    let hist = res[0].hist.toNimSeq(float)
    let bins = res[0].binning.toNimSeq(float)

    let xFit = res[0].x_pl.toNimSeq(float)
    let yFit = res[0].y_pl.toNimSeq(float)

    let dfPython = seqsToDf({ "hist" : hist,
                              "bins" : bins,
                              "xFit" : xFit,
                              "yFit" : yFit })
    let dfNim = seqsToDf({ "hist" : resNim.hist,
                           "bins" : resNim.binning,
                           "xFit" : resNim.xFit,
                           "yFit" : resNim.yFit})
    let df = bind_rows(("Nim", dfNim), ("Python", dfPython), id = "Lang")
    ggplot(df.filter(fn {`Lang` == "Nim"}), aes("bins", "hist")) +
      geom_histogram(stat = "identity", position = "identity") +
      geom_line(data = df, aes = aes(x = "xFit", y = "yFit", color = "Lang")) +
      ggsave(&"fe_spec_charge_run_{runNumber}_chip_{chipNumber}_both.pdf")

    # call energy calib
    let resEnergy = fitEnergyCalib(resNim, isPixel = false)
    plotFeEnergyCalib(resEnergy, 241, isPixel = false)
    # roughly compare energy caibration fits
    # fit parameters for fe spec are "very" different
    for i in 0 ..< resNim.pRes.len:
      check abs(resNim.pRes[i] - res[0].popt[i].to(float)) < 5.0
    # location of peaks is pretty good though, so energy calib yields satisfying
    # results
    for i in 0 ..< resEnergy.pRes.len:
      check abs(resEnergy.pRes[i] - res[1].popt[i].to(float)) < 0.1
