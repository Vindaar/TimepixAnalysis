import nimhdf5, ggplotnim
import std / [strutils, os, sequtils, sets, strformat]
import ingrid / [tos_helpers, ingrid_types]
import ingrid / calibration / [calib_fitting, calib_plotting]
import ingrid / calibration

let UseTeX = getEnv("USE_TEX", "false").parseBool

template toEDF*(data: seq[float], isCumSum = false): untyped =
  ## Computes the EDF of binned data
  var dataCdf = data
  if not isCumSum:
    seqmath.cumsum(dataCdf)
  let integral = dataCdf[^1]
  let baseline = min(data) # 0.0
  dataCdf.mapIt((it - baseline) / (integral - baseline))

import numericalnim / interpolate
import arraymancer
proc plotROC(dfB, dfC: DataFrame, outpath, tSuffix, suffix: string) =
  # 1. compute cumulative sum from each type of data that is binned in the same way
  # 2. plot cumsum, (1 - cumsum)
  when false:
    proc toInterp(df: DataFrame): InterpolatorType[float] =
      let data = df["riseTime", float].toSeq1D.sorted
      let edf = toEdf(data)
      ggplot(toDf(data, edf), aes("data", "edf")) +
        geom_line() +
        ggsave("/tmp/test_edf.pdf")
      result = newLinear1D(data, edf)
    let interpS = toInterp(dfC)
    let interpB = toInterp(dfB)
  proc doit(df: DataFrame) =
    let data = df["riseTime", float]
    let xs = linspace(data.min, data.max, 1000)
    let kde = kde(data)
  proc eff(data: seq[float], val: float, isBackground: bool): float =
    let cutIdx = data.lowerBound(val)
    result = cutIdx.float / data.len.float
    if isBackground:
      result = 1.0 - result

  let dataB = dfB["riseTime", float].toSeq1D.sorted
  let dataC = dfC["riseTime", float].toSeq1D.sorted
  var xs = newSeq[float]()
  var ysC = newSeq[float]()
  var ysB = newSeq[float]()
  var ts = newSeq[string]()
  for i in 0 ..< 200: # rise time
    xs.add i.float
    ysC.add dataC.eff(i.float, isBackground = false)
    ysB.add dataB.eff(i.float, isBackground = true)

  let sigEff = "signalEff"
  let backRej = "backSup"
  let df = toDf(xs, ysC, ysB)
    .rename(f{sigEff <- "ysC"},
            f{backRej <- "ysB"})
  ggplot(df, aes(sigEff, backRej)) +
    geom_line() +
    ggtitle("ROC curve of FADC rise time cut (only upper), ⁵⁵Fe vs. background in $#" % tSuffix) +
    xlab("Signal efficiency [%]") + ylab("Background suppression [%]") +
    themeLatex(fWidth = 0.9, width = 600, baseTheme = singlePlot, useTeX = UseTeX) +
    ggsave(outpath / ("fadc_rise_time_roc_curve_$#.pdf" % suffix), width = 800, height = 480)

  let dfG = df.gather([sigEff, backRej], "Type", "ys")
  ggplot(dfG, aes("xs", "ys", color = "Type")) +
    geom_line() +
    xlab("Rise time [ns]") + ylab("Signal efficiency / background suppression [%]") +
    themeLatex(fWidth = 0.9, width = 600, baseTheme = singlePlot, useTeX = UseTeX) +
    ggtitle("FADC rise time cut (only upper) efficiency in $#" % tSuffix) +
    ggsave(outpath / ("fadc_rise_time_efficiencies_$#.pdf" % suffix), width = 800, height = 480)

proc plotFallTimeRiseTime(df: DataFrame, outpath, tSuffix, suffix: string, isCdl, energyDep: bool, riseTimeHigh: float) =
  ## Given a full run of FADC data, create the
  ## Note: it may be sensible to compute a truncated mean instead
  # local copy filtered to maximum allowed rise time
  let df = df.filter(f{`riseTime` <= riseTimeHigh})

  proc plotDset(dset: string) =

    let title = if energyDep:
                   &"Comparison of FADC {dset} in photo- and escape peak data of $#" % tSuffix
                elif isCdl:
                   &"Comparison of FADC {dset} in CDL data by target/filter"
                else:
                   &"FADC signal {dset} in ⁵⁵Fe vs background data in $#" % tSuffix

    proc genOutfile(outpath, dset, suffix: string, isKde, isCdl, energyDep: bool): string =
      result = outpath / &"fadc_{dset}"
      if isKde:
        result = result & "_kde"
      if energyDep:
        result = result & "_energy_dep"
      elif isCdl:
        result = result & "_cdl"
      else:
        result = result & "_signal_vs_background"
      result = result & "_" & suffix & ".pdf"

    for (tup, subDf) in groups(group_by(df, "Type")):
      echo "============================== ", dset, " =============================="
      echo "Type: ", tup
      echo "Percentiles:"
      echo "\t 1-th: ", subDf[dset, float].percentile(1)
      echo "\t 5-th: ", subDf[dset, float].percentile(5)
      echo "\t50-th: ", subDf[dset, float].percentile(50)
      echo "\t mean: ", subDf[dset, float].mean
      echo "\t80-th: ", subDf[dset, float].percentile(80)
      echo "\t95-th: ", subDf[dset, float].percentile(95)
      echo "\t99-th: ", subDf[dset, float].percentile(99)
    df.writeCsv("/tmp/fadc_data_$#.csv" % suffix)
    #let df = df.filter(f{`Type` == "Cu-Ni-15kV"})
    ## Histogram of dataset
    ggplot(df, aes(dset, fill = "Type")) +
      geom_histogram(position = "identity", bins = 100, hdKind = hdOutline, alpha = 0.7) +
      ggtitle(title) +
      xlab(dset & " [ns]") +
      themeLatex(fWidth = 0.9, width = 600, baseTheme = singlePlot, useTeX = UseTeX) +
      ggsave(genOutfile(outpath, dset, suffix, false, isCdl, energyDep),
             width = 800, height = 480)
    ## KDE of dataset
    ggplot(df, aes(dset, fill = "Type")) +
      geom_density(normalize = true, alpha = 0.7, adjust = 2.0) +
      ggtitle(title) +
      xlab(dset & " [ns]") +
      themeLatex(fWidth = 0.9, width = 600, baseTheme = singlePlot, useTeX = UseTeX) +
      ggsave(genOutfile(outpath, dset, suffix, true, isCdl, energyDep),
             width = 800, height = 480)
    let dfLow = df.filter(f{`riseTime` < 200})
    ## Histogram of dataset of all entries < 200
    ggplot(dfLow, aes(dset, fill = "Type")) +
      geom_histogram(position = "identity", bins = 100, hdKind = hdOutline, alpha = 0.7) +
      ggtitle(title) +
      xlab(dset & " [ns]") +
      themeLatex(fWidth = 0.9, width = 600, baseTheme = singlePlot, useTeX = UseTeX) +
      ggsave(genOutfile(outpath, dset, "less_200_rise_" & suffix, false, isCdl, energyDep),
             width = 800, height = 480)
    ## KDE of dataset of all entries < 200
    ggplot(dfLow, aes(dset, fill = "Type")) +
      geom_density(normalize = true, alpha = 0.7, adjust = 2.0) +
      ggtitle(title) +
      xlab(dset & " [ns]") +
      themeLatex(fWidth = 0.9, width = 600, baseTheme = singlePlot, useTeX = UseTeX) +
      ggsave(genOutfile(outpath, dset, "less_200_rise_" & suffix, true, isCdl, energyDep),
             width = 800, height = 480)

    if isCdl:
      ## Ridgeline KDE of the CDL data
      let xrayRef = getXrayRefTable()
      var labelOrder = initTable[Value, int]()
      for idx, el in xrayRef:
        labelOrder[%~ el] = idx
      ggplot(dfLow, aes(dset, fill = "Type")) +
        ggridges("Type", overlap = 1.5, labelOrder = labelOrder) +
        geom_density(normalize = true, alpha = 0.7, adjust = 2.0, color = "black") +
        ggtitle(title) +
        themeLatex(fWidth = 0.9, width = 600, baseTheme = singlePlot, useTeX = UseTeX) +
        ggsave(genOutfile(outpath, dset, "ridgeline_less_200_rise_" & suffix, true, isCdl, energyDep),
               width = 800, height = 480)
    ## KDE of the different FADC settings
    ggplot(df, aes(dset, fill = "Settings")) +
      geom_density(normalize = true, alpha = 0.7, adjust = 2.0, color = "black") +
      ggtitle(dset & " of different FADC settings used") +
      themeLatex(fWidth = 0.9, width = 600, baseTheme = singlePlot, useTeX = UseTeX) +
      ggsave(genOutfile(outpath, dset, "different_fadc_amp_settings_" & suffix, true, isCdl, energyDep),
             width = 800, height = 480)
    ## KDE of the different runs in data
    ggplot(df, aes(dset, fill = factor("runNumber"))) +
      geom_density(normalize = true, alpha = 0.7, adjust = 2.0, color = "black") +
      ggtitle(dset & " of different runs") +
      themeLatex(fWidth = 0.9, width = 600, baseTheme = singlePlot, useTeX = UseTeX) +
      ggsave(genOutfile(outpath, dset, "different_runs_" & suffix, true, isCdl, energyDep),
             width = 800, height = 480)

    ## Mean values of the dataset
    let dfG = df.group_by("runNumber").summarize(f{float: dset << truncMean(col(dset).toSeq1D, 0.05)})
    let other = if dset == "riseTime": "fallTime" else: "riseTime"
    ggplot(dfG, aes(runNumber, dset, color = other)) +
      geom_point() +
      ggtitle(&"Mean {dset} in input data by run number, colored by {other}.") +
      themeLatex(fWidth = 0.9, width = 600, baseTheme = singlePlot, useTeX = UseTeX) +
      ggsave(outpath / &"fadc_mean_{dset}_$#.pdf" % suffix,
             width = 800, height = 480)

  plotDset("fallTime")
  plotDset("riseTime")


proc read(fname, typ: string, eLow, eHigh: float,
          isCdl: bool): DataFrame =
  var h5f = H5open(fname, "r")
  let fileInfo = h5f.getFileInfo()

  var peakPos = newSeq[float]()
  result = newDataFrame()
  for run in fileInfo.runs:
    if recoBase() & $run / "fadc" notin h5f: continue # skip runs that were without FADC
    var df = h5f.readRunDsets(
      run,
      fadcDsets = @["eventNumber",
                    "baseline",
                    "riseStart",
                    "riseTime",
                    "fallStop",
                    "fallTime",
                    "minvals",
                    "noisy",
                    "argMinval"]
    )
    if typ != "back": # if looking at any kind of X-ray data perform rough cuts
      let xrayRefCuts = getXrayCleaningCuts()
      let runGrp = h5f[(recoBase() & $run).grp_str]
      let tfKind = if not isCdl: tfMnCr12
                   else: runGrp.attrs["tfKind", string].parseEnum[:TargetFilterKind]()
      if isCdl:
        df["Type"] = $tfKind
      let cut = xrayRefCuts[$tfKind]
      let grp = h5f[(recoBase() & $run / "chip_3").grp_str]
      let passIdx = cutOnProperties(
        h5f,
        grp,
        crSilver, # try cutting to silver
        (toDset(igRmsTransverse), cut.minRms, cut.maxRms),
        (toDset(igEccentricity), 0.0, cut.maxEccentricity),
        (toDset(igLength), 0.0, cut.maxLength),
        (toDset(igHits), cut.minPix, Inf),
        (toDset(igEnergyFromCharge), eLow, eHigh)
      )

      let dfChip = h5f.readRunDsets(run, chipDsets = some((chip: 3, dsets: @["eventNumber"])))
      let allEvNums = dfChip["eventNumber", int]
      let evNums = passIdx.mapIt(allEvNums[it]).toSet
      # filter to allowed events & remove any noisy events
      df = df.filter(f{int: `eventNumber` in evNums and `noisy`.int < 1})
    df["runNumber"] = run

    df["Settings"] = "Setting " & $(@[75.5, 100.5, 120.5].lowerBound(run.float))
    result.add df
  if not isCdl:
    result["Type"] = typ
  echo result

proc main(calib: string, year: int,
          back = "",
          energyLow = 0.0, energyHigh = Inf,
          outpath = "Figs/statusAndProgress/FADC", # local!
          isCdl = false,
          riseTimeHigh = Inf) =
  let energyDep = if back.len > 0: false else: true
  if energyDep or not isCdl:
    let is2017 = year == 2017
    let is2018 = year == 2018
    if not is2017 and not is2018:
      raise newException(IOError, "The input file is neither clearly a 2017 nor 2018 calibration file!")

    let yearToRun = if is2017: 2 else: 3
    let tSuffix = "Run-$#" % $yearToRun
    let suffix = "run$#" % $yearToRun

    var df = newDataFrame()
    if energyDep:
      df.add read(calib, "escape", 2.5, 3.5, isCdl = false)
      df.add read(calib, "photo", 5.5, 6.5, isCdl = false)
    else:
      let dfC = read(calib, "⁵⁵Fe", energyLow, energyHigh, isCdl = false)
      let dfB = read(back, "back", 5.5, 6.5, isCdl = false)
      plotROC(dfB, dfC, outpath, tSuffix, suffix)
      df.add dfC
      df.add dfB

    plotFallTimeRiseTime(df, outpath, tSuffix, suffix, isCdl, energyDep, riseTimeHigh)
  else:
    let df = read(calib, "CDL", 0.0, Inf, isCdl = true)
    plotFallTimeRiseTime(df, outpath, "CDL", "CDL", isCdl, energyDep, riseTimeHigh)
when isMainModule:
  import cligen
  dispatch main
