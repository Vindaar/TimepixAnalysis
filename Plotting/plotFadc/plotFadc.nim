import nimhdf5, ggplotnim
import std / [strutils, os, sequtils, sets, strformat]
import ingrid / [tos_helpers, ingrid_types]
import ingrid / calibration / [calib_fitting, calib_plotting]
import ingrid / calibration

proc plotFallTimeRiseTime(df: DataFrame, suffix: string, isCdl: bool) =
  ## Given a full run of FADC data, create the
  ## Note: it may be sensible to compute a truncated mean instead
  proc plotDset(dset: string) =

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
    ggplot(df, aes(dset, fill = "Type")) +
      geom_histogram(position = "identity", bins = 100, hdKind = hdOutline, alpha = 0.7) +
      ggtitle(&"Comparison of FADC signal {dset} in ⁵⁵Fe vs background data in $#" % suffix) +
      ggsave(&"Figs/statusAndProgress/FADC/fadc_{dset}_energy_dep_$#.pdf" % suffix)
    ggplot(df, aes(dset, fill = "Type")) +
      geom_density(normalize = true, alpha = 0.7, adjust = 2.0) +
      ggtitle(&"Comparison of FADC signal {dset} in ⁵⁵Fe vs background data in $#" % suffix) +
      ggsave(&"Figs/statusAndProgress/FADC/fadc_{dset}_kde_energy_dep_$#.pdf" % suffix)
    let df = df.filter(f{`riseTime` < 200})
    ggplot(df, aes(dset, fill = "Type")) +
      geom_histogram(position = "identity", bins = 100, hdKind = hdOutline, alpha = 0.7) +
      ggtitle(&"Comparison of FADC signal {dset} in ⁵⁵Fe vs background data in $#" % suffix) +
      ggsave(&"Figs/statusAndProgress/FADC/fadc_{dset}_energy_dep_less_200_rise_$#.pdf" % suffix)
    ggplot(df, aes(dset, fill = "Type")) +
      geom_density(normalize = true, alpha = 0.7, adjust = 2.0) +
      ggtitle(&"Comparison of FADC signal {dset} in ⁵⁵Fe vs background data in $#" % suffix) +
      ggsave(&"Figs/statusAndProgress/FADC/fadc_{dset}_kde_energy_dep_less_200_rise_$#.pdf" % suffix)

    if isCdl:
      let xrayRef = getXrayRefTable()
      var labelOrder = initTable[Value, int]()
      for idx, el in xrayRef:
        labelOrder[%~ el] = idx
      ggplot(df, aes(dset, fill = "Type")) +
        ggridges("Type", overlap = 1.5, labelOrder = labelOrder) +
        geom_density(normalize = true, alpha = 0.7, adjust = 2.0, color = "black") +
        ggtitle(&"Comparison of FADC signal {dset} in ⁵⁵Fe vs background data in $#" % suffix) +
        ggsave(&"Figs/statusAndProgress/FADC/fadc_{dset}_ridgeline_kde_energy_dep_less_200_rise_$#.pdf" % suffix)

    ggplot(df, aes(dset, fill = "Settings")) +
      geom_density(normalize = true, alpha = 0.7, adjust = 2.0, color = "black") +
      ggtitle(dset & " of different FADC settings used") +
      ggsave(&"Figs/statusAndProgress/FADC/fadc_{dset}_kde_different_fadc_ampb_settings_$#.pdf" % suffix)

    ggplot(df, aes(dset, fill = factor("runNumber"))) +
      geom_density(normalize = true, alpha = 0.7, adjust = 2.0, color = "black") +
      ggtitle(dset & " of different runs") +
      ggsave(&"Figs/statusAndProgress/FADC/fadc_{dset}_kde_different_runs_$#.pdf" % suffix)


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
    let xrayRefCuts = getXrayCleaningCuts()
    let runGrp = h5f[(recoBase() & $run).grp_str]
    let tfKind = if not isCdl: tfMnCr12
                 else: runGrp.attrs["tfKind", string].parseEnum[:TargetFilterKind]()
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
    if isCdl:
      df["Type"] = $tfKind

    df["Settings"] = "Setting " & $(@[80, 101, 121].lowerBound(run))
    result.add df
  if not isCdl:
    result["Type"] = typ
  echo result

proc main(fname: string, year: int,
          energyLow = 0.0, energyHigh = Inf,
          isCdl = false) =
  if not isCdl:
    var df = newDataFrame()
    df.add read(fname, "escape", 2.5, 3.5, isCdl = false)
    df.add read(fname, "photo", 5.5, 6.5, isCdl = false)

    let is2017 = year == 2017
    let is2018 = year == 2018
    if not is2017 and not is2018:
      raise newException(IOError, "The input file is neither clearly a 2017 nor 2018 calibration file!")

    let yearToRun = if is2017: 2 else: 3
    let suffix = "run$#" % $yearToRun
    plotFallTimeRiseTime(df, suffix, isCdl)
  else:
    let df = read(fname, "", 0.0, Inf, isCdl = true)
    plotFallTimeRiseTime(df, "CDL", isCdl)
when isMainModule:
  import cligen
  dispatch main
