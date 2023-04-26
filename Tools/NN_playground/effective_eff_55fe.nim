import std / [os, sequtils, tables, stats, strformat, strutils, options, random]
import ingrid / [tos_helpers, ingrid_types]
import pkg / [nimhdf5, datamancer, ggplotnim]
import ./io_helpers
import ./nn_predict
import std / sets
import ingrid / fake_event_generator

import unchained


## XXX: a bit annoying that this is here...
const CdlFile = "/home/basti/CastData/data/CDL_2019/calibration-cdl-2018.h5"
const RmsCleaningCut = 1.5

proc readGasGains(fnames: seq[string]): Table[int, float] =
  ## Returns the mean gas gains of every run
  result = initTable[int, float]()
  for c in fnames:
    if c.len == 0: continue
    withH5(c, "r"):
      for (run, chip, grp) in chipGroups(h5f):
        if chip != 3: continue
        var gains = newSeq[float]()
        let chp = h5f[grp.grp_str]
        for slice in iterGainSlicesFromDset(h5f, chp):
          gains.add slice.G
        result[run] = gains.mean()

proc linear(p: Tensor[float], x: float): float =
  result = p[0] * x + p[1]

import pkg / numericalnim except linspace

proc resultStr(pRes, pErr: Tensor[float], χ²dof: float): string =
  result.add &"χ²/dof = {χ²dof:.4f}\n"
  for i in 0 ..< pRes.size.int:
    result.add &"p[{i}] = {pRes[i]:.4e} ± {pErr[i]:.4e}\n"
  result.strip()

proc echoResult(pRes, pErr: Tensor[float], χ²dof: float) =
  echo resultStr(pRes, pErr, χ²dof)

proc high[T](x: Tensor[T]): int = x.size.int - 1

const CutVal = "CutVal" # the actual cut value required to get `ε`!
const Err = 0.025

proc toFitDf(params: Tensor[float], xmin, xmax: float): DataFrame =
  let xFit = seqmath.linspace(xmin, xmax, 10)
  let yFit = xFit.mapIt(linear(params, it))
  result = toDf({ "Gain" : xFit,
                  CutVal : yFit })

proc fitLine(df: DataFrame, cutTab: CutValueInterpolator, gainTab: Table[string, float]): (DataFrame, seq[float], string) =
  # 1. first filter data to only Run-2 (or Run-3) data
  let dataMin = df["Gain", float].min
  let df = df.filter(f{`Dataset` != "CDL"})
  if df.len == 0: return (newDataFrame(), newSeq[float](), "") ## Return early, input is CDL data
  # 2. get the mean CDL gas gain
  let tfKind = parseEnum[TargetFilterKind](df["Target"].unique.item(string))
  let energy = toXrayLineEnergy(tfKind)

  let cdlGain = gainTab[$tfKind]
  # 3. get the cut value corresponding to this target kind
  let cutVal = cutTab[energy]
  let cdlGainErr = 0.005

  echo "Fitting with: gain = ", cdlGain, " and cutVal = ", cutVal

  let cuts = df[CutVal, float]
  let cutsErr = df[CutVal, float].map_inline(x * Err)
  let xT = concat(df["Gain", float], [cdlGain].toTensor, axis = 0)
  let yT = concat(cuts, [cutVal].toTensor, axis = 0)
  let ey = concat(cutsErr, [cdlGainErr].toTensor, axis = 0)
  #let xT = df["Gain", float]
  #let yT = cuts
  #let ey = cutsErr
  let m = (yT[yT.high] - yT[0]) / (xT[xT.high] - xT[0])
  let b = 15.0
  let params = @[m, b].toTensor
  echo "fitting"
  let pRes = levmarq(linear, params,
                     xT, yT, yError = ey)
  let pErr = paramUncertainties(pRes, linear, xT, yT, ey)
  #let pRes = params
  #let pErr = params

  let xFit = numericalnim.linspace(dataMin, xT.max, 10)
  let yFit = xFit.mapIt(linear(pRes, it)).toTensor
  let dof = xT.len - params.len

  let yAtFit = xT.map_inline(linear(pRes, x))
  let χ²dof = chi2(yT, yAtFit, ey) / dof.float
  let resStr = resultStr(pRes, pErr, χ²dof)

  echo resStr
  result = (toFitDf(pRes, dataMin, xT.max), pRes.toSeq1D, resStr)

proc getSlope(params: seq[float], energy: float, cutTab: CutValueInterpolator): float =
  # 6. compute the effective cut value given the slope and y intercept
  let slopeScale = cutTab[5.9] / cutTab[energy] # reference energy. Scale slope according to the
                                                # ratio of the cut values that are present to "compress" or
                                                # "extend" the whole fit
  #echo "SCALING SLOPE: ", slopeScale
  let castToCdl = 4.8 / cutTab[5.9]
  let energyRatio = energy / 5.9
  result = params[0] #/ 0.6060 #/ energyRatio #(pow(slopeScale, energyRatio)) #(slopeScale*castToCdl)   #slopeScale)

proc yIntercept(params: seq[float],
                energy: float,
                cutTab: CutValueInterpolator,
                gainTab: Table[int, float],
                cdlGainTab: Table[string, float]): float =
  let cutVal = cutTab[energy]
  echo "CUT VALUE AT ENERGY: ", energy, " IS ", cutVal
  # 3. get the gas gain corresponding to the CDL data giving rise to the cut value
  let gainCdl = cdlGainTab[energy.toRefDset()]
  # 4. compute the y intercept of the slope for this cut value & energy
  let slope = getSlope(params, energy, cutTab)
  result = cutVal - gainCdl * slope

proc effCut(params: seq[float],
            yIntercept, energy: float, run: int,
            cutTab: CutValueInterpolator,
            gainTab: Table[int, float],
            cdlGainTab: Table[string, float]): float =
  let cutVal = cutTab[energy]
  # 3. get the gas gain corresponding to the CDL data giving rise to the cut value
  let gainCdl = cdlGainTab[energy.toRefDset()]
  # 5. get the gas gain of this run
  let gain = gainTab[run]
  # 6. compute the effective cut value given the slope and y intercept
  let slope = getSlope(params, energy, cutTab)
  result = gain * slope + yIntercept
  #let effCutVal = gain * params[0] + params[1] #yIntercept
  echo "Y INTERCEPT: ", yIntercept, " vs originalYInter ", params[1], " gainCDL: ", gainCdl, " effCutVal ", result, " compared: ", cutVal

proc predictCut(model: string, df: DataFrame, cutVal: float): (int, seq[float]) =
  ## Performs the prediction of the given input data and energy and applies the
  ## cut to the data. Returns the number of clusters that pass the cuts!
  let pred = predict(model, df)
  if pred.len < 100: return (-1, @[])
  # 2. cut based on local prediction.
  ## XXX: ideally we would also look at the energy once we have `nkInterpolated`!
  var kept = 0
  for x in pred:
    if x >= cutVal:
      inc kept
  result = (kept, pred)

proc predictCut(model: string, df: DataFrame, energy: float, cutTab: CutValueInterpolator): (int, seq[float]) =
  ## Performs the prediction of the given input data and energy and applies the
  ## cut to the data. Returns the number of clusters that pass the cuts!
  result = predictCut(model, df, cutTab[energy])

const RunNumbers = @[239, 307]
const Dataset = @["Run-2", "Run-3", "CDL"]

proc printKeptInfo(kept: int, run: int, target, dataType: string, dfLen: int, cutVal: float, pred: seq[float], ε: float) =
  echo "Run: ", run, " for target: ", target, " of data type ", dataType
  let effectiveEff = kept.float / dfLen.float
  let perc = pred.percentile((100 - (ε * 100.0)).round.int)
  echo "Keeping : ", kept, " of ", dfLen, " = ", effectiveEff, " at cutVal: ", cutVal, " vs desired percentile cut: ", perc
  echo ""

proc printEffStats(s: string, data: seq[float]) =
  echo s
  echo "\tmean = ", data.mean
  echo "\tstd = ", data.standardDeviation

proc evaluateFit(model: string, df: DataFrame,
                 cdlGainTab: Table[string, float],
                 gainTab: Table[int, float],
                 cutTab: CutValueInterpolator,
                 params: seq[float],
                 ε: float) = #slope: float) =
  ## Let's evaluate the effective efficiencies we get if we adjust the cut value
  ## using the linear fit slope.
  var effTab = initTable[(string, float), seq[float]]()
  for (tup, subDf) in groups(df.group_by(["DataType", "runNumber"])):
    let dataType = tup[0][1].toStr
    let target = subDf["Target"].unique().item(string)
    let run = tup[1][1].toInt
    let energy = toXrayLineEnergy(target)

    # calculate the effective cut value to use based on the gain fit
    let yIntercept = yIntercept(params, energy, cutTab, gainTab, cdlGainTab)
    let effCutVal = effCut(params, yIntercept, energy, run, cutTab, gainTab, cdlGainTab)

    let (kept, pred) = predictCut(model, subDf, effCutVal)
    if kept < 0: continue # too little data, skip this run
    # 3. count each
    printKeptInfo(kept, run, target, dataType, subDf.len, effCutVal, pred, ε)
    let effectiveEff = kept.float / subDf.len.float
    if Dataset[RunNumbers.lowerBound(run)] != "CDL":
      echo "Run number: ", run
      if (dataType, energy) notin effTab:
        effTab[(dataType, energy)] = newSeq[float]()
      effTab[(dataType, energy)].add effectiveEff
  for k, v in effTab:
    printEffStats($k, v)

proc readCdlGains(cdlFile: string): Table[string, float] =
  ## Returns a table that contains each target/filter kind and maps it to the
  ## mean gas gain of all
  let gainTab = readGasGains(@[cdlFile])
  var tab = initTable[string, seq[float]]()
  var weights = initTable[string, seq[float]]()
  withH5(cdlFile, "r"):
    for grp in items(h5f, recoGroupGrpStr().string, depth = 1):
      let tfKind = grp.attrs["tfKind", string]
      let run = grp.attrs["runNumber", int]
      let weight = h5f[(grp.name / "chip_3/eventNumber").dset_str].shape[0]
      if tfKind notin tab:
        tab[tfKind] = @[]
        weights[tfKind] = @[]
      tab[tfKind].add gainTab[run]
      weights[tfKind].add weight.float

  proc mean(x, w: seq[float]): float =
    for i in 0 ..< x.len:
      result += x[i] * w[i]
    result /= (w.sum())

  for k, v in tab:
    let w = weights[k]
    result[k] = v.mean(w)

#proc readCdlData(cdlFile: string, dsets: seq[InGridDsetKind], tfKind: TargetFilterKind): DataFrame =
#  ## Returns a table that contains each target/filter kind and maps it to the
#  ## mean gas gain of all
#  let gainTab = readGasGains(@[cdlFile])
#  result = newDataFrame()
#  withH5(cdlFile, "r"):
#    for grp in items(h5f, recoGroupGrpStr().string, depth = 1):
#      let thisKind = grp.attrs["tfKind", string]
#      let run = grp.attrs["runNumber", int]
#      if thisKind != tfKind: continue
#      for dset in dsets:
#        result[

proc plotCdlCutEffs(cutTab: CutValueInterpolator, cdlGainTab: Table[string, float]) =
  let energies = getXrayFluorescenceLines()
  var df = newDataFrame()
  for E in energies:
    let tf = E.toRefDset()
    df.add toDf({ "E" : E,
                  "Eff" : cutTab[E],
                  "Gain" : cdlGainTab[tf] })
  ggplot(df, aes("E", "Eff", color = "Gain")) +
    geom_point() +
    ggsave("/tmp/effective_cut_values.pdf")

  ggplot(df, aes("E", "Gain", color = "Eff")) +
    geom_point() +
    ggsave("/tmp/effective_cut_values_gain.pdf")


  proc to2750(cutVal, gain: float): float =
    #let params = @[-9.1225e-04, 8.0926e+00]
    let slope = -8.8244e-04
    result = cutVal + slope * (2500 - gain)

  df = df.mutate(f{"CutValGainNorm" ~ to2750(`Eff`, `Gain`)})
  echo df.pretty(-1)

  ggplot(df, aes("E", "CutValGainNorm", color = "Eff")) +
    geom_point() +
    ggsave("/tmp/effective_cut_values_eff_div_gain.pdf")

proc handleFakeData(ctx: LikelihoodContext, rnd: var Rand, fname: string, typ: string,
                    fakeDesc: FakeDesc,
                    run = -1,
                    nmc = 1000): DataFrame =
  let h5f = H5open(fname, "r")
  ## Override `nFake`
  var fakeDesc = fakeDesc
  fakeDesc.nFake = nmc
  var data = generateFakeData(ctx, rnd, h5f,
                              fakeDesc = fakeDesc,
                              run = run,
                              useCache = true)
    .cutXrayCleaning(fakeDesc.tfKind)

  result = data
  result["eventNumber"] = 0
  result["Type"] = $dtSignal
  result["DataType"] = typ
  result["Target"] = $fakeDesc.tfKind
  result["isFake?"] = true
  result.drop("likelihood")
  discard h5f.close()

from ginger import transparent
proc plotDatasets(df: DataFrame, plotPath: string) =
  let dsets = df.getKeys().filterIt(it notin ["Type", "Target", "eventNumber", "DataType", "runNumber", "isFake?", "File"])
  for dset in dsets:
    let d = dset
    echo "Dataset: ", d #, " for DF ", df
    let dfL = df.filter(f{float -> bool: idx(d) > percentile(col(d), 1) and idx(d) < percentile(col(d), 99)})
    if dfL.len > 0 and dfL[d].unique.len > 100: ## XXX: fix me!
      #echo dfL, " and ", dfL[d]
      ggplot(dfL, aes(dset, color = "DataType")) +
        geom_histogram(bins = 100, hdKind = hdOutline, fillColor = transparent, position = "identity", density = true, lineWidth = 1.0) +
        ggsave(&"{plotPath}/dsetPlots/{dset}_comparison.pdf")

      ggplot(dfL, aes(dset, color = "DataType")) +
        geom_density(fillColor = transparent, normalize = true, size = 1.0) +
        ggsave(&"{plotPath}/dsetPlots/{dset}_kde_comparison.pdf")

      #let labTab = { %~ "3.0" : 0,
      #               %~ "FakeRemove3.0" : 1,
      #               %~ "Fake3.0" : 2,
      #               %~ "CDL3.0" : 3,
      #               %~ "5.9" : 4,
      #               %~ "FakeFixed5.9" : 5,
      #               %~ "Fake5.9" : 6,
      #               %~ "CDL5.9" : 7 }.toTable()
      ggplot(dfL, aes(dset, fill = "DataType")) +
        ggridges("DataType", overlap = 1.5) + #, labelOrder = labTab) +
        geom_histogram(bins = 300, hdKind = hdOutline, position = "identity", density = true, color = "black", lineWidth = 1.0) +
        ggsave(&"{plotPath}/dsetPlots/{dset}_ridgeline_comparison.pdf")

      ggplot(dfL, aes(dset, fill = "DataType")) +
        ggridges("DataType", overlap = 1.5) + #, labelOrder = labTab) +
        geom_density(normalize = true, color = "black", size = 1.0) +
        ggsave(&"{plotPath}/dsetPlots/{dset}_ridgeline_kde_comparison.pdf")
    else:
      echo "[INFO]: Skipping dataset ", d, " as no elements left after removing quantile 1 and 99. Input column likely constant."


proc analyzeIntermediateEvents(model: string, df: DataFrame, simCutVal, realCutVal: float, plotPath: string) =
  ## Generates some plots of only those events in the 5.9 keV data that are between the
  ## simulated cut value and the real data based one.
  var df = df
  df["pred"] = predict(model, df)
  df = df.mutate(f{float -> string:
    "DataType" ~ (
      if `pred` >= min(simCutVal, realCutVal) and `pred` <= max(simCutVal, realCutVal):
        "intermediate"
      elif `pred` >= max(simCutVal, realCutVal):
        "both"
      else:
        "neither"
    )
  })
    .filter(f{`DataType` != "neither"})
  plotDatasets(df, plotPath / "intermediate")
  let evs = df.filter(f{`DataType` == "intermediate"})["eventNumber"]
  echo "Event numbers that are intermediate: ", evs
  echo "As arguments: ", evs.toTensor(int).toSeq1D.mapIt("--events " & $it).join(" ")

proc evaluateEffectiveEfficiency(model: string, df: DataFrame,
                                 cutTab: CutValueInterpolator,
                                 gainTab: Table[int, float],
                                 ε: float): DataFrame =
  ## Evaluate the effective efficiency of the data on a run by run basis as well as
  ## the mean of all runs (each run assumed the same weight)
  result = newDataFrame()
  #var effTab = { 2.98 : newSeq[float](), 5.89: newSeq[float]() }.toTable
  var effTab = initTable[(string, float), seq[float]]()
  for (tup, subDf) in groups(df.group_by(["DataType", "runNumber"])):
    let dataType = tup[0][1].toStr
    let target = subDf["Target"].unique().item(string)
    let run = tup[1][1].toInt

    let energy = toXrayLineEnergy(target)
    let (kept, pred) = predictCut(model, subDf, energy, cutTab)
    if kept < 0: continue # too little data, skip this run

    printKeptInfo(kept, run, target, dataType, subDf.len, cutTab[energy], pred, ε)
    let thisCutVal = pred.percentile(100 - (ε * 100.0).round.int)
    let Idx = RunNumbers.lowerBound(run)
    let effectiveEff = kept.float / subDf.len.float
    result.add toDf({ "Eff" : effectiveEff,
                      "DataType" : dataType,
                      "Run" : run,
                      "Gain" : gainTab[run],
                      CutVal : thisCutVal,
                      "Target" : target,
                      "isFake?" : subDf["isFake?"].unique().item(bool),
                      "Dataset" : Dataset[Idx],
    })
    if Dataset[RunNumbers.lowerBound(run)] != "CDL":
      if (dataType, energy) notin effTab:
        effTab[(dataType, energy)] = newSeq[float]()
      effTab[(dataType, energy)].add effectiveEff

  for k, v in effTab:
    printEffStats($k, v)
  #printEffStats("3.0 keV", effTab[2.98])
  #printEffStats("5.9 keV", effTab[5.89])

proc plotEfficiencyVsGain(df: DataFrame, plotPath, suffix: string) =
  ## Plot the efficiency vs the gas gain for the *real* 55Fe CAST data and the CDL data.
  echo df
  let df = df.filter(f{idx("isFake?") == false})
  let dfLong = df.mutate(f{float: "GainRel" ~ `Gain` / max(col("Gain"))})
    .gather(["Eff", "GainRel"], "Data", "Value")
  echo dfLong.pretty(-1)

  ggplot(dfLong, aes("Run", "Value", color = "Data", shape = "Target")) +
    geom_point() +
    ggtitle("NN cut efficiency for 55Fe data and normalized gas gain for all runs") +
    ggsave(&"{plotPath}/nn_cut_efficiency_and_gas_gains_{suffix}.pdf")

  ggplot(df, aes("Gain", "Eff", color = "Run")) +
    facet_wrap("Target") +
    geom_point() +
    ggtitle("NN cut efficiency for 55Fe data against gas gain for all runs") +
    ggsave(&"{plotPath}/nn_cut_efficiency_vs_gas_gains_{suffix}.pdf", width = 800, height = 480)

  ggplot(df, aes("Gain", "Eff", color = "Run")) +
    facet_wrap("Target") +
    geom_point() +
    ggtitle("NN cut efficiency for 55Fe data against gas gain for all runs") +
    ggsave(&"{plotPath}/nn_cut_efficiency_vs_gas_gains_{suffix}.pdf", width = 800, height = 480)

  ggplot(df, aes("Gain", CutVal, color = "Run")) +
    facet_wrap("Target") +
    geom_point() +
    ggtitle("NN cut value for 55Fe data against gas gain for all runs") +
    ggsave(&"{plotPath}/nn_cut_value_vs_gas_gains_{suffix}.pdf", width = 800, height = 480)

proc plotAllVsGasGain(df: DataFrame, plotPath: string, suffix: string) =
  ## now plot all cut values together
  ggplot(df, aes("Gain", CutVal, color = "DataType")) +
    facet_wrap("DataType") +
    geom_point() +
    ggsave(&"{plotPath}/all_gas_gain_vs_cutval{suffix}.pdf")

proc readRealCalib(f: string, typ: string, eLow, eHigh: float,
                   tfKind: TargetFilterKind,
                   validDsets: set[InGridDsetKind] = {}): DataFrame =
  if validDsets.card == 0:
    result = readCalibData(f, typ, eLow, eHigh, tfKind = some(tfKind))
  else:
    result = readCalibData(f, typ, eLow, eHigh, tfKind = some(tfKind), validDsets = validDsets)
  result = result.rename(f{"DataType" <- "CalibType"})
  result["isFake?"] = false
  ## Overwrite the `Target` field in the read data, because that is based on calling `toRefDset`
  ## for each cluster. Some may just slightly cross a boundary causing the rest of our code
  ## that uses the target as a factor to break
  result["Target"] = $tfKind

proc readCdlData(cdlFile: string,
                 tfKind: TargetFilterKind,
                 eLow, eHigh: float,
                 typeSuffix = "_CDL"): DataFrame =
  ## Reads the data fro the given tfKind filtering to the desired energy range
  const dsets = ValidReadDsets + {igEnergyFromCdlFit} - {igLikelihood}
  proc filterAfter(df: DataFrame, frm, to: float): DataFrame =
    let E = igEnergyFromCdlFit.toDset()
    #df
    df.filter(f{float -> bool: idx(E) >= frm and idx(E) <= to})
      #.drop([E])
    #result = newDataFrame()
    #for (tup, subDf) in groups(df.group_by("runNumber")):
    #  result.add subDf.filter(f{float -> bool: idx(E) >= frm and idx(E) <= to})
  let energy = toXrayLineEnergy(tfKind)
  result = filterAfter(
    readRealCalib(cdlFile, $energy & $typeSuffix,
                  NegInf, Inf, # deactivates energy cut in `cutXrayCleaning`
                  tfKind = tfKind,
                  validDsets = dsets),
    eLow, eHigh
  )

  if true: # tfKind == tfCuEpic0_9:
    let dfP = result.gather(["energyFromCdlFit", "energyFromCharge"], "Ds", "Energy")
      .filter(f{float: `Energy` > percentile(col("Energy"), 1) and `Energy` < percentile(col("Energy"), 99) })
    #let dfP = result.filter(f{float -> bool: `energyFromCharge` > percentile(col("energyFromCharge"), 1) and `energyFromCharge` < percentile(col("energyFromCharge"), 99) })
    #ggplot(dfP, aes("energyFromCharge", fill = "Ds")) +
    ggplot(dfP, aes("Energy", fill = "Ds")) +
      facet_wrap("runNumber") +
      geom_histogram(bins = 100, alpha = 0.5, hdKind = hdOutline, position = "identity") +
      ggtitle("Energy CDL vs charge by run for " & $tfKind) +
      ggsave(&"/home/basti/Sync/energy_for_data_tfKind_{tfKind}.pdf", width = 1200, height = 800)
  result = result.drop([igEnergyFromCdlFit.toDset()])

proc computeLMH(df: DataFrame, dset: string): DataFrame =
  let dfS = df.arrange(dset)
  let data = dfS[dset, float]
  let nnCut = dfS["NN", float]
  let perc33 = (data.len * 0.33).round.int
  result = toDf( {
    "NN" : @[nnCut[0        ..< perc33].mean,
             nnCut[perc33   ..< 2*perc33].mean,
             nnCut[2*perc33 ..< data.high].mean],
    dset : @[data[0        ..< perc33].mean,
              data[perc33   ..< 2*perc33].mean,
              data[2*perc33 ..< data.high].mean],
    "Perc" : @["low", "mid", "high"]
    }
  )

proc studyRmsTransverseGasGain(df, dfEff: DataFrame, model, plotPath, suffix: string) =
  var rmsT = newSeq[float]()
  var dfLoc = newDataFrame()
  for (tup, subDf) in groups(df.group_by(["DataType", "runNumber"])):
    let rmsT = subDf["rmsTransverse", float].percentile(98)
    let typ = tup[0][1].toStr
    let run = tup[1][1].toInt
    dfLoc.add toDf({ "DataType" : typ, "Run" : run, "rmsT" : rmsT})
  var dfS = dfEff.arrange(["DataType", "Run"])
  let dfLS = dfLoc.arrange(["DataType", "Run"])
  dfS["rmsT"] = dfLS["rmsT"]
  ggplot(dfS, aes("Gain", "rmsT", color = CutVal, shape = "DataType")) +
    geom_point() +
    ggsave(&"{plotPath}/rmsTransverse_vs_gain_{suffix}.pdf")
  ggplot(dfS, aes("rmsT", CutVal, color = "Gain", shape = "DataType")) +
    geom_point() +
    ggsave(&"{plotPath}/rmsTransverse_vs_cutVal_{suffix}.pdf")

  let dsets = df.getKeys().filterIt(it notin ["Type", "Target", "eventNumber", "DataType", "runNumber", "isFake?"])

  var df = df
  # get prediction for every cluster
  df["NN"] = predict(model, df)
  for dset in dsets:
    var dfSum = newDataFrame()
    when true: # not very useful, way too much data. Well for SGD we can see many parameters "banana shaped"
      ggplot(df, aes(dset, "NN", color = "DataType")) +
        geom_point(size = 1.0) +
        ggtitle(dset & " versus NN cut value") +
        ggsave(&"{plotPath}/nn_cut_scatter_{dset}_{suffix}.png", width = 1200, height = 800)
    for (tup, subDf) in groups(df.group_by(["DataType", "runNumber"])):
      let typ = tup[0][1].toStr
      let run = tup[1][1].toInt
      var dfLoc = computeLMH(subDf, dset)
      dfLoc["DataType"] = typ
      dfLoc["Run"] = run
      dfSum.add dfLoc
    ggplot(dfSum, aes(dset, "NN", color = "Perc", shape = "DataType")) +
      facet_wrap("DataType") +
      geom_point(size = 1.0) +
      ggtitle("Lower, mid and upper 33% quantiles of property: " & $dset) +
      ggsave(&"{plotPath}/nn_cut_vs_quantiled_scatter_{dset}_{suffix}.pdf", width = 1200, height = 800)

proc evaluateEffectiveEfficiencyByFakeRunCutVal(
  ctx: LikelihoodContext, rnd: var Rand, model: string, fnames: seq[string], cdlFile: string,
  ε: float,
  gainTab: Table[int, float],
  plotPath: string,
  run = -1
     ) =
  var df = newDataFrame()
  template dfFile(expr, file: untyped): untyped =
    var dfLoc = expr
    dfLoc["File"] = file
    dfLoc

  for c in fnames:
    # 3 keV
    df.add dfFile(readRealCalib(c, "3.0", 2.5, 3.5, tfAgAg6), c)
    # 5.9 keV
    df.add dfFile(readRealCalib(c, "5.9", 4.9, 6.9, tfMnCr12), c)

  if cdlFile.len > 0:
    let energies = getXrayFluorescenceLines()
    for i, E in energies:
      # if E > 6.5: continue
      let bins = concat(@[0.0], getEnergyBinning())
      echo "Ctuting energy to ", bins[i], " and ", bins[i+1]
      df.add dfFile(readCdlData(cdlFile, E.toRefTfKind(),
                                #0.4, 0.7,
                                bins[i], bins[i+1],
                                typeSuffix = ""),
                    cdlFile)

    #df.add dfFile(readCdlData(cdlFile, tfAgAg6, 2.5, 3.5), cdlFile)
    #df.add dfFile(readCdlData(cdlFile, tfMnCr12, 5.0, 7.0), cdlFile)
  #var cutTab = initTable[(string, int)]()
  var dfEff = newDataFrame()
  for (tup, subDf) in groups(df.group_by(["DataType", "runNumber"])):

    let typ = tup[0][1].toStr
    let runNumber = tup[1][1].toInt

    if run > 0 and run != runNumber: continue
    echo "run ", run, " raw events: ", subDf

    let target = parseEnum[TargetFilterKind](subDf["Target"].unique.item(string))
    # generate fake data for this data
    let file = subDf["File"].unique.item(string)
    let energy = toXrayLineEnergy(target).keV

    let energyDset = if file == cdlFile: igEnergyFromCdlFit else: igEnergyFromCharge
    let fakeDesc = if typ == "3.0": FakeDesc(tfKind: target,
                                             energyDset: energyDset,
                                             kind: fkGainDiffusion,
                                             λ: 2.2.cm)
                   else: FakeDesc(tfKind: target,
                                  energyDset: energyDset,
                                  kind: fkGainDiffusion)

    let dfFake = ctx.handleFakeData(
      rnd,
      file,
      "Fake" & $target,
      fakeDesc,
      run = runNumber,
      nmc = 5000
    )

    ## plots for the real + fake data for each case
    plotDatasets(bind_rows([("Real", subDf), ("Fake", dfFake)]), plotPath = plotPath / typ & "_" & $runNumber)

    let runCutVal = model.determineCutValue(dfFake, ε)
    # now predict output of `subDf`
    let (kept, pred) = predictCut(model, subDf, runCutVal)
    let eff = (kept.float / subDf.len.float)

    ## XXX: filter out events that are in range between cut of simulated data for 80% and
    ## cut of 80% for real data.
    ## -> write proc similar to `predictCut` that also determines cut value of real data.
    ##  Then iterate all events and extract the indices that are within the two efficiencies.
    let realCutVal = model.determineCutValue(subDf, ε)
    analyzeIntermediateEvents(model, subDf, runCutVal, realCutVal, plotPath)

    echo "Fake data for ", typ, " at run ", runNumber, " and energy ", energy, " target ", target, " has a cut value ", runCutVal, " and effective eff ", eff

    if eff < 0.5 and target == tfMnCr12:
      echo "Efficiency is less than 50% in run: ", runNumber, " of ", typ, " tar ", target, " and ", energy, " w/ ", fakeDesc
      #if true: quit()
    dfEff.add toDf({ "DataType" : typ,
                     "Run" : runNumber,
                     "Gain" : gainTab[runNumber],
                     "Eff" : eff })
  ggplot(dfEff, aes("Run", "Eff", color = "Gain", shape = "DataType")) +
    geom_point() +
    ggsave(&"{plotPath}/efficiency_based_on_fake_data_per_run_cut_val.pdf")

proc main(fnames: seq[string], model: string, ε: float,
          cdlFile: string = "",
          plotPath: string = "",
          evaluateFit = false,
          plotDatasets = false,
          run = -1) =
  ## The `cdlFile` currently must be the `CDL_2019_Reco.h5` file!
  let mlpDesc = initMLPDesc(model, plotPath)
  let gainTab = readGasGains(concat(fnames, @[cdlFile]))

  let ctx = initLikelihoodContext(CdlFile,
                                  year = yr2018,
                                  energyDset = igEnergyFromCharge,
                                  region = crSilver,
                                  timepix = Timepix1,
                                  morphKind = mkLinear) # morphing to plot interpolation

  #let cutTab = calcCutValueTab(ctx) #NeuralNetCutValueTab(model, nkLocal, ε)
  #let cutTab = calcNeuralNetCutValueTab(model, nkLocal, ε)
  let cdlGainTab = readCdlGains(cdlFile)

  # 1. plot the cut efficiency vs the gain data
  #plotCdlCutEffs(cutTab, cdlGainTab)


  ## XXX: this should maybe be replaced by a function call that also performs
  ## data reading? Because the question is what do we want to include in the plotting?
  #if plotDatasets:
  #  plotDatasets(df, mlpDesc.plotPath)
  #  return

  var rnd = initRand(1337)

  ## For each run, extract the diffusion, then generate fake data for that run,
  ## use it to compute a bespoke effective cut value for each run and then apply
  ## that to the real data and see what efficiencies we end up with!
  ctx.evaluateEffectiveEfficiencyByFakeRunCutVal(rnd, mlpDesc.path, fnames, cdlFile, ε, gainTab, mlpDesc.plotPath, run)
  if true: quit()
  when false:
    # 2. evaluate effective efficiency of the real 55Fe CAST data
    var dfReal = newDataFrame()
    block:
      var df = newDataFrame()
      for c in fnames:
        # 3 keV
        df.add readRealCalib(c, "3.0", 2.5, 3.5, tfAgAg6)
        # 5.9 keV
        df.add readRealCalib(c, "5.9", 5.0, 7.0, tfMnCr12)
      if cdlFile.len > 0:
        df.add readCdlData(cdlFile, tfAgAg6, 2.5, 3.5)
        df.add readCdlData(cdlFile, tfMnCr12, 5.0, 7.0)
        df.add readCdlData(cdlFile, tfCEpic0_6, 0.0, 0.4)
      # now evaluate and plot
      let dfEff = evaluateEffectiveEfficiency(model, df, cutTab, gainTab, ε)
      # 3. plot (apparent) dependency of gas gain vs the cut value
      plotEfficiencyVsGain(dfEff, mlpDesc.plotPath, "RealAndCdl")

      dfReal = df

      studyRmsTransverseGasGain(df, dfEff, mlpDesc.path, mlpDesc.plotPath, "RealAndCdl")
      dfReal.drop("NN")

    # 2b.
    #block:
    when false:
      var df = dfReal.filter(f{`DataType` != "3.0"}) # everything but escape data
      for c in fnames:
        let σTs = linspace(500.0, 700.0, 8)
        for σT in σTs:
          df.add ctx.handleFakeData(rnd, c, $σT, FakeDesc(tfKind: tfMnCr12, kind: fkDiffusionFromData, σT: σT))
      # now evaluate and plot
      let dfEff = evaluateEffectiveEfficiency(model, df, cutTab, gainTab, ε)
      # 3. plot (apparent) dependency of gas gain vs the cut value
      plotEfficiencyVsGain(dfEff, mlpDesc.plotPath, "Diffusion")
      studyRmsTransverseGasGain(df, dfEff, mlpDesc.path, mlpDesc.plotPath, "Diffusion")

    # 4d. behavior of all CDL targets vs fake data of all of them in gain and their distributions
    #block:
    when false:
      var df = newDataFrame()
      let energies = getXrayFluorescenceLines()
      for c in fnames:
        for i, E in energies:
          df.add ctx.handleFakeData(rnd, c, $E, FakeDesc(tfKind: E.toRefTfKind(), kind: fkDiffusionFromData))
          # add CDL data of this tfKind
          let bins = concat(@[0.0], getEnergyBinning())
          df.add readCdlData(cdlFile, E.toRefTfKind(),
                             bins[i], bins[i+1],
                             typeSuffix = "")

      # now evaluate and plot
      let dfEff = evaluateEffectiveEfficiency(model, df, cutTab, gainTab, ε)
      ## 3. plot (apparent) dependency of gas gain vs the cut value
      plotEfficiencyVsGain(dfEff, mlpDesc.plotPath, "CDL_Energies")
      plotAllVsGasGain(dfEff, mlpDesc.plotPath, "CDL_Energies")
      studyRmsTransverseGasGain(df, dfEff, mlpDesc.path, mlpDesc.plotPath, "CDL_Energies")

      ## XXX: Add CDL target/filter kind and read all CDL data and show it in the same
      ## plot as NN cut value for fake vs real CDL data.
      ## Also run only with 2017 or only 2018

    # 4. study fake data generation:
    # 4a. 'equivalent' fake events for 5.9 keV, 3.0 keV & escape photons (i.e. *not* real 3 keV data)
    #     together in the same gas gain plot?
    block:
      var df = dfReal
      for c in fnames:
        # add 3 keV fake data
        # like real 3 keV X-rays
        df.add ctx.handleFakeData(rnd, c, "3.0_Fake", FakeDesc(tfKind: tfAgAg6, kind: fkDiffusionFromData))
        # like escape photon events, i.e. absorption length of 5.9 keV as those are the origin
        df.add ctx.handleFakeData(rnd, c, "3.0_Like5.9", FakeDesc(tfKind: tfAgAg6, kind: fkDiffusionFromData, λ: 2.2.cm))
        # add 5.9 keV fake data
        df.add ctx.handleFakeData(rnd, c, "5.9_Fake", FakeDesc(tfKind: tfMnCr12, kind: fkDiffusionFromData))

        #df.add ctx.handleFakeData(rnd, c, "1.5_Fake", 1.5, FakeDesc(kind: fkDiffusionFromData))

        #df.add ctx.handleFakeData(rnd, c, "0.25_Fake", 0.254, FakeDesc(kind: fkDiffusionFromData))


      df.add ctx.handleFakeData(rnd, cdlFile, "5.9_CDL_Fake", FakeDesc(tfKind: tfMnCr12, kind: fkDiffusionFromData))
      df.add ctx.handleFakeData(rnd, cdlFile, "3.0_CDL_Fake", FakeDesc(tfKind: tfAgAg6, kind: fkDiffusionFromData))

      echo df
      # now evaluate and plot
      let dfEff = evaluateEffectiveEfficiency(model, df, cutTab, gainTab, ε)
      ## 3. plot (apparent) dependency of gas gain vs the cut value
      plotEfficiencyVsGain(dfEff, mlpDesc.plotPath, "FakeLikeReal")
      plotAllVsGasGain(dfEff, mlpDesc.plotPath, "FakeLikeReal")
      studyRmsTransverseGasGain(df, dfEff, mlpDesc.path, mlpDesc.plotPath, "FakeLikeReal")
      if true: quit()



    # 4b. behavior of NN cut value on absorption length by generating same energy data at different absorption
    #     lengths in same plot
    block:
      ## XXX: Need dataset plots for these!!!!!!
      var df = newDataFrame()
      for c in fnames:
        ## add 3 keV fake data
        ## like real 3 keV X-rays
        #df.add ctx.handleFakeData(rnd, c, "3.0_Fake", FakeDesc(tfKind: tfAgAg, kind: fkDiffusionFromData))
        ## like escape photon events, i.e. absorption length of 5.9 keV as those are the origin
        #df.add ctx.handleFakeData(rnd, c, "3.0_Like5.9", FakeDesc(tfKind: tfAgAg6, kind: fkDiffusionFromData, λ: 2.2))
        # add 5.9 keV fake data

        let absLengths = linspace(0.0, 4.0, 8)
        for l in absLengths:
          df.add ctx.handleFakeData(rnd, c, $l, FakeDesc(tfKind: tfMnCr12, kind: fkDiffusionFromData, λ: l.cm))

      # now evaluate and plot
      let dfEff = evaluateEffectiveEfficiency(model, df, cutTab, gainTab, ε)
      ## 3. plot (apparent) dependency of gas gain vs the cut value
      #plotEfficiencyVsGain(dfEff, mlpDesc.plotPath)
      plotAllVsGasGain(dfEff, mlpDesc.plotPath, "AbsLenDep")

    # 4c. behavior of NN cut value on diffusion coefficient by generating same energy data at different
    #     diffusion coefficients in the same plot
    block:
      var df = newDataFrame()
      for c in fnames:
        ## add 3 keV fake data
        ## like real 3 keV X-rays
        #df.add ctx.handleFakeData(rnd, c, "3.0_Fake", FakeDesc(tfKind: tfAgAg6, kind: fkDiffusionFromData))
        ## like escape photon events, i.e. absorption length of 5.9 keV as those are the origin
        #df.add ctx.handleFakeData(rnd, c, "3.0_Like5.9", FakeDesc(tfKind: tfAgAg6, kind: fkDiffusionFromData, λ: 2.2))
        # add 5.9 keV fake data

        let σTs = linspace(450.0, 700.0, 8)
        for σT in σTs:
          df.add ctx.handleFakeData(rnd, c, $σT, FakeDesc(tfKind: tfMnCr12, kind: fkDiffusionFromData, σT: σT))

      # now evaluate and plot
      let dfEff = evaluateEffectiveEfficiency(model, df, cutTab, gainTab, ε)
      ## 3. plot (apparent) dependency of gas gain vs the cut value
      #plotEfficiencyVsGain(dfEff, mlpDesc.plotPath)
      plotAllVsGasGain(dfEff, mlpDesc.plotPath, "Diffusions")


    #var df = newDataFrame()
    #for c in fnames:
    #  # 3 keV
    #  df.add readRealCalib(c, "0_3.0", 2.5, 3.5, tfAgAg6)
    #  df.add ctx.handleFakeData(rnd, c, "1_FakeRemove3.0", FakeDesc(tfKind: tfAgAg6, kind: fkRemovePixels))
    #  df.add ctx.handleFakeData(rnd, c, "2_Fake3.0", FakeDesc(tfKind: tfAgAg6, kind: fkDiffusionFromData)) # , λ: 3.3))
    #  # 5.9 keV
    #  df.add readRealCalib(c, "4_5.9", 5.0, 7.0, tfMnCr12)
    #  df.add ctx.handleFakeData(rnd, c, "5_FakeFixed5.9", FakeDesc(tfKind: tfMnCr12, kind: fkDiffusionFromDataFixedLoc, zDrift: 1.5))
    #  df.add ctx.handleFakeData(rnd, c, "6_Fake5.9", FakeDesc(tfKind: tfMnCr12, kind: fkDiffusionFromData)) # , λ: 2.0))
    #
    #  # 1 keV
    #  df.add ctx.handleFakeData(rnd, c, "8_Fake1.0", FakeDesc(tfKind: tfCuEpic2, kind: fkDiffusionFromData))#, λ: 1.0 / 7.0))
    #  df.add ctx.handleFakeData(rnd, c, "9_Fake0.7", 0.70, FakeDesc(tfKind: tfCuEpic2, kind: fkDiffusionFromData))#, λ: 1.0 / 7.0))



    #proc plotWithFit(df, dfFit: DataFrame, typ: string, text = "") =
    #  ggplot(df, aes("Gain", CutVal)) +
    #    geom_point(aes = aes(color = "Run")) +
    #    geom_errorbar(aes = aes(color = "Run",
    #                            yMin = f{idx(CutVal) - idx(CutVal) * Err},
    #                            yMax = f{idx(CutVal) + idx(CutVal) * Err})) +
    #    geom_line(data = dfFit, aes = aes("Gain", CutVal)) +
    #    annotate(text, 0.75, 0.25, font = font(10.0, family = "monospace")) +
    #    ggtitle("NN cut efficiency for 55Fe data against gas gain for " & $typ) +
    #    #ggsave(&"{mlpDesc.plotPath}/nn_cut_efficiency_vs_gas_gains.pdf", width = 800, height = 480)
    #    ggsave(&"/tmp/gain_cut_val_fit_{typ}.pdf", width = 800, height = 480)
    #
    #var params59: seq[float]
    #for (tup, subDf) in groups(dfEff.group_by("DataType")):
    #  let typ = tup[0][1].toStr
    #  let (dfFit, params, text) = fitLine(subDf, cutTab, cdlGainTab)
    #  if dfFit.len == 0: continue
    #  if typ == "5.9": params59 = params
    #  plotWithFit(subDf, dfFit, typ & "_direct_fit", text)
    #
    ## now plot all cut values together
    #ggplot(dfEff, aes("Gain", CutVal, color = "DataType")) +
    #  geom_point() +
    #  ggsave("/t/all_gas_gain_vs_cutval.pdf")
    #
    #if evaluateFit:
    #  # Parameter for Run-2
    #  evaluateFit(model, df, cdlGainTab, gainTab, cutTab, params59, ε) # -2.9794e-03)#  -1.3563e-03) #  ± 6.5452e-09
    #
    #  for (tup, subDf) in groups(dfEff.group_by("DataType")):
    #    # now recreate plot of 3 keV with 5.9 based parameters (except y intercept!)
    #    let typ = tup[0][1].toStr
    #    let energy = toXrayLineEnergy(subDf["Target"].unique().item(string))
    #    let yIntercept = yIntercept(params59, energy, cutTab, gainTab, cdlGainTab)
    #    let slope = getSlope(params59, energy, cutTab)
    #    let paramsNew = @[slope, yIntercept].toTensor
    #    plotWithFit(subDf,
    #                toFitDf(paramsNew,
    #                        dfEff["Gain", float].min,
    #                        dfEff["Gain", float].max),
    #                &"{typ}_with_adjusted_yInter")
    #
    #echo "3.0 keV CAST data cut value mean = ", dfEff.filter(f{`DataType` == "3.0"})["CutVal", float].mean

when isMainModule:
  import cligen
  dispatch main
