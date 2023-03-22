import std / [os, sequtils, tables, stats, strformat]
import ingrid / [tos_helpers, ingrid_types]
from ingrid / calibration import iterGainSlicesFromDset
import pkg / [nimhdf5, datamancer, ggplotnim]
import ./io_helpers
import ./nn_predict

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

proc main(fnames: seq[string], model: string, ε: float,
          cdlFile: string = "",
          plotPath: string = "") =
  ## The `cdlFile` currently must be the `CDL_2019_Reco.h5` file!
  let mlpDesc = initMLPDesc(model, plotPath)

  var df = newDataFrame()
  for c in fnames:
    df.add readCalibData(c, "3.0", 2.5, 3.5)
    df.add readCalibData(c, "5.9", 5.0, 7.0)

  if cdlFile.len > 0:
    ## will read other runs too, but dominated by MnCr. Well, but I don't care about
    ## amount of data!
    #df.add readCalibData(cdlFile, "MnCr", 5.0, 7.0)
    df.add readCalibData(cdlFile, "3.0", 2.5, 3.5)
    df.add readCalibData(cdlFile, "5.9", 5.0, 7.0)


  let gainTab = readGasGains(concat(fnames, @[cdlFile]))

  let cutTab = calcNeuralNetCutValueTab(model, nkLocal, ε)
  #let (model, device) = loadModelMakeDevice(model)
  var dfEff = newDataFrame()
  for (tup, subDf) in groups(df.group_by(["CalibType", "runNumber"])):
    let target = tup[0][1].toStr
    let run = tup[1][1].toInt
    # 1. predict all data
    let pred = predict(model, subDf)
    # 2. cut based on local prediction.
    ## XXX: ideally we would also look at the energy once we have `nkInterpolated`!
    let energy = if target == "3.0": 3.0 else: 5.9
    #let energy = if target == "AgAg": 3.0 else: 5.9
    var kept = 0
    for x in pred:
      if x >= cutTab[energy]:
        inc kept
    # 3. count each
    echo "Run: ", run, " for target: ", target
    echo "Keeping : ", kept, " of ", subDf.len, " = ", kept.float / subDf.len.float
    echo ""
    dfEff.add toDf({ "Eff" : kept.float / subDf.len.float,
                     "Type" : target,
                     "Run" : run,
                     "Gain" : gainTab[run] })

  let dfEffLong = dfEff.mutate(f{float: "GainRel" ~ `Gain` / max(col("Gain"))})
    .gather(["Eff", "GainRel"], "Data", "Value")
  echo dfEffLong.pretty(-1)
  ggplot(dfEffLong, aes("Run", "Value", color = "Data", shape = "Type")) +
    geom_point() +
    ggtitle("NN cut efficiency for 55Fe data and normalized gas gain for all runs") +
    ggsave(&"{mlpDesc.plotPath}/nn_cut_efficiency_and_gas_gains.pdf")

  ggplot(dfEff, aes("Gain", "Eff", color = "Run")) +
    facet_wrap("Type") +
    geom_point() +
    ggtitle("NN cut efficiency for 55Fe data against gas gain for all runs") +
    ggsave(&"{mlpDesc.plotPath}/nn_cut_efficiency_vs_gas_gains.pdf", width = 800, height = 480)


when isMainModule:
  import cligen
  dispatch main
