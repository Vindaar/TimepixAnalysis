import std / [os, sequtils]
import ingrid / [tos_helpers, ingrid_types]
import pkg / [nimhdf5, datamancer]
import ./io_helpers
import ./nn_predict

proc main(fname: string, model: string, ε: float) =
  var df = newDataFrame()
  df.add readCalibData(fname, "escape", 2.5, 3.5)
  df.add readCalibData(fname, "photo", 5.5, 6.5)

  let cutTab = calcNeuralNetCutValueTab(model, nkLocal, ε)
  #let (model, device) = loadModelMakeDevice(model)
  for (tup, subDf) in groups(df.group_by(["Type", "runNumber"])):
    let target = tup[0][1].toStr
    let run = tup[1][1].toInt
    # 1. predict all data
    let pred = predict(model, subDf)
    # 2. cut based on local prediction.
    ## XXX: ideally we would also look at the energy once we have `nkInterpolated`!
    let energy = if target == "escape": 3.0 else: 5.9
    var kept = 0
    for x in pred:
      if x >= cutTab[energy]:
        inc kept
    # 3. count each
    echo "Run: ", run, " for target: ", target
    echo "Keeping : ", kept, " of ", subDf.len, " = ", kept.float / subDf.len.float
    echo ""

when isMainModule:
  import cligen
  dispatch main
