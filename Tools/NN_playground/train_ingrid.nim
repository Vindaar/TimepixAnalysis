import flambeau / [flambeau_nn, tensors]
import std / [strformat, os, strutils, sequtils, random, algorithm, options]

import ingrid / [tos_helpers, ingrid_types]
import pkg / [datamancer, unchained, nimhdf5]

# have to include the type definitions
include ./nn_types
import ./io_helpers
include ./nn_cuts

{.experimental: "views".}

var ClampOutput = 50.0

## XXX: a bit annoying that this is here...
const CdlFile = "/home/basti/CastData/data/CDL_2019/calibration-cdl-2018.h5"

proc generateTrainTest(df: var DataFrame):
                      ((RawTensor, RawTensor), (RawTensor, RawTensor)) {.noInit.} =
  # - generate
  # - generate random numbers from 0 to df length
  # - add as column and sort by it

  ## XXX: this mixes signal and background events!
  ## XXX: Replace this `"x" in df` check by something sane!
  if "x" in df: # raw data
    const Num = 1000
    let nTrain = Num #df.len div 2
    let nTest = Num # df.len - nTrain
    var data = rawtensors.zeros(nTrain * 256 * 256).reshape([nTrain.int64, 256, 256].asTorchView())
    var target = rawtensors.zeros(nTrain * 2).reshape(sizes = [nTrain.int64, 2].asTorchView())
    var dataTest = rawtensors.zeros(nTest * 256 * 256).reshape(sizes = [nTest.int64, 256, 256].asTorchView())
    var targetTest = rawtensors.zeros(nTest * 2).reshape(sizes = [nTest.int64, 2].asTorchView())
    var i = 0
    for (tup, subDf) in groups(df.group_by(["eventNumber", "Type"])):
      let ev = tup[0][1].toInt
      # echo "Event number ", ev, " at index ", i
      let xs = subDf["x", int]
      let ys = subDf["y", int]
      let isSignal = tup[1][1].toStr == $dtSignal
      for j in 0 ..< xs.size:
        let x = xs[j]
        let y = ys[j]
        if i < nTrain:
          data[i, x, y] = 1.0
          target[i, _] = if isSignal: [1, 0].toRawTensorFromScalar
                         else: [0, 1].toRawTensorFromScalar
        else:
          dataTest[i - nTrain, x, y] = 1.0
          targetTest[i - nTrain, _] = if isSignal: [1, 0].toRawTensorFromScalar
                                      else: [0, 1].toRawTensorFromScalar
      inc i
      if i >= nTrain + nTest: break
    result = (train: (data, target), test: (dataTest, targetTest))
  else:
    df["Idx"] = toSeq(0 .. df.high).shuff()
    df = df.arrange("Idx")

    let dfTrain = df[0 .. df.high div 2]
    let dfTest = df[df.high div 2 + 1 .. df.high]

    ## Add another output each to get the target
    ## thus return
    ## train: (input, target)
    ## test: (input, target)
    ## tuples each.
    echo "SIZES\n\n"
    echo dfTrain.len
    echo dfTest.len
    result = (train: dfTrain.toInputTensor, test: dfTest.toInputTensor)

proc plotLikelihoodDist(df: DataFrame) =
  echo df
  let df = df.mutate(f{float -> float: "likelihood" ~ (if classify(idx("likelihood")) == fcInf: 50.0
                                                       else: idx("likelihood"))})
  ggplot(df, aes("likelihood", fill = "Type")) +
    geom_histogram(bins = 100, position = "identity", alpha = some(0.5), hdKind = hdOutline) +
    scale_x_continuous() +
    ggsave("/tmp/likelihood.pdf")

proc plotTraining(predictions: seq[float], targets: seq[int],
                  outfile = "/tmp/test.pdf") =
  #echo "INPUT ", predictions.len
  #echo "TARG ", targets.len
  let dfPlt = toDf(predictions, targets)
    .mutate(f{"isSignal" ~ `targets` == 1})
    .filter(f{`predictions` > -ClampOutput and `predictions` < ClampOutput})
  #dfPlt.showBrowser()
  #echo "Number of signals: ", dfPlt.filter(f{`isSignal` == true})
  #echo "Number of backs: ", dfPlt.filter(f{`isSignal` == false})
  #if true: quit()
  try:
    ggplot(dfPlt, aes("predictions", fill = "isSignal")) +
      geom_histogram(bins = 100, position = "identity", alpha = some(0.5), hdKind = hdOutline) +
      scale_x_continuous() +
      ggsave(outfile)
  except:
    discard

proc plotType(epoch: int, data, testData: seq[float], typ: string, outfile: string) =
  if data.len < 2: return
  var df = toDf(data, testData)
  df["Epoch"] = linspace(0, epoch, data.len)
  df = df
    .gather(["data", "testData"], "Type", typ)
  try:
    ggplot(df, aes("Epoch", typ, color = "Type")) +
      geom_line() +
      scale_y_log10() +
      ggsave(outfile)
  except:
    discard

proc rocCurve(df: DataFrame,
              suffix = "", plotPath = "/tmp",
              yLow = 0.96) =
  ## plots the ROC curve of the predictions vs the targets
  var plt = if "Type" in df:
              ggplot(df, aes("sigEff", "backRej", color = "Type"))
            else:
              ggplot(df, aes("sigEff", "backRej"))
  try:
    plt +
      geom_line() +
      ylim(yLow, 1.0) +
      ggsave(&"{plotPath}/roc_curve{suffix}.pdf")
  except:
    discard

proc rocCurve(predictions: seq[float], targets: seq[int],
              suffix = "", plotPath = "/tmp") =
  let dfRoc = calcRocCurve(predictions, targets)
  rocCurve(dfRoc, suffix, plotPath)

proc logLValues(df: DataFrame): (seq[float], seq[int]) =
  let logl = df["likelihood", float].map_inline:
    if classify(x) == fcInf:
      50.0
    else: x
  let targets = df["Type", string].map_inline:
    if x == $dtBack: 0
    else: 1
  result = (logL.toSeq1D, targets.toSeq1D)

proc plotRocCurve(dfLnL, dfMLP: DataFrame, suffix = "_likelihood", plotPath = "/tmp") =
  ## plots the ROC curve of the predictions vs the targets
  let (logl, targets) = logLValues(dfLnL)
  let (preds, targetsMlp) = logLValues(dfMlp) ## Note: dfMlp must have its dataset named `likelihood`! # .rename(f{"likelihood" <- "preds"}))
  let dfRoc = calcRocCurve(logl, targets, preds, targetsMlp)
  rocCurve(dfRoc, suffix, plotPath)

proc prepareDataframe(fname: string, run: int, readRaw: bool): DataFrame =
  let dfCdl = prepareCdl(readRaw)
  let dfBack = prepareBackground(fname, run, readRaw) # .drop(["centerX", "centerY"])
  echo dfCdl
  echo dfBack
  result = newDataFrame()
  result.add dfCdl
  result.add dfBack

proc prepareMixedDataframe(calib, back: seq[string], readRaw: bool,
                           subsetPerRun: int): DataFrame =
  let dfCdl = prepareCdl(readRaw)
  var dfFe = newDataFrame()
  for c in calib:
    dfFe.add readCalibData(c, "escape", 2.5, 3.5, subsetPerRun div 2)
    dfFe.add readCalibData(c, "photo", 5.0, 7.0, subsetPerRun div 2)
  var dfBack = newDataFrame()
  for b in back:
    dfBack.add prepareAllBackground(b, readRaw, subsetPerRun = subsetPerRun * 6) # .drop(["centerX", "centerY"])
  echo "CDL: ", dfCdl
  echo "55Fe: ", dfFe
  echo "Back: ", dfBack
  result = newDataFrame()
  result.add dfCdl
  result.add dfFe.drop(["runNumber", "CalibType"])
  result.add dfBack.drop(["runNumber"])

template printTensorInfo(arg: untyped): untyped =
  echo astToStr(arg), ": ", typeof(arg), " on device ", arg.get_device(), " is cuda ", arg.is_cuda()

proc genCheckpointName(s: string, epoch: int, loss, acc: float): string =
  ## Generate a name for a model checkpoint given a base name (including file extension)
  ## and the epoch and test loss & accuracy
  result = s
  result.removeSuffix(".pt")
  result.add &"checkpoint_epoch_{epoch}_loss_{loss:.4f}_acc_{acc:.4f}.pt"

proc train(model: AnyModel, optimizer: var Optimizer,
           input, target: RawTensor,
           testInput, testTarget: RawTensor,
           device: Device,
           readRaw: bool,
           modelOutpath: string,
           plotPath = "/tmp/") =
  let dataset_size = input.size(0)
  var toPlot = false
  var accuracies = newSeq[float]()
  var testAccuracies = newSeq[float]()
  var losses = newSeq[float]()
  var testLosses = newSeq[float]()

  ## XXX: Implement snapshotting of the model if better loss in test than before!

  const PlotEvery = 5000
  for epoch in 0 .. 100000:
    var correct = 0
    if epoch mod 50 == 0:
      echo "Epoch is:" & $epoch
    if epoch mod PlotEvery == 0:
      toPlot = true
    var predictions = newSeqOfCap[float](dataset_size)
    var targets = newSeqOfCap[int](dataset_size)
    ## XXX: Adjust the upper end similar to in `test` to get all data!
    var sumLoss = 0.0
    var count = 0
    for batch_id in 0 ..< (dataset_size.float / bsz.float).ceil.int:
      # Reset gradients.
      optimizer.zero_grad()

      # minibatch offset in the Tensor
      ## TODO: generalize the batching and make sure to take `all` elements! (currently skip last X)
      # minibatch offset in the Tensor
      let offset = batch_id * bsz
      let stop = min(offset + bsz, dataset_size)
      var x: RawTensor
      if not readRaw:
        x = input[offset ..< stop, _ ].to(device)
      else:
        x = input[offset ..< stop, _, _ ].unsqueeze(1).to(device)
        #x = input[offset ..< offset + bsz, _, _ ].unsqueeze(1)
        echo "INPUT TENSOR SHAPE ", x.sizes

      let target = target[offset ..< stop].to(device)
      #printTensorInfo(target)
      # Running input through the network
      let output = model.forward(x)
      #echo "computed forward"
      #printTensorInfo(output)
      let pred = output.argmax(1)
      #echo "pred ", pred
      #printTensorInfo(pred)

      if toPlot:
        # take 0th column
        predictions.add output[_, 0].toNimSeq[:float]
        targets.add target[_, 0].toNimSeq[:int]
        correct += pred.eq(target.argmax(1)).sum().item(int)
      # Computing the loss
      #var loss = l1_loss(output, target)
      var loss = sigmoid_cross_entropy(output, target)
      # Compute the gradient (i.e. contribution of each parameter to the loss)
      loss.backward()
      #echo "did backward"
      # Correct the weights now that we have the gradient information
      optimizer.step()
      #echo "stepped"
      sumLoss += loss.item(float)
      inc count

    if toPlot:
      let train_acc = correct.float / dataset_size.float64() # loss.item(float)
      let averageLoss = sumLoss / count.float
      echo &"\nTrain set: Average loss: {averageLoss:.4f} " &
           &"| Accuracy: {correct.float64() / dataset_size.float64():.3f}"
      let (testAccuracy, testLoss, tOut, tTarget) = model.test(testInput, testTarget, device)
      # Accuracies
      accuracies.add train_acc
      testAccuracies.add testAccuracy
      # Losses
      losses.add(min(1.0, averageLoss))
      testLosses.add(min(1.0, testLoss))

      # now save the model as a checkpoint
      model.save(genCheckpointName(modelOutpath, epoch, testLoss, testAccuracy))

    ## create output plot
    if toPlot:
      plotTraining(predictions, targets, outfile = &"{plotPath}/test_training.pdf")
      plotType(epoch, accuracies, testAccuracies, "Accuracy", outfile = &"{plotPath}/accuracy.pdf")
      plotType(epoch, losses, testLosses, "Loss", outfile = &"{plotPath}/loss.pdf")
      let preds = predictions.mapIt(clamp(it, -ClampOutput, ClampOutput))
      rocCurve(preds, targets, plotPath = plotPath)
      toPlot = false

proc test(model: AnyModel,
          input, target: RawTensor,
          device: Device,
          plotOutfile = "/tmp/test_validation.pdf",
          neuron = 0): (float, float, seq[float], seq[int]) =
  ## returns the predictions / targets
  let dataset_size = input.size(0)
  var correct = 0
  var predictions = newSeqOfCap[float](dataset_size)
  var targets = newSeqOfCap[int](dataset_size)
  var sumLoss = 0.0
  var count = 0
  no_grad_mode:
    for batch_id in 0 ..< (dataset_size.float / bsz.float).ceil.int:
      # minibatch offset in the Tensor
      let offset = batch_id * bsz
      let stop = min(offset + bsz, dataset_size)
      let x = input[offset ..< stop, _ ].to(device)
      let target = target[offset ..< stop].to(device)
      # Running input through the network
      let output = model.forward(x)
      # get the larger prediction along axis 1 (the example axis)
      let pred = output.argmax(1)
      # take 0th column
      predictions.add output[_, neuron].toNimSeq[:float]
      targets.add target[_, 0].toNimSeq[:int]
      correct += pred.eq(target.argmax(1)).sum().item(int)
      # Computing the loss
      #let loss = l1_loss(output, target)
      let loss = sigmoid_cross_entropy(output, target)
      sumLoss += loss.item(float)
      inc count

  let testAcc = correct.float / dataset_size.float64()
  let testLoss = sumLoss / count.float
  echo &"\nTest set: Average loss: {testLoss:.4f} " &
       &"| Accuracy: {testAcc:.4f}"

  ## create output plot
  plotTraining(predictions, targets, outfile = plotOutfile)
  # will fail probably...
  # doAssert target == targets
  let preds = predictions.mapIt(clamp(it, -ClampOutput, ClampOutput))
  result = (testAcc, testLoss, preds, targets)

proc predict(model: AnyModel,
             input, target: RawTensor,
             device: Device,
             cutVal: float,
             neuron: int = 0,
             plotPath = "/tmp/"): seq[int] =
  ## returns the predictions / targets
  let dataset_size = input.size(0)
  var correct = 0
  var predictions = newSeqOfCap[float](dataset_size)
  var targets = newSeqOfCap[int](dataset_size)
  var sumLoss = 0.0
  var count = 0
  no_grad_mode:
    for batch_id in 0 ..< (dataset_size.float / bsz.float).ceil.int:
      # minibatch offset in the Tensor
      let offset = batch_id * bsz
      let stop = min(offset + bsz, dataset_size)
      let x = input[offset ..< stop, _ ].to(device)
      let target = target[offset ..< stop].to(device)
      # Running input through the network
      let output = model.forward(x)
      # get the larger prediction along axis 1 (the example axis)
      let pred = output.argmax(1)
      let predSeq = pred.toNimSeq[:int]
      # take 0th column
      predictions.add output[_, neuron].toNimSeq[:float]
      targets.add target[_, 0].toNimSeq[:int]
      correct += pred.eq(target.argmax(1)).sum().item(int)
      # Computing the loss
      let loss = sigmoid_cross_entropy(output, target)
      sumLoss += loss.item(float)
      inc count
      for i in 0 ..< min(bsz, stop - offset):
        if predSeq[i] == 0:
          result.add (offset + i).int # else add the index of this event that looks like signal

  let predAcc = correct.float / dataset_size.float64()
  let predLoss = sumLoss / count.float
  echo &"\nPred set: Average loss: {predLoss:.4f} " &
       &"| Accuracy: {predAcc:.4f}"

  ## create output plot
  #plotTraining(predictions, targets, outfile = plotOutfile)
  # will fail probably...
  # doAssert target == targets
  #let preds = predictions.mapIt(clamp(it, -50.0, 50.0))
  #result = preds
  when false:
    let dataset_size = input.size(0)
    var correct = 0
    var predictions = newSeq[float]()
    no_grad_mode:
      ## XXX: Adjust the upper end similar to in `test` to get all data!
      for batch_id in 0 ..< (dataset_size.float / bsz.float).ceil.int:
        # minibatch offset in the Tensor
        let offset = batch_id * bsz
        let stop = min(offset + bsz, dataset_size)
        let x = input[offset ..< stop, _ ].to(device)
        let target = target[offset ..< stop].to(device)
        # Running input through the network
        let output = model.forward(x)
        # get the larger prediction along axis 1 (the example axis)
        let pred = output.argmax(1).toNimSeq[:float]()
        # take 0th column
        #echo pred # [i], " vs ", output[i, _]
        predictions.add output[_, 0].toNimSeq[:float]
        for i in 0 ..< min(bsz, stop - offset):
          #if target[i, _].argmax() == output[i, _].argmax():
          #  #echo "Is correct! ", offset + i # in this case it is indeed classified as background
          #  discard
          #if pred[i].item(int) == 0: # target is [1, 0] for signal. i.e. argmax for signal should be `0` for signal
          #if output[i, _].argmax(0).item(int) == 0: # target is [1, 0] for signal. i.e. argmax for signal should be `0` for signal
          #  result.add offset + i
          #else:
          if output[i, 0].item(float) > cutVal:
            #echo "is false! ", offset + i
            result.add (offset + i).int # else add the index of this event that looks like signal
        correct += pred.eq(target.argmax(1)).sum().item(int)
    let test_loss = correct.float / dataset_size.float64()
    echo &"\nPredict set: Average loss: {test_loss:.4f} " &
         &"| Accuracy: {correct.float64() / dataset_size.float64():.3f}"

    let df = toDf(predictions)
      .filter(f{`predictions` > -ClampOutput})
    ggplot(df, aes("predictions")) +
      geom_histogram(bins = 100) +
      ggsave(&"{plotPath}/predictions_all_background.pdf")

## XXX: remove the plotting, scaling etc logic and rather output the
## prediction data correctly. or rather don't do any prediction in this tool
## per se, move that to `likelihood` and use it as a regular `fkMlp` / `fkCNN` method!
proc scaleDset(data: Column, totalTime, factor: float): Column =
  ## scales the data in `data` according to the area of the gold region,
  ## total time and bin width. The data has to be pre binned of course!
  let area = pow(0.95 - 0.45, 2) # area of gold region!
  const bin_width = 0.2 # 0.392
  const shutter_open = 1.0 # Does not play any role, because totalDuration takes
                           # this into account already!
  let scale = factor / (totalTime * shutter_open * area * bin_width) #* 1e5
  result = toColumn data.toTensor(float).map_inline(x * scale)

proc histogram(df: DataFrame): DataFrame =
  ## Calculates the histogam of the energy data in the `df` and returns
  ## a histogram of the binned data
  ## TODO: allow to do this by combining different `File` values
  let (hist, bins) = histogram(df["energies"].toTensor(float).toSeq1D,
                               range = (0.0, 20.0), bins = 100)
  result = toDf({ "energies" : bins, "count" : concat(hist, @[0]) })

proc plotBackground(data: seq[float], totalTime: Hour, plotPath = "/tmp/") =
  let dfE = toDf({"energies" : data})
    .filter(f{`energies` < 20.0})
  #dfE.showBrowser()
  var dfH = histogram(dfE)
  let t = if totalTime > 0.0.Hour: totalTime else: 2144.67.Hour ## Correct time for Run-2
  dfH["Rate"] = dfH["count"].scaleDset(t.to(Second).float, 1e5)
  echo dfH
  #ggplot(dfE, aes("energies")) +
  #  geom_histogram(bins = 100) +
  #  ggsave("/tmp/simple_background_rate.pdf")
  ggplot(dfH, aes("energies", "Rate")) +
    geom_histogram(stat = "identity", hdKind = hdOutline) +
    xlab("Energy [keV]") + ylab("Rate [10⁻⁵ keV⁻¹•cm⁻²•s⁻¹]") +
    ggsave(&"{plotPath}/simple_background_rate.pdf")

proc targetSpecificRoc(dfLnL, dfMlp: DataFrame, plotPath = "/tmp") =
  ## computes the CDL target specific ROC curves for logL and MLP

  # 1. split the input df by target and perform regular ROC calculations
  let dfRoc = bind_rows([("LogL", dfLnL), ("MLP", dfMlp)],
                        "Method")
  var df = newDataFrame()
  for tup, subDf in groups(dfRoc.group_by(["Method", "Target"])):
    let meth = tup[0][1].toStr # lnL or MLP
    let target = tup[1][1].toStr
    let (output, targets) = subDf.logLValues()
    let verbose = meth == "LogL" and target == "C-EPIC-0.6kV"
    var dfRoc = calcRocCurve(output, targets, verbose)
    dfRoc["Target"] = target
    dfRoc["Method"] = meth
    df.add dfRoc
  ggplot(df, aes("sigEff", "backRej", color = "Target", shape = "Method")) +
    geom_line() +
    ylim(0.5, 1.0) +
    ggsave(&"{plotPath}/roc_curve_combined_split_by_target.pdf")

proc determineCdlEfficiency(model: AnyModel, device: Device, ε: float, readRaw: bool,
                            global = true, plotPath = "/tmp/") =
  ##
  let dfCdl = prepareCdl(readRaw)
  let cutValGlobal = determineCutValue(model, device, ε, readRaw)
  let cutValLocal = determineLocalCutValue(model, device, ε, readRaw)
  # for reference determine cut value based on full data (expect to get out `ε`)
  proc getEff(model: AnyModel, df: DataFrame, cutVal: float, outfile: string): float =
    let (cdlInput, cdlTarget) = df.toInputTensor()
    let (_, _, cdlPredict, predTargets) = model.test(cdlInput, cdlTarget, device,
                                                     plotOutfile = outfile)
    result = determineEff(cdlPredict.sorted, cutVal)
  let totalEff = model.getEff(dfCdl, cutValGlobal, &"{plotPath}/cdl_prediction_all_data.pdf")
  echo "Total global efficiency = ", totalEff
  for (tup, subDf) in groups(dfCdl.group_by("Target")):
    let target = tup[0][1].toStr
    #echo "Current target: ", target, " data: ", subDf
    let cut = if global: cutValGlobal
              else: cutValLocal[target]
    let effectiveEff = model.getEff(subDf, cut, &"{plotPath}/cdl_prediction_{target}.pdf")
    let suffix = if global: "global" else: "local"
    echo "Target ", suffix, " ", target, " eff = ", effectiveEff

proc predictBackground(model: AnyModel, device: Device, fname: string, ε: float, totalTime: Hour,
                       readRaw: bool, plotPath: string) =
  # classify
  let cutVal = determineCutValue(model, device, ε, readRaw)
  ## XXX: make sure the cut value stuff & the sign of the predictions is correct everywhere! We flipped the
  ## sign initially to compare with logL!
  let dfAll = prepareAllBackground(fname, readRaw)
  let (allInput, allTarget) = dfAll.toInputTensor()
  let passedInds = model.predict(allInput, allTarget, device, cutVal)
  echo "p inds len ", passedInds.len, " compared to all ", dfAll.len

  # get all energies of these passing events
  let energiesAll = dfAll["energyFromCharge", float]
  var energies = newSeq[float](passedInds.len)
  for i, idx in passedInds:
    energies[i] = energiesAll[idx]
  plotBackground(energies, totalTime, plotPath)

proc predictCalib(model: AnyModel, device: Device, fname: string, ε: float) =
  # classify
  var df = newDataFrame()
  df.add readCalibData(fname, "escape", 2.5, 3.5)
  df.add readCalibData(fname, "photo", 5.5, 6.5)

  #let (model, device) = loadModelMakeDevice(model)
  for (tup, subDf) in groups(df.group_by(["Type", "runNumber"])):
    let (inp, target) = subDf.toInputTensor()
    let passedInds = model.predict(inp, target, device, 0.0)
    echo "p inds len ", passedInds.len, " compared to all ", subDf.len, " Efficiency: ", passedInds.len.float / subDf.len.float

proc readAllData(model: AnyModel, device: Device,
                 calib, back: seq[string],
                 subsetPerRun: int,
                 calcLogL: bool,
                 mlpOutput: bool): (DataFrame, DataFrame, DataFrame) =
  var dfCdl = prepareCdl(false)
    #.filter(f{`Target` == "Mn-Cr-12kV"}) #<- to only look at Photo peak equivalent. Better look at ridge line!
  var dfFe = newDataFrame()
  for c in calib:
    dfFe.add readCalibData(c, "escape", 2.5, 3.5) #, subsetPerRun div 2)
    dfFe.add readCalibData(c, "photo", 5.0, 7.0) #, subsetPerRun div 2)
  var dfBack = newDataFrame()
  for b in back:
    dfBack.add prepareAllBackground(b, false) # , subsetPerRun = subsetPerRun * 6) # .drop(["centerX", "centerY"])

  if calcLogL: # also calculate the likelihood of every event
    ## Set up the likelihood context to compute lnL values!
    let ctx = initLikelihoodContext(CdlFile,
                                    year = yr2018,
                                    energyDset = igEnergyFromCharge,
                                    region = crSilver,
                                    timepix = Timepix1,
                                    morphKind = mkLinear) # morphing to plot interpolation

    template lnL(ctx, df: untyped): untyped =
      df[igLikelihood.toDset()] = ctx.calcLikelihood(df).mapIt(clamp(it, -ClampOutput, ClampOutput))
    ctx.lnL(dfCdl)
    ctx.lnL(dfFe)
    ctx.lnL(dfBack)
  if mlpOutput:
    proc modelOutput(model: AnyModel, device: Device, df: var DataFrame, neuron: int) =
      df["Neuron_" & $neuron] = model.forward(df.toTorchTensor(), device, neuron)
        .mapIt(clamp(it, -ClampOutput, ClampOutput))
    model.modelOutput(device, dfCdl, neuron = 0)
    model.modelOutput(device, dfCdl, neuron = 1)
    model.modelOutput(device, dfFe, neuron = 0)
    model.modelOutput(device, dfFe, neuron = 1)
    model.modelOutput(device, dfBack, neuron = 0)
    model.modelOutput(device, dfBack, neuron = 1)

  result = (dfCdl, dfFe, dfBack)

proc predictAll(model: AnyModel, device: Device,
                calib, back: seq[string],
                subsetPerRun: int,
                plotPath: string) =
  let (dfCdl, dfFe, dfBack) = model.readAllData(device, calib, back,
                                                subsetPerRun,
                                                calcLogL = true,
                                                mlpOutput = true)

  # for each get prediction, neuron 0 and neuron 1
  template p(dfIn, neuron, typ: untyped): untyped =
    block:
      var df = toDf({"preds" : dfIn["Neuron_" & $neuron, float]})
      df["Type"] = typ
      df["Neuron"] = neuron
      df

  block AllPredictions:
    var dfAll = newDataFrame()
    dfAll.add p(dfCdl, 0, "CDL")
    dfAll.add p(dfCdl, 1, "CDL")
    dfAll.add p(dfFe, 0, "55Fe")
    dfAll.add p(dfFe, 1, "55Fe")
    dfAll.add p(dfBack, 0, "back")
    dfAll.add p(dfBack, 1, "back")
    ggplot(dfAll, aes("preds", fill = "Type")) +
      facet_wrap("Neuron") +
      geom_histogram(bins = 100, hdKind = hdOutline, alpha = 0.5, position = "identity") +
      ggsave(&"{plotPath}/all_predictions.pdf", width = 1000, height = 500)

  block AllByTypeRidge:
    template lnL(dfIn, typ: untyped): untyped =
      toDf({ "lnLs" : dfIn["likelihood", float],
             "Type" : typ })

    template byIt(df: DataFrame, by: string, useLnL = false): DataFrame =
      var res = newDataFrame()
      for (tup, subDf) in groups(df.group_by(by)):
        let typ = tup[0][1].toStr
        if not useLnL:
          res.add p(subDf, 0, typ)
        else:
          res.add lnL(subDf, typ)
      res

    # predictions for MLP
    var dfMlp = newDataFrame()
    dfMlp.add byIt(dfCdl, "Target")
    dfMlp.add byIt(dfFe, "CalibType")
    dfMlp.add p(dfBack, 0, "back")
    # likelihood values
    var dfLnL = newDataFrame()
    dfLnL.add byIt(dfCdl, "Target", true)
    dfLnL.add byIt(dfFe, "CalibType", true)
    dfLnL.add lnL(dfBack, "back")

    var labelOrder = initTable[Value, int]()
    labelOrder[%~ "back"] = 0
    labelOrder[%~ "escape"] = 1
    labelOrder[%~ "photo"] = 2
    const xrayRef = getXrayRefTable()
    for idx, el in xrayRef:
      labelOrder[%~ el] = idx + 3

    proc plotRidges(df: DataFrame, x, plotPath, suffix: string) =
      ggplot(df, aes(x, fill = "Type")) +
        ggridges("Type", overlap = 1.75, labelOrder = labelOrder) +
        geom_histogram(bins = 300, color = "black", lineWidth = 0.75,
                       hdKind = hdOutline, position = "identity",
                       density = true) +
        ggsave(&"{plotPath}/all_predictions_by_type_ridgeline{suffix}.pdf", width = 800, height = 500)
    plotRidges(dfMlp, "preds", plotPath, "_mlp")
    # remove anything in the clamp region to not blow up our y scale with an 'overfilled' bin
    plotRidges(dfLnL.filter(f{`lnLs` < 49.0}), "lnLs", plotPath, "_lnL")

proc generateRocCurves(model: AnyModel, device: Device, calib, back: seq[string],
                       subsetPerRun: int, plotPath: string) =
  let (dfCdl, dfFe, dfBack) = model.readAllData(device, calib, back,
                                                subsetPerRun,
                                                calcLogL = true,
                                                mlpOutput = true)
  proc asDf(dfIn: DataFrame, typ, dset: string): DataFrame =
    toDf({ "likelihood" : (if dset == "Neuron_0": dfIn[dset, float].map_inline(-x)
                           else: dfIn[dset, float]),
           "Target" : dfIn["Target", string],
           "Type" : typ })
  template buildLong(dfCdl, dfFe, dfBack, dset: untyped): untyped =
    # need the `likelihood` dataset and the data type
    var df = newDataFrame()
    df.add asDf(dfCdl, $dtSignal, dset)
    df.add asDf(dfFe,  $dtSignal, dset)
    df.add asDf(dfBack,  $dtBack, dset)
    df
  let dfLnL = buildLong(dfCdl, dfFe, dfBack, "likelihood")
  let dfMlp = buildLong(dfCdl, dfFe, dfBack, "Neuron_0")
  plotRocCurve(dfLnL, dfMlp, "_lnL_vs_mlp", plotPath)
  # target specific roc curves
  targetSpecificRoc(dfLnL, dfMlp, plotPath)

proc initDesc(datasets: seq[InGridDsetKind], numHidden: int, modelOutpath, plotPath: string): MLPDesc =
  if numHidden == 0:
    raise newException(ValueError, "Please provide a number of neurons on the hidden layer. 0 is invalid.")
  let dsets = if datasets.len == 0: CurrentDsets else: datasets
  CurrentDsets = dsets
  let plotPath = if plotPath.len == 0: "/tmp/" else: plotPath
  result = initMLPDesc(dsets, numHidden, modelOutpath, plotPath)
  # potentially create the output path
  discard existsOrCreateDir(result.plotPath)
  discard existsOrCreateDir(result.path.parentDir)
  result.toH5(result.path.parentDir / MLPDescName)

proc initDesc(model, plotPath: string): MLPDesc =
  result = initMLPDesc(model, plotPath)
  CurrentDsets = result.datasets

#Converter toModule[T: AnyModel](m: T): Module = Module(m)
#
## XXX: inside of a generic proc (what we would normally do) the `parameters` call
## breaks! Nim doesn'- understand that the MLP / ConvNet type can be converted
## to a `Module`!
template trainModel(Typ: typedesc,
                    device: Device,
                    calib: seq[string],
                    back: seq[string],
                    datasets: seq[InGridDsetKind],
                    numHidden: int,
                    modelOutpath = "/tmp/trained_model.pt",
                    subsetPerRun = 1000,
                    plotPath = ""
                   ): untyped {.dirty.} =
  ## If we are training, construct a type appropriate to
  let desc = initDesc(datasets, numHidden, modelOutpath, plotPath)
  echo "MLPDesc: ", desc
  var model = Typ.init(desc)
  model.to(device)
  when Typ is MLP:
    const readRaw = false
  else:
    const readRaw = true
  echo "Reading data"
  #var df = prepareDataframe(fname, run, readRaw = readRaw)
  var df = prepareMixedDataframe(calib, back, readRaw = readRaw, subsetPerRun = subsetPerRun)
  # get training & test dataset
  echo "Splitting data into train & test set"
  let (trainTup, testTup) = generateTrainTest(df)
  let (trainIn, trainTarg) = trainTup
  let (testIn, testTarg) = testTup
  # echo trainIn.to(kCPU)[0, _, _] # ConvNet test?

  # check if model already exists as trained file
  if not fileExists(modelOutpath):
    #if true:
    #  echo "fle does not exist"
    #  quit()
    # Learning loop
    # Stochastic Gradient Descent
    var optimizer = SGD.init(
      model.parameters(),
      SGDOptions.init(0.005).momentum(0.2) # .weight_decay(0.001)
      #learning_rate = 0.005
    )
    # Adam optimizer
    #var optimizer = Adam.init(
    #  model.parameters(),
    #  AdamOptions.init(0.005)
    #)
    echo "Training model"
    model.train(optimizer,
                trainIn.to(kFloat32).to(device),
                trainTarg.to(kFloat32).to(device),
                testIn.to(kFloat32).to(device),
                testTarg.to(kFloat32).to(device),
                device,
                readRaw,
                desc.path,
                desc.plotPath)
    model.save(modelOutpath)
  else:
    # load the model
    model.load(modelOutpath)
  # perform validation
  let (acc, loss, testPredict, testTargets) = model.test(testIn, testTarg, device,
                                                         plotOutfile = desc.plotPath / "test_validation.pdf")
  echo "Test loss after training: ", loss, " with accuracy ", acc

template predictModel(Typ: typedesc,
                      device: Device,
                      calib: seq[string],
                      back: seq[string],
                      modelPath: string,
                      ε = 0.8, # signal efficiency for background rate prediction
                      totalTime = -1.0.Hour, # total background rate time in hours. Normally read from input file
                      rocCurve = false,
                      subsetPerRun = 1000,
                      plotPath = ""
                     ): untyped {.dirty.} =
  doAssert fileExists(modelPath), "When using the `--predict` option, the trained model must exist!"
  ## If we are training, construct a type appropriate to
  let desc = initDesc(modelPath, plotPath)
  echo "MLPDesc: ", desc
  var model = Typ.init(desc)
  model.to(device)
  when Typ is MLP:
    const readRaw = false
  else:
    const readRaw = true

  # load the model
  model.load(desc.path)
  model.predictBackground(device, back[0], ε, totalTime, readRaw, desc.plotPath)
  model.predictCalib(device, calib[0], ε)

  model.predictAll(device, calib, back, subsetPerRun, desc.plotPath)

  model.determineCdlEfficiency(device, ε, readRaw, plotPath = desc.plotPath)
  model.determineCdlEfficiency(device, ε, readRaw, global = false, plotPath = desc.plotPath)

  model.generateRocCurves(device, calib, back, subsetPerRun, desc.plotPath)

proc writeMLPDesc(model, plotPath: string,
                  numHidden: int,
                  datasets: seq[InGridDsetKind]) =
  ## This helper can be used to write an `MLPDesc` for a given model in case it doesn't have
  ## this file for some reason (too old, deleted, etc.)
  ## You are responsible for making sure the input data is correct of course!
  let desc = initDesc(datasets, numHidden, model, plotPath)
  desc.toH5(model.parentDir / MLPDescName)

proc main(calib, back: seq[string] = @[],
          datasets: seq[InGridDsetKind] = @[],
          ε = 0.8, # signal efficiency for background rate prediction
          totalTime = -1.0.Hour, # total background rate time in hours. Normally read from input file
          rocCurve = false,
          predict = false,
          model = "MLP", # MLP or ConvNet #model = mkMLP ## parsing an enum here causes weird CT error in cligen :/
          modelOutpath = "/tmp/trained_model.pt",
          numHidden = 0,
          subsetPerRun = 1000,
          plotPath = "",
          clampOutput = 50.0,
          printDefaultDatasets = false,
          writeMLPDesc = false) =
  # 1. set up the model
  if printDefaultDatasets:
    echo "Total default datasets: ", CurrentDsets.len
    for d in CurrentDsets:
      echo d
    return
  if writeMLPDesc:
    let datasets = if datasets.len == 0: CurrentDsets else: datasets
    writeMLPDesc(modelOutpath, plotPath, numHidden, datasets)
    return

  ClampOutput = clampOutput

  Torch.manual_seed(1)
  var device_type: DeviceKind
  if Torch.cuda_is_available():
    echo "CUDA available! Training on GPU."
    device_type = kCuda
  else:
    echo "Training on CPU."
    device_type = kCPU
  let device = Device.init(device_type)

  let mKind = parseEnum[ModelKind](model)
  if mKind == mkMLP:
    if not predict:
      MLP.trainModel(device,
                     calib, back,
                     datasets, numHidden,
                     modelOutpath,
                     subsetPerRun,
                     plotPath)
    else:
      MLP.predictModel(device,
                       calib, back,
                       modelOutpath,
                       ε, totalTime,
                       rocCurve,
                       subsetPerRun,
                       plotPath)
  #else:
  #  ConvNet.trainModel(fname, device, run, ε, totalTime, rocCurve, predict)

when isMainModule:
  import cligen/argcvt
  proc argParse(dst: var Hour, dfl: Hour,
                a: var ArgcvtParams): bool =
    var val = a.val
    doAssert val.endsWith(".Hour"), "Please give total time in the form `<value>.Hour`"
    val.removeSuffix(".Hour")
    dst = val.parseFloat.Hour
    result = true
  proc argParse(dst: var seq[InGridDsetKind], dfl: seq[InGridDsetKind],
                a: var ArgcvtParams): bool =
    var val = a.val
    try:
      dst.add parseEnum[InGridDsetKind](val)
      result = true
    except ValueError:
      raise newException(Exception, "Invalid dataset given: " & $val)
  proc argHelp*(dfl: Hour; a: var ArgcvtParams): seq[string] =
    result = @[ a.argKeys, "<value>.Hour", $dfl ]
  import cligen
  dispatch main
