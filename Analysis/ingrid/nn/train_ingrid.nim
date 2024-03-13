import flambeau / [flambeau_nn, tensors]
import std / [strformat, os, strutils, sequtils, random, algorithm, options, macros]

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

let UseTeX = getEnv("USE_TEX", "false").parseBool
let Width = getEnv("WIDTH", "600").parseFloat
let Height = getEnv("Height", "420").parseFloat

proc thL(fWidth: float, width: float,
         baseTheme: (proc(): Theme),
         height = -1.0, ratio = -1.0,
         textWidth = 458.29268, # 455.24411
        ): Theme =
  if UseTeX:
    let texOptions = toTeXOptions(UseTeX, onlyTikZ = false,
                                  standalone = true,
                                  texTemplate = "", caption = "", label = "", placement = "htbp")
    result = themeLatex(fWidth, width, baseTheme, height, ratio, textWidth,
                        useTeX = UseTeX, texOptions = texOptions)
  else:
    result = Theme()

proc generateTrainTest(df: DataFrame, rnd: var Rand):
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
    var df = df # mutable copy
    df["Idx"] = rnd.shuff(toSeq(0 .. df.high))
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

proc plotPredictions(predictions: seq[float], targets: seq[int],
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
      xlab("Predictions") + ylab("Count") +
      scale_x_continuous() +
      margin(right = 3.5, left = 3.5) +
      thL(fWidth = 0.5, width = Width, baseTheme = sideBySide) +
      ggsave(outfile, width = 600, height = 420)
  except:
    discard

  # now a log10 version
  var dfL = newDataFrame()
  var maxH = 0.0
  for (tup, subDf) in groups(dfPlt.group_by("isSignal")):
    let (hist, bins) = histogram(subDf["predictions", float].toSeq1D, bins = 100)
    let dfH = toDf({"count" : hist.extend(0), "predictions" : bins, "isSignal" : tup[0][1].toBool})
    dfL.add dfH
    maxH = max(maxH, hist.max.float)
  try:
    ggplot(dfL, aes("predictions", "count", fill = "isSignal")) +
      geom_histogram(stat = "identity", position = "identity", alpha = some(0.5), hdKind = hdOutline) +
      scale_x_continuous() +
      xlab("Predictions") + ylab("Count") +
      scale_y_log10() +
      ylim(0.0, log10(maxH)) +
      margin(right = 3.5, left = 3.5) +
      thL(fWidth = 0.5, width = Width, baseTheme = sideBySide) +
      ggsave(outfile.replace(".pdf", "_log10.pdf"), width = 600, height = 420)
  except:
    discard

proc plotType(epochs: seq[int], data, testData: seq[float], typ: string,
              outfile: string) =
  ## XXX: this is broken for `continueAfterEpochs`!
  if data.len < 2: return
  let df = toDf({"Epoch" : epochs, data, testData})
    .gather(["data", "testData"], "Type", typ)
  try:
    ggplot(df, aes("Epoch", typ, color = "Type")) +
      geom_line() +
      scale_y_log10() +
      scale_x_continuous(breaks = 6) +
      margin(right = 5.0, left = 3.0) +
      thL(fWidth = 0.5, width = Width, baseTheme = sideBySide) +
      ggsave(outfile, width = 600, height = 420)
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

proc predictions(df: DataFrame): (seq[float], seq[int]) =
  ## The `Prediction` column contains either `"likelihood"` values or
  ## `"Neuron_0"` values. As likelihood may be saturated, we clamp to 50.
  let pred = df["Prediction", float].map_inline:
    if classify(x) == fcInf:
      50.0
    else: x
  let targets = df["Type", string].map_inline:
    if x == $dtBack: 0
    else: 1
  result = (pred.toSeq1D, targets.toSeq1D)

proc plotRocCurve(dfLnL, dfMLP: DataFrame, suffix = "_likelihood", plotPath = "/tmp") =
  ## plots the ROC curve of the predictions vs the targets
  let (logL, targets) = predictions(dfLnL)
  let (preds, targetsMlp) = predictions(dfMlp) ## Note: dfMlp must have its dataset named `likelihood`! # .rename(f{"likelihood" <- "preds"}))
  let dfRoc = calcRocCurve(logL, targets, preds, targetsMlp)
  rocCurve(dfRoc, suffix, plotPath)

import ./simulate_xrays
import ingridDatabase / databaseRead
proc prepareSimDataframe(mlpDesc: MLPDesc, readRaw: bool, rnd: var Rand): DataFrame =
  ## Generates a dataframe to train on that contains real background events and
  ## simulated X-rays. The calibration input files are only used to determine the
  ## boundaries of the gas gain and diffusion in which to generate events in!
  # for now we just define hardcoded bounds
  var dfSim = newDataFrame()
  if mlpDesc.simFiles.len == 0:
    let gains = @[2400.0, 4000.0]
    let diffusion = @[550.0, 700.0] # μm/√cm, will be converted to `mm/√cm` when converted to DF
    # the theta parameter describing the Pólya also needs to be varied somehow
    var dfSim = newDataFrame()
    for c in mlpDesc.calibFiles:
      withH5(c, "r"):
        let calibInfo = h5f.initCalibInfo()
        dfSim.add rnd.generateFakeEventsDf(calibInfo, gains, diffusion, mlpDesc.nFake)
  else:
    dfSim.add readSimData(mlpDesc.simFiles)

  var dfBack = newDataFrame()
  for b in mlpDesc.backFiles:
    ## NOTE: For all training runs before 2023/10/30 we multiplied `subsetPerRun` by 6.
    ## As we never adjusted `subsetPerRun` by hand, to reproduce old trainings, use
    ## `--subsetPerRun 6000`!
    dfBack.add prepareAllBackground(b, readRaw, subsetPerRun = mlpDesc.subsetPerRun,
                                    region = mlpDesc.backgroundRegion,
                                    chips = mlpDesc.backgroundChips)
  echo "Simulated: ", dfSim
  echo "Back: ", dfBack
  result = newDataFrame()
  result.add dfSim.drop(["runNumber", "eventNumber", "Target", "likelihood"])
  result.add dfBack.drop(["runNumber", "eventNumber", "Target"])

proc prepareMixedDataframe(mlpDesc: MLPDesc, readRaw: bool): DataFrame =
  let dfCdl = prepareCdl(readRaw)
  var dfFe = newDataFrame()
  let spr = mlpDesc.subsetPerRun
  for c in mlpDesc.calibFiles:
    dfFe.add readCalibData(c, "escape", 2.5, 3.5, spr div 2)
    dfFe.add readCalibData(c, "photo", 5.0, 7.0, spr div 2)
  var dfBack = newDataFrame()
  for b in mlpDesc.backFiles:
    dfBack.add prepareAllBackground(b, readRaw, subsetPerRun = spr, chips = mlpDesc.backgroundChips) # .drop(["centerX", "centerY"])
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

proc test(model: AnyModel,
          input, target: RawTensor,
          device: Device,
          desc: MLPDesc,
          plotOutfile = "/tmp/test_validation.pdf",
          neuron = 0,
          toPlot = true): (float, float, seq[float], seq[int])

proc calcLoss(output, target: RawTensor, desc: MLPDesc): RawTensor =
  ## Applies the correct loss function based on the MLPDesc
  case desc.lossFunction
  of lfL1loss              : result = l1_loss(output, target)
  of lfMSEloss             : result = mse_loss(output, target)
  of lfSigmoidCrossEntropy : result = sigmoid_cross_entropy(output, target)

proc train(model: AnyModel, optimizer: var Optimizer,
           input, target: RawTensor,
           testInput, testTarget: RawTensor,
           device: Device,
           readRaw: bool,
           desc: MLPDesc,
           epochs: int,
           continueAfterEpoch = 0) =
  let dataset_size = input.size(0)
  var toPlot = false

  var mlpDesc = desc # local mutable copy to store losses, accuracy etc in
  let plotPath = mlpDesc.plotPath

  let start = continueAfterEpoch
  let stop = start + epochs
  for epoch in start .. stop:
    var correct = 0
    if epoch mod 50 == 0:
      echo "Epoch is:" & $epoch
    if epoch mod desc.plotEvery == 0:
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
      let output = model.forward(desc, x)
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
      var loss = calcLoss(output, target, desc)
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
      let (testAccuracy, testLoss, tOut, tTarget) = model.test(testInput, testTarget, device, desc, toPlot = false)
      # epoch
      mlpDesc.epochs.add epoch
      # Accuracies
      mlpDesc.accuracies.add train_acc
      mlpDesc.testAccuracies.add testAccuracy
      # Losses
      mlpDesc.losses.add(min(1.0, averageLoss))
      mlpDesc.testLosses.add(min(1.0, testLoss))

      # now save the model as a checkpoint
      model.save(genCheckpointName(mlpDesc.path, epoch, testLoss, testAccuracy))
      # now write the mlpDesc again
      mlpDesc.toH5(mlpDesc.path.parentDir / MLPDescName)

      ## generate the plots
      plotPredictions(predictions, targets, outfile = &"{plotPath}/training_output.pdf")
      plotPredictions(tOut, tTarget, outfile = &"{plotPath}/validation_output.pdf")
      plotType(mlpDesc.epochs, mlpDesc.accuracies, mlpDesc.testAccuracies, "Accuracy", outfile = &"{plotPath}/accuracy.pdf")
      plotType(mlpDesc.epochs, mlpDesc.losses, mlpDesc.testLosses, "Loss", outfile = &"{plotPath}/loss.pdf")
      let preds = predictions.mapIt(clamp(it, -ClampOutput, ClampOutput))
      rocCurve(preds, targets, plotPath = plotPath)
      toPlot = false

proc test(model: AnyModel,
          input, target: RawTensor,
          device: Device,
          desc: MLPDesc,
          plotOutfile = "/tmp/test_validation.pdf",
          neuron = 0,
          toPlot = true): (float, float, seq[float], seq[int]) =
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
      let output = model.forward(desc, x)
      # get the larger prediction along axis 1 (the example axis)
      let pred = output.argmax(1)
      # take 0th column
      predictions.add output[_, neuron].toNimSeq[:float]
      targets.add target[_, 0].toNimSeq[:int]
      correct += pred.eq(target.argmax(1)).sum().item(int)
      # Computing the loss
      var loss = calcLoss(output, target, desc)
      sumLoss += loss.item(float)
      inc count

  let testAcc = correct.float / dataset_size.float64()
  let testLoss = sumLoss / count.float
  echo &"\nTest set: Average loss: {testLoss:.4f} " &
       &"| Accuracy: {testAcc:.4f}"

  ## create output plot
  if toPlot:
    plotPredictions(predictions, targets, outfile = plotOutfile)
  # will fail probably...
  # doAssert target == targets
  let preds = predictions.mapIt(clamp(it, -ClampOutput, ClampOutput))
  result = (testAcc, testLoss, preds, targets)

proc predict(model: AnyModel,
             input, target: RawTensor,
             device: Device,
             desc: MLPDesc,
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
      let output = model.forward(desc, x)
      # get the larger prediction along axis 1 (the example axis)
      let pred = output.argmax(1)
      let predSeq = pred.toNimSeq[:int]
      # take 0th column
      predictions.add output[_, neuron].toNimSeq[:float]
      targets.add target[_, 0].toNimSeq[:int]
      correct += pred.eq(target.argmax(1)).sum().item(int)
      # Computing the loss
      var loss = calcLoss(output, target, desc)
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
  #plotPredictions(predictions, targets, outfile = plotOutfile)
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
        let output = model.forward(desc, x)
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
    let (output, targets) = subDf.predictions()
    let verbose = meth == "LogL" and target == "C-EPIC-0.6kV"
    var dfRoc = calcRocCurve(output, targets, verbose)
    dfRoc["Target"] = target
    dfRoc["Method"] = meth
    df.add dfRoc
  ggplot(df, aes("sigEff", "backRej", color = "Target", shape = "Method")) +
    geom_line() +
    ylim(0.5, 1.0) +
    xlab("Signal efficiency") + ylab("Background rejection") +
    thL(fWidth = 0.9, width = Width, baseTheme = singlePlot) +
    ggsave(&"{plotPath}/roc_curve_combined_split_by_target.pdf")

proc determineCdlEfficiency(model: AnyModel, device: Device, desc: MLPDesc, ε: float, readRaw: bool,
                            global = true, plotPath = "/tmp/") =
  ##
  let dfCdl = prepareCdl(readRaw)
  let cutValGlobal = determineCutValue(model, device, desc, ε, readRaw)
  let cutValLocal = determineLocalCutValue(model, device, desc, ε, readRaw)
  # for reference determine cut value based on full data (expect to get out `ε`)
  proc getEff(model: AnyModel, df: DataFrame, cutVal: float, outfile: string): float =
    let (cdlInput, cdlTarget) = df.toInputTensor()
    let (_, _, cdlPredict, predTargets) = model.test(cdlInput, cdlTarget, device,
                                                     desc,
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
    echo "Target ", suffix, " ", target, " cutValue = ", cut, " eff = ", effectiveEff

proc predictBackground(model: AnyModel, device: Device, desc: MLPDesc, fname: string, ε: float, totalTime: Hour,
                       readRaw: bool, plotPath: string) =
  # classify
  let cutVal = determineCutValue(model, device, desc, ε, readRaw)
  ## XXX: make sure the cut value stuff & the sign of the predictions is correct everywhere! We flipped the
  ## sign initially to compare with logL!
  let dfAll = prepareAllBackground(fname, readRaw)
  let (allInput, allTarget) = dfAll.toInputTensor()
  let passedInds = model.predict(allInput, allTarget, device, desc, cutVal)
  echo "p inds len ", passedInds.len, " compared to all ", dfAll.len

  # get all energies of these passing events
  let energiesAll = dfAll["energyFromCharge", float]
  var energies = newSeq[float](passedInds.len)
  for i, idx in passedInds:
    energies[i] = energiesAll[idx]
  plotBackground(energies, totalTime, plotPath)

proc predictCalib(model: AnyModel, device: Device, desc: MLPDesc, fname: string, ε: float) =
  # classify
  var df = newDataFrame()
  df.add readCalibData(fname, "escape", 2.5, 3.5)
  df.add readCalibData(fname, "photo", 5.5, 6.5)

  #let (model, device) = loadModelMakeDevice(model)
  for (tup, subDf) in groups(df.group_by(["Type", "runNumber"])):
    let (inp, target) = subDf.toInputTensor()
    let passedInds = model.predict(inp, target, device, desc, 0.0)
    echo "p inds len ", passedInds.len, " compared to all ", subDf.len, " Efficiency: ", passedInds.len.float / subDf.len.float

proc readAllData(model: AnyModel, device: Device,
                 desc: MLPDesc,
                 calcLogL: bool,
                 mlpOutput: bool): (DataFrame, DataFrame, DataFrame) =
  var dfCdl = prepareCdl(false)
    #.filter(f{`Target` == "Mn-Cr-12kV"}) #<- to only look at Photo peak equivalent. Better look at ridge line!
  var dfFe = newDataFrame()
  for c in desc.calibFiles:
    dfFe.add readCalibData(c, "escape", 2.5, 3.5) #, subsetPerRun div 2)
    dfFe.add readCalibData(c, "photo", 5.0, 7.0) #, subsetPerRun div 2)
  var dfBack = newDataFrame()
  for b in desc.backFiles:
    dfBack.add prepareAllBackground(b, false) # , subsetPerRun = subsetPerRun * 6) # .drop(["centerX", "centerY"])

  if calcLogL: # also calculate the likelihood of every event
    ## Set up the likelihood context to compute lnL values!
    let ctx = initLikelihoodContext(CdlFile,
                                    year = yr2018,
                                    energyDset = igEnergyFromCharge,
                                    region = crSilver,
                                    timepix = Timepix1,
                                    morphKind = mkLinear,
                                    useLnLCut = true) # morphing to plot interpolation

    template lnL(ctx, df: untyped): untyped =
      df[igLikelihood.toDset()] = ctx.calcLikelihood(df).mapIt(clamp(it, -ClampOutput, ClampOutput))
    ctx.lnL(dfCdl)
    ctx.lnL(dfFe)
    ctx.lnL(dfBack)
  if mlpOutput:
    proc modelOutput(model: AnyModel, device: Device, df: var DataFrame, neuron: int) =
      df["Neuron_" & $neuron] = model.forward(df.toTorchTensor(), device, desc, neuron)
        .mapIt(clamp(it, -ClampOutput, ClampOutput))
    model.modelOutput(device, dfCdl, neuron = 0)
    model.modelOutput(device, dfCdl, neuron = 1)
    model.modelOutput(device, dfFe, neuron = 0)
    model.modelOutput(device, dfFe, neuron = 1)
    model.modelOutput(device, dfBack, neuron = 0)
    model.modelOutput(device, dfBack, neuron = 1)

  result = (dfCdl, dfFe, dfBack)

proc predictAll(model: AnyModel, device: Device,
                desc: MLPDesc) =
  let (dfCdl, dfFe, dfBack) = model.readAllData(device, desc,
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
      ggsave(&"{desc.plotPath}/all_predictions.pdf", width = 1000, height = 500)

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
    plotRidges(dfMlp, "preds", desc.plotPath, "_mlp")
    # remove anything in the clamp region to not blow up our y scale with an 'overfilled' bin
    plotRidges(dfLnL.filter(f{`lnLs` < 49.0}), "lnLs", desc.plotPath, "_lnL")

proc generateRocCurves(model: AnyModel, device: Device, desc: MLPDesc) =
  let (dfCdl, dfFe, dfBack) = model.readAllData(device, desc,
                                                calcLogL = true,
                                                mlpOutput = true)
  proc asDf(dfIn: DataFrame, typ, dset: string): DataFrame =
    toDf({ "Prediction" : (if dset == "Neuron_0": dfIn[dset, float].map_inline(-x)
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
  plotRocCurve(dfLnL, dfMlp, "_lnL_vs_mlp", desc.plotPath)
  # target specific roc curves
  targetSpecificRoc(dfLnL, dfMlp, desc.plotPath)

proc modelDirFile(model, outpath: string): (string, string) =
  if model.len == 0 and outpath.len == 0:
    raise newException(ValueError, "`model` and `modelOutpath` cannot be empty both.")
  elif model.len > 0 and outpath.len == 0: ## user gave a model, use that
    result = (model, model.parentDir)
  elif model.len > 0 and outpath.len > 0:
    result = (model, outpath)
  else:
    result = ("", outpath)

proc initDesc(calib, back, sim: seq[string], # data
              model, modelOutpath, plotPath: string,
              datasets: seq[string], # datasets used as input neurons
              numHidden: seq[int], # number of neurons on each hidden layer
              activation: ActivationFunction,
              outputActivation: OutputActivation,
              lossFunction: LossFunction,
              optimizer: OptimizerKind,
              learningRate: float,
              subsetPerRun: int,
              simulatedData: bool,
              rngSeed: int,
              backgroundRegion: ChipRegion,
              nFake: int,
              plotEvery: int,
              backgroundChips: set[uint8]): MLPDesc =
  if numHidden.len == 0:
    raise newException(ValueError, "Please provide a number of neurons for the hidden layers.")
  # 1. initialize the MLPDesc from the given parameters
  let plotPath = if plotPath.len == 0: "/tmp/" else: plotPath
  # 2. check if `modelOutpath` is a file or a directory
  let (modelFile, modelDir) = modelDirFile(model, modelOutpath)

  result = initMLPDesc(calib, back, sim, datasets,
                       modelFile, modelDir, plotPath,
                       numHidden,
                       activation, outputActivation, lossFunction, optimizer,
                       learningRate,
                       subsetPerRun,
                       simulatedData,
                       rngSeed,
                       backgroundRegion,
                       nFake,
                       plotEvery,
                       backgroundChips)

  # 3. check if such a file already exists to possibly merge it or just return that
  let outfile = result.modelDir / MLPDescName
  # potentially create the output paths
  discard existsOrCreateDir(result.plotPath)
  discard existsOrCreateDir(result.modelDir)
  if not fileExists(outfile):
    ## Existing file with same name. Possibly an older version of it?
    ## Try deserializing as a V1 MLPDesc
    let (version, newestFile) = findNewestFile(result.modelDir)
    if newestFile.len > 0:
      ## NOTE: for now we only have one older version to worry about. In the future we might
      ## need to generalize it based on the actual version.
      let descV1 = deserializeH5[MLPDescV1](result.modelDir / newestFile)
      if descV1.numHidden != 0: # means it really was V1. Therefore copy over
        ## As a V1, just copy over the training related fields
        macro copyFields(fs: varargs[untyped]): untyped =
          result = newStmtList()
          for f in fs:
            result.add quote do:
              result.`f` = descV1.`f`
        copyFields(epochs, accuracies, testAccuracies, losses, testLosses)
        if result.datasets.len == 0:
          result.datasets = descV1.datasets
          result.numInputs = result.datasets.len
        else:
          doAssert descV1.datasets == result.datasets, "Datasets in existing file and input don't match!"
    # no file exists at all
    # write back the either modified new version MLPDesc object or initial file
    result.toH5(outfile)
  else:
    ## is actually V2 or higher, just return it! This branch is
    # Note: most input parameters are ignored in this case!
    let input = result
    result = deserializeH5[MLPDesc](outfile)
    result.inputModel = result.inputModel.strip(chars = Whitespace + {'\0'})
    # Update all fields user wishes to change (that are supported!)
    var anySet = false
    template setAny(field: untyped): untyped =
      if input.field != result.field:
        result.field = input.field
        anySet = true

    if input.inputModel.len > 0:
      result.inputModel = input.inputModel
      anySet = true
    if input.learningRate != result.learningRate:
      # user wishes to change learning rate, add last + epoch to rates
      doAssert result.epochs.len > 0, "This should not happen. We did not train before?"
      result.pastLearningRates.add (lr: result.learningRate, toEpoch: result.epochs[^1])
      result.learningRate = input.learningRate
      anySet = true
    setAny(path)
    setAny(plotPath)
    setAny(modelDir)
    setAny(path)
    setAny(version)
    setAny(plotEvery)
    if anySet: # If anything changed due to input, write new H5 file
      result.toH5(result.modelDir / MLPDescName)

  # update the global datasets!
  CurrentDsets = result.datasets

## XXX: inside of a generic proc (what we would normally do) the `parameters` call
## breaks! Nim doesn't understand that the MLP / ConvNet type can be converted
## to a `Module`!
proc trainModel[T](_: typedesc[T],
                   device: Device,
                   mlpDesc: MLPDesc,
                   epochs = 100_000,
                   skipTraining = false
                  ) =
  ## If we are training, construct a type appropriate to
  var model = T.init(mlpDesc)
  model.to(device)
  if mlpDesc.inputModel.len > 0:
    echo "Loading model file: `", mlpDesc.inputModel, "`"
    model.load(mlpDesc.inputModel)
  elif skipTraining:
    raise newException(ValueError, "`skipTraining` only allowed if a trained input model is handed explicitly " &
      "via the `--model` argument.")
  when T is MLP:
    const readRaw = false
  else:
    const readRaw = true
  echo "Reading data"
  var rnd = initRand(mlpDesc.rngSeed)
  var df = if mlpDesc.simulatedData:
             echo "Using simulated data."
             prepareSimDataframe(mlpDesc, readRaw = readRaw, rnd = rnd)
           else:
             prepareMixedDataframe(mlpDesc, readRaw = readRaw)
  # get training & test dataset
  echo "Splitting data into train & test set"
  let (trainTup, testTup) = generateTrainTest(df, rnd)
  let (trainIn, trainTarg) = trainTup
  let (testIn, testTarg) = testTup

  ## Get a mutable MLPDesc so that we can set the `nTrain/Test/...` fields
  var mlpDesc = mlpDesc
  proc assignNumbers(desc: var MLPDesc, df: DataFrame, nTrain, nTest: int) =
    let dfS = df.filter(f{`Type` == $dtSignal})
    desc.nFakeTotal = dfS.len
    desc.nBack = df.len - dfS.len
    desc.nTrain = nTrain
    desc.nTest = nTest
  mlpDesc.assignNumbers(df, trainIn.sizes[0], testIn.sizes[0])


  ## Continue training from the last epoch if any training already done
  let continueAfterEpoch = if mlpDesc.epochs.len > 0: mlpDesc.epochs[^1] else: 0

  let lr = mlpDesc.learningRate
  if not skipTraining:
    template callTrain(arg: typed): untyped =
      model.train(arg,
                  trainIn.to(kFloat32).to(device),
                  trainTarg.to(kFloat32).to(device),
                  testIn.to(kFloat32).to(device),
                  testTarg.to(kFloat32).to(device),
                  device,
                  readRaw,
                  mlpDesc,
                  epochs,
                  continueAfterEpoch)
    withOptim(model, mlpDesc): # injects `optimizer` variable of correct type into body
      callTrain(optimizer)
    model.save(mlpDesc.path)
  else: # perform validation
    # reproduce the accuracy and loss plots
    plotType(mlpDesc.epochs, mlpDesc.accuracies, mlpDesc.testAccuracies, "Accuracy", outfile = &"{mlpDesc.plotPath}/accuracy.pdf")
    plotType(mlpDesc.epochs, mlpDesc.losses, mlpDesc.testLosses, "Loss", outfile = &"{mlpDesc.plotPath}/loss.pdf")
    # reproduce the training validation plots by using `test` procedure & overriding the plot names
    let (acc, loss, trainPredict, trainTargets) = model.test(trainIn, trainTarg, device, mlpDesc,
                                                             plotOutfile = mlpDesc.plotPath / "train_validation.pdf")
    echo "Train loss after training: ", loss, " with accuracy ", acc
    let preds = trainPredict.mapIt(clamp(it, -ClampOutput, ClampOutput))
    rocCurve(preds, trainTargets, plotPath = mlpDesc.plotPath)

  let (acc, loss, testPredict, testTargets) = model.test(testIn, testTarg, device, mlpDesc,
                                                         plotOutfile = mlpDesc.plotPath / "test_validation.pdf")
  echo "Test loss after training: ", loss, " with accuracy ", acc

template predictModel(Typ: typedesc,
                      device: Device,
                      desc: MLPDesc,
                      ε = 0.8, # signal efficiency for background rate prediction
                      totalTime = -1.0.Hour, # total background rate time in hours. Normally read from input file
                      rocCurve = false,
                     ): untyped {.dirty.} =
  var model = Typ.init(desc)
  model.to(device)
  when Typ is MLP:
    const readRaw = false
  else:
    const readRaw = true

  # load the model
  echo "Loading model file: ", desc.inputModel
  model.load(desc.inputModel)
  model.predictBackground(device, desc, desc.backFiles[0], ε, totalTime, readRaw, desc.plotPath)
  model.predictCalib(device, desc, desc.calibFiles[0], ε)

  model.predictAll(device, desc)

  model.determineCdlEfficiency(device, desc, ε, readRaw, plotPath = desc.plotPath)
  model.determineCdlEfficiency(device, desc, ε, readRaw, global = false, plotPath = desc.plotPath)

  model.generateRocCurves(device, desc)

proc main(calib, back, simFiles: seq[string] = @[],
          datasets: seq[string] = @[],
          ε = 0.8, # signal efficiency for background rate prediction
          totalTime = -1.0.Hour, # total background rate time in hours. Normally read from input file
          rocCurve = false,
          predict = false,
          model = "", ## Optional path to an existing model
          modelKind = "MLP", # MLP or ConvNet #model = mkMLP ## parsing an enum here causes weird CT error in cligen :/
          modelOutpath = "", ## Output path where to store the model. Can be empty if `model` given.
          numHidden: seq[int] = @[], ## number of neurons on the hidden layers. One number per layer.
          activation: ActivationFunction = afReLU,
          outputActivation: OutputActivation = ofLinear,
          lossFunction: LossFunction = lfSigmoidCrossEntropy,
          optimizer: OptimizerKind = opSGD,
          learningRate = Inf,
          subsetPerRun = 1000, ## If multiple `backgroundChips` given, read `subsetPerRun` for *each chip!
          backgroundChips = none(set[uint8]),
          plotPath = "",
          clampOutput = 50.0,
          printDefaultDatasets = false,
          simulatedData = false,
          rngSeed = 1337,
          backgroundRegion = crAll,
          nFake = 100_000,
          epochs = 100_000,
          plotEvery = 5000,
          skipTraining = false # If given only read a model and perform test validation
         ) =
  # 1. set up the model
  if printDefaultDatasets:
    echo "Total default datasets: ", CurrentDsets.len
    for d in CurrentDsets:
      echo d
    return

  let backgroundChips = if backgroundChips.isSome: backgroundChips.get else: {3'u8}
  let desc = initDesc(calib, back, simFiles, model, modelOutpath, plotPath,
                      datasets,
                      numHidden,
                      activation, outputActivation, lossFunction, optimizer,
                      learningRate,
                      subsetPerRun,
                      simulatedData,
                      rngSeed,
                      backgroundRegion,
                      nFake,
                      plotEvery,
                      backgroundChips)

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
  let mKind = parseEnum[ModelKind](modelKind)
  if mKind == mkMLP:
    if not predict:
      MLP.trainModel(device,
                     desc,
                     epochs,
                     skipTraining)
    else:
      MLP.predictModel(device,
                       desc,
                       ε, totalTime,
                       rocCurve)
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

  proc argParse[T: enum](dst: var T, dfl: T, a: var ArgcvtParams): bool =
    var val = a.val
    try:
      dst = parseEnum[T](val)
      result = true
    except ValueError:
      raise newException(Exception, "Invalid enum value: " & $val)
  proc argHelp*(dfl: Hour; a: var ArgcvtParams): seq[string] =
    result = @[ a.argKeys, "<value>.Hour", $dfl ]

  ## Allow `Option` inputs
  proc argParse[T](dst: var Option[T], dfl: Option[T],
                   a: var ArgcvtParams): bool =
    ## If this is being called it means we have _something_, so just
    ## dispath to type `T`
    var res: T
    var dfl: T
    result = argParse(res, dfl, a)
    # merge sets
    when T is set:
      if dst.isNone:
        dst = some(res)
      else:
        var dest = dst.get
        dest.incl res
        dst = some(dest)
    else:
      dst = some(res)

  proc argHelp*[T](dfl: Option[T]; a: var ArgcvtParams): seq[string] =
    result = @[ a.argKeys, $T, "none" ]

  import cligen
  dispatch main
