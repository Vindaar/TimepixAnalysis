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

import ingrid / [fake_event_generator, gas_physics]
from pkg / xrayAttenuation import FluorescenceLine
proc generateFakeEvents(rnd: var Rand,
                        ctx: LikelihoodContext,
                        calibInfo: CalibInfo,
                        gains: seq[float],
                        diffusion: seq[float],
                        nFake = 100_000): DataFrame =
  #let fakeEvs = generateAndReconstruct(rnd, 1000,
  #result = toDf(
  var count = 0
  var fakeEvs = newSeqOfCap[FakeEvent](nFake)
  # Note: tfKind and nFake irrelevant for us!
  var fakeDesc = FakeDesc(kind: fkGainDiffusion, gasMixture: initCASTGasMixture())
  while count < nFake:
    if count mod 5000 == 0:
      echo "Generated ", count, " events."
    # 1. sample an energy
    let energy = rnd.rand(0.1 .. 10.0).keV
    let lines = @[FluorescenceLine(name: "Fake", energy: energy, intensity: 1.0)]
    # 2. sample a gas gain
    let G = rnd.gauss(mu = (gains[1] + gains[0]) / 2.0, sigma = (gains[1] - gains[0]) / 4.0)
    let gain = GainInfo(N: 100_000.0, G: G, theta: rnd.rand(0.4 .. 2.4))
    # 3. sample a diffusion
    let σT = rnd.gauss(mu = 660.0, sigma = (diffusion[1] - diffusion[0] / 4.0))
    #let σT = rnd.rand(diffusion[0] .. diffusion[1])
    # update σT field of FakeDesc
    fakeDesc.σT = σT
    # now generate and add
    let fakeEv = rnd.generateAndReconstruct(fakeDesc, lines, gain, calibInfo, energy, ctx)
    if not fakeEv.valid:
      continue
    fakeEvs.add fakeEv
    inc count
  result = fakeToDf( fakeEvs )
  result["Type"] = $dtSignal

import ingridDatabase / databaseRead
proc prepareSimDataframe(mlpDesc: MLPDesc, readRaw: bool, rnd: var Rand): DataFrame =
  ## Generates a dataframe to train on that contains real background events and
  ## simulated X-rays. The calibration input files are only used to determine the
  ## boundaries of the gas gain and diffusion in which to generate events in!
  # for now we just define hardcoded bounds
  let gains = @[2400.0, 4000.0]
  let diffusion = @[550.0, 700.0] # μm/√cm, will be converted to `mm/√cm` when converted to DF
  # the theta parameter describing the Pólya also needs to be varied somehow
  let ctx = initLikelihoodContext(CdlFile,
                                  year = yr2018,
                                  energyDset = igEnergyFromCharge,
                                  region = crSilver,
                                  timepix = Timepix1,
                                  morphKind = mkLinear) # morphing to plot interpolation
  var dfSim = newDataFrame()
  for c in mlpDesc.calibFiles:
    withH5(c, "r"):
      let calibInfo = h5f.initCalibInfo()
      dfSim.add rnd.generateFakeEvents(ctx, calibInfo, gains, diffusion, mlpDesc.nFake)

  var dfBack = newDataFrame()
  for b in mlpDesc.backFiles:
    dfBack.add prepareAllBackground(b, readRaw, subsetPerRun = mlpDesc.subsetPerRun * 6,
                                    region = mlpDesc.backgroundRegion) # .drop(["centerX", "centerY"])
  echo "Simulated: ", dfSim
  echo "Back: ", dfBack
  result = newDataFrame()
  result.add dfSim.drop(["runNumber", "likelihood"])
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
    dfBack.add prepareAllBackground(b, readRaw, subsetPerRun = spr * 6) # .drop(["centerX", "centerY"])
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
           continueAfterEpoch = 0) =
  let dataset_size = input.size(0)
  var toPlot = false

  var mlpDesc = desc # local mutable copy to store losses, accuracy etc in
  let plotPath = mlpDesc.plotPath

  const PlotEvery = 5000
  let start = continueAfterEpoch
  let stop = start + 100000
  for epoch in start .. stop:
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
      plotTraining(predictions, targets, outfile = &"{plotPath}/training_output.pdf")
      plotTraining(tOut, tTarget, outfile = &"{plotPath}/validation_output.pdf")
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
    plotTraining(predictions, targets, outfile = plotOutfile)
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
    echo "Target ", suffix, " ", target, " eff = ", effectiveEff

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
                                    morphKind = mkLinear) # morphing to plot interpolation

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
  plotRocCurve(dfLnL, dfMlp, "_lnL_vs_mlp", desc.plotPath)
  # target specific roc curves
  targetSpecificRoc(dfLnL, dfMlp, desc.plotPath)

proc initDesc(calib, back: seq[string], # data
              modelOutpath, plotPath: string,
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
              nFake: int): MLPDesc =
  if numHidden.len == 0:
    raise newException(ValueError, "Please provide a number of neurons for the hidden layers.")
  # 1. initialize the MLPDesc from the given parameters
  let plotPath = if plotPath.len == 0: "/tmp/" else: plotPath
  result = initMLPDesc(calib, back, datasets,
                       modelOutpath, plotPath,
                       numHidden,
                       activation, outputActivation, lossFunction, optimizer,
                       learningRate,
                       subsetPerRun,
                       simulatedData,
                       rngSeed,
                       backgroundRegion,
                       nFake)

  # 2. check if such a file already exists to possibly merge it or just return that
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
    ## is actually V2, just return it! This branch is
    # Note: input parameters are ignored in this case!
    result = deserializeH5[MLPDesc](outfile)
  # update the global datasets!
  CurrentDsets = result.datasets

## XXX: inside of a generic proc (what we would normally do) the `parameters` call
## breaks! Nim doesn't understand that the MLP / ConvNet type can be converted
## to a `Module`!
proc trainModel[T](_: typedesc[T],
                   device: Device,
                   mlpDesc: MLPDesc,
                   continueAfterEpoch = -1
                  ) =
  ## If we are training, construct a type appropriate to
  var model = T.init(mlpDesc)
  model.to(device)
  if continueAfterEpoch > 0:
    model.load(mlpDesc.path)
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
  # check if model already exists as trained file
  let lr = mlpDesc.learningRate
  if not fileExists(mlpDesc.path) or continueAfterEpoch > 0:
    template callTrain(arg: typed): untyped =
      model.train(arg,
                  trainIn.to(kFloat32).to(device),
                  trainTarg.to(kFloat32).to(device),
                  testIn.to(kFloat32).to(device),
                  testTarg.to(kFloat32).to(device),
                  device,
                  readRaw,
                  mlpDesc,
                  continueAfterEpoch)
    withOptim(model, mlpDesc): # injects `optimizer` variable of correct type into body
      callTrain(optimizer)
    model.save(mlpDesc.path)
  else:
    # load the model
    model.load(mlpDesc.path)
  # perform validation
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
  model.load(desc.path)
  model.predictBackground(device, desc, desc.backFiles[0], ε, totalTime, readRaw, desc.plotPath)
  model.predictCalib(device, desc, desc.calibFiles[0], ε)

  model.predictAll(device, desc)

  model.determineCdlEfficiency(device, desc, ε, readRaw, plotPath = desc.plotPath)
  model.determineCdlEfficiency(device, desc, ε, readRaw, global = false, plotPath = desc.plotPath)

  model.generateRocCurves(device, desc)

proc main(calib, back: seq[string] = @[],
          datasets: seq[string] = @[],
          ε = 0.8, # signal efficiency for background rate prediction
          totalTime = -1.0.Hour, # total background rate time in hours. Normally read from input file
          rocCurve = false,
          predict = false,
          model = "MLP", # MLP or ConvNet #model = mkMLP ## parsing an enum here causes weird CT error in cligen :/
          modelOutpath = "/tmp/trained_model.pt",
          numHidden: seq[int] = @[], ## number of neurons on the hidden layers. One number per layer.
          activation: ActivationFunction = afReLU,
          outputActivation: OutputActivation = ofLinear,
          lossFunction: LossFunction = lfSigmoidCrossEntropy,
          optimizer: OptimizerKind = opSGD,
          learningRate = Inf,
          subsetPerRun = 1000,
          plotPath = "",
          clampOutput = 50.0,
          printDefaultDatasets = false,
          simulatedData = false,
          continueAfterEpoch = -1,
          rngSeed = 1337,
          backgroundRegion = crAll,
          nFake = 100_000) =
  # 1. set up the model
  if printDefaultDatasets:
    echo "Total default datasets: ", CurrentDsets.len
    for d in CurrentDsets:
      echo d
    return

  let desc = initDesc(calib, back, modelOutpath, plotPath,
                      datasets,
                      numHidden,
                      activation, outputActivation, lossFunction, optimizer,
                      learningRate,
                      subsetPerRun,
                      simulatedData,
                      rngSeed,
                      backgroundRegion,
                      nFake)

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
                     desc,
                     continueAfterEpoch)
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
  import cligen
  dispatch main
