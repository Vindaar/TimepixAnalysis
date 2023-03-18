import flambeau / [flambeau_nn, tensors]
import std / [strformat, os, strutils, sequtils, random, algorithm, options]

import ingrid / [tos_helpers, ingrid_types]
import pkg / [datamancer, unchained, nimhdf5]

# have to include the type definitions
include ./nn_types, ./io_helpers, ./nn_cuts

{.experimental: "views".}

## The batch size we use!
const bsz = 8192 # batch size

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
    .filter(f{`predictions` > -50.0 and `predictions` < 50.0})
  #dfPlt.showBrowser()
  #echo "Number of signals: ", dfPlt.filter(f{`isSignal` == true})
  #echo "Number of backs: ", dfPlt.filter(f{`isSignal` == false})
  #if true: quit()
  ggplot(dfPlt, aes("predictions", fill = "isSignal")) +
    geom_histogram(bins = 100, position = "identity", alpha = some(0.5), hdKind = hdOutline) +
    scale_x_continuous() +
    ggsave(outfile)

proc rocCurve(predictions: seq[float], targets: seq[int],
              suffix = "") =
  ## plots the ROC curve of the predictions vs the targets
  let dfRoc = calcRocCurve(predictions, targets)
  ggplot(dfRoc, aes("sigEff", "backRej")) +
    geom_line() +
    ggsave("/tmp/roc_curve" & suffix & ".pdf")

proc logLValues(df: DataFrame): (seq[float], seq[int]) =
  let logl = df["likelihood", float].map_inline:
    if classify(x) == fcInf:
      50.0
    else: x
  let targets = df["Type", string].map_inline:
    if x == $dtBack: 0
    else: 1
  result = (logL.toSeq1D, targets.toSeq1D)

proc plotLogLRocCurve(df: DataFrame) =
  ## plots the ROC curve of the predictions vs the targets
  let (logl, targets) = logLValues(df)
  rocCurve(logl, targets, "_likelihood")

proc prepareDataframe(fname: string, run: int, readRaw: bool): DataFrame =
  let dfCdl = prepareCdl(readRaw)
  let dfBack = prepareBackground(fname, run, readRaw) # .drop(["centerX", "centerY"])
  echo dfCdl
  echo dfBack
  result = newDataFrame()
  result.add dfCdl
  result.add dfBack
  # create likelihood plot
  #result.plotLikelihoodDist()
  #result.plotLogLRocCurve()

template printTensorInfo(arg: untyped): untyped =
  echo astToStr(arg), ": ", typeof(arg), " on device ", arg.get_device(), " is cuda ", arg.is_cuda()

proc train(model: AnyModel, optimizer: var Optimizer,
           input, target: RawTensor,
           device: Device,
           readRaw: bool) =
  let dataset_size = input.size(0)
  var toPlot = false
  for epoch in 0 .. 100000:
    var correct = 0
    if epoch mod 50 == 0:
      echo "Epoch is:" & $epoch
    if epoch mod 5000 == 0:
      toPlot = true
    var predictions = newSeqOfCap[float](dataset_size)
    var targets = newSeqOfCap[int](dataset_size)
    ## XXX: Adjust the upper end similar to in `test` to get all data!
    for batch_id in 0 ..< dataset_size div bsz:
      # Reset gradients.
      optimizer.zero_grad()

      # minibatch offset in the Tensor
      ## TODO: generalize the batching and make sure to take `all` elements! (currently skip last X)
      let offset = batch_id * bsz
      var x: RawTensor
      if not readRaw:
        x = input[offset ..< offset + bsz, _ ]
      else:
        x = input[offset ..< offset + bsz, _, _ ].unsqueeze(1)
        echo "INPUT TENSOR SHAPE ", x.sizes

      let target = target[offset ..< offset + bsz]
      #printTensorInfo(target)
      # Running input through the network
      let output = model.forward(x)
      #printTensorInfo(output)
      let pred = output.argmax(1)
      #printTensorInfo(pred)

      if toPlot:
        # take 0th column
        predictions.add output[_, 0].toNimSeq[:float]
        targets.add target[_, 0].toNimSeq[:int]
        correct += pred.eq(target.argmax(1)).sum().item(int)
      # Computing the loss
      var loss = sigmoid_cross_entropy(output, target)
      # Compute the gradient (i.e. contribution of each parameter to the loss)
      loss.backward()
      # Correct the weights now that we have the gradient information
      optimizer.step()

    if toPlot:
      let train_loss = correct.float / dataset_size.float64() # loss.item(float)
      echo &"\nTrain set: Average loss: {train_loss:.4f} " &
           &"| Accuracy: {correct.float64() / dataset_size.float64():.3f}"

    ## create output plot
    if toPlot:
      plotTraining(predictions, targets, outfile = "/tmp/test_training.pdf")
      let preds = predictions.mapIt(clamp(it, -50.0, 50.0))
      rocCurve(preds, targets)
      toPlot = false

proc test(model: AnyModel,
          input, target: RawTensor,
          device: Device,
          plotOutfile = "/tmp/test_validation.pdf"): (seq[float], seq[int]) =
  ## returns the predictions / targets
  let dataset_size = input.size(0)
  var correct = 0
  var predictions = newSeqOfCap[float](dataset_size)
  var targets = newSeqOfCap[int](dataset_size)
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
      predictions.add output[_, 0].toNimSeq[:float]
      targets.add target[_, 0].toNimSeq[:int]
      correct += pred.eq(target.argmax(1)).sum().item(int)
      # Computing the loss
      # var loss = sigmoid_cross_entropy(output, target)

  let test_loss = correct.float / dataset_size.float64()
  echo &"\nTest set: Average loss: {test_loss:.4f} " &
       &"| Accuracy: {correct.float64() / dataset_size.float64():.3f}"

  ## create output plot
  plotTraining(predictions, targets, outfile = plotOutfile)
  # will fail probably...
  # doAssert target == targets
  let preds = predictions.mapIt(clamp(it, -50.0, 50.0))
  result = (preds, targets)

proc predict(model: AnyModel,
             input, target: RawTensor,
             device: Device,
             cutVal: float): seq[int] =
  ## returns the predictions / targets
  let dataset_size = input.size(0)
  var correct = 0
  var predictions = newSeq[float]()
  no_grad_mode:
    ## XXX: Adjust the upper end similar to in `test` to get all data!
    for batch_id in 0 ..< dataset_size div bsz:
      # minibatch offset in the Tensor
      let offset = batch_id * bsz
      let x = input[offset ..< offset + bsz, _ ].to(device)
      let target = target[offset ..< offset + bsz].to(device)
      # Running input through the network
      let output = model.forward(x)
      # get the larger prediction along axis 1 (the example axis)
      let pred = output.argmax(1)
      # take 0th column
      predictions.add output[_, 0].toNimSeq[:float]
      for i in 0 ..< bsz:
        #if target[i, _].argmax() == output[i, _].argmax():
        #  #echo "Is correct! ", offset + i # in this case it is indeed classified as background
        #  discard
        #else:
        if output[i, 0].item(float) > cutVal:
          #echo "is false! ", offset + i
          result.add (offset + i).int # else add the index of this event that looks like signal
      correct += pred.eq(target.argmax(1)).sum().item(int)
  let test_loss = correct.float / dataset_size.float64()
  echo &"\nPredict set: Average loss: {test_loss:.4f} " &
       &"| Accuracy: {correct.float64() / dataset_size.float64():.3f}"

  let df = toDf(predictions)
    .filter(f{`predictions` > -50.0})
  ggplot(df, aes("predictions")) +
    geom_histogram(bins = 100) +
    ggsave("/tmp/predictions_all_background.pdf")

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

proc plotBackground(data: seq[float], totalTime: Hour) =
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
    ggsave("/tmp/simple_background_rate.pdf")

proc targetSpecificRoc(model: AnyModel, df: DataFrame, device: Device) =
  ## computes the CDL target specific ROC curves for logL and MLP

  # 1. split the input df by target and perform regular ROC calculations
  var dfSplit = newDataFrame()
  for tup, subDf in groups(df.group_by("Target")):
    let bin = tup[0][1].toStr
    let (logL, logLTargets) = subDf.logLValues()
    let (input, target) = subDf.toInputTensor()
    let (predict, mlpTargets) = model.test(input, target, device)
    # combined ROC
    var dfLogLRoc = calcRocCurve(logL, logLTargets)
    let dfMLPRoc = calcRocCurve(predict, mlpTargets)
    var dfRoc = bind_rows([("LogL", dfLogLRoc), ("MLP", dfMLPRoc)],
                          "Type")
    dfRoc["Target"] = bin
    dfSplit.add dfRoc

  ggplot(dfSplit, aes("sigEff", "backRej", color = "Target", shape = "Type")) +
    geom_line() +
    ylim(0.5, 1.0) +
    ggsave("/tmp/roc_curve_combined_split_by_target.pdf")

proc determineCdlEfficiency(model: AnyModel, device: Device, ε: float, readRaw: bool) =
  ##
  let dfCdl = prepareCdl(readRaw)
  let cutVal = determineCutValue(model, device, ε, readRaw)
  # for reference determine cut value based on full data (expect to get out `ε`)
  proc getEff(model: AnyModel, df: DataFrame, outfile: string): float =
    let (cdlInput, cdlTarget) = df.toInputTensor()
    let (cdlPredict, predTargets) = model.test(cdlInput, cdlTarget, device,
                                               plotOutfile = outfile)
    result = determineEff(cdlPredict.sorted, cutVal)
  let totalEff = model.getEff(dfCdl, &"/tmp/cdl_prediction_all_data.pdf")
  echo "Total efficiency = ", totalEff
  for (tup, subDf) in groups(dfCdl.group_by("Target")):
    let target = tup[0][1].toStr
    #echo "Current target: ", target, " data: ", subDf
    let effectiveEff = model.getEff(subDf, &"/tmp/cdl_prediction_{target}.pdf")
    echo "Target: ", target, " eff = ", effectiveEff

proc predictBackground(model: AnyModel, fname: string, ε: float, totalTime: Hour,
                       device: Device,
                       readRaw: bool) =
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
  plotBackground(energies, totalTime)

#converter toModule[T: AnyModel](m: T): Module = Module(m)
#
## XXX: inside of a generic proc (what we would normally do) the `parameters` call
## breaks! Nim doesn'- understand that the MLP / ConvNet type can be converted
## to a `Module`!
template trainModel(Typ: typedesc,
                    fname: string,
                    device: Device,
                    run = 186, # default background run to use. If none use a mix of all
                    ε = 0.8, # signal efficiency for background rate prediction
                    totalTime = -1.0.Hour, # total background rate time in hours. Normally read from input file
                    rocCurve = false,
                    predict = false,
                    modelOutpath = "/tmp/trained_model.pt"
                   ): untyped {.dirty.} =
  var model = Typ.init()
  model.to(device)
  when Typ is MLP:
    const readRaw = false
  else:
    const readRaw = true
  if not predict: # all data that needs CDL & specific run data (i.e. training & test dataset)
    var df = prepareDataframe(fname, run, readRaw = readRaw)
    # get training & test dataset
    let (trainTup, testTup) = generateTrainTest(df)
    let (trainIn, trainTarg) = trainTup
    let (testIn, testTarg) = testTup

    # check if model already exists as trained file
    if not fileExists(modelOutpath):
      # Learning loop
      # Stochastic Gradient Descent
      var optimizer = SGD.init(
        model.parameters(),
        SGDOptions.init(0.005).momentum(0.2)
        #learning_rate = 0.005
      )
      model.train(optimizer,
                  trainIn.to(kFloat32).to(device),
                  trainTarg.to(kFloat32).to(device),
                  device,
                  readRaw)
      model.save(modelOutpath)
    else:
      # load the model
      model.load(modelOutpath)
    # perform validation
    let (testPredict, testTargets) = model.test(testIn, testTarg, device)

    # combined ROC
    if rocCurve: ## XXX: this should not really live here either
      let (logL, logLTargets) = df.logLValues()
      var dfLogLRoc = calcRocCurve(logL, logLTargets)
      let dfMLPRoc = calcRocCurve(testPredict, testTargets)
      let dfRoc = bind_rows([("LogL", dfLogLRoc), ("MLP", dfMLPRoc)],
                            "Type")
      ggplot(dfRoc, aes("sigEff", "backRej", color = "Type")) +
        geom_line() +
        ggsave("/tmp/roc_curve_combined.pdf")

      # target specific roc curves
      targetSpecificRoc(model, df, device)

    ## XXX: also not really serving a purpose here
    # now predict background rate
    model.predictBackground(fname, ε, totalTime, device, readRaw)

    model.determineCdlEfficiency(device, ε, readRaw)
  else:
    doAssert fileExists(modelOutpath), "When using the `--predict` option, the trained model must exist!"
    # load the model
    model.load(modelOutpath)
    model.predictBackground(fname, ε, totalTime, device, readRaw)

    model.determineCdlEfficiency(device, ε, readRaw)

proc main(fname: string, run = 186, # default background run to use
          ε = 0.8, # signal efficiency for background rate prediction
          totalTime = -1.0.Hour, # total background rate time in hours. Normally read from input file
          rocCurve = false,
          predict = false,
          model = "MLP", # MLP or ConvNet #model = mkMLP ## parsing an enum here causes weird CT error in cligen :/
          modelOutpath = "/tmp/trained_model.pt") =
  # 1. set up the model
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
    MLP.trainModel(fname, device, run, ε, totalTime, rocCurve, predict, modelOutpath)
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
  proc argHelp*(dfl: Hour; a: var ArgcvtParams): seq[string] =
    result = @[ a.argKeys, "<value>.Hour", $dfl ]
  import cligen
  dispatch main
