import flambeau/[flambeau_nn]
import flambeau / tensors
import strformat


import nimhdf5
import ingrid / [tos_helpers, ingrid_types]
import os, strutils, sequtils, random, algorithm, options, cligen
import datamancer

{.experimental: "views".}

let bsz = 100 # batch size

# We will build the following network:
# Input --> Linear(out_features = 12) --> relu --> Linear(out_features = 1)

defModule:
  type
    XorNet* = object of Module
      hidden* = Linear(12, 500)
      #hidden2* = Linear(100, 100)
      classifier* = Linear(500, 2)
      #conv2* = Conv2d(10, 20, 5)
      #conv2_drop* = Dropout2d()
      #fc1* = Linear(320, 50)
      #fc2* = Linear(50, 10)

proc forward(net: XorNet, x: RawTensor): RawTensor =
  #var x = net.hidden2.forward(net.hidden.forward(x).relu()).relu()
  var x = net.hidden.forward(x).relu()
  return net.classifier.forward(x).squeeze(1)

proc readDset(h5f: H5File, grpName: string, igKind: InGridDsetKind): seq[float64] =
  result = h5f.readAs(grp_name / igKind.toDset(fkTpa), float64)

let validReadDSets = XrayReferenceDsets - { igNumClusters,
                                           igFractionInHalfRadius,
                                           igRadiusDivRmsTrans,
                                           igRadius, igBalance,
                                           igLengthDivRadius,
                                           igTotalCharge } + {
                                           igCenterX, igCenterY }

let validDsets = validReadDSets - { igLikelihood, igCenterX, igCenterY, igHits, igEnergyFromCharge }

proc readDsets(h5f: H5File, grpName: string): DataFrame =
  result = newDataFrame()
  for dst in XrayReferenceDsets:
    if dst in validReadDsets:
      result[dst.toDset(fkTpa)] = h5f.readDset(grpName, dst)
  result["eventNumber"] = h5f.readAs(grp_name / "eventNumber", int)

proc buildLogLHist*(h5f: H5file, dset: string): seq[bool] =
  ## returns the a boolean to either keep an event or throw it out.
  ## `true` if we keep it
  var grp_name = cdlPrefix($yr2018) & dset
  # create global vars for xray and normal cuts table to avoid having
  # to recreate them each time
  let xrayCutsTab = getXrayCleaningCuts()
  let cutsTab = getEnergyBinMinMaxVals2018()
  var frameworkKind = fkTpa
  if "FrameworkKind" in h5f.attrs:
    frameworkKind = parseEnum[FrameworkKind](h5f.attrs["FrameworkKind", string])
  # open h5 file using template
  let
    energyStr = igEnergyFromCharge.toDset(frameworkKind)
    centerXStr = igCenterX.toDset(frameworkKind)
    centerYStr = igCenterY.toDset(frameworkKind)
    eccStr = igEccentricity.toDset(frameworkKind)
    lengthStr = igLength.toDset(frameworkKind)
    chargeStr = igTotalCharge.toDset(frameworkKind)
    rmsTransStr = igRmsTransverse.toDset(frameworkKind)
    npixStr = igHits.toDset(frameworkKind)

  let
    energy = h5f.readAs(grp_name / energyStr, float64)
    centerX = h5f.readAs(grp_name / centerXStr, float64)
    centerY = h5f.readAs(grp_name / centerYStr, float64)
    ecc = h5f.readAs(grp_name / eccStr, float64)
    length = h5f.readAs(grp_name / lengthStr, float64)
    charge = h5f.readAs(grp_name / chargeStr, float64)
    rmsTrans = h5f.readAs(grp_name / rmsTransStr, float64)
    npix = h5f.readAs(grp_name / npixStr, float64)
    # get the cut values for this dataset
    cuts = cutsTab[dset]
    xrayCuts = xrayCutsTab[dset]
  result = newSeq[bool](energy.len)
  for i in 0 .. energy.high:
    let
      # first apply Xray cuts (see C. Krieger PhD Appendix B & C)
      regionCut  = inRegion(centerX[i], centerY[i], crSilver)
      xRmsCut    = rmsTrans[i] >= xrayCuts.minRms and rmsTrans[i] <= xrayCuts.maxRms
      xLengthCut = length[i]   <= xrayCuts.maxLength
      xEccCut    = ecc[i]      <= xrayCuts.maxEccentricity
      # then apply reference cuts
      chargeCut  = charge[i]   > cuts.minCharge and charge[i]   < cuts.maxCharge
      rmsCut     = rmsTrans[i] > cuts.minRms    and rmsTrans[i] < cuts.maxRms
      lengthCut  = length[i]   < cuts.maxLength
      pixelCut   = npix[i]     > cuts.minPix
    # add event to likelihood if all cuts passed
    if allIt([regionCut, xRmsCut, xLengthCut, xEccCut, chargeCut, rmsCut, lengthCut, pixelCut], it):
      # take it
      result[i] = true

proc shuff(x: seq[int]): seq[int] =
  result = x
  result.shuffle()

proc prepareCDL(): DataFrame =
  var h5f = H5file("/home/basti/CastData/data/CDL_2019/calibration-cdl-2018.h5", "r")
  let tb = getXrayRefTable()
  var df = newDataFrame()
  for k, bin in tb:
    var dfLoc = h5f.readDsets(cdlPrefix($yr2018) & bin)
    let pass = h5f.buildLogLHist(bin)
    doAssert pass.len == dfLoc.len
    dfLoc["pass?"] = pass
    dfLoc["Target"] = bin
    # remove all that don't pass
    dfLoc = dfLoc.filter(f{bool: idx("pass?") == true})
    df.add dfLoc
  discard h5f.close()
  df["Type"] = "signal"
  result = df
  result["Idx"] = toSeq(0 ..< result.len).shuff()
  result = result.arrange("Idx")
  #result = result[0 ..< 55000]
  result.drop("Idx")
  echo result

proc prepareBackground(h5f: H5File, run: int): DataFrame =
  let path = "/reconstruction/run_$#/chip_3" % $run
  let dsets = toSeq(validReadDsets).mapIt(it.toDset(fkTpa))
  let evNumDset = "eventNumber"
  let grp = h5f[path.grp_str]
  for dset in dsets:
    result[dset] = h5f.readAs(grp.name / dset, float)
  result["eventNumber"] = h5f.readAs(grp.name / "eventNumber", int)
  result["pass?"] = true # does not matter!
  #result["runNumber
  result["Type"] = "back"

  let energyBins = getEnergyBinning()
  let targetTab = getXrayRefTable()
  ## TODO: the following is broken? CHECK!
  # result = result.mutate(f{"Target" ~ targetTab[energyBins.lowerBound(idx("energyFromCharge"))]})
  result = result.mutate(f{"Target" ~ `energyFromCharge`.toRefDset}) # ["energyFromCharge", float].toRawSeq.mapIt(targetTab[energyBins.lowerBound(it)])

proc prepareBackground(fname: string, run: int): DataFrame =
  var h5f = H5open(fname, "r")
  result = h5f.prepareBackground(run)
  discard h5f.close()

proc prepareAllBackground(fname: string): DataFrame =
  var h5f = H5open(fname, "r")
  for run, grp in runs(h5f):
    var df = h5f.prepareBackground(run)
    df["runNumber"] = run
    result.add df

  # filter gold region
  result = result.filter(f{float -> bool: inRegion(`centerX`, `centerY`, crGold)})
  discard h5f.close()

proc toInputTensor(df: DataFrame): (RawTensor, RawTensor) {.noInit.} =
  ## Converts an appropriate data frame to a tuple of input / target tensors
  let cols = validDsets.card
  var input = rawtensors.zeros(df.len * cols).reshape(sizes = [df.len.int64, cols].asTorchView())
  for i, c in validDsets.toSeq.mapIt(it.toDset(fkTpa)).sorted:
    let xp = fromBlob[float](cast[pointer](df[c, float].unsafe_raw_offset()), df.len).convertRawTensor()
    echo xp.sizes(), " weird ", c
    input[_, i] = xp
  var target = rawtensors.zeros(df.len * 2).reshape([df.len.int64, 2].asTorchView())
  let typ = df["Type", string]
  for i in 0 ..< typ.size:
    if typ[i] == "signal":
      target[i, _] = [1, 0].toRawTensorFromScalar #.toTensor.convertRawTensor()
    else:
      target[i, _] = [0, 1].toRawTensorFromScalar #toTensor.convertRawTensor()
  result = (input, target)

proc generateTrainTest(df: var DataFrame):
                      ((RawTensor, RawTensor), (RawTensor, RawTensor)) {.noInit.} =
  # - generate
  # - generate random numbers from 0 to df length
  # - add as column and sort by it
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

proc toNimSeq[T](t: RawTensor): seq[T] =
  doAssert t.sizes().len == 1
  result = newSeq[T](t.size(0))
  for i in 0 ..< result.len:
    result[i] = t[i].item(T)

proc plotLikelihoodDist(df: DataFrame) =
  echo df
  let df = df.mutate(f{float -> float: "likelihood" ~ (if classify(idx("likelihood")) == fcInf: 50.0
                                                       else: idx("likelihood"))})
  ggplot(df, aes("likelihood", fill = "Type")) +
    geom_histogram(bins = 100, position = "identity", alpha = some(0.5), hdKind = hdOutline) +
    scale_x_continuous() +
    ggsave("/tmp/likelihood.pdf")

proc plotTraining(predictions: seq[float], targets: seq[int]) =
  echo "INPUT ", predictions.len
  echo "TARG ", targets.len
  let dfPlt = seqsToDf(predictions, targets)
    .mutate(f{"isSignal" ~ `targets` == 1})
    .filter(f{`predictions` > -50.0 and `predictions` < 50.0})
  #dfPlt.showBrowser()
  echo "Number of signals: ", dfPlt.filter(f{`isSignal` == true})
  echo "Number of backs: ", dfPlt.filter(f{`isSignal` == false})
  #if true: quit()
  ggplot(dfPlt, aes("predictions", fill = "isSignal")) +
    geom_histogram(bins = 100, position = "identity", alpha = some(0.5), hdKind = hdOutline) +
    scale_x_continuous() +
    ggsave("/tmp/test.pdf")


proc determineEff(pred: seq[float], cutVal: float,
                  isBackground = true): float =
  ## returns the efficiency given the sorted (!) predictions, a
  ## cut value `cutVal` and whether it's background or signal
  let cutIdx = pred.lowerBound(cutVal)
  result = cutIdx.float / pred.len.float
  if isBackground:
    result = (1.0 - result)

proc calcSigEffBackRej(df: DataFrame, bins: seq[float],
                       isBackground = true): DataFrame =
  ## returns the signal eff and backround rej for all bins of the
  ## given data frame, split by the `bins` column (that is CDL classes)
  let vals = df.arrange("predictions")["predictions", float]
  var effs = newSeqOfCap[float](bins.len)
  for l in bins:
    let eff = determineEff(vals.toRawSeq, l, isBackground = isBackground)
    effs.add eff
  result = seqsToDf({ "eff" : effs,
                      "cutVals" : bins })

proc calcRocCurve(predictions: seq[float], targets: seq[int]): DataFrame =
  # now use both to determine signal and background efficiencies
  # essentially have to generate some binning in `logL` we deem appropriate,
  let dfPlt = seqsToDf(predictions, targets)
    .mutate(f{"isSignal" ~ `targets` == 1})
  const nBins = 1000
  let bins = linspace(predictions.min, predictions.max, nBins)
  let dfSignal = dfPlt.filter(f{`isSignal` == true})
  let dfBackground = dfPlt.filter(f{`isSignal` == false})

  let sigEffDf = calcSigEffBackRej(dfSignal, bins, isBackground = false)
    .rename(f{"sigEff" <- "eff"})
  let backRejDf = calcSigEffBackRej(dfBackground, bins, isBackground = true)
    .rename(f{"backRej" <- "eff"})
  result = innerJoin(sigEffDf, backRejDf, by = "cutVals")

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
    if x == "back": 0
    else: 1
  result = (logL.toRawSeq, targets.toRawSeq)

proc plotLogLRocCurve(df: DataFrame) =
  ## plots the ROC curve of the predictions vs the targets
  let (logl, targets) = logLValues(df)
  rocCurve(logl, targets, "_likelihood")

proc train(model: XorNet, optimizer: var Optimizer,
           input, target: RawTensor,
           device: Device) =
  let dataset_size = input.size(0)
  echo "Starting size ", dataset_size
  var toPlot = false
  for epoch in 0 .. 1000:
    var correct = 0
    echo "Epoch is:" & $epoch
    if epoch mod 25 == 0:
      toPlot = true
    var predictions = newSeqOfCap[float](dataset_size)
    var targets = newSeqOfCap[int](dataset_size)
    for batch_id in 0 ..< dataset_size div bsz:
      # Reset gradients.
      optimizer.zero_grad()

      # minibatch offset in the Tensor
      ## TODO: generalize the batching and make sure to take `all` elements! (currently skip last X)
      let offset = batch_id * bsz
      let x = input[offset ..< offset + bsz, _ ]
      let target = target[offset ..< offset + bsz]
      # Running input through the network
      let output = model.forward(x)
      let pred = output.argmax(1)

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
      plotTraining(predictions, targets)
      let preds = predictions.mapIt(clamp(-it, -50.0, 50.0))
      rocCurve(preds, targets)
      toPlot = false

proc test(model: XorNet,
          input, target: RawTensor,
          device: Device): (seq[float], seq[int]) =
  ## returns the predictions / targets
  let dataset_size = input.size(0)
  var correct = 0
  var predictions = newSeqOfCap[float](dataset_size)
  var targets = newSeqOfCap[int](dataset_size)
  no_grad_mode:
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
      targets.add target[_, 0].toNimSeq[:int]
      correct += pred.eq(target.argmax(1)).sum().item(int)
      # Computing the loss
      # var loss = sigmoid_cross_entropy(output, target)

  let test_loss = correct.float / dataset_size.float64()
  echo &"\nTest set: Average loss: {test_loss:.4f} " &
       &"| Accuracy: {correct.float64() / dataset_size.float64():.3f}"

  ## create output plot
  plotTraining(predictions, targets)
  # will fail probably...
  # doAssert target == targets
  let preds = predictions.mapIt(clamp(-it, -50.0, 50.0))
  result = (preds, targets)

proc predict(model: XorNet,
             input, target: RawTensor,
             device: Device,
             cutVal: float): seq[int] =
  ## returns the predictions / targets
  let dataset_size = input.size(0)
  var correct = 0
  var predictions = newSeq[float]()
  no_grad_mode:
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

  let df = seqsToDf(predictions)
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
  let (hist, bins) = histogram(df["energies"].toTensor(float).toRawSeq,
                               range = (0.0, 20.0), bins = 100)
  result = seqsToDf({ "energies" : bins, "count" : concat(hist, @[0]) })

proc plotBackground(data: seq[float]) =
  let dfE = seqsToDf({"energies" : data})
    .filter(f{`energies` < 20.0})
  var dfH = histogram(dfE)
  dfH["Rate"] = dfH["count"].scaleDset((2401 * 3600).float, 1e5)
  echo dfH
  #ggplot(dfE, aes("energies")) +
  #  geom_histogram(bins = 100) +
  #  ggsave("/tmp/simple_background_rate.pdf")
  ggplot(dfH, aes("energies", "Rate")) +
    geom_histogram(stat = "identity", hdKind = hdOutline) +
    ggsave("/tmp/simple_background_rate.pdf")

proc targetSpecificRoc(model: XorNet, df: DataFrame, device: Device) =
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

proc main(fname: string) =

  let dfCdl = prepareCdl()
  let dfBack = prepareBackground(fname, 186).drop(["centerX", "centerY"])
  echo dfCdl
  echo dfBack
  #if true: quit()
  var df = newDataFrame()
  df.add dfCdl
  df.add dfBack
  # create likelihood plot
  df.plotLikelihoodDist()
  let (logL, logLTargets) = df.logLValues()
  df.plotLogLRocCurve()
  #df.drop(igLikelihood.toDset(fkTpa))
  let (trainTup, testTup) = generateTrainTest(df)
  let (trainIn, trainTarg) = trainTup
  let (testIn, testTarg) = testTup
  Torch.manual_seed(1)
  var device_type: DeviceKind
  if Torch.cuda_is_available():
    echo "CUDA available! Training on GPU."
    device_type = kCuda
  else:
    echo "Training on CPU."
    device_type = kCPU
  let device = Device.init(device_type)

  var model = XorNet.init()
  model.to(device)

  # Stochastic Gradient Descent
  var optimizer = SGD.init(
    model.parameters(),
    SGDOptions.init(0.005).momentum(0.2)
    #learning_rate = 0.005
  )

  # Learning loop
  model.train(optimizer, trainIn.to(kFloat32).to(device), trainTarg.to(kFloat32).to(device), device)
  let (testPredict, testTargets) = model.test(testIn, testTarg, device)

  # combined ROC
  var dfLogLRoc = calcRocCurve(logL, logLTargets)
  let dfMLPRoc = calcRocCurve(testPredict, testTargets)
  let dfRoc = bind_rows([("LogL", dfLogLRoc), ("MLP", dfMLPRoc)],
                        "Type")
  ggplot(dfRoc, aes("sigEff", "backRej", color = "Type")) +
    geom_line() +
    ggsave("/tmp/roc_curve_combined.pdf")

  # target specific roc curves
  targetSpecificRoc(model, df, device)

  # classify

  # determine cut value at 80%
  echo dfMLPRoc
  let dfme = dfMLPRoc.filter(f{`sigEff` >= 0.8}).head(20) # assumes sorted by sig eff (which it is)
  let cutVal = dfme["cutVals", float][0]
  echo "Cut value: ", cutVal
  echo dfMLPRoc.filter(f{`sigEff` >= 0.8}).tail(20)

  let dfAll = prepareAllBackground(fname)
  let (allInput, allTarget) = dfAll.toInputTensor()
  let passedInds = model.predict(allInput, allTarget, device, -cutVal)
  echo "p inds len ", passedInds.len, " compared to all "
  echo dfAll


  # get all energies of these passing events
  let energiesAll = dfAll["energyFromCharge", float]
  var energies = newSeq[float](passedInds.len)
  for i, idx in passedInds:
    energies[i] = energiesAll[idx]
  plotBackground(energies)
when isMainModule:
  dispatch main
