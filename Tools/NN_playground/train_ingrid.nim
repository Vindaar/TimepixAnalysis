import flambeau/[flambeau_nn]
import flambeau / tensors
import strformat


import nimhdf5
import ingrid / [tos_helpers, ingrid_types]
import os, strutils, sequtils, random, algorithm, options, cligen
import datamancer, unchained

{.experimental: "views".}

let bsz = 8192 # batch size

# We will build the following network:
# Input --> Linear(out_features = 12) --> relu --> Linear(out_features = 1)

defModule:
  type
    MLP* = object of Module
      hidden* = Linear(12, 500)
      #hidden2* = Linear(100, 100)
      #hidden2* = Linear(5000, 5000)
      classifier* = Linear(500, 2)
      #conv2* = Conv2d(10, 20, 5)
      #conv2_drop* = Dropout2d()
      #fc1* = Linear(320, 50)
      #fc2* = Linear(50, 10)

proc forward(net: MLP, x: RawTensor): RawTensor =
  #var x = net.hidden2.forward(net.hidden.forward(x).relu()).relu()
  var x = net.hidden.forward(x).relu()
  #x = net.hidden2.forward(x).relu()
  return net.classifier.forward(x).squeeze(1)

## XXX: Defining two models in a single file is currently broken. When trying to use it
## the nim compiler assigns the wrong destructor to the second one (reusing the one
## from the first type). This breaks it. To use it, we currently need to replace
## the logic instead (putting this one first). Once things work we can try to ask Hugo &
## check if submoduling individual models helps.
defModule:
  type
    ConvNet* = object of Module
      conv1* = Conv2d(1, 50, 15)
      conv2* = Conv2d(50, 70, 15)
      conv3* = Conv2d(70, 100, 15)
      #lin1* = Linear(100 * 15 * 15, 1000)
      lin1* = Linear(3610, 1000)
      lin2* = Linear(1000, 50)
      lin3* = Linear(50, 2)

proc forward(net: ConvNet, x: RawTensor): RawTensor =
  var x = net.conv1.forward(x).relu().max_pool2d([2, 2])
  x = net.conv2.forward(x).relu().max_pool2d([2, 2])
  x = net.conv3.forward(x).relu().max_pool2d([2, 2])
  x = net.lin1.forward(x).relu()
  x = net.lin2.forward(x).relu()
  x = net.lin3.forward(x).relu()
  return x

type
  ModelKind = enum
    mkMLP = "MLP"
    mkCNN = "ConvNet"

  AnyModel = MLP | ConvNet


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

proc readRaw(h5f: H5File, grpName: string, idxs: seq[int] = @[]): DataFrame =
  result = newDataFrame()
  let
    xs = h5f[grpName / "x", special_type(uint8), uint8]
    ys = h5f[grpName / "y", special_type(uint8), uint8]
    ev = h5f.readAs(grp_name / "eventNumber", int)
  doAssert xs.len == ev.len
  var xsAll = newSeqOfCap[int](xs.len * 100)
  var ysAll = newSeqOfCap[int](xs.len * 100)
  var evAll = newSeqOfCap[int](xs.len * 100)
  for i in idxs:
    for j in 0 ..< xs[i].len:
      xsAll.add xs[i][j].int
      ysAll.add ys[i][j].int
      evAll.add ev[i]
  result = toDf({"x" : xsAll, "y" : ysAll, "eventNumber" : evAll})


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

proc prepareCDL(readRaw: bool,
                cdlPath = "/home/basti/CastData/data/CDL_2019/calibration-cdl-2018.h5"
               ): DataFrame =
  var h5f = H5file(cdlPath, "r")
  let tb = getXrayRefTable()
  var df = newDataFrame()
  for k, bin in tb:
    var dfLoc = newDataFrame()
    if not readRaw:
      dfLoc = h5f.readDsets(cdlPrefix($yr2018) & bin)
      let pass = h5f.buildLogLHist(bin)
      doAssert pass.len == dfLoc.len
      dfLoc["pass?"] = pass
    else:
      let pass = h5f.buildLogLHist(bin)
      var idxs = newSeqOfCap[int](pass.len)
      for i, p in pass:
        if p:
          idxs.add i
      dfLoc = h5f.readRaw(cdlPrefix($yr2018) & bin, idxs)
      dfLoc["pass?"] = true
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

proc prepareBackground(h5f: H5File, run: int, readRaw: bool): DataFrame =
  let path = "/reconstruction/run_$#/chip_3" % $run
  let dsets = toSeq(validReadDsets).mapIt(it.toDset(fkTpa))
  let evNumDset = "eventNumber"
  let grp = h5f[path.grp_str]
  if not readRaw:
    for dset in dsets:
      result[dset] = h5f.readAs(grp.name / dset, float)
    result["eventNumber"] = h5f.readAs(grp.name / "eventNumber", int)
  else:
    result = h5f.readRaw(grp.name)
  result["pass?"] = true # does not matter!
  #result["runNumber
  result["Type"] = "back"

  let energyBins = getEnergyBinning()
  let targetTab = getXrayRefTable()
  ## TODO: the following is broken? CHECK!
  # result = result.mutate(f{"Target" ~ targetTab[energyBins.lowerBound(idx("energyFromCharge"))]})
  if not readRaw:
    result = result.mutate(f{"Target" ~ `energyFromCharge`.toRefDset}) # ["energyFromCharge", float].toSeq1D.mapIt(targetTab[energyBins.lowerBound(it)])

proc prepareBackground(fname: string, run: int, readRaw: bool): DataFrame =
  var h5f = H5open(fname, "r")
  result = h5f.prepareBackground(run, readRaw)
  discard h5f.close()

proc prepareAllBackground(fname: string, readRaw: bool): DataFrame =
  var h5f = H5open(fname, "r")
  for run, grp in runs(h5f):
    var df = h5f.prepareBackground(run, readRaw)
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
      let isSignal = tup[1][1].toStr == "signal"
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

proc plotTraining(predictions: seq[float], targets: seq[int],
                  outfile = "/tmp/test.pdf") =
  echo "INPUT ", predictions.len
  echo "TARG ", targets.len
  let dfPlt = toDf(predictions, targets)
    .mutate(f{"isSignal" ~ `targets` == 1})
    .filter(f{`predictions` > -50.0 and `predictions` < 50.0})
  #dfPlt.showBrowser()
  echo "Number of signals: ", dfPlt.filter(f{`isSignal` == true})
  echo "Number of backs: ", dfPlt.filter(f{`isSignal` == false})
  #if true: quit()
  ggplot(dfPlt, aes("predictions", fill = "isSignal")) +
    geom_histogram(bins = 100, position = "identity", alpha = some(0.5), hdKind = hdOutline) +
    scale_x_continuous() +
    ggsave(outfile)

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
    let eff = determineEff(vals.toSeq1D, l, isBackground = isBackground)
    effs.add eff
  result = toDf({ "eff" : effs,
                  "cutVals" : bins })

proc calcRocCurve(predictions: seq[float], targets: seq[int]): DataFrame =
  # now use both to determine signal and background efficiencies
  # essentially have to generate some binning in `logL` we deem appropriate,
  let dfPlt = toDf(predictions, targets)
    .mutate(f{"isSignal" ~ `targets` == 1})
  const nBins = 1000
  let bins = linspace(predictions.min, predictions.max, nBins)
  let dfSignal = dfPlt.filter(f{`isSignal` == true})
  let dfBackground = dfPlt.filter(f{`isSignal` == false})

  var
    sigEffDf = newDataFrame()
    backRejDf = newDataFrame()
  if dfSignal.len > 0:
    sigEffDf = calcSigEffBackRej(dfSignal, bins, isBackground = false)
      .rename(f{"sigEff" <- "eff"})
  if dfBackground.len > 0:
    backRejDf = calcSigEffBackRej(dfBackground, bins, isBackground = true)
      .rename(f{"backRej" <- "eff"})
  if sigEffDf.len > 0 and backRejDf.len > 0:
    result = innerJoin(sigEffDf, backRejDf, by = "cutVals")
  elif sigEffDf.len > 0:
    result = sigEffDf
  elif backRejDf.len > 0:
    result = backRejDf
  else:
    doAssert false, "Both signal and background dataframes are empty!"

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
  echo "Starting size ", dataset_size
  var toPlot = false
  for epoch in 0 .. 100000:
    var correct = 0
    if epoch mod 50 == 0:
      echo "Epoch is:" & $epoch
    if epoch mod 5000 == 0:
      toPlot = true
    var predictions = newSeqOfCap[float](dataset_size)
    var targets = newSeqOfCap[int](dataset_size)
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
  let t = if totalTime > 0.0.Hour: totalTime else: 2401.0.Hour
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

proc predictBackground(model: AnyModel, fname: string, ε: float, totalTime: Hour,
                       device: Device,
                       readRaw: bool) =
  # classify
  # determine cut value at 80%
  let dfCdl = prepareCdl(readRaw)
  let (cdlInput, cdlTarget) = dfCdl.toInputTensor()
  let (cdlPredict, predTargets) = model.test(cdlInput, cdlTarget, device,
                                             plotOutfile = "/tmp/cdl_prediction.pdf")
  # Note: `predTargets` & `cdlTarget` contain same data. May not have same order though!
  let dfMLPRoc = calcRocCurve(cdlPredict, predTargets)
  dfMlpRoc.showBrowser()

  let dfme = dfMLPRoc.filter(f{`sigEff` >= ε}).head(20) # assumes sorted by sig eff (which it is)
  let cutVal = dfme["cutVals", float][0]
  echo "Cut value: ", cutVal
  dfMLPRoc.filter(f{`sigEff` >= ε}).showBrowser()

  ## XXX: make sure the cut value stuff & the sign of the predictions is correct everywhere! We flipped the
  ## sign initially to compare with logL!
  let dfAll = prepareAllBackground(fname, readRaw)
  let (allInput, allTarget) = dfAll.toInputTensor()
  let passedInds = model.predict(allInput, allTarget, device, cutVal)
  echo "p inds len ", passedInds.len, " compared to all "
  echo dfAll

  # get all energies of these passing events
  let energiesAll = dfAll["energyFromCharge", float]
  var energies = newSeq[float](passedInds.len)
  for i, idx in passedInds:
    energies[i] = energiesAll[idx]
  plotBackground(energies, totalTime)

converter toModule[T: AnyModel](m: T): Module = Module(m)

proc trainModel[T: AnyModel](_: typedesc[T],
                             fname: string,
                             device: Device,
                             run = 186, # default background run to use. If none use a mix of all
                             ε = 0.8, # signal efficiency for background rate prediction
                             totalTime = -1.0.Hour, # total background rate time in hours. Normally read from input file
                             rocCurve = false,
                             predict = false
                            ) =
  var model = T.init()
  model.to(device)
  when T is MLP:
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
    if not fileExists("/tmp/trained_model.pt"):
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
      model.save("/tmp/trained_model.pt")
    else:
      # load the model
      model.load("/tmp/trained_model.pt")
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
  else:
    doAssert fileExists("/tmp/trained_model.pt"), "When using the `--predict` option, the trained model must exist!"
    # load the model
    model.load("/tmp/trained_model.pt")
    model.predictBackground(fname, ε, totalTime, device, readRaw)

proc main(fname: string, run = 186, # default background run to use
          ε = 0.8, # signal efficiency for background rate prediction
          totalTime = -1.0.Hour, # total background rate time in hours. Normally read from input file
          rocCurve = false,
          predict = false,
          model = "MLP") = # MLP or ConvNet
          #model = mkMLP) = # if true, only predict background rate from all data in input file

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
    MLP.trainModel(fname, device, run, ε, totalTime, rocCurve, predict)
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
  dispatch main
