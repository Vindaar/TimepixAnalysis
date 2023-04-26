import flambeau / [flambeau_nn, tensors]
import std / [strformat, os, strutils, sequtils, random, algorithm, options, macros]

import ingrid / [tos_helpers, ingrid_types]
import pkg / [datamancer, unchained, nimhdf5]

# have to include the type definitions
include ./nn_types
import ./io_helpers
include ./nn_cuts

{.experimental: "views".}

proc generateTrainTest(df: var DataFrame):
                      ((RawTensor, RawTensor), (RawTensor, RawTensor)) {.noInit.} =
  df["Idx"] = toSeq(0 .. df.high).shuff()
  df = df.arrange("Idx")
  let dfTrain = df[0 .. df.high div 2]
  let dfTest = df[df.high div 2 + 1 .. df.high]
  result = (train: dfTrain.toInputTensor, test: dfTest.toInputTensor)

proc train(model: AnyModel, #optimizer: var Optimizer,#GenericOptimizer,
           input, target: RawTensor,
           testInput, testTarget: RawTensor,
           device: Device,
           readRaw: bool,
           desc: MLPDesc,
           continueAfterEpoch = 0) =
  # initialize the optimizer
  var optimizer = initGenericOptim(model, desc)
  let dataset_size = input.size(0)
  var toPlot = false

  var mlpDesc = desc # local mutable copy to store losses, accuracy etc in
  let plotPath = mlpDesc.plotPath

  const PlotEvery = 5000
  let start = continueAfterEpoch
  let stop = start + 100000
  for epoch in start .. stop:
    for batch_id in 0 ..< (dataset_size.float / bsz.float).ceil.int:
      # Reset gradients.
      optimizer.zero_grad()

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
              rngSeed: int): MLPDesc =
  if numHidden.len == 0:
    raise newException(ValueError, "Please provide a number of neurons for the hidden layers.")
  # 1. initialize the MLPDesc from the given parameters
  let dsets = if datasets.len == 0: CurrentDsets else: datasets
  let plotPath = if plotPath.len == 0: "/tmp/" else: plotPath
  result = initMLPDesc(calib, back, dsets,
                       modelOutpath, plotPath,
                       numHidden,
                       activation, outputActivation, lossFunction, optimizer,
                       learningRate,
                       subsetPerRun,
                       simulatedData,
                       rngSeed)

  # 2. check if such a file already exists to possibly merge it or just return that
  let outfile = result.path.parentDir / MLPDescName
  if fileExists(outfile):
    ## Existing file with same name. Possibly an older version of it?
    ## Try deserializing as a V1 MLPDesc
    let descV1 = deserializeH5[MLPDescV1](outfile)
    if descV1.numHidden != 0: # means it really was V1. Therefore copy over
      ## As a V1, just copy over the training related fields
      macro copyFields(fs: varargs[untyped]): untyped =
        result = newStmtList()
        for f in fs:
          result.add quote do:
            result.`f` = descV1.`f`
      copyFields(epochs, accuracies, testAccuracies, losses, testLosses)
      doAssert descV1.datasets == result.datasets, "Datasets in existing file and input don't match!"
      # write back the now modified new version MLPDesc object
      result.toH5(outfile)
    else:
      ## is actually V2, just return it! This branch is
      # Note: input parameters are ignored in this case!
      result = deserializeH5[MLPDesc](outfile)
  else:
    # potentially create the output path, serialize the object
    discard existsOrCreateDir(result.plotPath)
    discard existsOrCreateDir(result.path.parentDir)
    result.toH5(outfile)
  # update the global datasets!
  CurrentDsets = result.datasets

## XXX: inside of a generic proc (what we would normally do) the `parameters` call
## breaks! Nim doesn't understand that the MLP / ConvNet type can be converted
## to a `Module`!
proc trainModel[T](Typ: typedesc[T],
                   device: Device,
                   mlpDesc: MLPDesc,
                   continueAfterEpoch = -1
                  ) = # : untyped {.dirty.} =
  ## If we are training, construct a type appropriate to
  var model = Typ.init(mlpDesc)
  model.to(device)
  if continueAfterEpoch > 0:
    model.load(mlpDesc.path)
  when Typ is MLP:
    const readRaw = false
  else:
    const readRaw = true
  echo "Reading data"
  var df = newDataFrame()
  # get training & test dataset
  echo "Splitting data into train & test set"
  let (trainTup, testTup) = generateTrainTest(df)
  let (trainIn, trainTarg) = trainTup
  let (testIn, testTarg) = testTup
  # check if model already exists as trained file
  let lr = mlpDesc.learningRate
  if not fileExists(mlpDesc.path) or continueAfterEpoch > 0:
    model.train(
                trainIn.to(kFloat32).to(device),
                trainTarg.to(kFloat32).to(device),
                testIn.to(kFloat32).to(device),
                testTarg.to(kFloat32).to(device),
                device,
                readRaw,
                mlpDesc,
                continueAfterEpoch)
    model.save(mlpDesc.path)

proc main(calib, back: seq[string] = @[],
          datasets: seq[string] = @[],
          ε = 0.8, # signal efficiency for background rate prediction
          rocCurve = false,
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
          simulatedData = false,
          continueAfterEpoch = -1,
          rngSeed = 1337) =

  let desc = initDesc(calib, back, modelOutpath, plotPath,
                      datasets,
                      numHidden,
                      activation, outputActivation, lossFunction, optimizer,
                      learningRate,
                      subsetPerRun,
                      simulatedData,
                      rngSeed)
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
    MLP.trainModel(device,
                   desc,
                   continueAfterEpoch)
  #else:
  #  ConvNet.trainModel(fname, device, run, ε, totalTime, rocCurve, predict)

when isMainModule:
  import cligen/argcvt
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
