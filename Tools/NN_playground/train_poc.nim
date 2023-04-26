import flambeau / [flambeau_nn, tensors]
import std / [strformat, os, strutils, sequtils, random, algorithm, options, macros]

import pkg / [cppstl]

# have to include the type definitions
#include ./nn_types
#import ./io_helpers
#include ./nn_cuts

type
  MLPImpl* {.pure, header: "mlp_impl.hpp", importcpp: "MLPImpl".} = object of Module
    hidden*: Linear
    hidden2*: Linear
    classifier*: Linear

  MLP* = CppSharedPtr[MLPImpl]

  ActivationFunction* = enum
    afReLU, afTanh, afELU, afGeLU # , ...?

  OutputActivation* = enum
    ofLinear, ofSigmoid, ofTanh # , ...?

  LossFunction* = enum
    lfSigmoidCrossEntropy, lfMLEloss, lfL1Loss # , ...?

  OptimizerKind* = enum
    opNone, opSGD, opAdam, opAdamW # , opAdaGrad, opAdaBoost ?

  GenericOptimizer* = object
    case kind: OptimizerKind
    of opNone: discard
    of opSGD: sgd: CppSharedPtr[SGD]
    of opAdam: adam: CppSharedPtr[Adam]
    of opAdamW: adamW: CppSharedPtr[AdamW]

  ## A helper object that describes the layers of an MLP
  ## The number of input neurons and neurons on the hidden layer.
  ## This is serialized as an H5 file next to the trained network checkpoints.
  ## On loading a network this file is parsed first and then used to initialize
  ## the correct size of a network.
  ## In addition it contains the datasets that are used for the input.
  MLPDesc* = object
    path*: string # model path to the checkpoint files including the default model name!
    plotPath*: string # path in which plots are placed
    calibFiles*: seq[string] ## Path to the calibration files
    backFiles*: seq[string] ## Path to the background data files
    simulatedData*: bool
    numInputs*: int
    numHidden*: seq[int]
    numLayers*: int
    learningRate*: float
    datasets*: seq[string] # Not `InGridDsetKind` to support arbitrary new columns
    subsetPerRun*: int
    rngSeed*: int
    #
    activationFunction*: ActivationFunction
    outputActivation*: OutputActivation
    lossFunction*: LossFunction
    optimizer*: OptimizerKind
    # fields that store training information
    epochs*: seq[int] ## epochs at which plots and checkpoints are generated
    accuracies*: seq[float]
    testAccuracies*: seq[float]
    losses*: seq[float]
    testLosses*: seq[float]

  ModelKind* = enum
    mkMLP = "MLP"
    mkCNN = "ConvNet"

  ## Placeholder for `ConvNet` above as currently defining both is problematic
  ConvNet* = object

  AnyModel* = MLP | ConvNet

## The batch size we use!
const bsz = 8192 # batch size

proc init*(T: type MLP): MLP =
  result = make_shared(MLPImpl)
  result.hidden = result.register_module("hidden_module", init(Linear, 13, 500))
  result.hidden2 = result.register_module("hidden2_module", init(Linear, 13, 500))
  result.classifier = result.register_module("classifier_module",
      init(Linear, 500, 2))

#func init*(T: type MLPImpl): T {.constructor, importcpp: "MLPImpl::MLPImpl(Linear(13, 500); Linear(500, 2))".}

proc init*(T: type MLP, numInput: int, numLayers: int, numHidden: seq[int], numOutput = 2): MLP =
  result = make_shared(MLPImpl)
  if numLayers != numHidden.len:
    raise newException(ValueError, "Number of layers does not match number of neurons given for each " &
      "hidden layer! Layers: " & $numLayers & " and neurons per layer: " & $numHidden)
  if numLayers > 2:
    raise newException(ValueError, "Only up to 2 hidden layers supported with the `MLPImpl` type.")
  result.hidden = result.register_module("hidden_module", init(Linear, numInput, numHidden[0]))
  if numLayers == 1: # single hidden layer MLP
    result.classifier = result.register_module("classifier_module", init(Linear, numHidden[0], numOutput))
  else: # dual hidden layer MLP
    result.hidden2 = result.register_module("hidden2_module", init(Linear, numHidden[0], numHidden[1]))
    result.classifier = result.register_module("classifier_module", init(Linear, numHidden[1], numOutput))

proc init*(T: type MLP, desc: MLPDesc): MLP =
  result = MLP.init(desc.numInputs, desc.numLayers, desc.numHidden)

proc initMLPDesc*(calib, back, datasets: seq[string],
                  modelPath: string, plotPath: string,
                  numHidden: seq[int],
                  activation: ActivationFunction,
                  outputActivation: OutputActivation,
                  lossFunction: LossFunction,
                  optimizer: OptimizerKind,
                  learningRate: float,
                  subsetPerRun: int,
                  simulatedData: bool,
                  rngSeed: int): MLPDesc =
  result = MLPDesc(calibFiles: calib,
                   backFiles: back,
                   datasets: datasets,
                   path: modelPath, plotPath: plotPath,
                   numInputs: datasets.len,
                   numHidden: numHidden,
                   numLayers: numHidden.len,
                   activationFunction: activation,
                   outputActivation: outputActivation,
                   lossFunction: lossFunction,
                   optimizer: optimizer,
                   learningRate: learningRate,
                   subsetPerRun: subsetPerRun,
                   simulatedData: simulatedData,
                   rngSeed: rngSeed)

proc initGenericOptim*(model: AnyModel, mlpDesc: MLPDesc): GenericOptimizer = # {.noInit.} =
  let lr = mlpDesc.learningRate
  result = GenericOptimizer(kind: mlpDesc.optimizer) #opNone)
  case mlpDesc.optimizer
  of opNone: doAssert false
  of opSGD:
    var sgd = makeShared(SGD)
    var sgdD = sgd.deref
    sgdD = SGD.init(
      model.deref.parameters(),
      SGDOptions.init(lr).momentum(0.2) # .weight_decay(0.001)
    )
    #result = GenericOptimizer(kind: opSGD, sgd: sgd)
    result.sgd = sgd
  of opAdam:
    var adam = makeShared(Adam)
    var adamD = adam.deref
    adamD = Adam.init(
      model.deref.parameters(),
      AdamOptions.init(lr)
    )
    #result = GenericOptimizer(kind: opAdam, adam: adam)
    result.adam = adam
  of opAdamW:
    var adamW = makeShared(AdamW)
    var adamWD = adamW.deref
    adamWD = AdamW.init(
      model.deref.parameters(),
      AdamWOptions.init(lr)
    )
    #result = GenericOptimizer(kind: opAdamW, adamW: adamW)
    result.adamW = adamW

#template initGenericOptim*(model: AnyModel, mlpDesc: MLPDesc): untyped =
#  let lr = mlpDesc.learningRate
#  var optim: Optimizer = Optimizer()
#  case mlpDesc.optimizer
#  of opSGD:
#    #var sgd = makeShared(SGD)
#    #var sgdD = sgd.deref
#    var sgdD = SGD.init(
#      model.deref.parameters(),
#      SGDOptions.init(lr).momentum(0.2) # .weight_decay(0.001)
#    )
#    #result = GenericOptimizer(kind: opSGD, sgd: sgd)
#    optim = sgdD
#  of opAdam:
#    #var adam = makeShared(Adam)
#    #var adamD = adam.deref
#    var adamD = Adam.init(
#      model.deref.parameters(),
#      AdamOptions.init(lr)
#    )
#    #result = GenericOptimizer(kind: opAdam, adam: adam)
#    optim = adamD
#  of opAdamW:
#    #var adamW = makeShared(AdamW)
#    #var adamWD = adamW.deref
#    var adamWD = AdamW.init(
#      model.deref.parameters(),
#      AdamWOptions.init(lr)
#    )
#    #result = GenericOptimizer(kind: opAdamW, adamW: adamW)
#    optim = adamWD
#  optim

proc `=destroy`*(go: var GenericOptimizer) = discard
#proc `=sink`*(go: var GenericOptimizer, b: GenericOptimizer) = {.error: "Not available".}

proc zero_grad*(go: var GenericOptimizer) =
  case go.kind
  of opNone: doAssert false
  of opSGD: go.sgd.deref.zero_grad()
  of opAdam: go.adam.deref.zero_grad()
  of opAdamW: go.adamW.deref.zero_grad()

proc step*(go: var GenericOptimizer) =
  case go.kind
  of opNone: doAssert false
  of opSGD: go.sgd.deref.step()
  of opAdam: go.adam.deref.step()
  of opAdamW: go.adamW.deref.step()

{.experimental: "views".}

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
    for batch_id in 0 ..< (dataset_size.float / bsz.float).int:
      # Reset gradients.
      optimizer.zero_grad()

proc initDesc(calib, back: seq[string], # data
              modelOutpath, plotPath: string,
              numHidden: seq[int], # number of neurons on each hidden layer
              activation: ActivationFunction,
              outputActivation: OutputActivation,
              lossFunction: LossFunction,
              optimizer: OptimizerKind,
              learningRate: float,
              subsetPerRun: int,
              simulatedData: bool,
              rngSeed: int): MLPDesc =
  result = initMLPDesc(calib, back, @[],
                       modelOutpath, plotPath,
                       numHidden,
                       activation, outputActivation, lossFunction, optimizer,
                       learningRate,
                       subsetPerRun,
                       simulatedData,
                       rngSeed)

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
  # get training & test dataset
  echo "Splitting data into train & test set"
  var
    trainIn, trainTarg, testIn, testTarg: RawTensor
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
  proc argParse[T: enum](dst: var T, dfl: T, a: var ArgcvtParams): bool =
    var val = a.val
    try:
      dst = parseEnum[T](val)
      result = true
    except ValueError:
      raise newException(Exception, "Invalid enum value: " & $val)
  import cligen
  dispatch main
