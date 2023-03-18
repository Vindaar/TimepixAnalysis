import std / [os, sequtils]
import pkg / [datamancer, nimhdf5]
import flambeau/[flambeau_nn, tensors]
import ingrid / [tos_helpers, ingrid_types]

{.experimental: "views".}

# have to include the type definitions
include ./nn_types
import ./io_helpers
include ./nn_cuts

proc readEvent(h5f: H5File, run, idx: int, path = recoBase()): DataFrame =
  let grp = path & $run & "/chip_3" # XXX: use center chip
  result = newDataFrame()
  for dsetName in ValidDsets:
    let dset = h5f[(grp / (dsetName).toDset).dset_str]
    result[dsetName.toDset] = dset.readAs(@[idx], float)
  result["Type"] = $dtBack

proc toTorchTensor(df: DataFrame): RawTensor {.noInit.} =
  ## Converts an appropriate data frame to a 2D RawTensor of the data to be fed to the network
  let cols = ValidDsets.card
  var input = rawtensors.zeros(df.len * cols).reshape(sizes = [df.len.int64, cols].asTorchView())
  for i, c in ValidDsets.toSeq.mapIt(it.toDset(fkTpa)).sorted:
    let xp = fromBlob[float](cast[pointer](df[c, float].unsafe_raw_offset()), df.len).convertRawTensor()
    input[_, i] = xp
  result = input

# if this is ``imported``, we'll create a global classifier that we use for prediction.
# 1. set up the model
template loadModelMakeDevice*(modelPath: string): untyped {.dirty.} =
  var device_type: DeviceKind
  if Torch.cuda_is_available():
    #echo "CUDA available! Training on GPU."
    device_type = kCuda
  else:
    #echo "Training on CPU."
    device_type = kCPU
  let device = Device.init(device_type)
  var model = MLP.init()
  Torch.manual_seed(1)
  model.to(device)
  model.load(modelPath)

proc predict*(h5f: H5File, modelPath: string, run, idx: int): float =
  ## Returns the prediction of the (globally declared!) network for the given run & event index, assuming the
  ## center chip, number 3
  loadModelMakeDevice(modelPath)

  let inp = toTorchTensor(h5f.readEvent(run, idx))
  echo "INP ", inp
  result = model.predictSingle(inp, device)

import ingrid / ingrid_types
proc predict*(h5f: H5File, modelPath: string, grp: string): seq[float] =
  ## Returns the prediction of the (globally declared!) network for the given run & event index, assuming the
  ## center chip, number 3
  ##
  ## Currently this proc implies reading the model and copying it to the GPU each call!
  if modelPath.len == 0:
    raise newException(ValueError, "Cannot predict events without a path to a trained model!")

  loadModelMakeDevice(modelPath)
  # 1. read all required data
  let data = h5f.readValidDsets(grp)
  # 2. convert to torch tensor
  let inp = toTorchTensor(data)
  # 3. forward pass through network
  result = model.predict(inp, device)

proc predict*(model: AnyModel, device: Device, df: DataFrame): seq[float] =
  ## Returns the prediction of the (globally declared!) network for the given data, which
  ## must be a dataframe containing all required datasets!
  let inp = toTorchTensor(df)
  # 3. forward pass through network
  result = model.predict(inp, device)

proc predict*(modelPath: string, df: DataFrame): seq[float] =
  ## Returns the prediction of the (globally declared!) network for the given data, which
  ## must be a dataframe containing all required datasets!
  loadModelMakeDevice(modelPath)
  let inp = toTorchTensor(df)
  # 3. forward pass through network
  result = model.predict(inp, device)

proc predict*(ctx: LikelihoodContext, h5f: H5File, grp: string): seq[float] =
  if not ctx.useNeuralNetworkCut: return
  result = predict(h5f, ctx.nnModelPath, grp)

proc calcNeuralNetCutValueTab*(modelPath: string, cutKind: NeuralNetCutKind, ε: float): CutValueInterpolator =
  if modelPath.len == 0:
    raise newException(ValueError, "Cannot predict events without a path to a trained model!")

  loadModelMakeDevice(modelPath)

  result = initCutValueInterpolator(cutKind)
  case cutKind
  of nkGlobal: result.cut = determineCutValue(model, device, ε, readRaw = false)
  of nkLocal: result.nnCutTab = determineLocalCutValue(model, device, ε, readRaw = false)
  else:
    doAssert false, "Not supported yet!"

proc calcNeuralNetCutValueTab*(ctx: LikelihoodContext): CutValueInterpolator =
  if not ctx.useNeuralNetworkCut: return
  result = calcNeuralNetCutValueTab(ctx.nnModelPath, ctx.nnCutKind, ctx.nnSignalEff)
