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
  for dsetName in CurrentDsets:
    let dset = h5f[(grp / (dsetName)).dset_str]
    result[dsetName] = dset.readAs(@[idx], float)
  result["Type"] = $dtBack

proc initDesc*(model: string): MLPDesc =
  result = initMLPDesc(model, "")
  CurrentDsets = result.datasets

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
  let desc = initDesc(modelPath)
  var model = MLP.init(desc)
  Torch.manual_seed(1)
  model.to(device)
  model.load(modelPath)

proc predict*(h5f: H5File, modelPath: string, run, idx: int): float =
  ## Returns the prediction of the (globally declared!) network for the given run & event index, assuming the
  ## center chip, number 3
  loadModelMakeDevice(modelPath)

  let inp = toTorchTensor(h5f.readEvent(run, idx))
  echo "INP ", inp
  result = model.predictSingle(inp, device, desc)

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
  let data = h5f.readValidDsets(grp, filterNan = false)
  # 2. convert to torch tensor
  let inp = toTorchTensor(data)
  # 3. forward pass through network
  result = model.forward(inp, device, desc)

proc predict*(model: AnyModel, device: Device, desc: MLPDesc, df: DataFrame): seq[float] =
  ## Returns the prediction of the (globally declared!) network for the given data, which
  ## must be a dataframe containing all required datasets!
  let inp = toTorchTensor(df)
  # 3. forward pass through network
  result = model.forward(inp, device, desc)

proc predict*(modelPath: string, df: DataFrame): seq[float] =
  ## Returns the prediction of the (globally declared!) network for the given data, which
  ## must be a dataframe containing all required datasets!
  loadModelMakeDevice(modelPath)
  let inp = toTorchTensor(df)
  # 3. forward pass through network
  result = model.forward(inp, device, desc)

proc predict*(ctx: LikelihoodContext, h5f: H5File, grp: string): seq[float] =
  if not ctx.vetoCfg.useNeuralNetworkCut: return
  result = predict(h5f, ctx.vetoCfg.nnModelPath, grp)

proc determineCutValue*(modelPath: string, df: DataFrame, ε: float): float =
  ## Determines the cut value for the given model & df at `ε`
  loadModelMakeDevice(modelPath)
  result = model.determineCutValue(device, desc, df, ε)

proc calcNeuralNetCutValueTab*(modelPath: string, cutKind: NeuralNetCutKind, ε: float): CutValueInterpolator =
  if modelPath.len == 0:
    raise newException(ValueError, "Cannot predict events without a path to a trained model!")

  result = initCutValueInterpolator(cutKind)
  case cutKind
  of nkGlobal:
    loadModelMakeDevice(modelPath)
    result.cut = determineCutValue(model, device, desc, ε, readRaw = false)
  of nkLocal:
    loadModelMakeDevice(modelPath)
    result.nnCutTab = determineLocalCutValue(model, device, desc, ε, readRaw = false)
  of nkRunBasedLocal:
    echo "[INFO]: nkRunBasedLocal requires run based determination of the local cut values ",
     "using the rmsTransverse data for a reference of the diffusion."
  else:
    doAssert false, "Not supported yet!"

proc calcNeuralNetCutValueTab*(ctx: LikelihoodContext): CutValueInterpolator =
  if not ctx.vetoCfg.useNeuralNetworkCut: return
  result = calcNeuralNetCutValueTab(ctx.vetoCfg.nnModelPath, ctx.vetoCfg.nnCutKind, ctx.vetoCfg.nnSignalEff)

import ingrid / fake_event_generator
from pkg / unchained import FemtoFarad
from std / random import Rand
const CacheTabFile = "/dev/shm/cacheTab_runLocalCutVals.h5"
type
  TabKey = (int, string, float)
  #         ^-- run number
  #              ^-- sha1 hash of the NN model `.pt` file
  #                      ^-- target signal efficiency
  TabVal = seq[(string, float)]
  #             ^-- CDL target
  #                     ^-- MLP cut value
  CacheTabTyp = Table[TabKey, TabVal]
var CacheTab =
  if fileExists(CacheTabFile):
    tryDeserializeH5[CacheTabTyp](CacheTabFile)
  else:
    initTable[TabKey, TabVal]()
proc fileAvailable(run: int, modelHash: string, ε: float): bool =
  if (run, modelHash, ε) in CacheTab:
    result = true
  else:
    # try rereading & updating file
    if fileExists(CacheTabFile):
      let tab = tryDeserializeH5[CacheTabTyp](CacheTabFile)
      # merge `tab` and `CacheTab`
      for k, v in tab:
        CacheTab[k] = v # overwrite possible existing keys in table
      # write merged table
      CacheTab.tryToH5(CacheTabFile)
    result = (run, modelHash, ε) in CacheTab # still not in: not available

import std / sha1

proc fromSeq(s: TabVal): CutValueInterpolator =
  result = initCutValueInterpolator(nkRunBasedLocal)
  for x in s:
    let (target, cutVal) = x
    result.nnCutTab[target] = cutVal

proc toSeq(cutTab: CutValueInterpolator): TabVal =
  result = newSeq[(string, float)]()
  for target, cutVal in cutTab.nnCutTab.pairs:
    result.add (target, cutVal)

proc calcLocalNNCutValueTab*(ctx: LikelihoodContext,
                             rnd: var Rand,
                             h5f: H5File,
                             runType: RunTypeKind,
                             run, chipNumber: int,
                             chipName: string,
                             capacitance: FemtoFarad
                            ): CutValueInterpolator =
  let model = ctx.vetoCfg.nnModelPath
  let modelHash = $(model.readFile.secureHash)
  let ε = ctx.vetoCfg.nnSignalEff
  if fileAvailable(run, modelHash, ε):
    let data = CacheTab[(run, modelHash, ε)]
    result = data.fromSeq()
  else:
    loadModelMakeDevice(model)
    var dfFake = newDataFrame()
    for tf in TargetFilterKind:
      let fakeDesc = FakeDesc(nFake: 5000,
                              tfKind: tf,
                              kind: fkGainDiffusion)
      var dfLoc = generateRunFakeData(rnd, h5f, run, chipNumber, chipName, capacitance, fakeDesc, runType, DataFrame, some(ctx))
      dfLoc["Target"] = $tf
      dfFake.add dfLoc

    # now use fake data to determine cuts
    result = initCutValueInterpolator(nkRunBasedLocal)
    result.nnCutTab = determineRunLocalCutValue(model, device, desc, dfFake, ε)

    CacheTab[(run, modelHash, ε)] = result.toSeq()
    CacheTab.tryToH5(CacheTabFile)

proc main(calib, back: seq[string] = @[],
          model: string,
          subsetPerRun = 1000) =
  discard

when isMainModule:
  import cligen
  dispatch main
