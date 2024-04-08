import std / [os, sequtils, strutils]
import pkg / [datamancer, nimhdf5]
import flambeau/[flambeau_nn, flambeau_raw]
import arraymancer / tensor
import ingrid / [tos_helpers, ingrid_types]

import latexdsl_nochecks

{.experimental: "views".}

# have to include the type definitions
include ./nn_types
import ./io_helpers
include ./nn_cuts

let tmpl = readFile("template.tex")

proc initDesc*(model: string): MLPDesc =
  result = initMLPDesc(model, "")
  CurrentDsets = result.datasets

proc forward(model: AnyModel, device: Device, desc: MLPDesc, df: DataFrame): (float, float) =
  doAssert df.len == 1
  let inp = toTorchTensor(df)
  # 3. forward pass through network
  result[0] = model.forward(inp, device, desc, neuron = 0)[0]
  result[1] = model.forward(inp, device, desc, neuron = 1)[0]

proc toRow(v: seq[float], addNewline: bool): string =
  result = v.mapIt(&"{it:.4f}").join(" & ")
  if addNewline:
    result = result & r"\\" & "\n"

proc toCol(v: seq[float]): string =
  result = v.mapIt(&"{it:.4f}").join(r" \\ ")

proc toTeX(t: RawTensor): string =
  for row in 0 ..< t.sizes[0]:
    var rowData = newSeq[float](t.sizes[1])
    for j in 0 ..< t.sizes[1]:
      rowData[j] = t[row.int, j.int].item(float)
    let last = row == t.sizes[0] - 1
    result.add toRow(rowData, addNewline = not last)

proc toTeXVec(t: RawTensor): string =
  var rowData = newSeq[float](t.sizes[0])
  for row in 0 ..< t.sizes[0]:
    rowData[row] = t[row.int].item(float)
  result = toCol(rowData)

proc manualForward(x: Tensor[float], W_h1, b1, W_h2, b2, W_c, b_c: Tensor[float]): (float, float) =
  ## Manual forward pass
  doAssert @(x.shape) == @[14], "No x.sizes was " & $x.shape
  let l1 = tanh(W_h1 * x  + b_1)
  let l2 = tanh(W_h2 * l1 + b_2)
  let res = W_c * l2 + b_c
  doAssert res.rank == 1, "No, res.rank was " & $res.rank
  doAssert @(res.shape) == @[2], "No, res.shape was " & $res.shape
  result = (res[0], res[1])

proc toAmancer(t: RawTensor): Tensor[float] =
  let shape = @(t.sizes.asNimView()).mapIt(it.int)
  result = zeros[float](shape)
  doAssert t.sizes.len == 2, "No, t.sizes.len was " & $t.sizes.len
  for row in 0 ..< t.sizes[0]:
    for col in 0 ..< t.sizes[1]:
      result[row.int, col.int] = t[row.int, col.int].item(float)

proc toAmancerVec(t: RawTensor): Tensor[float] =
  let shape = @(t.sizes.asNimView()).mapIt(it.int)
  result = zeros[float](shape)
  doAssert t.sizes.len == 1, "No, t.sizes.len was " & $t.sizes.len
  for col in 0 ..< t.sizes[0]:
    result[col.int] = t[col.int].item(float)

proc writeTeX(W_h1, b1, W_h2, b2, W_c, bc: RawTensor): string =
  result = tmpl % [ W_h1.toTex(), b1.toTexVec(),
                    W_h2.toTex(), b2.toTexVec(),
                    W_c.toTex(), bc.toTexVec() ]


template loadModelMakeDevice*(modelPath: string): untyped {.dirty.} =
  var device_type: DeviceKind
  if Torch.cuda_is_available():
    echo "CUDA available! Training on GPU."
    device_type = kCuda
  else:
    echo "Training on CPU."
    device_type = kCPU
  let device = Device.init(device_type)
  let desc = initDesc(modelPath)
  var model = MLP.init(desc)
  Torch.manual_seed(1)
  model.to(device)
  echo "Loading model: ", modelPath
  model.load(modelPath, device)

import mlp_schematic
proc rawToATensor(t: RawTensor): Tensor[float] =
  echo "T sizes: ", t.sizes
  result = zeros[float](t.numel().int)
  let sh = [t.numel()]
  let tx = t.reshape(sh.asTorchView())
  for i in 0 ..< t.numel():
    result[i.int] = tx[i.int].item(float)
  result = result.reshape(@(t.sizes.asNimView).mapIt(it.int))
  echo "result shape: ", result.shape

proc generateMlpSchematic(desc: MLPDesc, mlp: MLP, file: string) =
  ## Plots a schematic of the given MLP
  doAssert desc.numHidden.len == 2
  let w1 = mlp.hidden.weight.rawToATensor
  let w2 = mlp.hidden2.weight.rawToATensor
  let w3 = mlp.classifier.weight.rawToATensor
  plotSchematic(desc.numInputs, desc.numHidden[0], desc.numHidden[1], 2, # 2 output neurons
                w1, w2, w3,
                desc.datasets,
                file)

proc main(modelPath: string, simFile: string = "", writeTeX = false,
          schematic = "") =
  loadModelMakeDevice(modelPath)
  echo "loading model success"

  let w1 = model.hidden.weight
  let b1 = model.hidden.bias

  let w2 = model.hidden2.weight
  let b2 = model.hidden2.bias

  let wc = model.classifier.weight
  let bc = model.classifier.bias

  block WriteTensors:
    echo "Biases:"
    echo "\t", b1
    echo "\t", b2
    echo "\t", bc

    echo "Sizes:"
    echo w1
    echo w2
    echo wc

  echo "For input datasets: ", CurrentDsets.sorted
  for i, d in CurrentDsets.sorted:
    echo "Neuron ", i+1, " : ", d

  # and now predict:
  if simFile.len > 0:
    const fname = "/tmp/sim_head30.csv"
    let df = if fileExists(fname):
               readCsv(fname)
             else:
               readSimData(@[simFile])
    if not fileExists(fname):
      df.head(30).writeCsv("/tmp/sim_head30.csv")
    for rowDf in rows(df.head(10)):
      let (pred0, pred1) = model.forward(device, desc, rowDf)
      echo "Prediction: ", (pred0, pred1)
      echo "Manual prediction: ", manualForward(rowDf.toTorchTensor.squeeze.toAmancerVec(),
                                                w1.toAmancer, b1.toAmancerVec,
                                                w2.toAmancer, b2.toAmancerVec,
                                                wc.toAmancer, bc.toAmancerVec)
      ## XXX: Note that there is a discrepancy between our manual multiplication and the libtorch
      ## result. I don't quite understand why, but honestly I don't care right now. The basic gist
      ## is correct.

  if writeTeX:
      let body = writeTeX(w1, b1, w2, b2, wc, bc)
      compile("/tmp/tex_mlp.tex", body, fullBody = true)

  if schematic.len > 0:
    generateMlpSchematic(desc, model, schematic)


when isMainModule:
  import cligen
  dispatch main
