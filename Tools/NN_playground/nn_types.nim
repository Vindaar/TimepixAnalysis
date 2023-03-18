## WARNING:
## This file needs to be included and cannot be imported. Otherwise the C++ emitted
## types are not available in the calling scope. :(
## Otherwise:
#[
CC: train_ingrid.nim
/home/basti/.cache/nim/train_ingrid_r/@mtrain_ingrid.nim.cpp:470:25: error: 'MLPImpl' was not declared in this scope
  470 | typedef std::shared_ptr<MLPImpl>  TY__2pAqwDbe669bAzOuOuBbWug;
      |                         ^~~~~~~
compilation terminated due to -Wfatal-errors.
Error: execution of an external compiler program 'g++ -c -std=gnu++17 -funsigned-char  -w -fmax-errors=3 -fpermissive -pthread -I/home/basti/CastData/ExternC
ode/flambeau/vendor/libtorch/include -I/home/basti/CastData/ExternCode/flambeau/vendor/libtorch/include/torch/csrc/api/include -Wfatal-errors -O3 -fno-strict
-aliasing -fno-ident -fno-math-errno   -I/home/basti/src/nim/nim_git_repo/lib -I/home/basti/CastData/ExternCode/TimepixAnalysis/Tools/NN_playground -o /home/
basti/.cache/nim/train_ingrid_r/@mtrain_ingrid.nim.cpp.o /home/basti/.cache/nim/train_ingrid_r/@mtrain_ingrid.nim.cpp' failed with exit code: 1

]#

# We will build the following network:
# Input --> Linear(out_features = 12) --> relu --> Linear(out_features = 1)

import flambeau / [flambeau_nn, tensors]
import ingrid / [tos_helpers, ingrid_types]

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

proc forward*(net: MLP, x: RawTensor): RawTensor =
  #var x = net.hidden2.forward(net.hidden.forward(x).relu()).relu()
  var x = net.hidden.forward(x).relu()
  #x = net.hidden2.forward(x).relu()
  return net.classifier.forward(x).squeeze(1)

## XXX: Defining two models in a single file is currently broken. When trying to use it
## the nim compiler assigns the wrong destructor to the second one (reusing the one
## from the first type). This breaks it. To use it, we currently need to replace
## the logic instead (putting this one first). Once things work we can try to ask Hugo &
## check if submoduling individual models helps.
#defModule:
#  type
#    ConvNet* = object of Module
#      conv1* = Conv2d(1, 50, 15)
#      conv2* = Conv2d(50, 70, 15)
#      conv3* = Conv2d(70, 100, 15)
#      #lin1* = Linear(100 * 15 * 15, 1000)
#      lin1* = Linear(3610, 1000)
#      lin2* = Linear(1000, 50)
#      lin3* = Linear(50, 2)
#
#proc forward(net: ConvNet, x: RawTensor): RawTensor =
#  var x = net.conv1.forward(x).relu().max_pool2d([2, 2])
#  x = net.conv2.forward(x).relu().max_pool2d([2, 2])
#  x = net.conv3.forward(x).relu().max_pool2d([2, 2])
#  x = net.lin1.forward(x).relu()
#  x = net.lin2.forward(x).relu()
#  x = net.lin3.forward(x).relu()
#  return x

type
  ModelKind* = enum
    mkMLP = "MLP"
    mkCNN = "ConvNet"

  ## Placeholder for `ConvNet` above as currently defining both is problematic
  ConvNet* = object

  AnyModel* = MLP | ConvNet

## The batch size we use!
const bsz = 8192 # batch size

proc predictSingle*(model: MLP, input: RawTensor, device: Device): float =
  # predict the output for the single input event
  no_grad_mode:
    # Running input through the network, get the 0th neuron output
    result = model.forward(input.to(device))[_, 0].item(float)

from ./io_helpers import toNimSeq
proc predict*(model: AnyModel,
              input: RawTensor,
              device: Device): seq[float] =
  ## Returns the predictions for all input data contained in `input`
  let dataset_size = input.size(0)
  result = newSeqOfCap[float](dataset_size)
  no_grad_mode:
    for batch_id in 0 ..< (dataset_size.float / bsz.float).ceil.int:
      # minibatch offset in the Tensor
      let offset = batch_id * bsz
      let stop = min(offset + bsz, dataset_size)
      let x = input[offset ..< stop, _ ].to(device)
      # Running input through the network
      let output = model.forward(x)
      result.add output[_, 0].toNimSeq[:float]
