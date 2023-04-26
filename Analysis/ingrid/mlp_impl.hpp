#include "/home/basti/CastData/ExternCode/flambeau/vendor/libtorch/include/torch/csrc/api/include/torch/torch.h"

struct MLPImpl: public torch::nn::Module {
    torch::nn::Linear hidden{nullptr};
    torch::nn::Linear hidden2{nullptr};
    torch::nn::Linear classifier{nullptr};
};
typedef std::shared_ptr<MLPImpl> MLP;
