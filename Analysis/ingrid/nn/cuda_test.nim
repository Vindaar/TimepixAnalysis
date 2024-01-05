import flambeau/[flambeau_nn]
import flambeau
when not defined(cuda):
  {.error: "Need to compile with `-d:cuda` to have CUDA support!".}
echo Torch.cuda_is_available()
echo Torch.deviceCount()
