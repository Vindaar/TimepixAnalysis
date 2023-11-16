import std / [tables, os]
import nimhdf5
import nimhdf5 / serialize_tables
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

if fileExists(CacheTabFile):
  echo "File : ", CacheTabFile, " exists."
  let tab = tryDeserializeH5[CacheTabTyp](CacheTabFile)
  echo "Removing saved file..."
  removeFile(CacheTabFile)
  echo "Writing the file back..."
  tab.tryToH5(CacheTabFile)
  echo "Done."
