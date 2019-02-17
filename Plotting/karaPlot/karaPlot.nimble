# Package

version       = "0.1.0"
author        = "Sebastian Schmidt"
description   = "Interactive plotting tool for Timepix HDF5 files"
license       = "MIT"

# Dependencies

requires "nim >= 0.19.9"
requires "jswebsockets"
requires "parsetoml"

import shell

task server, "Build the server":
  exec "nim c -d:release -d:H5_FUTURE --threads:on plotData.nim"

task client, "Build the client":
  exec "nim js plotDataClient.nim"
