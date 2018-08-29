# Package

version       = "0.2.2"
author        = "Sebastian Schmidt"
description   = "A selection of functions, which help during analysis etc. of InGrid related data files (created by TOS and other)"
license       = "MIT"
skipDirs      = @["out", "data"]
skipExt       = @["h5"]

# Dependencies

requires "nim >= 0.18.1"
requires "loopfusion >= 0.0.1"
requires "arraymancer >= 0.4.0"
requires "https://github.com/vindaar/seqmath#head"
requires "nimhdf5#head"
requires "mpfit#head"
requires "nlopt#head"
requires "plotly#head"
requires "zero_functional#head"
requires "helpers >= 0.2.0"
