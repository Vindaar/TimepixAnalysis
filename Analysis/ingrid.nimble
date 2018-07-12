# Package

version       = "0.2.1"
author        = "Sebastian Schmidt"
description   = "A selection of functions, which help during analysis etc. of InGrid related data files (created by TOS and other)"
license       = "MIT"
skipDirs      = @["out", "data"]

# Dependencies

requires "nim >= 0.18.1"
requires "loopfusion >= 0.0.1"
requires "arraymancer >= 0.4.0"
requires "https://github.com/vindaar/seqmath >= 0.1.2"
requires "nimhdf5 >= 0.2.9"
requires "mpfit >= 0.1.0"
requires "nlopt >= 0.3.0"
requires "plotly >= 0.1.0"
requires "zero_functional#head"
requires "helpers >= 0.2.0"

