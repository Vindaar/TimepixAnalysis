# Package

version       = "0.3.3"
author        = "Sebastian Schmidt"
description   = "A selection of functions, which help during analysis etc. of InGrid related data files (created by TOS and other)"
license       = "MIT"
skipDirs      = @["out", "data"]
skipExt       = @["h5"]

# Dependencies

requires "nim >= 0.19.0"
requires "loopfusion#head"
requires "arraymancer#head"
requires "https://github.com/vindaar/seqmath >= 0.1.3"
requires "ggplotnim >= 0.2.1"
requires "nimhdf5 >= 0.3.3"
requires "docopt#head"
requires "mpfit >= 0.1.1"
requires "nlopt >= 0.3.1"
requires "plotly >= 0.2.0"
requires "zero_functional#head"
requires "helpers >= 0.2.0"
requires "nimpy#head"
requires "parsetoml"
requires "karax"
requires "https://github.com/yglukhov/threadpools#head"
