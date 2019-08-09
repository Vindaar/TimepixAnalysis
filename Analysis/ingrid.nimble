# Package

version       = "0.3.2"
author        = "Sebastian Schmidt"
description   = "A selection of functions, which help during analysis etc. of InGrid related data files (created by TOS and other)"
license       = "MIT"
skipDirs      = @["out", "data"]
skipExt       = @["h5"]

# Dependencies

requires "nim >= 0.19.0"
requires "loopfusion#head"
requires "arraymancer#head"
requires "https://github.com/vindaar/seqmath#head"
requires "https://github.com/vindaar/ggplotnim#head"
requires "nimhdf5#head"
requires "docopt#head"
requires "mpfit#head"
requires "nlopt#head"
requires "plotly#head"
requires "zero_functional#head"
requires "helpers >= 0.2.0"
requires "nimpy#head"
requires "parsetoml"
requires "karax"
requires "https://github.com/yglukhov/threadpools#head"
