# Package

version       = "0.3.4"
author        = "Sebastian Schmidt"
description   = "A selection of functions, which help during analysis etc. of InGrid related data files (created by TOS and other)"
license       = "MIT"
skipDirs      = @["out", "data"]
skipExt       = @["h5"]

# Dependencies

requires "nim >= 1.4.0"
# internal module
requires "helpers >= 0.2.0"
# major dependencies
requires "arraymancer >= 0.7.20"
requires "https://github.com/vindaar/seqmath >= 0.1.15"
requires "ingridDatabase"
requires "ggplotnim >= 0.5.9"
requires "nimhdf5 >= 0.5.3"
requires "unchained >= 0.3.8"
requires "xrayAttenuation >= 0.1.2"
requires "datamancer#head"
# optimization dependencies
requires "mpfit >= 0.1.1"
requires "nlopt >= 0.3.1"
# UI/UX deps and optional deps
requires "docopt#head"
requires "cligen"
requires "parsetoml"
requires "karax"
# additional deps (these are on their way out)
requires "plotly >= 0.2.0"
requires "zero_functional#head"
requires "nimpy >= 0.2.0"
requires "https://github.com/yglukhov/threadpools#head"
requires "weave#head"
requires "flambeau#head"
requires "cppstl"
requires "numericalnim"
