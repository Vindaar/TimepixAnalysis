# Package

version       = "0.4.0"
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
requires "arraymancer >= 0.7.23"
requires "ingridDatabase"
# UI/UX deps and optional deps
requires "cligen"
requires "docopt#head" # mostly not in use anymore! Replaced byyy `cligen` pretty much everywhere
requires "adix"
requires "parsetoml"
requires "karax"
# additional deps (these are on their way out)
requires "plotly >= 0.2.0"
requires "zero_functional#head"
requires "nimpy >= 0.2.0"
requires "https://github.com/yglukhov/threadpools#head"
requires "weave#head"
# requires "https://github.com/SciNim/flambeau#head"
requires "cppstl"
requires "numericalnim"
requires "alea"

# These dependencies below here are some of the most important.
# They are pinned to specific versions for this git tag to have the exact used versions
# of the most relevant libraries in place for the (printed) thesis results
requires "ggplotnim == 0.7.0"
requires "ginger == 0.6.1"
requires "nimhdf5 == 0.5.12"
requires "unchained == 0.4.0"
requires "datamancer == 0.4.0"
requires "https://github.com/SciNim/xrayAttenuation == 0.3.1"
requires "measuremancer == 0.2.6"
# optimization dependencies
requires "mpfit == 0.2.0"
requires "nlopt == 0.3.2"
# additional packages
requires "https://github.com/vindaar/seqmath == 0.2.0"
requires "latexdsl == 0.2.0"
requires "shell == 0.5.2"
requires "https://github.com/Vindaar/orgtables"
requires "https://github.com/Vindaar/flatBuffers == 0.1.0"
