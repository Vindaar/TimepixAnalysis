# Package

version       = "0.4.1"
author        = "Sebastian Schmidt"
description   = "A selection of functions, which help during analysis etc. of InGrid related data files (created by TOS and other)"
license       = "MIT"
skipDirs      = @["out", "data"]
skipExt       = @["h5"]

# Dependencies

requires "nim >= 1.4.0"
# major dependencies
requires "arraymancer >= 0.7.32"
# UI/UX deps and optional deps
requires "cligen >= 1.7.4"
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
requires "https://github.com/Vindaar/nblosc"

# These dependencies below here are some of the most important.
# They are pinned to specific versions for this git tag to have the exact used versions
# of the most relevant libraries in place for the (printed) thesis results
requires "ggplotnim == 0.7.0"
requires "ginger == 0.6.1"
requires "nimhdf5 == 0.6.1"
requires "unchained >= 0.4.3"
requires "datamancer >= 0.4.0"
requires "xrayAttenuation == 0.4.3"
requires "measuremancer == 0.2.8"
# optimization dependencies
requires "mpfit == 0.2.0"
requires "nlopt == 0.3.2"
# additional packages
requires "https://github.com/vindaar/seqmath == 0.2.2"
requires "latexdsl == 0.2.0"
requires "shell == 0.6.0"
requires "orgtables"
requires "https://github.com/Vindaar/flatBuffers == 0.1.0"


import std / [strutils, sequtils, strformat]
task koch, "Build all binaries in TPA":
  proc compile(bin: string, flags: seq[string]) =
    let f = @flags.mapIt("-d:" & it).join(" ")
    exec &"nim c {f} {bin}"

  let bins = @[
    ("ingrid/parse_raw_tpx3", @["danger", "blosc"]),
    ("ingrid/raw_data_manipulation", @["danger", "blosc"]),
    ("ingrid/reconstruction", @["danger"]),
    ("ingrid/likelihood", @["danger"]),
    ("ingrid/runAnalysisChain", @["release"]),
    ("ingrid/fake_event_generator", @["danger"]),
    ("createAllLikelihoodCombinations", @[""])
  ]
  for (b, f) in bins:
    compile(b, f)
