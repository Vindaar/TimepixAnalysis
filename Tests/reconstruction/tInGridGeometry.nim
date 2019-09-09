import ingrid / ingrid_types
import ingrid / tos_helpers
import ingridDatabase / databaseDefinitions
import sequtils, strutils, os, algorithm
import data/ingridTestData
import unittest
import arraymancer

import ggplotnim
import seqmath

const pwd = currentSourcePath().parentDir
const dataPwd = pwd / "marlinTestEvents"
const jsonData = dataPwd / "marlinEvents.json"

import times
let plotSuffix = $getTime().toUnix & ".pdf"

#proc checkEqual(ingrid: ): bool =
  # checks whther the ingrid event is equal to our expectation
#  return true

func almostEqual(x, y: float, ep = 1e-5): bool =
  ## checks very roughly if the values match. Anything beyond
  ## 1e-5 should be of no issue for us
  result = (x - y) < ep

suite "InGrid geometry calculations":
  test "Reconstructing single cluster manually":
    # in the first test we read the InGrid data and perform the reconstruction
    # and cluster search manually. This is done for an event that shows up in the
    # in Christoph's root tree to compare where differences come from
    let fileContent = readFile(fnameVbackground).strip.splitLines
    let data = concat(@[fnameVbackground], fileContent)

  test "Comparison of Marlin events with TPA reco":
    ## given 3 events in data/marlinTestEvents
    ## They are *not*
