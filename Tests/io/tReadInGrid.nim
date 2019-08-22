import ingrid / ingrid_types
import ingrid / tos_helpers
import sequtils, strutils, os, algorithm
import data/ingridTestData
import unittest
import arraymancer

import ggplotnim
import seqmath

const pwd = currentSourcePath().parentDir
const fnameOldTOS = pwd / "data/oldTos/data000013_1_141548204.txt"
const fnameSRSTOS = pwd / "data/srsTos/run_006615_data_000041_160502_17-53-31.txt"
const fnameVirtex = pwd / "data/newTos/data000003.txt"
const fnameBackground = pwd / "data/newTos/data022059.txt"

import times
let plotSuffix = $getTime().toUnix & ".pdf"

#proc checkEqual(ingrid: ): bool =
  # checks whther the ingrid event is equal to our expectation
#  return true

func almostEqual(x, y: float, ep = 1e-5): bool =
  ## checks very roughly if the values match. Anything beyond
  ## 1e-5 should be of no issue for us
  result = (x - y) < ep

suite "InGrid data":
  test "Reading old TOS InGrid data":
    let fileContent = readFile(fnameOldTos).strip.splitLines
    let data = concat(@[fnameOldTos], fileContent)

    let ev: OldEvent = processOldEventScanf(data)[]

    proc checkHeader(header: Table[string, string]): bool =
      check header["numChips"] == $1
      check header["shutterMode"] == "verylong"
      check header["shutterTime"] == "13"
      check header["eventNumber"] == "13"
      check header["timestamp"] == "141548204"
      check header["pathName"] == fnameOldTos.parentDir
      check header.len == 6
      return true
    check checkHeader(ev.evHeader)

    check ev.nChips == 1
    # TODO: make sure to calculate length here already?
    check ev.length == 0.0 # undefined for old TOS data at this point
    check ev.chips.len == 1

    let oldTosPix = oldTosPixelStr.strip.splitLines.mapIt((x:  parseInt(it.split[0]).uint8,
                                                           y:  parseInt(it.split[1]).uint8,
                                                           ch: parseInt(it.split[2]).uint16))

    proc checkChip(ch: ChipEvent): bool =
      check ch.chip.name == "D3 W63"
      check ch.chip.number == 0
      check ch.pixels.len == 182
      check ch.pixels == oldTosPix
      return true
    check checkChip(ev.chips[0])

    # TODO: add test for `.dat` file
    let oldTosRunInfo = getOldRunInformation(fnameOldTos.parentDir, 475, rfOldTos)
    check oldTosRunInfo[0] == 144498
    check oldTosRunInfo[1] == 26332
    check oldTosRunInfo[2] == 1444997747
    check oldTosRunInfo[3] == 1445001347
    check oldTosRunInfo[4] == 3600
    # NOTE: parsing of the file works fine, but we can't check whether reading all data from the run
    # actually results in these values, since providing a whole test run for this purpose is
    # too much data. We need a better way to store such a run than using the git repository
    # One possible way to avoid polluting our git repo would be to create a separate test data repo
    # which we clone on travis, link to the test directory and ignore the link in gitignore
    # an alternative is to use git lfs for this purpose.


  test "Reading SRS TOS InGrid data":
    discard
  test "Reading current Virtex TOS InGrid data":
    discard