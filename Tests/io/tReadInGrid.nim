import ingrid / ingrid_types
import ingrid / tos_helpers
import ingridDatabase / databaseDefinitions
import sequtils, strutils, os, algorithm
import data/ingridTestData
import unittest
import arraymancer

import seqmath

const pwd = currentSourcePath().parentDir
const fnameOldTOS = pwd / "data/oldTos/data000013_1_141548204.txt"
const fnameSRSTOS = pwd / "data/srsTOS/run_006615_data_000041_160502_17-53-31.txt"
const fnameSrsRunInfoRunMode1 = pwd / "data/srsTOS/runMode1/"
const fnameSrsRunInfoRunMode0 = pwd / "data/srsTOS/runMode0/"
const fnameVCalib = pwd / "data/newTos/data000003.txt"
const fnameVBackground = pwd / "data/newTos/data022059.txt"

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
    let fileContent = readFile(fnameOldTos).strip
    let data = ProtoFile(name: fnameOldTos,
                         fileData: fileContent)

    let ev: OldEvent = processOldEventScanf(data)

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
    let fileContent = readFile(fnameSrsTos).strip
    let data = ProtoFile(name: fnameSrsTos,
                         fileData: fileContent)

    let ev: SrsEvent = processSrsEventScanf(data)

    proc checkHeader(header: Table[string, string]): bool =
      check header["dateTime"] == "2016-05-02T17:53:31+02:00"
      check header["timestamp"] == "1462204411"
      check header["eventNumber"] == "41"
      check header["runNumber"] == "6615"
      check header["numChips"] == "1"
      check header["pathName"] == fnameSrsTos.parentDir
      check header.len == 6
      return true
    check checkHeader(ev.evHeader)

    check ev.nChips == 1
    # TODO: make sure to calculate length here already?
    check ev.length == 0.0 # undefined for old TOS data at this point
    check ev.chips.len == 1

    let srsTosPix = srsTosPixelStr.strip.splitLines.mapIt((x:  parseInt(it.split[0]).uint8,
                                                           y:  parseInt(it.split[1]).uint8,
                                                           ch: parseInt(it.split[2]).uint16))

    proc checkChip(ch: ChipEvent): bool =
      check ch.chip.name == SrsDefaultChipName
      check ch.chip.number == 0
      check ch.pixels.len == 219
      check ch.pixels == srsTosPix
      return true
    check checkChip(ev.chips[0])

  test "Reading SRS TOS run information, run mode == 1":
    let srsTosRunInfo = parseSrsRunInfo(fnameSrsRunInfoRunMode1)
    check srsTosRunInfo["runMode"] == "1"
    check srsTosRunInfo["runTime"] == "0"
    check srsTosRunInfo["shutterMode"] == "long"
    check srsTosRunInfo["shutterTime"] == "2"
    check srsTosRunInfo["runTimeFrames"] == "200000"
    check srsTosRunInfo[SrsNoChipId] == SrsNoChipIdMsg
    check "numChips" notin srsTosRunInfo
    check SrsRunIncomplete notin srsTosRunInfo
    check srsTosRunInfo.len == 6

  test "Reading SRS TOS run information, run mode == 0":
    let srsTosRunInfo = parseSrsRunInfo(fnameSrsRunInfoRunMode0)
    check srsTosRunInfo["runMode"] == "0"
    check srsTosRunInfo["runTime"] == "1000"
    check srsTosRunInfo["shutterMode"] == "long"
    check srsTosRunInfo["shutterTime"] == "2"
    check srsTosRunInfo["runTimeFrames"] == "1"
    check srsTosRunInfo[SrsNoChipId] == SrsNoChipIdMsg
    check "numChips" notin srsTosRunInfo
    check SrsRunIncomplete notin srsTosRunInfo
    check srsTosRunInfo.len == 6

  test "Reading SRS TOS run information w/ chip headers":
    # NOTE: it seems like we don't have a run with the correct FEC headers available right now
    discard

  test "Reading current Virtex TOS InGrid data, calibration example":
    let fileContent = readFile(fnameVcalib).strip
    let data = ProtoFile(name: fnameVcalib,
                         fileData: fileContent)

    let ev: Event = processEventWithScanf(data)

    proc checkHeader(header: Table[string, string]): bool =
      check header["runNumber"] == "302"
      check header["runTime"] == "7200"
      check header["runTimeFrames"] == "0"
      check header["pathName"] == "data/runs/Run_302_181217-14-18" # the real path name
      check header["dateTime"] == "2018-12-17.14:18:09"
      check header["numChips"] == "7"
      check header["shutterTime"] == "30"
      check header["shutterMode"] == "verylong"
      check header["runMode"] == "0"
      check header["fastClock"] == "0"
      check header["externalTrigger"] == "0"
      check header["eventNumber"] == "3"
      check header["useHvFadc"] == "1"
      check header["fadcReadout"] == "1"
      check header["szint1ClockInt"] == "0"
      check header["szint2ClockInt"] == "3427"
      check header["fadcTriggerClock"] == "996321"
      check header["timestamp"] == "1545052689"
      check header.len == 18
      return true
    check checkHeader(ev.evHeader)

    check ev.nChips == 7
    # TODO: make sure to calculate length here already?
    check ev.length == 0.0 # undefined for old TOS data at this point
    check ev.chips.len == 7

    let vcalibTosPixCh3 = virtexCalibPixelStrCh3
      .strip
      .splitLines
      .mapIt((x:  parseInt(it.split[0]).uint8,
              y:  parseInt(it.split[1]).uint8,
              ch: parseInt(it.split[2]).uint16))

    let vcalibTosPixCh6 = virtexCalibPixelStrCh6
      .strip
      .splitLines
      .mapIt((x:  parseInt(it.split[0]).uint8,
              y:  parseInt(it.split[1]).uint8,
              ch: parseInt(it.split[2]).uint16))


    proc checkChip(ch: ChipEvent, chNum: int): bool =
      check ch.chip.name.replace(" ", "") == getSeptemHChip(chNum).replace(" ", "")
      check ch.chip.number == chNum
      if chNum == 3:
        check ch.pixels.len == 257
        check ch.pixels == vcalibTosPixCh3
      elif chNum == 6:
        check ch.pixels.len == 2
        check ch.pixels == vcalibTosPixCh6
      else:
        check ch.pixels.len == 0
      return true
    for i in 0 .. 6:
      check checkChip(ev.chips[i], i)

  test "Reading current Virtex TOS InGrid data, background example":
    let fileContent = readFile(fnameVbackground).strip
    let data = ProtoFile(name: fnameVbackground,
                         fileData: fileContent)

    let ev: Event = processEventWithScanf(data)

    proc checkHeader(header: Table[string, string]): bool =
      check header["runNumber"] == "76"
      check header["runTime"] == "0"
      check header["runTimeFrames"] == "1"
      check header["pathName"] == "data/runs/Run_76_171030-18-39" # the real path name
      check header["dateTime"] == "2017-10-31.09:01:17"
      check header["numChips"] == "7"
      check header["shutterTime"] == "32"
      check header["shutterMode"] == "verylong"
      check header["runMode"] == "0"
      check header["fastClock"] == "0"
      check header["externalTrigger"] == "0"
      check header["eventNumber"] == "22059"
      check header["useHvFadc"] == "1"
      check header["fadcReadout"] == "1"
      check header["szint1ClockInt"] == "4095"
      check header["szint2ClockInt"] == "0"
      check header["fadcTriggerClock"] == "83508713"
      check header["timestamp"] == "1509436877"
      check header.len == 18
      return true
    check checkHeader(ev.evHeader)

    check ev.nChips == 7
    # TODO: make sure to calculate length here already?
    check ev.length == 0.0 # undefined for old TOS data at this point
    check ev.chips.len == 7

    func parsePix(data: (int, string)): (int, seq[(uint8, uint8, uint16)]) =
      result[0] = data[0]
      result[1] = data[1]
      .strip
      .splitLines
      .mapIt((x:  parseInt(it.split[0]).uint8,
              y:  parseInt(it.split[1]).uint8,
              ch: parseInt(it.split[2]).uint16))

    proc checkChip(ch: ChipEvent, chNum: int): bool =
      check ch.chip.name.replace(" ", "") == getSeptemHChip(chNum).replace(" ", "")
      check ch.chip.number == chNum
      let data = parsePix(virtexBackgroundPixelStr[chNum])
      check ch.pixels.len == data[0]
      check ch.pixels == data[1]
      return true
    for i in 0 .. 6:
      check checkChip(ev.chips[i], i)
