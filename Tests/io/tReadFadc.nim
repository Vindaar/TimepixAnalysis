import ingrid / [ingrid_types, fadc_helpers]
import sequtils, strutils, os
import data/fadcTestData
import unittest

const fname = currentSourcePath().parentDir / "data/data000004.txt-fadc"


proc checkEqual(fadc: FadcFile): bool =
  # checks whther the FadcFile is equal to our expectation
  check fadc.nChannels == 0
  check fadc.channelMask == 15
  check fadc.postTrig == 80
  check fadc.preTrig == 15000
  check fadc.trigRec == 65
  check fadc.frequency == 2
  check fadc.samplingMode == 0
  check fadc.pedestalRun == false
  check fadc.eventNumber == 4
  check fadc.data == dataArray
  return true

suite "Fadc data":
  test "Reading raw FADC data":
    # check fails if input not stripped
    expect AssertionError:
      discard readFadcFile(concat(@[fname], readFile(fname).splitLines))
    let fdata = readFile(fname).strip.splitLines

    # test `readFadcFile`
    let namePlusData = concat(@[fname], fdata)
    let fadc = readFadcFile(namePlusData)

    # test readFadcFileMem
    let fadcMem = readFadcFileMem(fname)

    check fadc[] == fadcMem[]

    # check equal to expecations
    check fadc[].checkEqual
    check fadcMem[].checkEqual

  test "Fadc returns nil if file broken":
    # write a broken file to disk and read it. See it fails.
    var data = readFile(fname).strip.splitLines
    # create a broken file
    data = data[0 .. ^100]
    let brokenFname = fname.parentDir / "brokenFile.txt"
    writeFile(brokenFname, data.join("\n"))

    let brokenFadc = readFadcFileMem(brokenFname)
    check brokenFadc == nil
    removeFile(brokenFname)

  test "Fadc returns nil if cannot be opened":
    # try to open file that doesn't exist to trigger OSError
    let brokenFadc = readFadcFileMem("I_dont_exist.txt")
    check brokenFadc == nil

  test "Fadc returns nil if header is broken":
    # write a broken file to disk and read it. See it fails.
    var data = readFile(fname).strip.splitLines
    # create a broken file
    data[2] = "I don't belong here"
    let brokenFname = fname.parentDir / "brokenFile.txt"
    writeFile(brokenFname, data.join("\n"))

    let brokenFadc = readFadcFileMem(brokenFname)
    check brokenFadc == nil
    removeFile(brokenFname)

  test "Checking FADC data transformation":
    let fadcMem = readFadcFileMem(fname)
