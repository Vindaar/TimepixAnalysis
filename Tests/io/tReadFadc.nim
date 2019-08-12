import ingrid / [ingrid_types, fadc_helpers]
import sequtils, strutils, os, algorithm
import data/fadcTestData
import unittest
import arraymancer

import ggplotnim
import seqmath

const pwd = currentSourcePath().parentDir
const fname = pwd / "data/data000004.txt-fadc"

import times
let plotSuffix = $getTime().toUnix & ".pdf"

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
  check fadc.bitMode14 == false
  check fadc.data == fadcTestData.dataArray
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
    let pedestalRun = getPedestalRun()
    let fData = fadcFileToFadcData(fadcMem[], pedestalRun)

    check fData.data.shape == @[2560]

    # Now manually calculate the steps for the data conversion and at then
    # end also check whether this is actually what happens
    # according to the CAEN FADC manual p. 15, need to rotate by
    let fadcPedApplied = applyFadcPedestalRun(fadcMem.data, pedestalRun)
    let fadcPedManual = map(
      zip(fadcMem.data, pedestalRun),
      proc(val: (uint16, uint16)): float = float(val[0]) - float(val[1])
    )
    for i in 0 .. fadcPedApplied.high:
      check fadcPedApplied[i] == fadcPedManual[i]

    # extract indices of channel 0
    let idxCh0 = arange(3, 4*2560, 4)
    check idxCh0 == getCh0Indices()
    # set first two elements to 0, since they correspond to our "broken registers"
    var ch0data = fadcPedApplied[idxCh0]
    ch0data[0] = 0
    ch0data[1] = 0

    let nRoll = (fadcMem.trigRec - fadcMem.postTrig) * 20
    let rolled = rotatedLeft(ch0data, nRoll)

    # check the absolute minimum is shifted by exactly `nRoll` indices
    check rolled.argmin + nRoll == ch0data.argmin

    # compare with temporal correction proc
    let tempCorrected = performTemporalCorrection(ch0data, fadcMem.trigRec, fadcMem.postTrig)
    check tempCorrected == rolled

    # finally convert tick values to a voltage
    let convFactor = 1.0 / 2048.0
    let voltage = rolled.mapIt(it * convFactor)

    check voltage == fData.data.toRawSeq

    # compare the calculated spectrum with our expectation for this
    # event
    # NOTE: this comparison does not mean our calculation is bug free, since
    # the sequence we compare to was originally calculated  the same way.
    # It only guarantees that we catch possible errors / changes introduced
    # in the future
    func almostEqual(x, y: float, ep = 1e-5): bool =
      ## checks very roughly if the values match. Anything beyond
      ## 1e-5 should be of no issue for us
      result = (x - y) < ep
    for i, val in fdata.data:
      check almostEqual(val, convertedDataArray[i[0]])

    let df = seqsToDf({ "x" : toSeq(0 ..< 2560),
                        "data": fdata.data.toRawSeq })

    # Comparison has to be done by hand unfortunately
    let path = pwd / "plots/fadc_spectrum"
    ggplot(df, aes(x ~ data)) + geom_line() +
      geom_point(color = color(0.1, 0.1, 0.1, 0.1)) +
      ggsave(path & plotSuffix)

    # NOTE: apparently we cannot do this. Calling the code twice in a row
    # with the exact same input, yields a different result :/
    # check if the plot is exactly the same (high level sanity check)
    # check sameFile(path & plotSuffix, path & ".pdf")
