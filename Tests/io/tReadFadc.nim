import ingrid / [ingrid_types, fadc_helpers, fadc_analysis]
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

func almostEqual(x, y: float, ep = 1e-5): bool =
  ## checks very roughly if the values match. Anything beyond
  ## 1e-5 should be of no issue for us
  result = (x - y) < ep


suite "Fadc data":
  test "Reading raw FADC data":
    # check fails if input not stripped
    let pf = ProtoFile(name: fname, fileData: readFile(fname))
    # test `readFadcFile`
    let fadc = readFadcFile(pf)

    let fadcDirect = readFadcFile(fname)

    # test readFadcFileMem
    let fadcMem = readFadcFileMem(fname)

    check fadc == fadcMem
    check fadc == fadcDirect

    # check equal to expecations
    check fadc.checkEqual
    check fadcMem.checkEqual
    check fadcDirect.checkEqual

  test "Fadc isValid false if file broken":
    # write a broken file to disk and read it. See it fails.
    var data = readFile(fname).strip.splitLines
    # create a broken file
    data = data[0 .. ^100]
    let brokenFname = fname.parentDir / "brokenFile.txt"
    writeFile(brokenFname, data.join("\n"))

    let brokenFadc = readFadcFileMem(brokenFname)
    check not brokenFadc.isValid
    removeFile(brokenFname)

  test "Fadc isValid false if cannot be opened":
    # try to open file that doesn't exist to trigger OSError
    let brokenFadc = readFadcFileMem("I_dont_exist.txt")
    check not brokenFadc.isValid

  test "Fadc isValid false if header is broken":
    # write a broken file to disk and read it. See it fails.
    var data = readFile(fname).strip.splitLines
    # create a broken file
    data[2] = "I don't belong here"
    let brokenFname = fname.parentDir / "brokenFile.txt"
    writeFile(brokenFname, data.join("\n"))

    let brokenFadc = readFadcFileMem(brokenFname)
    check not brokenFadc.isValid
    removeFile(brokenFname)

  test "Checking FADC data transformation":
    let fadcMem = readFadcFileMem(fname)
    let pedestalRun = getPedestalRun()
    let fData = fadcFileToFadcData(fadcMem, pedestalRun)

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
    var ch0data = fadcPedApplied[idxCh0].toSeq1D
    ch0data[0] = 0
    ch0data[1] = 0

    let nRoll = (fadcMem.trigRec - fadcMem.postTrig) * 20
    let rolled = rotatedLeft(ch0data, nRoll)

    # check the absolute minimum is shifted by exactly `nRoll` indices
    check rolled.argmin + nRoll == ch0data.argmin

    # compare with temporal correction proc
    let tempCorrected = performTemporalCorrection(ch0data.toTensor, fadcMem.trigRec, fadcMem.postTrig)
    check tempCorrected == rolled

    # finally convert tick values to a voltage
    let convFactor = 1.0 / 2048.0
    let voltage = rolled.mapIt(it * convFactor)
    echo "NRoll was ", nRoll
    check voltage == fData.data.toRawSeq

    # compare the calculated spectrum with our expectation for this
    # event
    # NOTE: this comparison does not mean our calculation is bug free, since
    # the sequence we compare to was originally calculated  the same way.
    # It only guarantees that we catch possible errors / changes introduced
    # in the future
    for i, val in fdata.data:
      check almostEqual(val, convertedDataArray[i[0]])

    let df = toDf({ "x" : toSeq(0 ..< 2560),
                        "data": fdata.data.toRawSeq })

    # Comparison has to be done by hand unfortunately
    let path = pwd / "plots/fadc_spectrum"
    ggplot(df, aes("x", "data")) + geom_line() +
      geom_point(color = color(0.1, 0.1, 0.1, 0.1)) +
      ggsave(path & plotSuffix)

    # NOTE: apparently we cannot do this. Calling the code twice in a row
    # with the exact same input, yields a different result :/
    # check if the plot is exactly the same (high level sanity check)
    # check sameFile(path & plotSuffix, path & ".pdf")

  test "Checking FADC property calculations":
    # starting from the working calculation of the spectrum reconstruction,
    # we now check all calculations that extract information from the spectrum
    let fadcMem = readFadcFileMem(fname)
    ## XXX: TO MAKE THIS TEST work again, we need a full run of FADC files for the
    ## pedestal calculation!
    let pedestalRun = getPedestalRun()
    let fData = fadcFileToFadcData(fadcMem, pedestalRun)
    let fadc = fData.data

    # - call `calcMinOfPulse`
    # - call `isFadcFileNoisy`
    #   - write tests for `findPeaks`
    # - call `calcRiseAndFallTime`
    #   - write tests for `findThresholdValue`
    let tup = calcRiseAndFallTime(fadc)
    echo tup
    let (baseline, xMin, riseStart, fallStop, riseTime, fallTime, skewness, noisy, minVal) = tup
    # again, the derived values correspond to the (at the time of writing the tests)
    # current results, which look reasonable to me, aside from the issues mentioned in
    # https://github.com/Vindaar/TimepixAnalysis/issues/41
    # Once #41 is fixed this test will be modified.
    # Also check out the plot below, which highlights the calculated values in the plot
    let exp = (baseline: -0.0006229694073016827, argMinval: 2260, riseStart: 2257, fallStop: 2264, riseTime: 2, fallTime: 2, skewness: -7.932674826190749, noisy: 0, minVal: -0.04831949869791666)
    check baseline.almostEqual(exp.baseLine)
    check xMin.int == exp.argMinval
    check riseStart.int == exp.riseStart
    check fallStop.int == exp.fallStop
    check riseTime == xMin - riseStart
    check riseTime.int == exp.riseTime
    check fallTime == fallStop - xMin
    check fallTime.int == exp.fallTime

    # now calculate manually and compare
    let xminMan = argmin(fadc, 0)[0]
    check xmin.int == xminMan
    let baseman = fadc.percentile(p = 50) + max(fadc) * 0.1
    check baseman == baseline

    const PercentileMean = 0.995 # 0.5% = 2560 * 0.005 = 12.8 registers around the minimum for the minimum val
    const OffsetToBaselineTop = 0.1 # 10 % below baseline due to noisy events
    const OffsetToBaselineBottom = 0.025 # 2.5 % below baseline seems reasonable
    let meanMinVal = calcMinOfPulse(fadc, PercentileMean)
    let offsetTop = abs(OffsetToBaselineTop * (meanMinVal - baseline)) # relative to the 'amplitude'
    let offsetBottom = abs(OffsetToBaselineBottom * (meanMinVal - baseline)) # relative to the 'amplitude'
    echo "mean min val: ", meanMinVal
    echo "offset top val: ", offsetTop
    echo "offset bottom val: ", offsetBottom
    #check meanMinVal ==

    let
      (riseSMan, riseStopMan) = findThresholdValue(fadc, xMin.int, meanMinVal + offsetBottom, baseline - offsetTop)
      (fallSMan, fallStartMan)  = findThresholdValue(fadc, xMin.int, meanMinVal + offsetBottom, baseline - offsetTop, left = false)
    check riseSMan == riseStart.int
    check fallSMan == fallStop.int

    let
      riseTMan = diffUnderModulo(xminMan, riseSMan, 2560)
      fallTMan = diffUnderModulo(xminMan, fallSMan, 2560)
    check riseTMan == riseTime.int
    check fallTMan == fallTime.int

    #let baselineY = toSeq(0 ..< 2560).mapIt(baseline.int)
    #let xminlineX = @[xmin, xmin] # one point for x of min, max
    #let xminlineY = linspace(fData.data.min, fData.data.max, xminlineX.len)

    let riseStartX = @[riseStart.int, riseStart.int]
    let fallStopX = @[fallStop.int, fallStop.int]


    # NOTE: regarding this plot. We could also use a combination of `gather` with
    # `dropNulls = true` and `bind_rows`, but this way it's easier by assigning the
    # individual data frames for the geoms
    let by = @[baseline.int, baseline.int]
    let rs = @[riseStart.int, riseStart.int]
    let fs = @[fallStop.int, fallStop.int]
    let rt = @[riseTime.int, riseTime.int]
    let ft = @[fallTime.int, fallTime.int]
    var df = toDf({ "baseline": by,
                    "riseStart" : rs,
                    "fallStop" : fs,
                    "riseTime" : rt,
                    "fallTime" : ft })
    # Comparison has to be done by hand unfortunately
    let path = pwd / "plots/fadc_spectrum_baseline"
    #ggplot(df, aes("x", "data")) + geom_line() +
    #  geom_point(color = color(0.1, 0.1, 0.1, 0.1)) +
    #  geom_line(aes = aes("x", "baseline"),
    #            color = "blue") +
    #  geom_line(data = df.filter(f{Value -> bool: isNull(`xminY`) == false}), aes = aes("xminX", "xminY"),
    #            color = "red") +
    #  geom_line(data = df.filter(f{Value -> bool: isNull(`xminY`) == false}), aes = aes("riseStart", "xminY"),
    #            color = "green") +
    #  geom_line(data = df.filter(f{Value -> bool: isNull(`xminY`) == false}), aes = aes("fallStop", "xminY"),
    #            color = "pink") +
    #  ggsave(path & plotSuffix)
    # TODO: compare the output of these plots! Can either be done the same way as we do
    # in gglotnim repo or via a SHA hash


    proc addAdditionalFields(df: var DataFrame, minVal: float) =
      let baseline = df["baseline", float][0]
      let riseStart = df["riseStart", float][0]
      let riseTime = df["riseTime", float][0]
      let fallStop = df["fallStop", float][0]
      let fallTime = df["fallTime", float][0]
      df["x"] = @[0, 2560]
      df["y"] = @[baseline, minVal]
      df["riseStop"] = @[riseStart + riseTime, riseStart + riseTime]
      df["fallStart"] = @[fallStop - fallTime, fallStop - fallTime]
      df["minEdges"] = @[riseStart + riseTime, fallStop - fallTime]
    let xFadc = toSeq(0 ..< 2560)
    let yFadc = fadc# [0,_].squeeze
    df.addAdditionalFields(yFadc.min)

    let dfFadc = toDf({"x" : xFadc, "y" : yFadc})

    let mLeft = getEnv("L_MARGIN", "2.5").parseFloat
    let mRight = getEnv("R_MARGIN", "1.0").parseFloat
    let mTop = getEnv("T_MARGIN", "2.0").parseFloat
    let mBottom = getEnv("B_MARGIN", "2.0").parseFloat

    echo dfFadc
    echo df

    ggplot(dfFadc, aes("x", "y")) +
      geom_line() +
      #geom_point(color = color(0.1, 0.1, 0.1, 0.1)) +
      #geom_line(data = df, aes = aes("x", "baseline"),
      #                 color = "blue") +
      #geom_line(data = df, aes = aes("argMinval", "y"),
      #                 color = "red") +
      #geom_line(data = df, aes = aes("riseStart", "y"),
      #                 color = "green") +
      #geom_line(data = df, aes = aes("riseStop", "y"),
      #                 color = "green", lineType = ltDashed) +
      #geom_line(data = df, aes = aes("fallStop", "y"),
      #                 color = "pink") +
      #geom_line(data = df, aes = aes("fallStart", "y"),
      #                 color = "pink", lineType = ltDashed) +
      #geom_line(data = df, aes = aes("minEdges", y = "minVal"),
      #                 color = "purple") +
      margin(left = mLeft, top = mTop, right = mRight, bottom = mBottom) +
      xlim(0.0, 2560.0) + # no need to go further!
      ggsave(path & ".pdf")
