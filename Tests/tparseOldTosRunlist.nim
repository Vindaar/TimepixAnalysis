import ingrid / [tos_helpers, ingrid_types]
import strformat, os, strutils
import unittest

suite "Parse old TOS CSV file (2014/15)":
  test "Correct number of runs found":
    const path = currentSourcePath().parentDir / "../resources/Runlist-CAST-D03-W0063.csv"

    let calib = parseOldTosRunlist(path, rtCalibration)
    let back = parseOldTosRunlist(path, rtBackground)
    let xray = parseOldTosRunlist(path, rtXrayFinger)

    echo &"Found {calib.card} calibration runs"
    echo &"Found {back.card} background runs"
    echo &"Found {xray.card} X-ray finger runs"
    check calib.card == 189
    check back.card == 224
    check xray.card == 4

    echo "Calibration runs: ", calib
    echo "Background runs: ", back
    echo "Xray-Finger runs: ", xray
