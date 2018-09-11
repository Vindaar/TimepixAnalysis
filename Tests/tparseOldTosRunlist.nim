import ingrid / [tos_helpers, ingrid_types]
import strformat

when isMainModule:
  const path = "/mnt/Daten/CAST/2014_15/D03-W0063/Runlist-CAST-D03-W0063.csv"

  let calib = parseOldTosRunlist(path, rtCalibration)
  let back = parseOldTosRunlist(path, rtBackground)
  let xray = parseOldTosRunlist(path, rtXrayFinger)

  echo &"Found {calib.card} calibration runs"
  echo &"Found {back.card} background runs"
  echo &"Found {xray.card} X-ray finger runs"

  echo "Calibration runs: ", calib
  echo "Background runs: ", back
  echo "Xray-Finger runs: ", xray
