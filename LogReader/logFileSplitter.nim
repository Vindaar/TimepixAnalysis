import ggplotnim, cligen
import strutils, os

import cast_log_reader

proc main(files: seq[string]) =
  ## given a list of input tracking files, splits out one huge CSV
  ## file with the following columns:
  ## timestamp, isTracking, isMoving
  ## and three additional CSV files with the data split by:
  ## 1. magnet at rest
  ## 2. magnet moving but no tracking
  ## 3. tracking
  var df = newDataFrame()
  for f in files:
    let tracking = read_tracking_logfile(f)
    let dfLoc = toDf({ "timestamp" : tracking.timestamps,
                           "isMoving" : tracking.isMoving,
                           "isTracking" : tracking.isTracking })
    df.add dfLoc

  df.writeCsv("/tmp/tracking_logfile_all_data.csv")

  block AtRest:
    let dfLoc = df.filter(f{bool -> bool: `isMoving` == false})
    df.writeCsv("/tmp/tracking_logfile_at_rest.csv")
  block NoTracking:
    let dfLoc = df.filter(f{bool -> bool: `isMoving` and `isTracking` == false})
    df.writeCsv("/tmp/tracking_logfile_is_moving_no_tracking.csv")
  block Tracking:
    let dfLoc = df.filter(f{bool -> bool: `isTracking` == true})
    df.writeCsv("/tmp/tracking_logfile_is_tracking.csv")


when isMainModule:
  dispatch main
