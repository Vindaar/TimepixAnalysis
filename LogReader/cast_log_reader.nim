## this program reads a set of log files from CAST (slow control and tracking)
## and appends the data to a HDF5 file, which was created using the raw_data_manipulation
## and reconstruction programs

import times
import strutils, strformat
import os
import ospaths
import docopt
import algorithm
import sequtils, future
import strutils
import hashes
import ggplotnim

when not defined(pure):
  import nimhdf5

import ingrid/tos_helpers

{.deadCodeElim: on.}

let doc = """
CAST log file reader. Reads log files (Slow Control or Tracking Logs)
and outputs information about the files.

Note: This tool can be compiled with the `-d:pure` flag to avoid pulling
in the HDF5 dependency.

Usage:
  cast_log_reader --sc <path> [options]
  cast_log_reader --tracking <path> [options]
  cast_log_reader --sc <sc_log> --tracking <tracking_log> [options]
  cast_log_reader <h5file> --delete


Options:
  --sc <sc_log>               Path to slow control log file
  --tracking <tracking_log>   Path to tracking log file
  --h5out <h5file>            Path to a H5 file, which contains InGrid
                              data for the tracking logs currently being
                              read. Used to store run start / end times.
  --delete                    Deletes all tracking related attributes in
                              given H5 file
  -h, --help                  Show this help
  --version                   Show version



Description:
  The user has to supply either the `--sc` flag or the `--tracking` flag
  in addition to a folder or single log file.

  If this tool is run with a single slow control and / or tracking log
  file (i.e. not a folder containing several log files), the tool only
  prints information from that file. In case of a tracking log this is
  the start and end date of a contained tracking (if any). For the slow
  control log currently no useful output is produced.

  If a path to tracking log files is given, and a H5 file is given via
  the `--h5out` option, which  InGrid data for the (some of) the times
  covered by the log files being read. The start and end times of a run
  will be added to the attributes to the TOS runs contained within.

  If only a single H5 file is given together with the `--delete` argument
  all tracking information will be deleted from the attributes in the file.
"""

type
  TrackingKind* = enum
    rkTracking, # if tracking took place in this log
    rkNoTracking # if no tracking took place in this log

  LogFileKind* = enum
    lkSlowControl, lkTracking
  # object storing all slow control log data
  SlowControlLog = object
    # store date of log file as string or Time?
    date: Time
    # time of day, stored as time duration starting from date
    # each log starts at 12am
    times: seq[Duration]
    # pressures relevant for InGrid
    pmm, p3, p3_ba: seq[float]
    # what's called `I_magnet` which is ``not`` the current, but the
    # magnetic field in tesla!
    B_magnet: seq[float]
    # environmental ambient temperature
    T_env: seq[float]
    # pointing positions as angles
    h_angle, v_angle: seq[float]
    # and encoder values
    h_encoder, v_encoder: seq[uint16]
    # MM gas alarm (if on, gas supply closed!)
    mm_gas: seq[bool]
  # object storing the tracking log data
  TrackingLog* = object
    # date the tracking log covers
    date*: Time
    case kind*: TrackingKind
    of rkNoTracking: discard
    of rkTracking:
      tracking_start*, tracking_stop*: Time
    # and the line based data
    timestamps*: seq[int]
    #speedH: seq[float]
    #speedV: seq[float]
    isMoving*: seq[bool]
    isTracking*: seq[bool]
    magB*: seq[float] # magnetic field
                      # apparently magB is not the real value in the tracking logs anymore!

proc newSlowControlLog(): SlowControlLog =
  result.date = fromUnix(0)
  result.times = @[]
  result.pmm = @[]
  result.p3 = @[]
  result.p3_ba = @[]
  result.B_magnet = @[]
  result.T_env = @[]
  result.h_angle = @[]
  result.v_angle = @[]
  result.h_encoder = @[]
  result.v_encoder = @[]
  result.mm_gas = @[]

proc hash(x: TrackingLog): Hash =
  ## proc to define a hash value for a tracking log object
  ## needed to store it in a hash Table
  var h: Hash = 0
  # can use hashing of Time
  h = h !& hash(x.date.toUnix)
  case x.kind
  of rkTracking:
    h = h !& hash(x.tracking_start.toUnix)
    h = h !& hash(x.tracking_stop.toUnix)
  else: discard
  result = !$h

proc `==`(x, y: TrackingLog): bool =
  ## equality of tracking logs, since variant types do not provide an equality
  ## operator by itself
  result = x.date == y.date
  result = result and (x.kind == y.kind)
  if x.kind == y.kind and x.kind == rkTracking:
    # if tracking true, compare start and end times
    result = result and (x.tracking_start == y.tracking_start)
    result = result and (x.tracking_stop == y.tracking_stop)

proc parse_time(time_str: string): Duration =
  ## proc to parse the time from a string hh:mm:ss returned as a
  ## time duration from midnight
  let t = time_str.split(":")
  # create the time interval based on the date given
  result = initDuration(hours = parseInt(t[0]), minutes = parseInt(t[1]), seconds = parseInt(t[2]))

proc is_magnet_moving(h_me, v_me: (int, int)): bool {.inline.} =
  ## determine whether magnet is moving horizontally and (!) vertically
  ## by checking previous and current motor encoder values
  result = if h_me[0] != h_me[1] and v_me[0] != v_me[1]: true else: false


when not defined(pure):
  proc map_log_to_run(logs: seq[TrackingLog], h5file: string): Table[TrackingLog, int] =
    ## maps a sequence of tracking logs to run numbers given in a H5 file
    ## inputs:
    ##   logs: seq[TrackingLog] = tracking logs to map
    ##   h5file: string = path to H5 file containing runs
    ## outputs:
    ##   Table[TrackingLog, int] = table mapping each tracking log to the
    ##     run number in which the tracking is contained
    ## throws:
    ##   HDF5LibraryError = if a call to the H5 library fails

    result = initTable[TrackingLog, int]()

    var h5f = H5open(h5file, "r")

    var run_times = initTable[int, (int, int)]()
    for grp in items(h5f, "/reconstruction", depth = 1):
      if "/reconstruction/run" in grp.name:
        # given a run, check its start and end time
        var tstamp = h5f[(grp.name & "/timestamp").dset_str]
        let
          t_start = tstamp[0, int]
          t_stop  = tstamp[tstamp.high, int]
        # for now need mutable attributes object to access
        var attrs = grp.attrs
        let run_num = attrs["runNumber", int]
        # add start and end time to table of run times
        run_times[run_num] = (t_start, t_stop)

    # now iterate over tracking logs and for each log determine in which run
    # it fits
    for log in logs:
      case log.kind
      of rkTracking:
        # now given start time, look for run for which it fits. Filter elements
        # where run covers tracking start - end
        proc filterRuns(runs: seq[int], log: TrackingLog): int =
          result = -1
          for r in runs:
            let
              t_start = int(log.tracking_start.toUnix)
              t_stop  = int(log.tracking_stop.toUnix)
            if run_times[r][0] < t_start and run_times[r][1] > t_stop:
              doAssert result < 0, "There are multiple trackings? Something is wrong. " & $result & " and " & $r
              result = r
        let run_num = filterRuns(toSeq(run_times.keys), log)
        if run_num > 0:
          result[log] = run_num
      else: discard
    discard h5f.close()

  proc deleteTrackingAttributes(h5file: string) =
    ## proc to delete all tracking related attributes in a H5 file
    withH5(h5file, "rw"):
      let baseName = if "/runs" in h5f: rawDataBase()
                     elif "/reconstruction" in h5f: recoBase()
                     else: raise newException(IOError, "Invalid input file " &
          $h5file & ". It contains neither `runs` nor `reconstruction` group!")
      for grp in items(h5f, baseName, depth = 1):
        if baseName / "run" in grp.name:
          # get number of tracking related attributes
          # for now need mutable attributes object to access
          var attrs = grp.attrs
          var num_tr = 0
          var mgrp = grp
          try:
            num_tr = attrs["num_trackings", int]
          except KeyError:
            # no trackings, continue
            continue
          for i in 0 ..< num_tr:
            var deleted = false
            let
              attr_start = "tracking_start_$#" % $i
              attr_stop  = "tracking_stop_$#" % $i
              attr_num   = "num_trackings"
            deleted = mgrp.deleteAttribute(attr_start)
            deleted = mgrp.deleteAttribute(attr_stop)
            deleted = mgrp.deleteAttribute(attr_num)
            if deleted == false:
              echo "Could not delete one of " &
                "$#, $# or $# in group $#" % [$attr_start, $attr_stop, $attr_num, $mgrp.name]

  proc write_tracking_h5(trck_tab: Table[TrackingLog, int], h5file: string) =
    ## proc to write the mapping of tracking logs to run numbers to the appropriate
    ## groups of the H5 file
    ## inputs:
    ##   trck_tab: Table[TrackingLog, int] = table mapping tracking logs to run numbers
    ##   h5file: string = h5 file in which to write the information
    ## throws:
    ##   HDF5LibraryError = if a call to the H5 library fails


    var runTab = initCountTable[int]()

    withH5(h5file, "rw"):
      for log, run_number in trck_tab:
        # increment run table of current run number
        inc(runTab, runNumber)
        let num = runTab[runNumber]

        let baseName = if "/runs" in h5f: rawDataBase()
                       elif "/reconstruction" in h5f: recoBase()
                       else: raise newException(IOError, "Invalid input file " &
            $h5file & ". It contains neither `runs` nor `reconstruction` group!")
        # take run number, build group name and add attributes
        let runName = baseName & $run_number
        var runGrp = h5f[runName.grp_str]
        # check if there is a number of trackings already

        if "num_trackings" notin runGrp.attrs:
          runGrp.attrs["num_trackings"] = 1
        else:
          runGrp.attrs["num_trackings"] = num

        # store trackings zero indexed
        let numIdx = num - 1

        # add tracking. Trackings will be zero indexed
        runGrp.attrs["tracking_start_$#" % $numIdx] = $log.tracking_start
        runGrp.attrs["tracking_stop_$#" % $numIdx] = $log.tracking_stop

proc sortTrackingLogs(tr_logs: seq[TrackingLog], order = SortOrder.Ascending): seq[TrackingLog] =
  ## proc to sort a sequence of tracking logs by date
  var
    gt = 0
    lt = 0
  # check the sorting order. Depending on case, we define the variables greater and
  # lesser than to either be 1 or 0
  case order
  of Ascending:
    gt = 1
    lt = 0
  of Descending:
    gt = 0
    lt = 1
  result = sorted(tr_logs) do (r, t: TrackingLog) -> int:
    let c = times.`<`(r.date, t.date)
    result = if c == true: lt else: gt

proc read_sc_logfile(filename: string): SlowControlLog =
  ## proc to read a slow control log file
  ## inputs:
  ##    filename: string = the filename of the log file
  ## outputs:
  ##    SlowControlLog = an object storing the data from a slow control
  ##    log file.
  # define the indices for each data column. See the ./sc_log_understand.org file
  const
    date_i = 0
    time_i = 1
    # for these: taking into account offset of 1
    pmm_i = 14
    p3_i = 17
    p3ba_i = 18
    Imag_i = 24 # is actually the magnetic field in ``T``
    tenv_i = 52
    mm_gas_i = 76
    hang_i = 103
    vang_i = 104
    hme_i = 106
    vme_i = 107

  result = newSlowControlLog()
  ## NOTE: instead of parsing manually, we could in principle also use
  ## ggplotnim's `read_csv` to parse it into a df!
  var parsedDate = false
  for line in lines(filename):
    # skip the first line (header)
    if line.startsWith("#"):
      continue
    elif line.strip.startsWith("DATE"):
      continue # skip actual header
    elif line.len == 0:
      # indicates we reached end of file, empty line
      break

    # else add the data to the SlowControlLog object
    let d = line.splitWhitespace
    if not parsedDate and not line.strip.startsWith("DATE"):
      # set date based on first row
      try:
        result.date = toTime(parse(d[date_i], "MM/dd/yyyy"))
        parsedDate = true
        if result.date.utc().year == 1970:
          echo "fuck this file: ", filename
          quit()
      except TimeParseError:
        echo "Failed to parse the following date line: ", d[date_i]
        echo "Full line: ", d
        echo "In file ", filename
        quit()
    if d.len > 107:
      ## NOTE: there are a number of lines in a few files with miss
      ## some colums. E.g.
      ## `SCDV4.00_2017-07-14.daq`
      ## at some points only has 87 instead of 108 columns
      result.times.add parse_time(d[time_i])
      result.pmm.add parseFloat(d[pmm_i])
      result.p3.add parseFloat(d[p3_i])
      result.p3_ba.add parseFloat(d[p3ba_i])
      result.B_magnet.add parseFloat(d[Imag_i])
      result.T_env.add parseFloat(d[tenv_i])
      result.h_angle.add parseFloat(d[hang_i])
      result.v_angle.add parseFloat(d[vang_i])
      result.h_encoder.add uint16(parseInt(d[hme_i]))
      result.v_encoder.add uint16(parseInt(d[vme_i]))
      let mm_gas_b = if d[mm_gas_i] == "0": false else: true
      result.mm_gas.add mm_gas_b

proc parseDateTime(date, time: string): Time =
  result = toTime(parse(date, "MM/dd/yy"))
  result += parseTime(time)

proc read_tracking_logfile*(filename: string): TrackingLog =
  ## reads a tracking log file and returns an object, which stores
  ## the relevant data
  const
    tracking_i = 0
    date_i = 6
    time_i = 7
    h_me = 9
    v_me = 10
    magB_i = 22

  var
    count = 0
    h_me_p = 0
    v_me_p = 0
    tracking_p = false
    tracking_start = fromUnix(0)
    tracking_stop = fromUnix(0)
    # helper bool, needed because some tracking logs start
    # before midnight
    date_set = false
  for line in lines(filename):
    if count < 2:
      # skip header
      inc count
      continue
    if line.len == 0:
      break
    let d = line.splitWhitespace
    if count > 1 and not date_set:
      # parse the date
      let date = toTime(parse(d[date_i], "MM/dd/yy"))
      if date.utc().year == 2003:
        echo "Log ", filename
        quit()

      # check whether this date is at the end of day from the
      # previous log file or already past midnight
      if date.utc().hour > 23:
        continue
      else:
        # set the date of the tracking log
        result.date = date
        date_set = true
    let
      h_me = int(parseFloat(d[h_me]))
      v_me = int(parseFloat(d[v_me]))
      timestamp = parse_time(d[time_i])
      tracking = if int(parseFloat(d[tracking_i])) == 1: true else: false
      magB = parseFloat(d[magB_i])
    # determine magnet movement and set old encoder values
    let move = is_magnet_moving((h_me, h_me_p), (v_me, v_me_p))
    h_me_p = h_me
    v_me_p = v_me
    if not tracking_p and tracking:
      tracking_start = result.date + timestamp
      tracking_p = not tracking_p
    elif tracking_p and not tracking:
      tracking_stop = result.date + timestamp
      tracking_p = not tracking_p

    # append seq data
    result.timestamps.add parseDateTime(d[date_i], d[time_i]).toUnix.int
    result.isMoving.add move
    result.isTracking.add tracking
    result.magB.add magB
    inc count

  # now set the tracking variant object depending on whether tracking took place
  # or not
  if tracking_start == tracking_stop:
    result = TrackingLog(kind: rkNoTracking,
                         date: result.date)
  else:
    result.tracking_start = tracking_start
    result.tracking_stop = tracking_stop

proc split_tracking_logs(logs: seq[TrackingLog]): (seq[TrackingLog], seq[TrackingLog]) =
  ## given a sequence of sorted tracking logs, splits them by `kind`, i.e.
  ## tracking or no tracking
  ## inputs:
  ##   logs: seq[TrackingLog] = seq of mixed TrackingLog kinds to be split by kind
  ## outputs:
  ##   (seq[TrackingLog], seq[TrackingLog]) =
  ##     0: TrackingLog.kind == rkTracking
  ##     1: TrackingLog.kind == rkNoTracking
  ##     both seqs are sorted by date
  # init result
  result[0] = @[]
  result[1] = @[]
  for log in logs:
    case log.kind
    of rkTracking: result[0].add log
    of rkNoTracking: result[1].add log

proc print_tracking_logs(logs: seq[TrackingLog], print_type: TrackingKind, sorted = true) =
  ## proc to pretty print a seq of TrackingLogs using org date format
  ## inputs:
  ##    logs: seq[TrackingLog] = seq of (potentially mixed) tracking logs to be
  ##      printed using org mode date
  ##    print_type: TrackingKind = the kind of logs to be printed. Either tracking
  ##      no tracking
  ##    sorted: bool = True, assumes `logs` is sorted already, if false sort internally
  ##      before printing

  var s_logs = logs
  if sorted == false:
    s_logs = sortTrackingLogs(s_logs)

  for log in s_logs:
    case log.kind
    of rkTracking:
      echo "<$#>    <$#>" % [formatAsOrgDate(log.tracking_start), formatAsOrgDate(log.tracking_stop)]
    of rkNoTracking:
      echo "<$#>" % formatAsOrgDate(log.date)

  case print_type
  of rkTracking:
    echo &"There are {s_logs.len} solar trackings found in the log file directory"
    var trackTime: Duration
    for log in logs:
      if log.kind == rkTracking:
        trackTime += (log.tracking_stop - log.tracking_start)
    echo &"The total time of all trackings: {trackTime.inHours()} h (exact: {tracktime})"
  of rkNoTracking:
    echo &"There are {s_logs.len} runs without solar tracking found in the log file directory"

  # compute total time magnet was on
  when false:
    ## NOTE: this is wrong for modern tracking logs!
    var sumB: int
    for log in logs:
      for i in 1 ..< log.timestamps.len:
        let diff = log.timestamps[i] - log.timestamps[i-1]
        if log.magB[i] > 8.0:
          sumB = sumB + diff
    echo &"Total time the magnet was on (> 1 T): {sumB.float / 60.0} h"

proc print_slow_control_logs(logs: seq[SlowControlLog]) =
  ## proc to pretty print useful information about the Slow Control data
  ## Mainly related to magnet activity.
  # compute total time magnet was on, here we do it based on difference between last and this
  # timestamp. Then add that diff if magnet is on now.
  var sumB: Duration
  var Bvals = newSeq[float]()
  var Tdiffs = newSeq[int]()
  let sortedLogs = logs.sortedByIt(it.date)
  var oldTime = sortedLogs[0].date + sortedLogs[0].times[0]
  for log in sortedLogs:
    for i in 1 ..< log.times.len:
      let curTime = log.date + log.times[i]
      let diff = (curTime - oldTime)
      Tdiffs.add(diff.inSeconds().int)
      if log.B_magnet[i] > 1.0: # and log.B_magnet[i-1] > 8.0:
        sumB = sumB + diff
      elif log.B_magnet[i] > 0.0 and log.B_magnet[i] < 1.0:
        Bvals.add log.B_magnet[i]
      oldTime = curTime

  let df = seqsToDf({"B" : Bvals})
  ggplot(df, aes("B")) + geom_histogram(bins = 100) + ggsave("/tmp/B_field_larger_0.pdf")
  let dfT = seqsToDf({"Tdiff" : Tdiffs})
    .filter(f{`Tdiff` > 0.1 and `Tdiff` < 500.0})
  echo dfT
  ggplot(dfT, aes("Tdiff")) + geom_histogram(bins = 100) +
    scale_y_continuous() + scale_x_continuous() +
    ggsave("/tmp/T_diffs.pdf")
  echo &"Total time the magnet was on (> 1 T): {sumB.inHours()} h"

proc read_tracking_log_folder(log_folder: string): seq[TrackingLog] =
  ## reads all log files from `log_folder` and returns a tuple of sorted `TrackingLog`
  ## objects, one set for logs w/ tracking, others without
  ## inputs:
  ##   log_folder: string = folder in which log files are stored
  ## outputs:
  ##   (seq[TrackingLog], seq[TrackingLog]) =
  ##     0: TrackingLog.kind == rkTracking
  ##     1: TrackingLog.kind == rkNoTracking
  ##     both seqs are sorted by date
  var tracking_logs: seq[TrackingLog] = @[]

  # for each file in log_folder we perform a check of which run corresponds
  # to the date of this tracking log
  for log_p in walkDir(log_folder):
    case log_p.kind
    of pcFile:
      let log = log_p.path
      # check whether actual log file (extension fits)
      let (dir, fn, ext) = splitFile(log)
      if ext == ".log":
        tracking_logs.add read_tracking_logfile(log)
      else:
        # skipping files other than log files
        continue
    else: discard

  result = sortTrackingLogs(tracking_logs)

proc read_sc_log_folder(log_folder: string) =
  ## reads all slow contro log files from `log_folder` and (at the moment)
  ## determines the time the magnet was turned on in total during the time
  ## covered in all the log files in the folder.
  # for each file in log_folder we perform a check of which run corresponds
  # to the date of this tracking log
  var dfDir = newDataFrame()
  echo log_folder
  var scLogs: seq[SlowControlLog]
  for log_p in walkDir(log_folder):
    case log_p.kind
    of pcFile:
      let log = log_p.path
      # check whether actual log file (extension fits)
      let (dir, fn, ext) = splitFile(log)
      echo "File ", log, " and ext ", ext
      if ext == ".daq":
        let sc = read_sc_logfile(log)
        var df = seqsToDf({ "Time / s" : sc.times.mapIt(it.inSeconds),
                            "B / T" : sc.B_magnet })
        df["Date"] = constantColumn(sc.date.toUnix, df.len)
        dfDir.add df
        scLogs.add sc
      else:
        # skipping files other than log files
        echo "Skipping file ", log_p
        continue
    else: discard

  ## NOTE: we assume that each row in the slow control file covers 1 minute,
  ## because normally the slow control software updates once a minute and writes
  ## at that interval. A more thorough student should make sure the time duration
  ## between each line is roughly 1 minute each!
  if dfDir.len > 0:
    dfDir = dfDir.arrange("Date", SortOrder.Ascending)
      .filter(f{c"Date" != 0})

    dfDir = dfDir.mutate(f{"Timestamp / s" ~ `Date` + c"Time / s"})
    ggplot(dfDir, aes("Timestamp / s", "B / T")) +
      geom_point() +
      ggsave("/tmp/B_against_time.png", width = 1920, height = 1080)

    let firstDate = fromUnix(dfDir["Date"][0, int] + dfDir["Time / s"][0, int])
    let lastDate = fromUnix(dfDir["Date"][dfDir.high, int] + dfDir["Time / s"][dfDir.high, int])
    echo firstDate
    echo lastDate
    let nRowsActive = dfDir.filter(f{c"B / T" > 8.0}).len
    echo &"Magnet was turned on for: {nRowsActive.float / 60.0} h between " &
         &"{$firstDate} and {$lastDate}."

  print_slow_control_logs(scLogs)

proc process_log_folder(folder: string, logKind: LogFileKind,  h5file = "") =
  case logKind
  of lkSlowControl:
    read_sc_log_folder(folder)
  of lkTracking:
    let
      tracking_logs = read_tracking_log_folder(folder)
      (trk, notrk) = split_tracking_logs(tracking_logs)
    echo "No tracking days : "
    print_tracking_logs(notrk, rkNoTracking)
    echo "Tracking days :"
    print_tracking_logs(trk, rkTracking)

    when not defined(pure):
      if h5file.len > 0:
        # given the H5 file, create a referential table connecting
        # the trackings with run numbers
        let trackmap = map_log_to_run(trk, h5file)

        # given mapping of tracking logs to run numbers, finally
        # add tracking information to H5 file
        write_tracking_h5(trackmap, h5file)

when isMainModule:
  # parse docopt string and determine
  # correct proc to call
  let args = docopt(doc)
  echo args

  let h5file = if $args["--h5out"] != "nil": $args["--h5out"] else: ""
  let scPath = $args["--sc"]
  let trackingPath = $args["--tracking"]
  if scPath.len > 0 and scPath != "nil":
    # check whether actual log file (extension fits)
    let (dir, fn, ext) = splitFile(scPath)
    if ext == ".daq":
      # single slow control file
      let sc = read_sc_logfile(scPath)
      echo "Parsed log file: ", sc
      echo "Lot's of useless output, huh?"
    else:
      doAssert ext.len == 0, "Invalid slow control with extension " & $ext &
        ". Instead we expect .daq"
      process_log_folder(scPath, lkSlowControl, h5file)
  elif trackingPath.len > 0 and trackingPath != "nil":
    # check whether actual log file (extension fits)
    let (dir, fn, ext) = splitFile(trackingPath)
    if ext == ".log":
      # single tracking file
      let t = read_tracking_logfile(trackingPath)
      echo "Parsed slow control file: ", t
      echo "Lot's of useless output, huh?"
    else:
      doAssert ext.len == 0, "Invalid tracking with extension " & $ext &
        ". Instead we expect .log"
      process_log_folder(trackingPath, lkTracking, h5file)
  elif $args["--delete"] != "nil":
    when not defined(pure):
      deleteTrackingAttributes($args["<h5file>"])
    else:
      discard
