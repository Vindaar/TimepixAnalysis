## this program reads a set of log files from CAST (slow control and tracking)
## and appends the data to a HDF5 file, which was created using the raw_data_manipulation
## and reconstruction programs

import times
import strutils
import os
import ospaths
import docopt
import algorithm
import sequtils, future
import strutils
import hashes

import ingrid/tos_helper_functions

import nimhdf5

{.deadCodeElim: on.}

let doc = """
CAST log file reader

Usage:
  cast_log_reader --sc <sc_log> [options]
  cast_log_reader --tracking <tracking_log> [options]
  cast_log_reader --sc <sc_log> --tracking <tracking_log> [options]
  cast_log_reader <folder> --h5out <h5file> [options]
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
  If this tool is run with a single slow control and / or tracking log
  file (i.e. not a folder containing several log files), the tool only
  prints information from that file. In case of a tracking log this is
  the start and end date of a contained tracking (if any). For the slow
  control log currently no useful output is produced.

  If a path to tracking log files is given, a H5 file is needed, which 
  contains InGrid data for the (some of) the times covered by the log
  files being read. The start and end times of a run will be added to
  the attributes to the TOS runs contained within.
"""

type
  TrackingKind = enum
    rkTracking, # if tracking took place in this log
    rkNoTracking # if no tracking took place in this log
  # object storing all slow control log data
  SlowControlLog = object
    # store date of log file as string or Time?
    date: Time
    # time of day, stored as time interval starting from date
    # each log starts at 12am 
    times: seq[TimeInterval]
    # pressures relevant for InGrid
    pmm, p3, p3_ba: seq[float]
    # magnet current
    # seems like magnet current is not actually part of the sc log file?
    # I_magnet: seq[float]
    # environmental ambient temperature
    T_env: seq[float]
    # pointing positions as angles 
    h_angle, v_angle: seq[float]
    # and encoder values
    h_encoder, v_encoder: seq[uint16]
    # MM gas alarm (if on, gas supply closed!)
    mm_gas: seq[bool]
  # object storing the tracking log data
  TrackingLog = object
    # date the tracking log covers
    date: Time
    case kind: TrackingKind
    of rkNoTracking: discard
    of rkTracking:
      tracking_start, tracking_stop: Time

proc newSlowControlLog(): SlowControlLog =
  result.date = fromSeconds(0)
  result.times = @[]
  result.pmm = @[]
  result.p3 = @[]
  result.p3_ba = @[]
  #result.I_magnet = @[]
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
  h = h !& hash(x.date.toSeconds)
  case x.kind
  of rkTracking:
    h = h !& hash(x.tracking_start.toSeconds) 
    h = h !& hash(x.tracking_stop.toSeconds)
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

proc parse_time(time_str: string): TimeInterval =
  ## proc to parse the time from a string hh:mm:ss returned as a
  ## time interval from midnight
  let t = time_str.split(":")
  # create the time interval based on the date given
  result = initInterval(hours = parseInt(t[0]), minutes = parseInt(t[1]), seconds = parseInt(t[2]))

proc is_magnet_moving(h_me, v_me: (int, int)): bool {.inline.} =
  ## determine whether magnet is moving horizontally and (!) vertically
  ## by checking previous and current motor encoder values
  result = if h_me[0] != h_me[1] and v_me[0] != v_me[1]: true else: false


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

  var h5f = H5File(h5file, "r")

  var run_times = initTable[int, (int, int)]()
  for grp in items(h5f, "/runs", depth = 1):
    if "/runs/run" in grp.name:
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
      let run_num = filter(toSeq(run_times.keys)) do (r: int) -> bool:
        let
          t_start = int(log.tracking_start.toSeconds)
          t_stop  = int(log.tracking_stop.toSeconds)
        if run_times[r][0] < t_start and run_times[r][1] > t_stop: true else: false
      if run_num.len == 1:
        result[log] = run_num[0]
    else: discard
      
  discard h5f.close()

proc deleteTrackingAttributes(h5file: string) =
  ## proc to delete all tracking related attributes in a H5 file
  let groups = @["/runs", "/reconstruction"]
  withH5(h5file, "rw"):
    for base_grp in groups:
      for grp in items(h5f, base_grp, depth = 1):
        if base_grp / "run" in grp.name:
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

  withH5(h5file, "rw"):
    for log, run_number in trck_tab:
      # take run number, build group name and add attributes
      let
        raw_name = rawDataBase() & $run_number
        reco_name = recoBase() & $run_number
      var
        raw_grp = h5f[raw_name.grp_str]
        reco_grp = h5f[reco_name.grp_str]
      # check if there is a number of trackings already
      var num = 0
      try:
        num = raw_grp.attrs["num_trackings", int]
        raw_grp.attrs["num_trackings"] = num + 1
        reco_grp.attrs["num_trackings"] = num + 1
      except KeyError:
        #if run_number == 124:
        #  echo "HOW DID we end up here again???"
        # in this case create key with num == 1
        raw_grp.attrs["num_trackings"] = 1
        reco_grp.attrs["num_trackings"] = 1
      # add tracking. Trackings will be zero indexed
      raw_grp.attrs["tracking_start_$#" % $num] = $log.tracking_start
      raw_grp.attrs["tracking_stop_$#" % $num] = $log.tracking_stop
      reco_grp.attrs["tracking_start_$#" % $num] = $log.tracking_start
      reco_grp.attrs["tracking_stop_$#" % $num] = $log.tracking_stop

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
  
  let file = readFile filename

  # define the indices for each data column. See the ./sc_log_understand.org file
  const
    date_i = 0
    time_i = 1
    # for these: taking into account offset of 1
    pmm_i = 14
    p3_i = 17
    p3ba_i = 18
    #Imag_i = 25
    tenv_i = 52
    mm_gas_i = 76
    hang_i = 103
    vang_i = 104
    hme_i = 106
    vme_i = 107

  result = newSlowControlLog()
    
  var count = 0
  for line in splitLines(file):
    if count == 0:
      # skip the first line (header)
      inc count
      continue
    elif line.len == 0:
      # indicates we reached end of file, empty line
      break
    # else add the data to the SlowControlLog object
    let d = line.splitWhitespace
    if count == 1:
      # set date based on first row
      result.date = toTime(parse(d[date_i], "MM/dd/yyyy"))
    result.times.add parse_time(d[time_i])
    result.pmm.add parseFloat(d[pmm_i])
    result.p3.add parseFloat(d[p3_i])
    result.p3_ba.add parseFloat(d[p3ba_i])
    #result.I_magnet.add parseFloat(d[Imag_i])
    result.T_env.add parseFloat(d[tenv_i])
    result.h_angle.add parseFloat(d[hang_i])
    result.v_angle.add parseFloat(d[vang_i])
    result.h_encoder.add uint16(parseInt(d[hme_i]))
    result.v_encoder.add uint16(parseInt(d[vme_i]))
    let mm_gas_b = if d[mm_gas_i] == "0": false else: true
    result.mm_gas.add mm_gas_b

proc read_tracking_logfile(filename: string): TrackingLog =
  ## reads a tracking log file and returns an object, which stores
  ## the relevant data
  let file = readFile filename

  # TrackingLog = object
  #   case kind: TrackingKind
  #   of rkNoTracking: discard
  #   of rkTracking:
  #     tracking_start, tracking_stop: Time
  const
    tracking_i = 0
    date_i = 6
    time_i = 7
    h_me = 9
    v_me = 10
  
  var
    count = 0
    h_me_p = 0
    v_me_p = 0
    tracking_p = false
    tracking_start = fromSeconds(0)
    tracking_stop = fromSeconds(0)
    # helper bool, needed because some tracking logs start
    # before midnight
    date_set = false
    
  for line in splitLines(file):
    if count < 2:
      # skip header
      inc count
      continue
    if line.len == 0:
      break
    let d = line.splitWhitespace
    if count > 1 and date_set == false:
      # parse the date
      let date = toTime(parse(d[date_i], "MM/dd/yy"))
      # check whether this date is at the end of day from the
      # previous log file or already past midnight
      if toTimeInterval(date).hours > 23:
        continue
      else:
        # set the date of the tracking log
        result.date = date
    let
      h_me = int(parseFloat(d[h_me]))
      v_me = int(parseFloat(d[v_me]))
      timestamp = parse_time(d[time_i])
      tracking = if int(parseFloat(d[tracking_i])) == 1: true else: false
    # determine magnet movement and set old encoder values
    let move = is_magnet_moving((h_me, h_me_p), (v_me, v_me_p))
    h_me_p = h_me
    v_me_p = v_me
    if tracking_p == false and tracking == true:
      tracking_start = result.date + timestamp
      tracking_p = true
    elif tracking_p == true and tracking == false:
      tracking_stop = result.date + timestamp
      tracking_p = false
    inc count

  # now set the tracking variant object depending on whether tracking took place
  # or not
  if tracking_start == tracking_stop:
    result.kind = rkNoTracking
  else:
    result.tracking_start = tracking_start
    result.tracking_stop = tracking_stop


proc single_file() =
  ## proc called, if a single file is being worked on
  let args = docopt(doc)

  let
    file_sc = $args["--sc"]
    file_tr = $args["--tracking"]

  if file_sc != "nil":
    let sc = read_sc_logfile(file_sc)
    echo sc

  if file_tr != "nil":
    let t = read_tracking_logfile(file_tr)
    echo t


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

proc read_log_folder(log_folder: string): seq[TrackingLog] =
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

proc process_log_folder() =

  let args = docopt(doc)

  let
    log_folder = $args["<folder>"]
    h5file = $args["--h5out"]
    tracking_logs = read_log_folder(log_folder)
    (trk, notrk) = split_tracking_logs(tracking_logs)

  echo "Tracking days :"
  print_tracking_logs(trk, rkTracking)
  echo "No tracking days : "
  print_tracking_logs(notrk, rkNoTracking)

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

  let folder = $args["<folder>"]
  if folder != "nil":
    process_log_folder()
  elif $args["--delete"] != "nil":
    deleteTrackingAttributes($args["<h5file>"])
  else:
    single_file()
