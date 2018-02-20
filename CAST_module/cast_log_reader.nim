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

import ingrid/tos_helper_functions

{.deadCodeElim: on.}

let doc = """
CAST log file reader

Usage:
  cast_log_reader --sc <sc_log> [options]
  cast_log_reader --tracking <tracking_log> [options]
  cast_log_reader --sc <sc_log> --tracking <tracking_log> [options]
  cast_log_reader <folder> --h5out <h5file> [options]


Options:
  --sc <sc_log>               Path to slow control log file
  --tracking <tracking_log>   Path to tracking log file
  --h5out <h5file>            Path to a H5 file, which contains InGrid 
                              data for the tracking logs currently being
                              read. Used to store run start / end times.
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
      tracking_start, tracking_end: Time

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

proc parse_time(time_str: string): TimeInterval =
  ## proc to parse the time from a string hh:mm:ss returned as a
  ## time interval from midnight
  let t = time_str.split(":")
  # create the time interval based on the date given
  result = initInterval(hours = parseInt(t[0]), minutes = parseInt(t[1]), seconds = parseInt(t[2]))

proc is_magnet_moving(h_me, v_me: (int, int)): bool =
  ## determine whether magnet is moving horizontally and (!) vertically
  ## by checking previous and current motor encoder values
  result = if h_me[0] != h_me[1] and v_me[0] != v_me[1]: true else: false
  
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
  #     tracking_start, tracking_end: Time
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
    tracking_end = fromSeconds(0)
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
      tracking_end = result.date + timestamp
      tracking_p = false
    inc count

  # now set the tracking variant object depending on whether tracking took place
  # or not
  if tracking_start == tracking_end:
    result.kind = rkNoTracking
  else:
    result.tracking_start = tracking_start
    result.tracking_end = tracking_end


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

proc log_folder() =
  ## proc called, if a folder of tracking logs is given
  let args = docopt(doc)

  let
    log_folder = $args["<folder>"]
    h5file = $args["--h5out"]

  var
    tracking_days: seq[TrackingLog] = @[]
    notracking_days: seq[TrackingLog] = @[]

  # for each file in log_folder we perform a check of which run corresponds
  # to the date of this tracking log
  for log_p in walkDir(log_folder):
    case log_p.kind
    of pcFile:
      let log = log_p.path
      # check whether actual log file (extension fits)
      let (dir, fn, ext) = splitFile(log)
      if ext == ".log":
        let tracking_log = read_tracking_logfile(log)
        case tracking_log.kind
        of rkTracking: tracking_days.add tracking_log
        else: notracking_days.add tracking_log
      else:
        # skipping files other than log files
        continue
    else: discard

  echo "Tracking days :"

  # sort the files by the time of the starting
  let
    s_tracking_days   = sortTrackingLogs(tracking_days)
    s_notracking_days = sortTrackingLogs(notracking_days)
  
  for t in s_tracking_days:
    echo formatAsOrgDate(t.tracking_start)
    
  echo "No tracking days : ", s_notracking_days.len
  #for t in s_notracking_days:
  #  echo formatAsOrgDate(t.date)
                     
  

when isMainModule:
  # parse docopt string and determine
  # correct proc to call
  let args = docopt(doc)
  echo args

  let folder = $args["<folder>"]
  if folder != "nil":
    log_folder()
  else:
    single_file()

  
  
