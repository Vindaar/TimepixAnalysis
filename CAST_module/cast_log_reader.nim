# this program reads a set of log files from CAST (slow control and tracking)
# and appends the data to a HDF5 file, which was created using the raw_data_manipulation
# and reconstruction programs

import times
import strutils
import os

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

proc is_magnet_moving(h_me, v_me: (int, int)): bool =
  ## determine whether magnet is moving horizontally and (!) vertically
  ## by checking previous and current motor encoder values
  result = if h_me[0] != h_me[1] and v_me[0] != v_me[1]: true else: false
  

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
    
  for line in splitLines(file):
    if count < 2:
      # skip header
      inc count
      continue
    if line.len == 0:
      break
    let d = line.splitWhitespace
    if count == 2:
      # set the date of the tracking log
      result.date = toTime(parse(d[date_i], "MM/dd/yy"))
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

when isMainModule:
  let file_sc = paramStr(1)
  let file_tr = paramStr(2)

  let sc = read_sc_logfile(file_sc)
  
  echo sc

  let t = read_tracking_logfile(file_tr)

  echo t

  
  
