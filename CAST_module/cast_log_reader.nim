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
    I_magnet: seq[float]
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
    case kind: TrackingKind
    of rkNoTracking: discard
    of rkTracking:
      tr_stamp_start, tr_stamp_end: int
      tracking_start, tracking_end: Time

proc newSlowControlLog(): SlowControlLog =
  result.date = fromSeconds(0)
  result.times = @[]
  result.pmm = @[]
  result.p3 = @[]
  result.p3_ba = @[]
  result.I_magnet = @[]
  result.T_env = @[]
  result.h_angle = @[]
  result.v_angle = @[]
  result.h_encoder = @[]
  result.v_encoder = @[]
  result.mm_gas = @[]


proc read_sc_logfile(filename: string): SlowControlLog =
  ## proc to read a slow control log file
  ## inputs:
  ##    filename: string = the filename of the log file
  ## outputs:
  ##    SlowControlLog = an object storing the data from a slow control
  ##    log file.
  
  let file = readFile(filename)

  const
    date_i = 0
    time_i = 1
    # for these: taking into account offset of 1
    pmm_i = 15 
    p3_i = 18
    p3ba_i = 19
    Imag_i = 25
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
      inc count
      continue
    let d = line.splitWhitespace
    result.date = toTime(parse(d[date_i], "MM/dd/yyyy"))
    let t = d[time_i].split(":")
    echo t
    result.times.add initInterval(hours = parseInt(t[0]), minutes = parseInt(t[1]), seconds = parseInt(t[2]))
    result.pmm.add parseFloat(d[pmm_i])
    result.p3.add parseFloat(d[p3_i])
    result.p3_ba.add parseFloat(d[p3ba_i])
    result.I_magnet.add parseFloat(d[Imag_i])
    result.T_env.add parseFloat(d[tenv_i])
    result.h_angle.add parseFloat(d[hang_i])
    result.v_angle.add parseFloat(d[vang_i])
    result.h_encoder.add uint16(parseInt(d[hme_i]))
    result.v_encoder.add uint16(parseInt(d[vme_i]))
    let mm_gas_b = if d[mm_gas_i] == "0": false else: true
    result.mm_gas.add mm_gas_b
    echo result
    quit()

when isMainModule:
  let file = paramStr(1)

  let t = read_sc_logfile(file)
  
  
    

  
  
