import ospaths
import os
import re
import sequtils, future
import strutils
import helper_functions
import algorithm

type
  # object to save FADC data from file into
  FadcFile* = object
    vals*: seq[float]
    posttrig*: int 
    trigrec*: int
    bit_mode14*: bool

  # object to store actual FADC data, which is
  # used (ch0 already extracted)
  FadcData* = object
    data*: seq[float]
    posttrig*: int
    trigrec*: int

proc getListOfFiles*(folder: string, regex = ""): seq[string] = 
  # returns a list of files from folder
  result = @[]
  for file in walkDirRec(folder):
    #if file.match re(regex):
    if match(file, re(regex)):
      result.add(file)

proc walkRunFolderAndGetFadcFiles*(folder: string): seq[string] = 
  # walks a run folder and returns a seq of FADC filename strings
  
  # let lf = system.listFiles(folder)
  result = getListOfFiles(folder, ".*-fadc")

proc convertFadcTicksToVoltage*[T](array: seq[T], bit_mode14: bool): seq[T] = 
  ## this function converts the channel arrays from FADC ticks to V, by 
  ## making use of the mode_register written to file.
  ## Mode register contains (3 bit register, see CAEN manual p.31):
  ##    bit 0: EN_VME_IRQ interruption tagging of VME bus?!
  ##    bit 1: 14BIT_MODE if set to 1, output uses 14 bit register, instead of 
  ##           backward compatible 12 bit
  ##    bit 2: AUTO_RESTART_ACQ if 1, automatic restart of acqusition at end of 
  ##           RAM readout
  result = @[]
  var conversion_factor: float = 1'f64
  if bit_mode14 == true:
    conversion_factor = 1 / 8192'f
  else:
    # should be 2048. instead of 4096 (cf. septemClasses.py)
    conversion_factor = 1 / 2048'f

  # calculate conversion using map and lambda proc macro
  result = map(array, (x: float) -> float => x * conversion_factor)

proc readFadcFile*(file: string): FadcFile = #seq[float] =
  result = FadcFile()
  var 
    vals: seq[float] = @[]
    posttrig, trigrec: int
    bit_mode14: bool
  for line in lines file:
    if "postrig" in line or "posttrig" in line:
      let line_spl = line.splitWhitespace
      posttrig = parseInt(line_spl[line_spl.high])
    elif "triggerrecord" in line:
      let line_spl = line.splitWhitespace
      trigrec  = parseInt(line_spl[line_spl.high])
    elif "sampling mode" in line:
      let line_spl = line.splitWhitespace
      let mode_register: int = parseInt(line_spl[line_spl.high])
      # now get bit 1 from mode_register by comparing with 0b010
      bit_mode14 = (mode_register and 0b010) == 0b010
    elif '#' notin line.string:
      vals.add(parseFloat(line))
  
  result.posttrig = posttrig
  result.trigrec = trigrec
  result.bit_mode14 = bit_mode14
  result.vals = vals

  return result


proc calcMinOfPulse*[T](array: seq[T], percentile: float): T = 
  # calculates the minimum of the input array (an FADC pulse) based on
  # a minimum percentile of the array 
  var filtered_array: seq[T] = @[]

  # first we're not interested in values close to zero (since we already applied 
  # the pedestal)
  # will filter out all elements, which are smaller (negative values in array)
  # then 5% of minimum of whole array
  let min: T = min(array)
  filtered_array = filter(array, (x: T) -> bool => x < 0.05 * min)
  # given resulting array, calculate percentile
  let n_elements = filtered_array.len
  sort(filtered_array, system.cmp[T])
  #echo filtered_array[0], filtered_array[filtered_array.high]
  #echo filtered_array[0..30]
  let ind: int = toInt(float(n_elements) * percentile)
  # echo n_elements
  # echo ind
  let threshold = filtered_array[n_elements - ind]
  filtered_array = filter(array, (x: T) -> bool => x < threshold)

  result = mean(filtered_array)

  
proc applyFadcPedestalRun*[T](fadc_vals, pedestal_run: seq[T]): seq[T] = 
  # applys the pedestal run given in the second argument to the first one
  # by zipping the two arrays and using map to subtract each element
  result = map(
    zip(fadc_vals, pedestal_run), 
    proc(val: (T, T)): T = val[0] - val[1]
  )

  

proc fadcFileToFadcData*[T](fadc_file: FadcFile, pedestal_run: seq[T]): FadcData =
  # this function converts an FadcFile object (read from a file) to
  # an FadcData object (extracted Ch0, applied pedestal run, converted
  # to volt)
  result = FadcData()

  # first apply the pedestal run
  var fadc_vals = applyFadcPedestalRun(fadc_file.vals, pedestal_run)
  
  # and cut out channel 3 (the one we take data with)
  let ch0_indices = arange(3, 4*2560, 4)
  let ch0_vals = getSubarrayByIndices(fadc_vals, ch0_indices)
  result.data = ch0_vals
  # set the two 'faulty' registers to 0
  result.data[0] = 0
  result.data[1] = 0

  # convert to volt
  result.data = convertFadcTicksToVoltage(result.data, fadc_file.bit_mode14)
  

proc getPedestalRun*(): seq[float] =
  # this convenience function returns the data array from
  # our local pedestal run

  const home = getHomeDir()
  const pedestal_file = joinPath(home, "CastData/Code/scripts/data/pedestalRuns/pedestalRun000042_1_182143774.txt-fadc")
  let pedestal = readFadcFile(pedestal_file)
  result = pedestal.vals

proc build_filename_from_event_number(number: string): string =
  # function receives event number as string and builds filename from it
  
  # first pad event number with 0
  let padded_number = align(number, 6, '0')
  # and concat strings
  let filename = join(["data", padded_number, ".txt"])
  return filename

proc buildListOfXrayFiles*(file: string): seq[string] = 
  # function reads the file (a p_y_given_x...) file from a classification 
  # done by a network, creates a list of files for X-ray like events
  # and returns a list of filenames
  var 
    event_list: seq[string] = @[]
  
  for line in lines file:
    if "#" notin line:
      let line_spl = line.split("\t")
      let P_sig: float = parseFloat(line_spl[1].strip())
      if P_sig > 0.5:
        # in this case more likely X-ray than background
        let event_number = line_spl[3].strip()
        let event_name = build_filename_from_event_number(event_number)
        event_list.add(event_name)
        
  return event_list
