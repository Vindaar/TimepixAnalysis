import ospaths
import os
import re
import sequtils, future
import strutils
import helper_functions
import threadpool
# import read list of files, to read FADC files in parallel
from tos_helper_functions import readListOfFiles, readMemFilesIntoBuffer
import ingrid_types
import algorithm

import arraymancer

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

proc convertFadcTicksToVoltage*(data: Tensor[float], bit_mode14: bool): Tensor[float] =
  ## equivalent proc to above with Tensor[T] instead of seq[T]
  ## this function converts the channel arrays from FADC ticks to V, by 
  ## making use of the mode_register written to file.
  ## Mode register contains (3 bit register, see CAEN manual p.31):
  ##    bit 0: EN_VME_IRQ interruption tagging of VME bus?!
  ##    bit 1: 14BIT_MODE if set to 1, output uses 14 bit register, instead of 
  ##           backward compatible 12 bit
  ##    bit 2: AUTO_RESTART_ACQ if 1, automatic restart of acqusition at end of 
  ##           RAM readout
  var conversion_factor: float = 1'f64
  if bit_mode14 == true:
    conversion_factor = 1 / 8192'f
  else:
    # should be 2048. instead of 4096 (cf. septemClasses.py)
    conversion_factor = 1 / 2048'f

  # calculate conversion using tensor map, first create mutable copy by assigning
  # to result
  result = data.map(x => x * conversion_factor)
  
proc readFadcFile*(file: seq[string]): ref FadcFile = #seq[float] =
  result = new FadcFile
  var
    # create a sequence with a cap size large enough to hold the whole file
    # speeds up the add, as the sequence does not have to be resized all the
    # time
    data = newSeqOfCap[float](10300)
    posttrig, trigrec, pretrig, n_channels, frequency, sampling_mode: int
    bit_mode14, pedestal_run: bool
    line_spl: seq[string]
  for line in file:
    if likely('#' notin line.string):
      # we add a likely statement, because almost all lines are data lines, hence without '#' 
      data.add(parseFloat(line))
    elif "nb of channels" in line:
      line_spl = line.splitWhitespace
      result.n_channels = parseInt(line_spl[line_spl.high])
    elif "channel mask" in line:
      line_spl = line.splitWhitespace
      result.channel_mask = parseInt(line_spl[line_spl.high])
    elif "postrig" in line or "posttrig" in line:
      line_spl = line.splitWhitespace
      result.posttrig = parseInt(line_spl[line_spl.high])
    elif "pretrig" in line:
      line_spl = line.splitWhitespace
      result.pretrig = parseInt(line_spl[line_spl.high])      
    elif "triggerrecord" in line:
      line_spl = line.splitWhitespace
      result.trigrec  = parseInt(line_spl[line_spl.high])
    elif "frequency" in line:
      line_spl = line.splitWhitespace
      result.frequency = parseInt(line_spl[line_spl.high])
    elif "sampling mode" in line:
      line_spl = line.splitWhitespace
      let mode_register = parseInt(line_spl[line_spl.high])
      # now get bit 1 from mode_register by comparing with 0b010
      result.bit_mode14 = (mode_register and 0b010) == 0b010
      result.sampling_mode = mode_register
    elif "pedestal run" in line:
      line_spl = line.splitWhitespace
      let p_run_flag = parseInt(line_spl[line_spl.high])
      result.pedestal_run = if p_run_flag == 0: false else: true

  # finally assign data sequence                 
  result.data = data

proc readFadcFile*(filename: string): ref FadcFile =
  # wrapper around readFadcFile(file: seq[string]), which first
  # reads all lines in the file before 
  let file = readFile(filename).strip.splitLines
  result = readFadcFile(file)

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

proc applyFadcPedestalRun*[T](fadc_data, pedestal_run: seq[T]): seq[T] = 
  # applys the pedestal run given in the second argument to the first one
  # by zipping the two arrays and using map to subtract each element
  result = map(
    zip(fadc_data, pedestal_run), 
    proc(val: (T, T)): T = val[0] - val[1]
  )

proc fadcFileToFadcData*[T](fadc_file: FadcFile, pedestal_run: seq[T]): FadcData =
  # this function converts an FadcFile object (read from a file) to
  # an FadcData object (extracted Ch0, applied pedestal run, converted
  # to volt)
  result = FadcData()

  # first apply the pedestal run
  # TODO: extend this to apply the closest pedestal run instead?
  var fadc_data = applyFadcPedestalRun(fadc_file.data, pedestal_run)
  
  # and cut out channel 3 (the one we take data with)
  let ch0_indices = arange(3, 4*2560, 4)
  let ch0_vals = fadc_data[ch0_indices]
  result.data = ch0_vals.toTensor
  # set the two 'faulty' registers to 0
  result.data[0] = 0
  result.data[1] = 0

  # convert to volt
  result.data = convertFadcTicksToVoltage(result.data, fadc_file.bit_mode14)

proc getFadcData*(filename: string): FadcData =
  # a convenience function, which performs all steps from reading an FADC
  # file and returns a calibrated FadcData object, only containing the
  # channel we're interested in
  # create the indices with a global pragma, which declares it as equivalent
  # to a static variable in C. It is only initialized once. Saves computation
  # on many subsequent calls.
  let ch0_indices {.global.} = arange(3, 4*2560, 4)
  # same for the pedestal run data
  const home = getHomeDir()  
  const pedestal_run = joinPath(home, "CastData/data/pedestalRuns/pedestalRun000042_1_182143774.txt-fadc")
  let pedestal_d {.global.} = readFadcFile(pedestal_run)
  
  let data = readFadcFile(filename)[]
  result = fadcFileToFadcData(data, pedestal_d.data)

proc getPedestalRun*(): seq[float] =
  # this convenience function returns the data array from
  # our local pedestal run

  const home = getHomeDir()
  const pedestal_file = joinPath(home, "CastData/Code/scripts/data/pedestalRuns/pedestalRun000042_1_182143774.txt-fadc")
  let pedestal = readFadcFile(pedestal_file)
  result = pedestal.data

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

# use experimental pragma to use parallel statement, which is contained in
# dirty template
{.experimental.}  
proc readListOfFadcFiles*(list_of_files: seq[string]): seq[FlowVar[ref FadcFile]] =
  # this procedure receives a list of files, reads them into memory (as a buffer)
  # and processes the content into a seq of ref FadcFile
  # inputs:
  #    list_of_files: seq[string] = a seq of fadc filenames, which are to be read in one go
  # outputs:
  #    seq[FlowVar[ref FadcFile]] = a seq of flow vars pointing to fadc events, since we read
  #                                 in parallel
  # the meat of the proc is in the readListOfFiles function. Here we simply tell it
  # what kind of datatype we are reading.
  result = readListOfFiles[FadcFile](list_of_files)

  
#import gnuplot
#import kissfft/kissfft
#import kissfft/binding
#proc plotFadcFile*(file: string) =
  # a very much work in progress function to plot an FADC file using gnuplot
  # and perform a simple FFT of the signal
#  discard

  # var fadc_file = readFadcFile(file)
  # echo fadc_file.data.len

  # let pedestal_run = getPedestalRun()

  # let fadc_data = fadcFileToFadcData(fadc_file, pedestal_run)
  # echo fadc_data.data.len

  # var
  #   kiss_fft = kissfft.newKissFFT(2560, false)
  #   # f_in: array[2560, binding.kiss_fft_cpx]
  #   # f_out: array[2560, binding.kiss_fft_cpx]
  #   f_in: array[2560, Complex]
  #   f_out: array[2560, Complex]

  # for i in 0..<fadc_data.data.len:
  #   f_in[i] = toComplex(fadc_data.data[i])

  # var r_in: seq[float] = @[]
  # var r_out: seq[float] = @[]  
  
  # let numbers = arange(0, 2560, 1)


  # echo f_out[0]
  # transform(kiss_fft, f_in, f_out)
  # echo f_out[0]

  
  # for i in 0..<f_in.len:
  #   r_in.add(f_in[i].r)
  #   r_out.add(f_out[i].r)

  # plot(numbers, r_in)
  # sleep(1000)
  # plot(numbers, r_out)

#Old code to plot also the peak locations of 
#   peaks in the FADC data
  # var peak_loc: seq[int] = @[]

  # for i in 0..<steps:
  #   let ind = i * lookahead
    
  #   let view = t[ind..(ind + lookahead - 1)]
    
  #   let min_ind = findArgOfLocalMin(view, ind)
    
  #   var min_range = min_ind - int(lookahead / 2)
  #   min_range = if min_range > 0: min_range else: 0
  #   var max_range = min_ind + int(lookahead / 2)
  #   max_range = if max_range < (t.size - 1): max_range else: t.size - 1
    
  #   let min_from_min = findArgOfLocalMin(t[min_range..max_range], min_range)
  #   if min_ind == min_from_min:
  #     peak_loc.add(min_ind)
      
  # #plotFadcFile(name)
  # echo "found peaks at"
  # echo peak_loc
  # echo "Variance of this file is ", variance(t)
  # echo "Mean of this file is ", mean(t)
  # echo "Now dropping all peaks, which are larger than the mean"
  # var peak_vals: seq[float] = @[]
  # var i = 0
  # let cut_value = mean(t) - std(t)
  # echo "Cut value is ", cut_value
  # while i < peak_loc.len:
  #   if t[peak_loc[i]] < cut_value:
  #     peak_vals.add(t[peak_loc[i]])
  #   else:
  #     echo "deleting ", peak_loc[i], " ", i
  #     peak_loc.delete(i)
  #     i -= 1
  #   i += 1

  # var peak_vals: seq[float] = @[]
  # for p in peak_loc:
  #   peak_vals.add(t[p])

  # let numbers = arange(0, 2560, 1)
  # plot(numbers, toRawSeq(t))
  # plot(peak_loc, peak_vals)
  #sleep(3000)


