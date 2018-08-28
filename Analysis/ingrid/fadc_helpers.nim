import ospaths
import os
import re
import sequtils, sugar
import strutils
import helpers/utils
import threadpool
# import read list of files, to read FADC files in parallel
from tos_helpers import readListOfFiles, readMemFilesIntoBuffer
import ingrid_types
import algorithm

import seqmath
import arraymancer

proc walkRunFolderAndGetFadcFiles*(folder: string): seq[string] =
  # walks a run folder and returns a seq of FADC filename strings

  # let lf = system.listFiles(folder)
  result = getListOfFiles(folder, ".*-fadc")

proc convertFadcTicksToVoltage*[T](array: seq[T], bitMode14: bool): seq[T] =
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
  if bitMode14 == true:
    conversion_factor = 1 / 8192'f
  else:
    # should be 2048. instead of 4096 (cf. septemClasses.py)
    conversion_factor = 1 / 2048'f

  # calculate conversion using map and lambda proc macro
  result = map(array, (x: float) -> float => x * conversion_factor)

proc convertFadcTicksToVoltage*(data: Tensor[float], bitMode14: bool): Tensor[float] =
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
  if bitMode14 == true:
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
    data = newSeqOfCap[uint16](10240)
    postTrig, trigRec, preTrig, n_channels, frequency, sampling_mode: int
    bitMode14, pedestalRun: bool
    line_spl: seq[string]
  # line 0 is the filename itself
  let filepath = file[0]

  for line in file[1 .. ^1]:
    if likely('#' notin line):
      # we add a likely statement, because almost all lines are data lines, hence without '#'
      data.add(uint16(parseInt(line)))
    elif "nb of channels" in line:
      line_spl = line.splitWhitespace
      result.nChannels = parseInt(line_spl[line_spl.high])
    elif "channel mask" in line:
      line_spl = line.splitWhitespace
      result.channel_mask = parseInt(line_spl[line_spl.high])
    elif "postrig" in line or "posttrig" in line:
      line_spl = line.splitWhitespace
      result.postTrig = parseInt(line_spl[line_spl.high])
    elif "pretrig" in line:
      line_spl = line.splitWhitespace
      result.preTrig = parseInt(line_spl[line_spl.high])
    elif "triggerrecord" in line:
      line_spl = line.splitWhitespace
      result.trigRec  = parseInt(line_spl[line_spl.high])
    elif "frequency" in line:
      line_spl = line.splitWhitespace
      result.frequency = parseInt(line_spl[line_spl.high])
    elif "sampling mode" in line:
      line_spl = line.splitWhitespace
      let mode_register = parseInt(line_spl[line_spl.high])
      # now get bit 1 from mode_register by comparing with 0b010
      result.bitMode14 = (mode_register and 0b010) == 0b010
      result.sampling_mode = mode_register
    elif "pedestal run" in line:
      line_spl = line.splitWhitespace
      let p_run_flag = parseInt(line_spl[line_spl.high])
      result.pedestalRun = if p_run_flag == 0: false else: true

  let evNumberRegex =  re".*data(\d{4,6})\.txt-fadc"
  var evNumChipNumStr: array[1, string]

  if match(filepath, evNumberRegex, evNumChipNumStr) == true:
    echo (evNumChipNumStr[0].parseInt)
    result.eventNumber = parseInt(evNumChipNumStr[0])
  else:
    echo "error. could not match file"

  # finally assign data sequence
  result.data = data

proc readFadcFile*(filename: string): ref FadcFile =
  # wrapper around readFadcFile(file: seq[string]), which first
  # reads all lines in the file before
  let file = readFile(filename).strip.splitLines
  result = readFadcFile(file)

proc calcMinOfPulse*(ar: Tensor[float], percentile: float): float =
  # calculates the minimum of the input ar (an FADC pulse) based on
  # a minimum percentile of the array
  let
    arg_min = argmin(ar)
    n_elements = int(float(ar.size) * (1'f - percentile) / 2'f)
    ind_min_r = arg_min - n_elements
    ind_min = if ind_min_r > 0: ind_min_r else: 0
    ind_max_r = arg_min + n_elements
    ind_max = if ind_max_r < ar.size: ind_max_r else: ar.size

  result = mean(ar[ind_min..ind_max])

proc calcMinOfPulseAlt*(array: Tensor[float], percentile: float): float =
  # calculates the minimum of the input array (an FADC pulse) based on
  # a minimum percentile of the array
  var filtered_array: seq[float] #Tensor[float]

  # first we're not interested in values close to zero (since we already applied
  # the pedestal)
  # will filter out all elements, which are smaller (negative values in array)
  # then 5% of minimum of whole array
  let `min` = min(array)
  # filtered_array = filter(array,(x: loat) -> bool => x < 0.05 * min)
  filtered_array = filterIt(array.toRawSeq, it < 0.05 * `min`)
  # given resulting array, calculate percentile
  let n_elements = filtered_array.len
  sort(filtered_array, system.cmp[float])
  echo "Filtered array 0 ", filtered_array[^10..^1]

  #echo filtered_array[0], filtered_array[filtered_array.high]
  #echo filtered_array[0..30]
  let ind: int = toInt(float(n_elements) * percentile)
  # echo n_elements
  # echo ind
  let threshold = filtered_array[n_elements - ind]
  #filtered_array = filter(array, (x: T) -> bool => x < threshold)
  filtered_array = filterIt(array.toRawSeq, it < threshold)
  #echo "Filtered array ", filtered_array
  echo "Min of array is ", `min`
  result = mean(filtered_array)

proc applyFadcPedestalRun*[T](fadc_data, pedestalRun: seq[T]): seq[float] =
  # applys the pedestal run given in the second argument to the first one
  # by zipping the two arrays and using map to subtract each element
  result = map(
    zip(fadc_data, pedestalRun),
    proc(val: (T, T)): float = float(val[0]) - float(val[1])
  )

proc getCh0Indices*(): seq[int] {.inline.} =
  # proc which simply returns the channel 0 indices
  # NOTE: this always creates a full sequence by using a loop
  # over the 2560 elements. So save a copy of this, if you need
  # it often!
  result = arange(3, 4*2560, 4)

proc performTemporalCorrection*[T](data: seq[T], trigRec, postTrig: int): seq[T] =
  ## performs the temporal correction of the FADC cyclic register
  ## see CAEN FADC manual p. 15
  ## It is done by rotating the data array according to
  ## .. code-block:
  ##   nRoll = (trigRec - postTrig) * 20
  let nRoll = (trigRec - postTrig) * 20
  # now simply roll
  result = rotatedLeft(data, nRoll)

proc fadcFileToFadcData*[T](fadc_file: FadcFile,
                            pedestalRun: seq[T],
                            ch0_indices: seq[int]): FadcData =
  # this function converts an FadcFile object (read from a file) to
  # an FadcData object (extracted Ch0, applied pedestal run, converted
  # to volt)
  result = FadcData()

  # first apply the pedestal run
  # TODO: extend this to apply the closest pedestal run instead?
  var fadc_data = applyFadcPedestalRun(fadc_file.data, pedestalRun)

  # and cut out channel 3 (the one we take data with)
  var ch0_vals = fadc_data[ch0_indices]
  # set the two 'faulty' registers to 0
  ch0_vals[0] = 0
  ch0_vals[1] = 0

  # now perform temporal correction
  let tempCorrected = performTemporalCorrection(ch0_vals, fadc_file.trigRec, fadc_file.postTrig)

  # assign result as tensor
  result.data = tempCorrected.toTensor.astype(float)

  # convert to volt
  result.data = convertFadcTicksToVoltage(result.data, fadc_file.bitMode14)

proc fadcFileToFadcData*[T](fadc_file: FadcFile, pedestalRun: seq[T]): FadcData =
  # proc which wraps above proc by first creating the indices needed for the
  # calculation
  let ch0_indices = getCh0Indices()
  result = fadcFileToFadcData(fadc_file, pedestalRun, ch0_indices)

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
  const pedestalRun = joinPath(home, "CastData/data/pedestalRuns/pedestalRun000042_1_182143774.txt-fadc")
  let pedestal_d {.global.} = readFadcFile(pedestalRun)

  let data = readFadcFile(filename)[]
  result = fadcFileToFadcData(data, pedestal_d.data)

proc getPedestalRun*(): seq[uint16] =
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




  ###################################################################################
  # The following procs all deal with the calculation of whether a given FADC event #
  # is noisy or not                                                                 #
  # NOTE: the general procs have been moved to helpers/utils module              #
  ###################################################################################

proc isFadcFileNoisy*(fadc: FadcData, n_dips: int): bool =
  # this procedure checks whether a given file name
  # is a noisy FADC file. Determined by the number of dips found in the
  # FADC signal.
  let peak_loc = findPeaks(fadc.data, 150)
  result = if len(peak_loc) >= n_dips: true else: false

proc isFadcFileNoisy*(fadc_data: Tensor[float], n_dips: int): int =
  # this procedure checks whether a given file name
  # is a noisy FADC file. Determined by the number of dips found in the
  # FADC signal.
  let peak_loc = findPeaks(fadc_data, 150)
  result = if len(peak_loc) >= n_dips: 1 else: 0

proc isFadcFileNoisy*(file: string, n_dips: int): bool =
  # overload of proc above, which first reads a file from disk,
  # performs conversion and then checks
  var fadc_file = readFadcFile(file)
  let pedestalRun = getPedestalRun()
  let fadc_data = fadcFileToFadcData(fadc_file, pedestalRun)
  result = isFadcFileNoisy(fadc_data, n_dips)

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
