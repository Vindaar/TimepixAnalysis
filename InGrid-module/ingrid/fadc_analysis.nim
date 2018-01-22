# this file contains procs, which deal with the analysis of the
# FADC data stored in a H5 file

import strutils, sequtils
import seqmath except shape
import docopt
import nimhdf5
import tables
import times
import future

import arraymancer

import tos_helper_functions
import helper_functions

let doc = """
InGrid FADC analysis tool

Usage:
  fadc_analysis <HDF5file> [options]
  fadc_analysis <HDF5file> --run_number <number> [options]
  fadc_analysis -h | --help
  fadc_analysis --version

Options:
  --run_number <number>   Only work on this run (currently only supported)
  -h --help               Show this help
  --version               Show version.
"""

# proc noiseAnalysis(h5f: var H5FileObj) =
#   ## proc which performs the analysis of the noise, i.e. ratio of active and dead
#   ## detector. Read the timestamps in the files, get event durations and in bin
#   ## ranges, which can be determined (~ 20 s to account for 1s accuracy in timestamp?)

#   for num, group in runs(h5f):
    
proc findThresholdValue[T](data: seq[seq[T]], x_min: seq[int], threshold: seq[T], left = true, positive = false): seq[int] =
  # left determines whether we start search left or right
  # positive sets the range of the data. postiive == false means we consider
  # dips instaed of peaks
  result = newSeq[int](x_min.len)
  for i, el in data:
    let thr     = threshold[i]
    var x = x_min[i]
    var count = 0
    while el[x] < thr and count < 2560:
      # positive == true would mean el[x] > thr
      if left == true:
        dec x
        # NOTE: the 2 here represents the fact that the FADC data is 0 for the first entry in the
        # array!
        if x < 2:
          x = el.len
      else:
        inc x
        if x > el.high:
          x = 2
      inc count
    # once we're out, we found our threshold in the data
    if count >= 2560:
      # set the threshold value to the start value, indicating, that
      # it could not be found
      result[i] = x_min[i]
    else:
      # in case of a good result..
      result[i] = x
  

proc reshape[T](s: seq[T], shape = int): seq[seq[T]] =
  # unfinished proc which reshapes the given sequence to a
  # higher order nested sequence. Currently hardcoded for FADC
  let dim2 = s.len
  var tmp: seq[T] = @[]
  echo dim2
  for i in 0 ..< dim2:
    if i == 0:
      result = @[]
    elif i mod 2560 == 0:
      result.add(tmp)
      tmp.setLen(0)
    tmp.add s[i]

proc diffUnderModulo[T](a, b: T, modulo: int): T {.inline.} =
  let
    d1 = abs(a - b)
    d2 = abs(modulo - abs(a - b))
  result = min(d1, d2)

template asType[T](s: seq[T], dtype: typedesc): untyped =
  mapIt(s, dtype(it))

proc calcRiseAndFallTimes*(h5f: var H5FileObj, run_number: int) =
  ## proc which reads the FADC data from the given file
  ## then performs the calculation of fall and rise tim
  ## starting from the index of the minimum

  let
    minvals_group = minvalsBasename(run_number)
    fadc_group = fadcDataBasename(run_number)
  var
    minvals = h5f[minvals_group.dset_str]
    fadc = h5f[fadc_group.dset_str]

  let t0 = epochTime()
  var f_data = fadc[float64].reshape
  let nevents = f_data.len

  # given the reshaped array, we can now compute the
  # fall and rise times

  # get location of minimum
  let x_min = mapIt(f_data, argmin(it.toTensor))

  let t1 = epochTime()
  echo "Getting all minima took $# seconds" % $(t1 - t0)

  # now determine baseline + 10 % of its maximum value
  let baseline = mapIt(f_data, it.percentile(p = 50) + system.max(it) * 0.1)
  let t2 = epochTime()
  echo "Baseline calc took $# seconds" % $(t2 - t1)

  # now determine rise and fall times
  let
    rise_start = findThresholdValue(f_data, x_min, baseline)
    fall_stop = findThresholdValue(f_data, x_min, baseline, left = false)

  # given this data, we could now in principle already calculate the fall and rise
  # times in nano seconds, but we're going to stick to registers for now
  let
    rise_times = mapIt(zip(x_min, rise_start), diffUnderModulo(it[0], it[1], 2560))
    fall_times = mapIt(zip(x_min, fall_stop), diffUnderModulo(it[1], it[0], 2560))

  # now write data back to h5file
  var
    base_dset = h5f.create_dataset(fadcBaselineBasename(run_number), nevents, float)
    x_min_dset = h5f.create_dataset(argMinvalBasename(run_number), nevents, uint16)
    rise_s_dset = h5f.create_dataset(riseStartBasename(run_number), nevents, uint16)
    fall_s_dset = h5f.create_dataset(fallStopBasename(run_number), nevents, uint16)
    rise_t_dset = h5f.create_dataset(riseTimeBasename(run_number), nevents, uint16)
    fall_t_dset = h5f.create_dataset(fallTimeBasename(run_number), nevents, uint16)

  # convert data to uint16
  let
    x_min_u = x_min.asType(uint16)
    rise_start_u = rise_start.asType(uint16)
    fall_stop_u = fall_stop.asType(uint16)
    rise_times_u = rise_times.asType(uint16)
    fall_times_u = fall_times.asType(uint16)

  # write the data
  let all = base_dset.all
  base_dset[all] = baseline
  x_min_dset[all] = x_min_u
  rise_s_dset[all] = rise_start_u
  fall_s_dset[all] = fall_stop_u
  rise_t_dset[all] = rise_times_u
  fall_t_dset[all] = fall_times_u

proc main() =
  let args = docopt(doc)
  echo args

  let
    h5file = $args["<HDF5file>"]
    run_number = $args["--run_number"]


  var h5f = H5File(h5file, "rw")
  calcRiseAndFallTimes(h5f, parseInt(run_number))
  echo "H5 library closed with ", h5f.close()

when isMainModule:
  main()
