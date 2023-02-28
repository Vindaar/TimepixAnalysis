# this file contains procs, which deal with the analysis of the
# FADC data stored in a H5 file

import std / [strutils, sequtils, strformat, heapqueue, os, tables, times, typetraits]
import nimhdf5

import arraymancer
import seqmath

import tos_helpers, ingrid_types, fadc_helpers
import helpers/utils


# for moving average statistics
import adix / stat

const doc = """
InGrid FADC analysis tool

Usage:
  fadc_analysis <HDF5file> [options]
  fadc_analysis <HDF5file> --run_number <number> [options]
  fadc_analysis <HDF5file> --noise_analysis <outfile> [options]
  fadc_analysis -h | --help
  fadc_analysis --version

Options:
  --run_number <number>       Only work on this run (currently only supported)
  --noise_analysis <outfile>  If set only perform analysis of FADC noise
  -h --help                   Show this help
  --version                   Show version.
"""

type
  FadcIntTime = enum
    fk50, fk100

proc getFadcIntTime(date: DateTime): FadcIntTime =
  ## given a `date` we return whether this corresponds to
  ## 50 ns or 100 ns FADC integration time
  # 29/11/17 change diff and amp
  # <2017-12-07 Thu 8:00>  int 50 ns -> 100 ns
  # <2017-12-08 Fri 17:50> int 100 ns -> 50 ns
  # <2017-12-15 Fri 10:20> int 50 ns -> 100 ns
  let
    t1 = initDateTime(7, mDec, 2017, 8, 0, 0)
    t2 = initDateTime(8, mDec, 2017, 17, 50, 0)
    t3 = initDateTime(15, mDec, 2017, 10, 20, 0)
  if date < t1:
    result = fk50
  elif date > t1 and date < t2:
    result = fk100
  elif date > t2 and date < t3:
    result = fk50
  elif date > t3:
    result = fk100

proc noiseAnalysis(h5f: H5File): tuple[f50, f100: Table[string, float64]] =
  ## proc which performs the analysis of the noise, i.e. ratio of active and dead
  ## detector. Read the timestamps in the files, gets start and end times of runs
  ## and compares total dead time for each run, depending on type of FADC setting
  ## NOTE: simpler approach to the problem compared to the (commented out )
  ## version above
  const date_syntax = getDateSyntax()

  # table to store (middle) dates of each run and the corresponding percentage
  # the detector was active. For 50 and 100 ns integration time
  var date_tab_50  = initTable[string, float64]()
  var date_tab_100 = initTable[string, float64]()

  for num, group in runs(h5f):
    # for each event we get the start and end time of the tracking
    # calculate the time difference between start and end (to get total
    # time) and sum up the time the detector was active using
    # eventDuration
    let grp = h5f[group.grp_str]
    var attrs = grp.attrs
    echo "Group is ", grp.name
    if "num_trackings" notin attrs:
      # if no tracking at all, skip this run
      continue
    let ntrackings = attrs["num_trackings", int]
    for i in 0 ..< ntrackings:
      let
        tr_start = attrs[&"tracking_start_{i}", string].parse(date_syntax)
        tr_stop  = attrs[&"tracking_stop_{i}", string].parse(date_syntax)
        tr_mean  = tr_start + initDuration(
          seconds = ((tr_stop - tr_start).inSeconds.float / 2).round.int
        )
      let
        durations = h5f[(group / "eventDuration").dset_str][float64]
        # get tracking indices of this tracking number
        tracking_inds = h5f.getTrackingEvents(grp, i, tracking = true)
        # get all durations for events within tracking from the tstamp indices
        durations_track = durations[tracking_inds]
        # sum all event durations within tracking
        duration_live = initDuration(seconds = durations_track.sum.int)
      let
        # calculate shutter open time
        ratio = duration_live.inSeconds.float / (tr_stop - tr_start).inSeconds.float
        fk_type = getFadcIntTime(tr_mean)
      # build the key for the table consisting of run number, FADC int time
      # and the date
      let key = &"{fk_type} \t {num} \t {$tr_mean}"
      # add to output tables, depending on which type it is
      case fk_type
      of fk50:
        date_tab_50[key] = ratio
      of fk100:
        date_tab_100[key] = ratio

      when not defined(release):
        # debugging..
        echo "Number of trackings ", ntrackings
        echo "Start ", tr_start
        echo "Stop  ", tr_stop
        echo "Middle time ", tr_mean
        echo "Tracking corresponds to ", getFadcIntTime(tr_mean)
        echo "Length of run ", tr_stop - tr_start
        echo "Total live time ", duration_live
        echo "Ratio is ", ratio

  if date_tab_50.len > 0:
    echo "Fk50 mean ratio ", toSeq(values(date_tab_50)).mean
  if date_tab_100.len > 0:
    echo "Fk100 mean ratio ", toSeq(values(date_tab_100)).mean
  result.f50 = date_tab_50
  result.f100 = date_tab_100

proc writeNoiseData(tab_tup: tuple[f50, f100: Table[string, float64]], outfile: string) =
  ## proc to write the ratio of noisy data from for both FADC integration
  ## times into a file
  let
    t50 = tab_tup.f50
    t100 = tab_tup.f100
  var outf = open(outfile, fmWrite)
  outf.write("# 50 ns integration time\n")
  outf.write("# fadcIntTime \t run_number \t run_date \t live_ratio\n")
  for k, v in t50:
    outf.write(&"{k} \t {v}\n")
  outf.write("# 100 ns integration time\n")
  for k, v in t100:
    outf.write(&"{k} \t {v}\n")

proc high[T](t: Tensor[T]): int =
  assert t.rank == 1
  result = t.size.int - 1

proc len[T](t: Tensor[T]): int =
  assert t.rank == 1
  result = t.size.int


## XXX: Instead of this, use a simple array that we index with a
## running (`mod` based) index! We just rotate around so that
## `idx - WindowSize` (mod) is always the oldest index!
type
  ## Helper data types to be aware of the current data in the MovingStat window.
  ## We use a queue to pop the correct values once they go out of the window.
  FifoElement = tuple[idx: int, val: float]
  FIFO = object
    idx: int # incrementing counter used for next element as the index
    windowSize: int
    data: HeapQueue[FifoElement]

proc `<`(a, b: FifoElement): bool = a.idx < b.idx
proc initFifo(windowSize: int): FIFO =
  result = FIFO(idx: 0,
                windowSize: windowSize,
                data: initHeapQueue[FifoElement]())

proc push(fifo: var FIFO, x: float) =
  fifo.data.push (fifo.idx, x)
  inc fifo.idx

proc pop(fifo: var FIFO): float =
  let removed = fifo.data.pop()
  result = removed.val

proc popIdx(fifo: var FIFO): int =
  let removed = fifo.data.pop()
  result = removed.idx

proc meanIdx(fifo: var FIFO): int =
  doAssert fifo.data.len > 0, "FIFO is empty! Cannot calculate a mean"
  var res = 0.0
  let num = fifo.data.len
  while fifo.data.len > 0:
    let idx = fifo.popIdx()
    res += idx.float
  result = (res / num.float).round.int

proc findThresholdValue*[T: seq | Tensor; U](data: T, x_min: int,
                                             minThreshold, threshold: U,
                                             left = true, positive = false): (int, int) =
  ## left determines whether we start search left or right
  ## positive sets the range of the data. postiive == false means we consider
  ## dips instead of peaks
  ##
  ## The current implementation is problematic as it is falling for outliers in the data.
  ## We will replace it by a version that calculates a truncated running mean.
  ##
  ## The `minThreshold` is the threshold we need to cross on the _minimum_ of the signal
  ## before we start counting the rise/fall time. `threshold` on the other hand is the
  ## point we need to cross to `stop` counting rise/fall time.
  ##
  ## Returns the register at which we cross the baseline threshold (`threshold`) and the
  ## register at which we cross the signal threshold (our start, `minThreshold`). Difference
  ## between the two is the rise/fall time.
  const WindowSize = 5

  var x = xMin # current FADC register we look at
  var count = 0 # number of registers we've looked at

  var removeQueue = initFifo(WindowSize)
  var stat = initMovingStat[float, int]()
  stat.push data[x]        # start with position of minimum
  removeQueue.push data[x] # in both

  var start = false # indicates whether we are above `minThreshold`
  var xMinThreshold = -1 # where we crossed the minimum threshold, i.e. `start` set to `true`
  while stat.mean() < threshold and count < 2560:
    if left:
      dec x
      if x <= 0:
        x = data.high
    else:
      inc x
      if x >= data.high:
        x = 0
    inc count
    # push current datapoint
    stat.push data[x]
    removeQueue.push data[x]
    if stat.n > WindowSize: # if window is full, time to pop elements
      let toRemove = removeQueue.pop()
      stat.pop toRemove

    if not start and stat.mean() > minThreshold: # start rise time from here!
      start = true
      xMinThreshold = x

  # once we're out, we found our threshold in the data
  if count >= 2560:
    # set the threshold value to the start value, indicating, that
    # it could not be found
    result = (x_min, x_min)
  else:
    # register is the mean index of the queue (we start index from 0,
    # increase for each add) add to the start index mod 2560
    var resIdx: int
    # remove the offset induced by where we start (`xMin`) and where we start
    # rise/fall time calc (`xMinThreshold`)
    let meanIdx = removeQueue.meanIdx() - (abs(xMin - xMinThreshold))
    if left:
      resIdx = xMinThreshold - meanIdx
      if resIdx < 0:
        resIdx += 2560 # push it to the right by one full window
    else:
      resIdx = meanIdx + xMinThreshold
    result = (resIdx mod 2560, xMinThreshold)

proc diffUnderModulo*[T](a, b: T, modulo: int): T {.inline.} =
  ## returns the difference between two values taking into account
  ## modulo a certain value
  let
    d1 = abs(a - b)
    d2 = abs(modulo - abs(a - b))
  result = min(d1, d2)

proc percIdx(q: float, len: int): int = (len.float * q).round.int
proc biasedTruncMean1D(t: Tensor[float], qLow, qHigh: float): float =
  let numH = t.size.int
  let plow = percIdx(qLow, numH)
  let phih = percIdx(qHigh, numH)
  let subSorted = t.sorted
  ## compute the biased truncated mean by slicing sorted data to lower and upper
  ## percentile index
  var red = 0.0
  for j in max(0, plow) ..< min(numH, phih): # loop manually as data is `uint16` to convert
    red += subSorted[j].float
  result = red / (phih - plow).float

proc skewness(t: Tensor[float]): float =
  ## Skewness from a tensor
  var ms = initMovingStat[float, int]()
  for i in 0 ..< t.size.int:
    ms.push t[i]
  result = ms.skewness()

proc calcRiseAndFallTime*(fadc: Tensor[float]): tuple[baseline: float,
                                                      argMinval,
                                                      riseStart,
                                                      fallStop,
                                                      riseTime,
                                                      fallTime: uint16,
                                                      skewness: float,
                                                      noisy: int,
                                                      minVal: float
                                                     ] =
  ## Calculates the baseline, minimum location, start of pulse rise, end of pulse fall,
  ## rise time and fall time for a ``single`` FADC spectrum
  ##
  ## WARNING: `calcRiseAndFallTime` return data *must* be in the order of the fields
  ## in `RecoFadc`!

  ## XXX: Make the parameters adjustable!
  # get location of minimum
  let xMin = argmin(fadc)

  let
    # we demand at least 4 dips, before we can consider an event as noisy
    n_dips = 4
    # the percentile considered for the calculation of the minimum
    min_percentile = 0.99 # this seems better ~ fine
  let noisy = fadc.isFadcFileNoisy(n_dips)
  let minVal = fadc.calcMinOfPulse(min_percentile)

  when false:
    let df = toDf({ "fadc" : fadc,
                    "idx" : toSeq(0 ..< 2560) })
    ggplot(df, aes("idx", "fadc")) + geom_line() + ggsave("/t/test.pdf")

  # determine baseline based on truncated mean. Remove good chunk of lower data to make sure we 'cut out'
  # the typical signal. 2560/0.3 = 768 channels is a reasonable width that covers most signals.
  let baseline = fadc.biasedTruncMean1D(0.30, 0.95)
  # now determine rise and fall times. The thresholds are defined as the point 5% below the
  # baseline (due to fluctuations). In addition we use `PercentileMean` to compute a tighter minimum
  # value and define a threshold of `OffsetToBaseline` from the minimum value which is the starting
  # point from where rise/fall times are measured!
  const PercentileMean = 0.995 # 0.5% = 2560 * 0.005 = 12.8 registers around the minimum for the minimum val
  const OffsetToBaseline = 0.025 # 2.5 % below baseline seems reasonable
  let meanMinVal = calcMinOfPulse(fadc, PercentileMean)
  let offset = abs(OffsetToBaseline * (meanMinVal - baseline)) # relative to the 'amplitude'

  let
    (riseStart, riseStop) = findThresholdValue(fadc, xMin, meanMinVal + offset, baseline - offset)
    (fallStop, fallStart)  = findThresholdValue(fadc, xMin, meanMinVal + offset, baseline - offset, left = false)

  # given this data, we could now in principle already calculate the fall and rise
  # times in nano seconds, but we're going to stick to registers for now
  let
    riseTime = diffUnderModulo(riseStop, riseStart, 2560)
    fallTime = diffUnderModulo(fallStart, fallStop, 2560)

  result = (baseline: baseline, argMinval: xMin.uint16, riseStart: riseStart.uint16, fallStop: fallStop.uint16,
            riseTime: riseTime.uint16, fallTime: fallTime.uint16,
            skewness: fadc.skewness(),
            noisy: noisy,
            minVal: minVal)

proc calcRiseAndFallTimes*(fadc: Tensor[float]): RecoFadc =
  ## Calculates the baseline, minimum location, start of pulse rise, end of pulse fall,
  ## rise time and fall time for ``all`` given FADC spectra
  let nSpectra = fadc.shape[0]
  for field, data in fieldPairs(result):
    data = newTensorUninit[get(genericParams(typeof data), 0)](nSpectra)
  # for i in `||`(0, fadc.high, ""):
  for i in 0 ..< nSpectra:
    ## WARNING: `calcRiseAndFallTime` return fields *must* have the same names as the
    ## fields in `RecoFadc`!
    let tup = calcRiseAndFallTime(fadc[i, _].squeeze)
    for field, data in fieldPairs(result):
      data[i] = getField(tup, field)

proc calcRiseAndFallTimes*(h5f: H5File, run_number: int) =
  ## proc which reads the FADC data from the given file
  ## then performs the calculation of fall and rise time
  ## starting from the index of the minimum

  if fadcDataBasename(runNumber) notin h5f:
    # if no FADC data available, do nothing
    return

  let
    fadc_group = fadcDataBasename(run_number)
  var
    fadc = h5f[fadc_group.dset_str]

  let fadcShape = fadc.shape
  var f_data = newTensorUninit[float](fadcShape)
  # use unsafe `read` to avoid seq->Tensor copy
  fadc.read(f_data.toUnsafeView())
  let nEvents = fadcShape[0]

  # given the reshaped array, we can now compute the
  # fall and rise times
  let t0 = epochTime()
  echo "Start fadc calc"
  let recoFadc = calcRiseAndFallTimes(f_data)
  echo "FADC minima calculations took: ", (epochTime() - t0)
  # now write data back to h5file
  for field, data in fieldPairs(recoFadc):
    type innerType = get(genericParams(typeof data), 0)
    if fadcRecoPath(runNumber) / "minvals" in h5f: # delete old school minvals naming
      discard h5f.delete(fadcRecoPath(runNumber) / "minvals")
    var dset = h5f.create_dataset(fadcRecoPath(runNumber) / field, nEvents, innerType)
    # write the data
    dset.unsafeWrite(data.toUnsafeView(), nEvents)

proc main(h5file: string, runNumber: int = -1,
          noise_analysis = false,
          outfile = "") =
  var h5f = H5open(h5file, "rw")
  if runNumber > 0:
    calcRiseAndFallTimes(h5f, run_number)
  elif noise_analysis:
    let t_tup = noiseAnalysis(h5f)
    if outfile.len > 0:
      writeNoiseData(t_tup, outfile)
  echo "H5 library closed with ", h5f.close()

when isMainModule:
  import cligen
  dispatch main
