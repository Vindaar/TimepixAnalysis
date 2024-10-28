import std / [tables, times, math]
import ../ingrid_types
import ./hdf5_utils
import arraymancer

proc tpx3EventHeader(tpx3: Tpx3RunConfig): (Table[string, string], Table[string, string]) =
  ## placeholder information for event header
  var e_header = initTable[string, string]()
  var c_header = initTable[string, string]()
  e_header["eventNumber"] = "0"
  # TODO: in this case simply use the "old default constants"
  e_header["shutterMode"] = "stream"
  e_header["shutterTime"] = "stream"
  # only support 1 chip for old storage format anyways
  e_header["numChips"] = $1
  e_header["timestamp"] = "0"
  c_header["numHits"] = "0"
  c_header["chipName"] = chipNameFromTpx3RunConfig(tpx3)
  # subtract 1, because in the past TOS was 1 indexed
  c_header["chipNumber"] = "0"
  e_header["pathName"] = ""
  e_header["dateTime"] = format(timeFromTpx3RunConfig(tpx3), "YYYY-MM-dd'.'HH:mm:ss")
  result = (e_header, c_header)

proc tpx3RunHeader(tpx3: Tpx3RunConfig, runNumber: int): Table[string, string] =
  result = initTable[string, string]()
  result["runNumber"] = $runNumber
  result["runTime"] = "0"
  result["runTimeFrames"] = "0"
  result["numChips"] = "1"
  result["shutterTime"] = "stream"
  result["runMode"] = "0"
  #result["fastClock"] = "0"
  #result["externalTrigger"] = "0"
  result["pathName"] = ""
  result["dateTime"] = format(timeFromTpx3RunConfig(tpx3), "YYYY-MM-dd'.'HH:mm:ss")
  result["shutterMode"] = "stream"

proc computeTpx3RunParameters*(data: seq[Tpx3Data], startIdx, clusterTimeCutoff, runNumber: int,
                               totCut: TotCut,
                               runConfig: Tpx3RunConfig): ProcessedRun =
  ## this procedure walks over the Timepix3 data and returns all data
  ## we can extract from it that fits into the `ProcessedRun`. This means
  ## that `ProcessedRun` is still incomplete after this!
  let (e_header, c_header) = tpx3EventHeader(runConfig)

  var
    clusters = newSeq[Event]()
    lengths = newSeq[float]()
    hits = newSeq[uint16]()
    tots = newSeq[uint16]()
    cluster = ChipEvent(version: Timepix3)
    ev = Event(isValid: true, chips: newSeq[ChipEvent](1), nChips: 1,
               evHeader: e_header)
    lastToa = 0 # int64.high
    clusterStart = 0 # start in ToA cycles
    startT = 0.0
    lastT = 0.0
    occ = zeros[int64]([1, 256, 256])
    eventIdx = startIdx
    numOverflows = 0
    # mutable local copy to assign to output
    totCut = totCut
  const overflow = 2^14

  for i, el in data:
    let toa = el.TOA_Combined.int
    ## 3 Cases decide whether we keep accumulating to this cluster
    ## (0.: this TOA is equal to last TOA)
    ## 1. this ToA value is *smaller* than the last + user defined cutoff (and larger than start)
    ## 2. this ToA value is *larger* than cluster start - user cutoff. (and smaller than last)
    ##    This implies ToA values are not sorted correctly, but the new ToA is still
    ##    within one cutoff of the starting time
    ## 3. this ToA value + 2^14 (counter max vaule) is *smaller* than last ToA + user defined cutoff
    ##    This is for overflows of the counter, such that we keep accumulating
    ##    as long as we're technically within one cutoff
    if not ( toa == lastToa or # 0.
      ( toa <= lastToa + clusterTimeCutoff and
        toa > clusterStart ) or  # 1.
      ( toa >= (clusterStart - clusterTimeCutoff) and
        toa < lastToa )): # 2.
      if false and cluster.pixels.len == 1:
        ## Debug output to help
        echo "looking at index ", i, " toa ", toa, " cluster ", lastToa
        echo "toa <= lastToa + clusterTimeCutoff ",  toa <= lastToa + clusterTimeCutoff
        echo "toa >= (clusterStart - clusterTimeCutoff) ", toa >= (clusterStart - clusterTimeCutoff)
        echo "(overflow + toa) <= (lastToa + clusterTimeCutoff) ", (overflow + toa) <= (lastToa)
      ev.length = lastT - startT
      if cluster.pixels.len > 0:
        ev.chips = @[cluster]
        ev.evHeader["eventNumber"] = $eventIdx
        ev.evHeader["numHits"] = $cluster.pixels.len
        ev.evHeader["timestamp"] = $(el.chunkStartTime.round.int)
        ev.evHeader["dateTime"] = format(fromUnixFloat(el.chunkStartTime.float), "YYYY-MM-dd'.'HH:mm:ss")
        clusters.add ev
        lengths.add ev.length
        hits.add cluster.pixels.len.uint16
        inc eventIdx
      cluster = ChipEvent(version: Timepix3, pixels: newSeqOfCap[Pix](400),
                          toa: newSeqOfCap[uint16](400),
                          toaCombined: newSeqOfCap[uint64](400))
      startT = el.chunk_start_time
      clusterStart = toa
      lastToa = 0 #int64.high
      numOverflows = 0

    # now apply ToT cut. Only if pixel passes ToTCut actually use it
    if el.ToT < totCut.low.uint16:
      inc totCut.rmLow
    elif el.ToT > totCut.high.uint16:
      inc totCut.rmHigh
    else:
      cluster.pixels.add (x: el.x, y: el.y, ch: el.TOT)
      if toa < clusterTimeCutoff and lastToa + clusterTimeCutoff > overflow:
        inc numOverflows
        # means we had an overflow. As it's a 14 bit counter, we have plenty of space
        # to correct the overflows
        # check `lastToa != int64.high` as we set it such when adding new
        echo "Overflow detected : ", el.TOA, " and ", lastToa
        if numOverflows > 3:
          echo "Current data ", cluster.toa
          echo "Maximum number of overflows detected! BAD BAD BAD"
          raise newException(IOError, "Input data has more than 3 overflows in a single cluster! " &
            "Our current data model cannot handle this.")
      ## XXX: note: this currently does not account for the fact that in theory maybe one pixel
      ## is overflowed, but then we receive a pixel that is *not* overflown yet? Order being wrong?
      cluster.toa.add (el.TOA + (numOverflows * overflow).uint16)
      cluster.toaCombined.add el.TOA_Combined
      tots.add el.TOT
      lastToa = toa
      lastT = el.chunk_start_time
      occ[0, el.y.int, el.x.int] += el.TOT.int64
  result.timepix = Timepix3
  result.events = clusters
  result.length = lengths
  result.tots = @[tots]
  result.hits = @[hits]
  result.occupancies = occ
  result.runHeader = tpx3RunHeader(runConfig, runNumber)
  result.totCut = totCut

  result.nChips = 1 ## TODO: allow multiple chips, find out where to best read from input file
  result.chips = @[(name: chipNameFromTpx3RunConfig(runConfig), number: 0)]
  ## XXX: generate a run number?. For now CL argument to `raw_data_manipulation`
  result.runNumber = runNumber
