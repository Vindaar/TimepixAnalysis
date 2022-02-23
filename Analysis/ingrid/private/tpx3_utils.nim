import std / [tables, times]
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

proc tpx3RunHeader(tpx3: Tpx3RunConfig): Table[string, string] =
  result = initTable[string, string]()
  result["runNumber"] = "0"
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

proc computeTpx3RunParameters*(data: seq[Tpx3Data], startIdx, clusterTimeCutoff: int,
                               runConfig: Tpx3RunConfig): ProcessedRun =
  ## this procedure walks over the Timepix3 data and returns all data
  ## we can extract from it that fits into the `ProcessedRun`. This means
  ## that `ProcessedRun` is still incomplete after this!
  let (e_header, c_header) = tpx3EventHeader(runConfig)

  var clusters = newSeq[Event]()
  var lengths = newSeq[float]()
  var hits = newSeq[uint16]()
  var tots = newSeq[uint16]()
  var cluster: ChipEvent
  var ev = Event(isValid: true, chips: newSeq[ChipEvent](1), nChips: 1,
                 evHeader: e_header)
  var clusterTime = 0
  var startT = 0.0
  var lastT = 0.0
  var occ = zeros[int64]([1, 256, 256])
  var eventIdx = 0
  for i, el in data:
    ## XXX: make sure we don't cut off things based on ToA overflowing!
    if el.TOA.int > clusterTime + clusterTimeCutoff:
      ev.length = lastT - startT
      if cluster.pixels.len > 0:
        ev.chips = @[cluster]
        ev.evHeader["eventNumber"] = $eventIdx
        clusters.add ev
        lengths.add ev.length
        hits.add cluster.pixels.len.uint16
        inc eventIdx
      startT = el.chunk_start_time
    cluster.pixels.add (x: el.x, y: el.y, ch: el.TOT)
    tots.add el.TOT
    clusterTime = el.TOA.int
    lastT = el.chunk_start_time
    occ[0, el.y.int, el.x.int] += el.TOT.int64
  result.timepix = Timepix3
  result.events = clusters
  result.length = lengths
  result.tots = @[tots]
  result.hits = @[hits]
  result.occupancies = occ
  result.runHeader = tpx3RunHeader(runConfig)

  result.nChips = 1 ## TODO: allow multiple chips, find out where to best read from input file
  result.chips = @[(name: chipNameFromTpx3RunConfig(runConfig), number: 0)]
  ## XXX: generate a run number?
  result.runNumber = 0
