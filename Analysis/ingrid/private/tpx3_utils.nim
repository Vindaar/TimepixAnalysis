import std / [tables]
import ../ingrid_types
import arraymancer

proc addTpx3Placeholder(): (Table[string, string], Table[string, string]) =
  ## placeholder information for event header
  var e_header = initTable[string, string]()
  var c_header = initTable[string, string]()
  e_header["eventNumber"] = "0"
  # TODO: in this case simply use the "old default constants"
  e_header["shutterMode"] = "0"
  e_header["shutterTime"] = "0"
  # only support 1 chip for old storage format anyways
  e_header["numChips"] = $1
  e_header["timestamp"] = "0"
  c_header["numHits"] = "0"
  c_header["chipName"] = ""
  # subtract 1, because in the past TOS was 1 indexed
  c_header["chipNumber"] = "0"
  e_header["pathName"] = ""
  e_header["dateTime"] = "2021-05-05.23:52:22"
  result = (e_header, c_header)

proc addTpx3RunPlaceholder(): Table[string, string] =
  result = initTable[string, string]()
  result["runNumber"] = "0"
  result["runTime"] = "0"
  result["runTimeFrames"] = "0"
  result["numChips"] = "1"
  result["shutterTime"] = "0"
  result["runMode"] = "0"
  result["fastClock"] = "0"
  result["externalTrigger"] = "0"
  result["pathName"] = ""
  result["dateTime"] = "2021-05-05.23:52:22"
  result["shutterMode"] = "stream"

proc computeTpx3RunParameters*(data: seq[Tpx3Data], clusterTimeCutoff: int): ProcessedRun =
  ## this procedure walks over the Timepix3 data and returns all data
  ## we can extract from it that fits into the `ProcessedRun`. This means
  ## that `ProcessedRun` is still incomplete after this!
  let (e_header, c_header) = addTpx3Placeholder()

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
  for i, el in data:
    if el.TOA.int > clusterTime + clusterTimeCutoff:
      ev.length = lastT - startT
      if cluster.pixels.len > 0:
        ev.chips = @[cluster]
        clusters.add ev
        lengths.add ev.length
        hits.add cluster.pixels.len.uint16
      cluster = ChipEvent(pixels: newSeqOfCap[Pix](400))
      startT = el.chunk_start_time
    cluster.pixels.add (x: el.x, y: el.y, ch: el.TOT)
    tots.add el.TOT
    clusterTime = el.TOA.int
    lastT = el.chunk_start_time
    occ[0, el.y.int, el.x.int] += el.TOT.int64
  result.events = clusters
  result.length = lengths
  result.tots = @[tots]
  result.hits = @[hits]
  result.occupancies = occ
  result.runHeader = addTpx3RunPlaceholder()
