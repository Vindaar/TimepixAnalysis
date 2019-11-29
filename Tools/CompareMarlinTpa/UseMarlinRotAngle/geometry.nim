import shell, os, macros, sequtils, tables
import ingrid / ingrid_types
import nimhdf5, arraymancer
import locks
## this is the monkey patched version of ingrid / geometry.nim
#var L = Lock()
#L.initLock()

const RunPath = "background-sc-not-only"

template argFilterIt(s, cond: untyped): untyped =
  ## returns indices where filter condition is true
  var idx = 0
  var res: seq[int]
  for it {.inject.} in s:
    if cond:
      res.add idx
    inc idx
  res

var h5fSTUFF {.global.} = H5file("../../../resources/background_splitted.2014+2015.h5", "r")
let grp = h5fStuff[(RunPath).grp_str]
let runs = h5fStuff[grp.name / "RunNumber", float32].mapIt(it.int)
let evs = h5fStuff[grp.name / "EventNumber", float32].mapIt(it.int)

# TODO: for now just read all X, Y, Charge data in global scope. If
# problematic from memory standpoint, move that to `readClusterData` and
# read by index
let specialU8 = special_type(uint8)
#let specialU16 = special_type(uint16)
let xVec = h5fStuff[grp.name / "XCoordinatesVector", specialU8, uint8]
let yVec = h5fStuff[grp.name / "YCoordinatesVector", specialU8, uint8]
let rotAngsMarlin = h5fStuff[grp.name / "RotationAngle", float32]

# Helper type to access correct data from Marlin data
type
  MarlinCluster = object
    data: Cluster[Pix]
    globalIdx: int
    rotAngle: float
  MarlinEvent = object
    eventNumber: int
    clusters: seq[MarlinCluster]
  MarlinRuns = object
    run: int
    events: seq[MarlinEvent]

proc splitBySame(ev: seq[int]): seq[(int, seq[int])] =
  ## splits the given sequence by the indices belonging to the same numbers
  result = newSeqOfCap[(int, seq[int])](ev.len)
  var old = -1
  var tmp = newSeq[int]()
  for i, x in ev:
    if old > 0 and old == x:
      tmp.add i
    else:
      if tmp.len > 0:
        result.add (ev[i - 1], tmp)
      tmp.setLen(0)
      tmp.add i
    old = x

proc readClusterData(idx: int): MarlinCluster =
  ## reads the cluster information for the cluster found at the
  ## global index ``idx``.
  ## For now we're only interested in the rotation angle
  ## NOTE: we have to use a fugly "hack" to get the real ToT values for the
  ## Marlin clusters because Marlin is such a stinker that it only provides the
  ## the actual charge values, for fucks sake, sigh. We either:
  ## - build a tensor for the given full event based on what we read from
  ##   the raw file
  ## - or we walk the different clusters and then search for each correct pixel
  ##   in the full pixel sequence we read before.
  ## Building the tensor is relatively
  ## slow, but access is near instant, searching requires no work before hand
  ## but is itself very slow
  ## The main problem is these things cannot be done here, since this code
  ## is run at global scope. That means we have to fill in the ToT values
  ## at a later time when we actually work on the data in TPA
  let x = xVec[idx]
  let y = yVec[idx]
  doAssert x.len == y.len
  result.data = newSeq[Pix](x.len)
  for i in 0 ..< x.len:
    # NOTE: we leave the charge empty! ToT values will be taken from raw files
    result.data[i] = (x: x[i], y: y[i], ch: 0'u16)

proc buildMarlinEvents(run: int, clusterIdx: seq[int]): seq[MarlinEvent] =
  ## reads all required information for the given event indices for run
  ## `run`
  # start by splitting the clusters by events
  let eventIdx = splitBySame(clusterIdx)
  for (eventNumber, indices) in eventIdx:
    # for each event now walk all clusters and read data
    var mEvent = MarlinEvent(eventNumber: eventNumber)
    for cluster in indices:
      mEvent.clusters.add readClusterData(cluster)
    result.add mEvent

proc mapRunsToEvents(runs, evs: seq[int]): Table[int, MarlinRuns] =
  # first split all runs into a bunch of seq[seq[int]] belonging
  # to the same run number
  let runSeqs = splitBySame(runs)
  # iterate all those runs and create a MarlinRuns object for each
  for (run, events) in runSeqs:
    # build the MarlinEvents
    let mEvents = buildMarlinEvents(run, events)
    result[run] =  MarlinRuns(run: run, events: mEvents)

let mappedRuns = mapRunsToEvents(runs, evs)

template findIt(s: typed, cond: untyped): untyped =
  ## finds the first appearence of `cond` in `s`, returns the element
  type retType = type(s[0])
  var res: retType
  for i, it {.inject.} in s:
    if cond:
      res = it
      break
  res

proc readMarlinEvent(evNum, run: int): MarlinEvent {.inline.} =
  {.gcsafe.}:
    let mRun = mappedRuns[run]
    result = mRun.events.findIt(it.eventNumber == evNum)

#proc injectMarlinRotAngles(): float =
template injectMarlinClusters() =
  # instead of using our normal cluster finder we read the clusters found by the
  # Marlin algorithm for event number `evNum` instead
  result = new RecoEvent[T]
  static: echo type(dat)
  result.event_number = dat.eventNumber
  result.chip_number = chip
  if dat[0].len > 0:
    {.gcsafe.}:
      echo "Event ", result.eventNumber
      let event = readMarlinEvent(result.eventNumber, run)
      # now build the tensor of the raw event
      let rawEvent = block:
        var tmp = newTensor[uint16](256, 256)
        for (x, y, ch) in dat[0]:
          tmp[y.int, x.int] = ch
        tmp
      result.cluster = newSeq[ClusterObject[T]](event.clusters.len)
      for i, cl in event.clusters:
        # now get the correct data from the raw event for the current cluster
        # and use it for reco
        var mcl = block:
          var tmp = cl.data
          for i, (x, y, ch) in cl.data:
            tmp[i] = (x: x, y: y, ch: rawEvent[y.int, x.int])
          tmp
        result.cluster[i] = recoCluster(mcl)

proc replaceBody*(procImpl: NimNode): NimNode =
  echo procImpl.repr
  let bod = procImpl.body
  var param = procImpl.params
  let name = postfix(procImpl.name, "*")
  #var res = nnkTemplateDef.newTree(name, newEmptyNode(),
  #                                 newEmptyNode(), param,
  #                                 newEmptyNode(), newEmptyNode())
  var res = procImpl
  let newBody = getAst(injectMarlinClusters())
  #res.add newBody
  res.body = newBody
  result = res
  echo result.repr

macro buildGeometry(): untyped =
  let ingridPath = shellVerbose:
    nimble path ingrid
  let geometryOriginal = shellVerbose:
    cat ($ingridPath[0])/"ingrid/private/geometry.nim"
  result = parseStmt(geometryOriginal[0])
  #echo result.repr
buildGeometry()
