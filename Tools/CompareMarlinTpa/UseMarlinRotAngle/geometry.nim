import shell, os, macros, sequtils, tables, algorithm, strutils
import ingrid / ingrid_types
import nimhdf5, arraymancer
import locks
import sugar # for collect

import ggplotnim
## this is the monkey patched version of ingrid / geometry.nim
#var L = Lock()
#L.initLock()

const BackgroundRunPath = "background-sc-not-only"
const CalibrationRunPrefix = "calibration-fe55-run$#"

#template argFilterIt(s, cond: untyped): untyped =
#  ## returns indices where filter condition is true
#  var idx = 0
#  var res: seq[int]
#  for it {.inject.} in s:
#    if cond:
#      res.add idx
#    inc idx
#  res

#var sortedRuns: seq[MarlinRuns]
#for run, mrun in mappedRuns:
#  sortedRuns.add mrun
#sortedRuns = sortedByIt(sortedRuns, it.run)

#template findIt(s: typed, cond: untyped): untyped =
#  ## finds the first appearence of `cond` in `s`, returns the element
#  type retType = type(s[0])
#  var res: retType
#  var found = false
#  for i, it {.inject.} in s:
#    if cond:
#      res = it
#      found = true
#      break
#  if not found:
#    echo "!!! coult not find "
#    quit()
#  res

#proc readMarlinEvent(evNum, run: int): MarlinEvent {.inline.} =
#  {.gcsafe.}:
#    when defined(hijackBackground):
#      let mRun = mappedRuns[run]
#    else:
#      let mRun = mapRunsToEvents(run)
#    result = mRun.events.findIt(it.eventNumber == evNum)
#    #echo "LAST EL ", mRun.events[mRun.events.high], "\n\n\n"


# Helper type to access correct data from Marlin data
type
  MarlinCluster = object
    data: Cluster[Pix]
    globalIdx: int
    rotAngle: float
  MarlinEvent = object
    numMismatches: int
    eventNumber: int
    clusters: seq[MarlinCluster]
  MarlinRuns = object
    run: int
    events: seq[MarlinEvent]

var
  xVec: seq[seq[uint8]]
  yVec: seq[seq[uint8]]
  evsMarlin: seq[int]
  #rotAngsMarlin: seq[float32]

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
  #echo "READING AT INDEX ", idx
  let x = xVec[idx]
  let y = yVec[idx]
  doAssert x.len == y.len
  doAssert x.len > 0
  result.data = newSeq[Pix](x.len)
  for i in 0 ..< x.len:
    # NOTE: we leave the charge empty! ToT values will be taken from raw files
    result.data[i] = (x: x[i], y: y[i], ch: 0'u16)
  doAssert result.data.len > 0

proc splitBySame(idx: seq[int]): seq[(int, seq[int])] =
  ## splits the given sequence by the indices belonging to the same numbers
  result = newSeqOfCap[(int, seq[int])](idx.len)
  var old = -1
  var globalIdx = 0
  var tmp = newSeq[int]()
  for i, x in idx:
    if old > 0 and old == x:
      tmp.add globalIdx
    else:
      if tmp.len > 0:
        result.add (idx[i - 1], tmp)
      tmp.setLen(0)
      tmp.add globalIdx
    old = x
    inc globalIdx
  result.add (idx[^1], tmp)

proc buildMarlinEvents(run: int, clusterIdx, events: seq[int]): seq[MarlinEvent] =
  ## reads all required information for the given event indices for run
  ## `run`
  # start by splitting the clusters by events
  #let eventIdx = splitBySame(clusterIdx)
  var oldEv = -1
  var mEvent: MarlinEvent
  #echo "for ", clusterIdx[0..15]
  #echo "!!! ", evsMarlin[clusterIdx[0] .. clusterIdx[15]]
  for cIdx in clusterIdx:
    # for each event now walk all clusters and read data
    let evNum = events[cIdx]
    if cIdx == clusterIdx[clusterIdx.high]:
      echo "Last event number: ", evNum, " at index ", cIdx, " For run ", run
    #if run == 245:
    #  echo "EVENT NUMBER ", evNum, " AND INDEX ", cIdx
    if evNum == oldEv:
      doAssert mEvent.clusters.len > 0
      mEvent.clusters.add readClusterData(cIdx)
      doAssert mEvent.clusters[^1].data.len > 0, " was " & $evNum & " at " & $cIDx
      if cIdx == clusterIdx.high:
        result.add mEvent
    else:
      #echo "MEVENT ", mEvent, " for evNum ", evNum, " w/ old ", oldEv
      if mEvent.clusters.len > 0:
        result.add mEvent
      mEvent = MarlinEvent(eventNumber: evNum)
      mEvent.clusters.add readClusterData(cIdx)
      doAssert mEvent.clusters[^1].data.len > 0, " was " & $evNum & " at " & $cIDx
      doAssert mEvent.clusters.len > 0

    oldEv = evNum


proc mapRunsToEvents(mruns: TableRef[int, MarlinRuns], runs: seq[int], evs: seq[int] = @[]) =
  # first split all runs into a bunch of seq[seq[int]] belonging
  # to the same run number
  doAssert runs.len == evs.len
  let runSeqs = splitBySame(runs)
  # iterate all those runs and create a MarlinRuns object for each
  for (run, eventForRunIdx) in runSeqs:
    # build the MarlinEvents
    echo "Building run ", run, " :) :) :)"
    let mEvents = buildMarlinEvents(run, eventForRunIdx, evs)
    # now combine mEvens with the run, if it already exists
    if run in mruns:
      # combine by inserting all elements of mEvents into existing at correct order
      var mr = mruns[run]
      let l1 = mr.events.mapIt(it.eventNumber).len
      for e in mEvents:
        mr.events.insert(e, mr.events.lowerBound(e,
          (e1, e2) => system.cmp(e1.eventNumber, e2.eventNumber)))
      let l2 = mr.events.mapIt(it.eventNumber).len
      echo "Len before ", l1, " after ", l2
      mruns[run] = mr
    else:
      echo "NOT IN IT! ", run
      mruns[run] = MarlinRuns(run: run, events: mEvents)

proc mapRunsToEvents(h5f: H5FileObj, run: int): TableRef[int, MarlinRuns] =
  # first split all runs into a bunch of seq[seq[int]] belonging
  # to the same run number
  result = newTable[int, MarlinRuns]()
  let evsMarlin = h5f[CalibrationRunPrefix % $run / "EventNumber", float32].mapIt(it.int)
  let eventForRunIdx = toSeq(0 ..< evsMarlin.len)
  let mEvents = buildMarlinEvents(run, eventForRunIdx, evsMarlin)
  result[run] = MarlinRuns(run: run, events: mEvents)

var mappedRuns: TableRef[int, MarlinRuns]

# TODO: for now just read all X, Y, Charge data in global scope. If
# problematic from memory standpoint, move that to `readClusterData` and
# read by index
let specialU8 = special_type(uint8)
#let specialU16 = special_type(uint16)

when defined(hijackBackground):
  var h5fMarlin {.global.} = H5file("/mnt/1TB/CAST/2014_15/oldBackground/background_splitted.2014+2015.h5", "r") #H5file("../../../resources/background_splitted.2014+2015.h5", "r")
  # iterate all groups and read data for each
  #let grp = h5fMarlin[(BackgroundRunPath).grp_str]
  mappedRuns = newTable[int, MarlinRuns]()
  for grp in items(h5fMarlin, depth = 1):
    echo "!!! ", grp
    let runsMarlin = h5fMarlin[grp.name / "RunNumber", float32].mapIt(it.int)
    evsMarlin = h5fMarlin[grp.name / "EventNumber", float32].mapIt(it.int)
    xVec = h5fMarlin[grp.name / "XCoordinatesVector", specialU8, uint8]
    yVec = h5fMarlin[grp.name / "YCoordinatesVector", specialU8, uint8]
    mapRunsToEvents(mappedRuns, runsMarlin, evsMarlin)

  # rotAngsMarlin = h5fMarlin[grp.name / "RotationAngle", float32]
elif defined(hijackCalibration):
  var h5fMarlin {.global.} = H5file("../../../resources/calibration-fe55-splitted.2014+2015.h5", "r")
  mappedRuns = newTable[int, MarlinRuns]()
else:
  static: error("Either define `-d:hijackBackground` or `-d:hijackCalibration`!")

var globalMarlinEvCounter = 0

var lastRun: int
var curRun: int

proc readEvByIdx(index, run: int): MarlinEvent {.inline.} =
  {.gcsafe.}:
    curRun = run
    if run in mappedRuns:
      if lastRun != curRun:
        globalMarlinEvCounter = 0
      #echo "Getting run ", run, "\n\n\n\n"
      let mRun = mappedRuns[run]
      #echo "Num index ", index, " w/ len ", mRun.events.len
      if index < mRun.events.len:
        result = mRun.events[index]
      else:
        echo "doesn't exist anymore! ", index
      lastRun = curRun
    else:
      # clear table to avoid running out of memory
      mappedRuns.clear()
      globalMarlinEvCounter = 0
      xVec = h5fMarlin[CalibrationRunPrefix % $run / "XCoordinatesVector", specialU8, uint8]
      yVec = h5fMarlin[CalibrationRunPrefix % $run / "YCoordinatesVector", specialU8, uint8]
      mappedRuns = mapRunsToEvents(h5fMarlin, run)
      result = readEvByIdx(globalMarlinEvCounter, run)

#iterator readMarlinEvent(run: int): iterator(): MarlinEvent =
#  {.gcsafe.}:
#    let mRun = mappedRuns[run]
#    result = iterator(): MarlinEvent =
#      for x in mRun.events:
#        yield x

var scratchTensor = newTensor[uint16](256, 256)
var scratchTab = initTable[(uint8, uint8), uint16]()

proc fillEvent(t: var Tensor[uint16], data: seq[Pix]) =
  ## fills the tensor with data building the full raw data event, from
  ## which we extract the ToT values
  #t.apply_inline(0)
  for (x, y, ch) in data:
    t[y.int, x.int] = ch


proc fillEvent(t: var TableRef[(uint8, uint8), uint16], data: seq[Pix]) =
  ## fills the tensor with data building the full raw data event, from
  ## which we extract the ToT values
  #t.clear()
  for (x, y, ch) in data:
    t[(y, x)] = ch

proc getClusterData(t: Tensor[uint16], cl: MarlinCluster): seq[Pix] =
  result = newSeq[Pix](cl.data.len)
  for i, (x, y, ch) in cl.data:
    result[i] = (x: x - 1, y: y, ch: t[y.int, x.int - 1])

proc getClusterData(t: TableRef[(uint8, uint8), uint16], cl: MarlinCluster): seq[Pix] =
  result = newSeq[Pix](cl.data.len)
  for i, (x, y, ch) in cl.data:
    result[i] = (x: x, y: y, ch: t[(y, x)])

let noisePixels = readFile("noisePixels.txt").strip.splitLines.filterIt(it.len > 0)
  .mapIt((it.splitWhitespace[0].parseInt.uint8,
          it.splitWhitespace[1].parseInt.uint8))
  .toHashSet()


proc handleInner(mdat: seq[Pix], eventNumber, run: int, recurse = true): MarlinEvent =
  {.gcsafe.}:
    #var mdat = dat[0].filterIt((it.x, it.y) != (167'u8, 200'u8))
    #var mdat = dat[0].filterIt(not ((it.x, it.y) == (167'u8, 200'u8) or
    #                                (it.x, it.y) == (111'u8, 103'u8) or
    #                                (it.x, it.y) == (240'u8, 129'u8) or
    #                                (it.x, it.y) == (241'u8, 129'u8) ) )
    #echo "MDAT!!!!! ", mdat
    var numMismatch {.global.}: int # for small mismatches
    var numBadEvents {.global.}: int # for big mismatches
    var numMismatchPixel {.global.}: int
    var numGoodReps {.global.}: int
    var consecutiveMismatches {.global.}: int

    if numBadEvents > 100:
      quit("COULDN'T FIND MATCH ANYMORE")

    result = readEvByIdx(globalMarlinEvCounter, run)
    #let event = readMarlinEvent(globalMarlinEvCounter, run)
    # now build the tensor of the raw event
    scratchTensor.fillEvent(mdat)
    #echo "RAW EVENT NUMBER ", dat[0]
    #echo "EVENT EVENT NUMBER ", event#.eventNumber
    #scratchTab.fillEvent(dat[0])
    let dfRaw = seqsToDf({ "x" : mdat.mapIt(it.x + 1),
                           "y" : mdat.mapIt(it.y),
                           "ch" : mdat.mapIt(it.ch) })
    var
      xM: seq[uint8]
      yM: seq[uint8]
      chM: seq[uint16]
    if result.clusters.len == 0:
      echo result
      echo "ev number ", eventNumber
      echo run
      return
      #quit()
    for i, cl in result.clusters:
      xM.add cl.data.mapIt(it.x)
      #echo "Len of xm ", xm.len, " compared ", dfRaw.len
      yM.add cl.data.mapIt(it.y)
      chM.add cl.data.mapIt(it.ch)
    #echo "----------------------------"

    #echo "At event ", eventNumber, " and ", globalMarlinEvCounter, " marlin event ", result.eventnumber
    #result.numMismatches = abs(xM.len - dfRaw.len)
    #let s1 = toHashSet(mdat.mapIt((it.x + 1, it.y)))
    #let s2 = toHashSet(toSeq(0 ..< xM.len).mapIt((xM[it], yM[it])))
    let s1 = toHashSet(mdat.mapIt((it.x + 1, it.y)))
    let s2 = toHashSet(toSeq(0 ..< xM.len).mapIt((xM[it], yM[it])))
    #echo s1
    #echo s2
    #echo symmetricdifference(s1, s2)
    #let dfRawMod = dfRaw.select("x", "y").arrange("x")
    #let dfMNoch = dfMarlin.select("x", "y").arrange("x")
    #result.numMismatches = setDiff(dfRawMod, dfMNoch).len
    result.numMismatches = symmetricdifference(s1, s2).card
    #result.numMismatches = abs(xM.len - dfRaw.len)
    if result.numMismatches > 5:
      echo "At event ", eventNumber, " and ", globalMarlinEvCounter, " marlin event ", result.eventnumber
      #let dfRawMod = dfRaw.select("x", "y").arrange("x")
      #let dfMNoch = dfMarlin.select("x", "y").arrange("x")
      #echo dfRawMod.pretty(-1)
      #echo dfMNoch.pretty(-1)
      ##echo xM.sorted
      ##echo mdat.mapIt(it.x + 1).sorted
      #let dfDiff = setDiff(dfRawMod, dfMNoch)
      #echo symmetricdifference(s1, s2)
      echo symmetricdifference(s1, s2).len
      #var noiseF = open("noisePixels.txt", fmAppend)
      #for row in dfDiff:
      #  noiseF.write($(row["x"].toInt - 1) & " " & $row["y"].toInt & "\n")
      #noiseF.close()
      let dfMarlin = seqsToDf({ "x" : xM,
                                "y" : yM,
                                "ch" : chM })

      let df = bind_rows([dfRaw, dfMarlin], id = "from")
      ggplot(df, aes = aes("x", "y", color = "from")) +
        geom_point(data = df.filter(f{"from" == "0"})) +
        geom_point(data = df.filter(f{"from" == "1"}), size = 1.5) +
        scale_x_continuous() +
        scale_y_continuous() +
        ggtitle("BAD dfRawlen: " & $dfRaw.len & " in dfMarlin: " & $dfMarlin.len &
          "of event: " & $result.eventNumber, & " MCount " & $globalMarlinEvCounter) +
        #ggsave("bad_event_" & $globalMarlinEvCounter & "_run_" & $run & ".pdf")
        ggsave("bad_event.pdf")

      #echo ".................."
      #echo "Mismatches was: ", numMismatch
      #echo "Mismatched pixels was: ", numMismatchPixel
      echo "Consecutive mismatches: ", consecutiveMismatches
      # skip this event without increasing global marlin counter
      if consecutiveMismatches > 1:
        # search ahead for N events
        var oldCount = globalMarlinEvCounter
        var oldMismatches = result.numMismatches
        var broken = false
        if recurse:
          for i in 0 ..< 20:
            echo "CHECKING FOR ", i, " w/ counter ", globalMarlinEvCounter
            inc globalMarlinEvCounter
            result = handleInner(mdat, eventNumber, run, false)
            if result.numMismatches <= 5:
              echo "FOUND GOOD REPLACEMENT! ", numGoodReps
              inc numGoodReps
              #inc globalMarlinEvCounter
              return result
          if result.numMismatches > 5:
            echo "BAD EVENT"
            broken = true
          if broken:
            # also search backwards
            globalMarlinEvCounter = oldCount
            for i in 0 ..< 10:
              echo "CHECKING backwards ", i, " w/ counter ", globalMarlinEvCounter
              dec globalMarlinEvCounter
              if globalMarlinEvCounter < 0:
                break
              result = handleInner(mdat, eventNumber, run, false)
              if result.numMismatches <= 5:
                echo "FOUND GOOD REPLACEMENT BACKWARDS", numGoodReps
                inc numGoodReps
                #inc globalMarlinEvCounter
                return result
        globalMarlinEvCounter = oldCount

        # did not find event
        if recurse:
          # so we know this was the first proc call
          inc numBadEvents
          echo "NUM BAD EVENTS ", numBadEvents
        #dec globalMarlinEvCounter, 2
        #consecutiveMismatches = 0
        #return handleInner(mdat, eventNumber, run)
      inc consecutiveMismatches
      return
      #quit()
      #echo df.pretty(-1)
    #elif xM.len != dfRaw.len:
    elif result.numMismatches > 0: # and implies result.numMismatches <= 5
      echo "Difference, but only of: ", result.numMismatches
      let dfMarlin = seqsToDf({ "x" : xM,
                                "y" : yM,
                                "ch" : chM })
      let df = bind_rows([dfRaw, dfMarlin], id = "from")
      ggplot(df, aes = aes("x", "y", color = "from")) +
        geom_point(data = df.filter(f{"from" == "0"})) +
        geom_point(data = df.filter(f{"from" == "1"}), size = 1.5) +
        scale_x_continuous() +
        scale_y_continuous() +
        ggtitle("SmallMismatch Number of elems in dfRaw: " & $dfRaw.len & " in dfMarlin: " & $dfMarlin.len) +
        ggsave("smallMismatch_event_" & $globalMarlinEvCounter & "_run_" & $run & ".pdf")
      inc numMismatch
      inc numMismatchPixel, abs(xm.len - dfRaw.len)
    else:
      consecutiveMismatches = 0

    #echo df
    #if df.len > 2:
      #ggplot(df, aes = aes("x", "y", color = "from")) +
      #  geom_point(data = df.filter(f{"from" == "0"})) +
      #  geom_point(data = df.filter(f{"from" == "1"}), size = 1.5) +
      #  scale_x_continuous() +
      #  scale_y_continuous() +
      #  ggtitle("Number of elems in dfRaw: " & $dfRaw.len & " in dfMarlin: " & $dfMarlin.len) +
      #  ggsave("event_" & $globalMarlinEvCounter & ".pdf")
      # discard

    inc globalMarlinEvCounter
    #if globalMarlinEvCounter == 100:
    #  quit()


template injectMarlinClusters() =
  ## This shit's slow af
  # instead of using our normal cluster finder we read the clusters found by the
  # Marlin algorithm for event number `evNum` instead
  result = new RecoEvent[Pix]
  static: echo type(dat)
  result.event_number = dat.eventNumber
  result.chip_number = chip
  var mdat: seq[Pix]
  {.gcsafe.}:
    for p in dat[0]:
      if (p.x, p.y) notin noisePixels and p.ch < 400:
        mdat.add p
  if mdat.len > 2:
    {.gcsafe.}:
      let event = handleInner(mdat, result.event_number, run, recurse = true)
      result.cluster = newSeq[ClusterObject[Pix]](event.clusters.len)
      for i, cl in event.clusters:
        # now get the correct data from the raw event for the current cluster
        # and use it for reco
        #let mcl = getClusterData(scratchTab, cl)
        let mcl = getClusterData(scratchTensor, cl)
        result.cluster[i] = recoCluster(mcl)

#template injectMarlinClusters() =
#  result = hijackWork(dat, chip, run)

proc replaceBody*(procImpl: NimNode): NimNode =
  echo procImpl.repr
  let bod = procImpl.body
  var param = procImpl.params
  let name = postfix(procImpl.name, "*")
  var res = procImpl
  let newBody = getAst(injectMarlinClusters())
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
