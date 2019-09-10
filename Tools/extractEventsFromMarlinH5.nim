
#[
A simple script to extract a few events from the
`TimepixAnalysis/resources/background_splitted.2014+2015.h5`
file to be able to compare them to the reconstructed raw data.

]#

import ingrid / ingrid_types
import nimhdf5, sequtils, seqmath, strutils, random, json, os
import options, optionsutils
export options except get, unsafeGet
export optionsutils

import ggplotnim

const path = currentSourcePath().parentDir / "../resources/background_splitted.2014+2015.h5"
const NumEvents = 3

template read(astype: untyped, name: string): untyped  =
  block:
    var dset = h5f[(group / name).dset_str]
    dset[ind, float32].asType

template readVlen(vlen, dtype: untyped, name: string): untyped  =
  block:
    var dset = h5f[(group / name).dset_str]
    dset[ind, vlen, dtype]

proc readCluster(h5f: var H5FileObj, group: string, ind, evNum: int): Option[ClusterObject[Pix]] =
  ## compares the event number of the cluster at index `evNum` with the event
  ## number at the index. If they match returns a `some[ClusterObject[Pix]]`, else
  ## returns a `none`.
  echo "Reading cluster at ", ind, " expecting ", evNum
  let readEvNum = read(int, "EventNumber")
  if evNum != readEvNum:
    return none[ClusterObject[Pix]]()
  # else we continue
  var cObj: ClusterObject[Pix]
  cObj.hits = read(int, "NumberOfPixels")
  cObj.centerX = read(float, "PositionX")
  cObj.centerY = read(float, "PositionY")
  # NOTE: total charge does not actually contain `charges`, but `ToT` values!
  cObj.sumTot = read(int, "TotalCharge")
  # to calculate energy from charge, have to perform charge calibration manually then
  cObj.energy = read(float, "EnergyFromCharge")
  # now geometry
  var cGeom: ClusterGeometry
  cGeom.rmsLongitudinal = read(float, "RmsLongitudinal")
  cGeom.rmsTransverse = read(float, "RmsTransverse")
  cGeom.eccentricity = read(float, "Excentricity")
  cGeom.rotationAngle = read(float, "RotationAngle")
  cGeom.skewnessLongitudinal = read(float, "SkewnessLongitudinal")
  cGeom.skewnessTransverse = read(float, "SkewnessTransverse")
  cGeom.kurtosisLongitudinal = read(float, "KurtosisLongitudinal")
  cGeom.kurtosisTransverse = read(float, "KurtosisTransverse")
  cGeom.length = read(float, "Length")
  cGeom.width = read(float, "Width")
  cGeom.fractionInTransverseRms = read(float, "FractionWithinRmsTransverse")
  cGeom.lengthDivRmsTrans = cGeom.length / cGeom.rmsTransverse

  cObj.geometry = cGeom

  # finally read data of pixels
  let vlenUint8 = special_type(uint8)
  let vlenUint16 = special_type(uint16)
  let xdata = readVlen(vlenUint8, uint8, "XCoordinatesVector")
  let ydata = readVlen(vlenUint8, uint8, "YCoordinatesVector")
  let totdata = readVlen(vlenUint16, uint16, "ChargeValuesVector")
  for i in 0 ..< xdata.len:
    cObj.data.add (x: xdata[i], y: ydata[i], ch: totData[i])
  doAssert cObj.data.len == cObj.hits
  result = some(cObj)

proc readRecoEvent(h5f: var H5FileObj, group: string, ind: int): Option[RecoEvent[Pix]] =
  ## reads event with number `ind` from the file. Assumes we read from the
  ## no tracking group
  ## Returns a `some` if the cluster at index `ind` matches the following conditions
  ## - 200 < hits < 400
  ## - 1.2 < eccentricitry < 1.5
  # the major "difficulty" lies in the mapping of RecoEvent fields to the correct
  # datasets. No way around hardcoding those

  # No, as a matter of fact the main issue stems from the fact that the indices are
  # not related to the event numbers. Instead we have to do the following:
  # a given index points to an event number and run number cluster
  # we *have* to check whether there exists more clusters of the same event numbers
  # In principle (since the clusters are sorted by event number) it should be enough
  # to check the indices above and below the given one for clusters from the same
  # event number. Should work easily if we split this proc into one that reads the cluster
  # object and another that sits on top returning the actual RecoEvent.

  # get event number from event number dataset at this index
  var res: RecoEvent[Pix]
  res.eventNumber = read(int, "EventNumber")
  res.chipNumber = 0
  # get dataset length
  let nElems = h5f[(group / "EventNumber").dset_str].shape[0]

  # check conditions of current cluster. If matched
  let hits = read(int, "NumberOfPixels")
  let ecc = read(float, "Excentricity")
  if hits < 200 or
     hits > 400 or
     ecc < 1.2 or
     ecc > 1.5:
    return none[RecoEvent[Pix]]()

  let runNumber = read(int, "RunNumber")
  echo "[INFO]: Run number of idx ", ind, " : ", runNumber, " and eventNumber ", res.eventNumber

  # now fill clusters
  withSome readCluster(h5f, group, ind, res.eventNumber):
    some cluster: res.cluster.add cluster
    none: echo "Error: by definition this cannot happen!"
  # now search to lower indices
  var idx = ind - 1
  template whileCond(cond, modVar: untyped): untyped =
    while cond:
      withSome readCluster(h5f, group, idx, res.eventNumber):
        some cluster:
          echo "INFO: Found another cluster at idx: ", idx, " started at ", ind
          res.cluster.add cluster
        none:
          # not part of same event anymore, jump out
          break
      modVar
  whileCond(idx > 0):
    dec idx
  idx = ind + 1
  whileCond(idx < nElems):
    inc idx

  result = some(res)

proc `%`(p: Pix): JsonNode =
  result = %* { "x" : % p.x,
                "y" : % p.y,
                "ch" : % p.ch }

proc bindToDf[T](df: var DataFrame, clusters: seq[ClusterObject[T]]) =
  var dfs = newSeq[DataFrame]()
  for i, cl in clusters:
    let ldf = seqsToDf({ "x" : cl.data.mapIt(it.x),
                         "y" : cl.data.mapIt(it.y),
                         "ch" :  cl.data.mapIt(it.ch)})
    dfs.add ldf
  df = bind_rows(dfs, id = "from")
  if "from" notin df:
    df["from"] = toVector(toSeq(0 ..< df.len).mapIt(%~ "-1"))


proc main =

  var h5f = H5file(path, "r")
  defer: discard h5f.close()
  var showPlots = false
  var numEvents = NumEvents
  if paramCount() > 0:
    echo "paramCount() ", paramCount()
    if paramStr(1) == "--plot":
      showPlots = true
      echo "Updating plots in recoed_cluster.pdf"
      numEvents = int.high

  # first generate some indices randomly
  const group = "background-sc-not-only"
  let nElems = h5f[(group / "EventNumber").dset_str].shape[0] - 1

  var r = initRand(299792458)
  var recos = newSeq[RecoEvent[Pix]]()
  while recos.len < NumEvents:
    let idx = r.rand(nElems)
    withSome readRecoEvent(h5f, group, idx):
      some reco:
        recos.add reco
        if showPlots:
          var df = DataFrame()
          df.bindToDf(reco.cluster)
          ggplot(df, aes("x", "y", color = "from")) +
            geom_point() +
            ggtitle("Event number: " & $reco.eventNumber & " with clusters: " & $reco.cluster.len) +
            ggsave("recoed_cluster.pdf")
          if reco.cluster.len > 1:
            echo "FOUND ", reco.cluster.len
            sleep(3000)
          else:
            sleep(100)
      none:
        echo "INFO: Skipping non matching index ", idx
        discard

  # given all read events, dump them to Json
  var outJson = newJArray()
  for i, r in recos:
    let rJson = %r
    ## assert conversion back and forth idempotent
    doAssert rJson.to(RecoEvent[Pix]) == recos[i]
    outJson.add rJson

  var outf = open("marlinEvents.json", fmWrite)
  defer: outf.close()
  outf.write(outJson.pretty)

when isMainModule:
  main()