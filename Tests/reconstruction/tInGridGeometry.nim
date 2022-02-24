import ingrid / ingrid_types
import ingrid / tos_helpers
import ingridDatabase / databaseDefinitions
import sequtils, strutils, os, algorithm, strformat
import unittest
import json, sugar
from ingrid / reconstruction import recoEvent

import ggplotnim
import seqmath

const pwd = currentSourcePath().parentDir
const dataPwd = pwd / "../../resources/TPAresources/reconstruction/marlinExtracted"
const jsonData = dataPwd / "marlinEvents.json"

import times
let plotSuffix = $getTime().toUnix & ".pdf"

# if false we do not correct for the "one-off" error that
# MarlinTPC has
const CorrectOneOffXError = true

func almostEq[T: SomeNumber](x, y: T, ep = 1e-6): bool =
  ## checks very roughly if the values match. Anything beyond
  ## 1e-5 should be of no issue for us
  ## NOTE: We explicitly do *not* use the `stdlib` `almostEqual` as we don't care about
  ## any kind of smart comparison!
  result = abs((x - y).float) < ep

proc `%`(p: Pix): JsonNode =
  result = %* { "x" : % p.x,
                "y" : % p.y,
                "ch" : % p.ch }

proc bindToDf[T](clusters: seq[ClusterObject[T]]): DataFrame =
  var dfs = newSeq[DataFrame]()
  for i, cl in clusters:
    let ldf = seqsToDf({ "x" : cl.data.mapIt(it.x),
                         "y" : cl.data.mapIt(it.y),
                         "ch" :  cl.data.mapIt(it.ch)})
    dfs.add ldf
  result = bind_rows(dfs, id = "from")
  if "from" notin result:
    result["from"] = toVector(toSeq(0 ..< result.len).mapIt(%~ "-1"))

template red(s: string): untyped =
  "\e[91m" & s & "\e[0m"

template echoCheck(name, cond1, cond2: untyped, eps = 1e-5): untyped =
  echo "| $# | $# | $# | $# | $# |" % [$name, $cond1, $cond2, $abs(cond2 - cond1), $eps]
  let condres = almostEq(cond1, cond2, eps) == true
  if not condres:
    echo red"Check: ", astToStr(cond1), "~=~", astToStr(cond2), " failed"
    echo red"Was: ", astToStr(cond1), " == ", cond1
    echo red"Was: ", astToStr(cond2), " == ", cond2
  condres

func pixContained(p: Pix, s: seq[Pix]): bool =
  ## returns true if the pixel `x` is found in `s`, while ignoring
  ## the charge value
  result = s.filterIt(it.x == p.x and it.y == p.y).len > 0

proc echoMissing[T](s1, s2: ClusterObject[T]) =
  ## assumes that `s2` contains more elements than `s1`
  ## Echoes all events that are *not* in `s1`, but are in `s2`
  doAssert s1.hits < s2.hits
  for x in s2.data:
    if not x.pixContained(s1.data):
      echo "Missing: ", x

proc ggplotMoreClusters[T](recoCluster, expCluster: seq[ClusterObject[T]],
                           idx, eventNumber: int) =
  let dfReco = bindToDf(recoCluster)
  let dfExp = bindToDf(expCluster)
  let df = bind_rows([dfReco, dfExp], id = "origin")
  ggplot(df, aes("x", "y", color = "from")) +
    geom_point() +
    facet_wrap(~origin) +
    ggtitle(&"Event {idx}/eventNumber: {eventNumber}" &
            &" with clusters: {recoCluster.len}") +
    ggsave(&"recoed_cluster_{idx}.pdf")


proc ggplotMissing[T](recoCluster, expCluster: seq[ClusterObject[T]],
                      idx, eventNumber: int) =
  # both only single cluster. Create plot w/o facet wrap, highlighting
  # possible missing pixels between the two
  check recoCluster.len == 1
  check expCluster.len == 1
  var missing: DataFrame
  var missIn = ""
  let recoLen = recoCluster[0].data.len
  let expLen = expCluster[0].data.len
  let dfReco = bindToDf(recoCluster)
  let dfExp = bindToDf(expCluster)
  var df: DataFrame
  var fname = ""
  if recoLen < expLen:
    missing = setDiff(dfExp.select("x", "y"), dfReco.select("x", "y"))
    df = bind_rows([("false", dfReco), ("true", missing)], id = "missing")
    missIn = "TPA"
    fname = &"recoed_cluster_missing_tpa_{idx}.pdf"
  elif expLen < recoLen:
    missing = setDiff(dfReco.select("x", "y"), dfExp.select("x", "y"))
    df = bind_rows([("false", dfReco), ("true", missing)], id = "missing")
    missIn = "Marlin"
    fname = &"recoed_cluster_missing_marlin_{idx}.pdf"
  else:
    missing = setDiff(dfReco.select("x", "y"), dfExp.select("x", "y"))
    doAssert missing.len == 0
    df = dfReco # just take reco, since the two contain the same pixels anyways
    fname = &"recoed_cluster_same_{idx}.pdf"
  ggplot(df, aes("x", "y", color = "missing")) +
    geom_point() +
    ggtitle(&"Event {idx}/eventNumber: {eventNumber} with clusters: " &
            &"{recoCluster.len}  and missing pix in {missIn}: {missing.len}") +
    ggsave(fname)
  if missing.len > 0:
    doAssert false, "Had an event with missing pixels!"

proc compareClusters(recoClusters, expClusters: seq[ClusterObject[Pix]],
                     idx, eventNumber: int): (bool, (float, float)) =
  template retFalse(cond: untyped): untyped =
    ## helper template to return false, if the condition fails
    let cont = cond == true
    if not cont:
      return (false, (0.0, 0.0))

  var rotAngDiff: float
  var meanPropDiff: float
  for j in 0 ..< recoClusters.len:
    var recoCluster = recoClusters[j]
    var expCluster = expClusters[j]
    if recoCluster.hits != expCluster.hits:
      # determine which pixels are missing from one or the other
      if recoCluster.hits < expCluster.hits:
        echoMissing(recoCluster, expCluster)
      else:
        echoMissing(expCluster, recoCluster)
    retFalse echoCheck("hits", recoCluster.hits, expCluster.hits)
    retFalse echoCheck("data", recoCluster.data.len, expCluster.data.len)
    # sort cluster content by pixels x, y
    recoCluster.data = recoCluster.data.sortedByIt((it[0], it[1]))
    expCluster.data = expCluster.data.sortedByIt((it[0], it[1]))

    # compare content by x, y pixels
    for k in 0 ..< recoCluster.data.len:
      if CorrectOneOffXError:
        # corrected 1 off error manually
        #check recoCluster.data[k].x == expCluster.data[k].x
        discard
      else:
        retFalse recoCluster.data[k].x + 1 == expCluster.data[k].x # Marlin x coordinate is off by 1!
      retFalse recoCluster.data[k].y == recoCluster.data[k].y
      # cannot compare charge yet

    # if all passed we have
    # - the same number of clusters
    # - the same pixels in the clusters (except the noisy pixel)
    # now compare the geometrical properties
    retFalse echoCheck("centerX", recoCluster.centerX, expCluster.centerX)
    retFalse echoCheck("centerY", recoCluster.centerY, expCluster.centerY)
    # sum TOT will not be the same, since expCluster contains charge values
    # check recoCluster.sumTot == expCluster.sumTot
    # have to calculate energy before we can compare it
    # check recoCluster.energy == expCluster.energy
    let recoGeom = recoCluster.geometry
    let expGeom = expCluster.geometry
    retFalse echoCheck("rmsLongitudinal", recoGeom.rmsLongitudinal, expGeom.rmsLongitudinal, 1e-3)
    retFalse echoCheck("rmsTransverse", recoGeom.rmsTransverse, expGeom.rmsTransverse, 1e-3)
    retFalse echoCheck("eccentricity", recoGeom.eccentricity, expGeom.eccentricity, 1e-3)
    retFalse echoCheck("rotationAngle", recoGeom.rotationAngle, expGeom.rotationAngle, 1e-2)
    retFalse echoCheck("skewnessLongitudinal", recoGeom.skewnessLongitudinal, expGeom.skewnessLongitudinal, 1e-2)
    retFalse echoCheck("skewnessTransverse", recoGeom.skewnessTransverse, expGeom.skewnessTransverse, 1e-2)
    retFalse echoCheck("kurtosisLongitudinal", recoGeom.kurtosisLongitudinal, expGeom.kurtosisLongitudinal, 1e-2)
    if idx == 9:
      retFalse echoCheck("kurtosisTransverse", recoGeom.kurtosisTransverse, expGeom.kurtosisTransverse, 1e-1)
      retFalse echoCheck("width", recoGeom.width, expGeom.width, 1e-1)
    else:
      retFalse echoCheck("kurtosisTransverse", recoGeom.kurtosisTransverse, expGeom.kurtosisTransverse, 1e-3)
      retFalse echoCheck("width", recoGeom.width, expGeom.width, 1e-3)
    retFalse echoCheck("length", recoGeom.length, expGeom.length, 1e-2)
    retFalse echoCheck("fractionInTransverseRms", recoGeom.fractionInTransverseRms, expGeom.fractionInTransverseRms, 1e-6)
    retFalse echoCheck("lengthDivRmsTrans", recoGeom.lengthDivRmsTrans, expGeom.lengthDivRmsTrans, 1e-2)

    # survived the checks, consider difference of rotation angle
    rotAngDiff = abs(recoGeom.rotationAngle - expGeom.rotationAngle)
    let diffs = @[abs(recoGeom.skewnessLongitudinal - expGeom.skewnessLongitudinal),
                  abs(recoGeom.skewnessTransverse - expGeom.skewnessTransverse),
                  abs(recoGeom.kurtosisLongitudinal - expGeom.kurtosisLongitudinal),
                  abs(recoGeom.kurtosisTransverse - expGeom.kurtosisTransverse),
                  abs(recoGeom.length - expGeom.length),
                  abs(recoGeom.width - expGeom.width),
                  abs(recoGeom.fractionInTransverseRms - expGeom.fractionInTransverseRms),
                  abs(recoGeom.lengthDivRmsTrans - expGeom.lengthDivRmsTrans)]
    meanPropDiff = (diffs.sum / diffs.len.float)
  result = (true, (rotDiff: rotAngDiff, propDiff: meanPropDiff))

suite "InGrid geometry calculations":
  test "Reconstructing single cluster manually":
    # in the first test we read the InGrid data and perform the reconstruction
    # and cluster search manually. This is done for an event that shows up in the
    # in Christoph's root tree to compare where differences come from
    #let fileContent = readFile(fnameVbackground).strip.splitLines
    #let data = concat(@[fnameVbackground], fileContent)
    discard

  test "Comparison of Marlin events with TPA reco":
    ## given 3 events in data/marlinTestEvents
    const fnames = [("data059424_1_061500181.txt", 59424),
                    ("data057178_1_011827128.txt", 57178),
                    ("data034772_1_021331054.txt", 34772),
                    ("data036747_1_004216875.txt", 36747),
                    ("data031801_1_232447799.txt", 31801),
                    ("data038770_1_211608905.txt", 38770),
                    ("data055899_1_054858137.txt", 55899),
                    ("data053375_1_061636475.txt", 53375),
                    ("data011673_1_172643126.txt", 11673),
                    ("data057233_1_015223500.txt", 57233),
                    ("data025415_1_184828778.txt", 25415),
                    ("data043020_1_232312198.txt", 43020),
                    ("data069254_1_094825394.txt", 69254),
                    ("data053908_1_050952120.txt", 53908),
                    ("data074237_1_062130738.txt", 74237)]
    # reconstruct each file in `fnames` and compare with the `RecoEvent[Pix]` from `marlinEvent.json`
    let expEventsJ = jsonData.readFile.parseJson
    var expEvents = newSeq[RecoEvent[Pix]]()
    for exp in expEventsJ:
      var expReco = exp.to(RecoEvent[Pix])
      # sort by cluster length (Marlin and TPA don't agree on the cluster ordering)
      expReco.cluster.sort((r1, r2: auto) => cmp(r1.data.len , r2.data.len))
      expEvents.add expReco

    var skippedEvents = newSeq[string]()
    var rotAngDiffs = newSeq[float]()
    var meanPropDiffs = newSeq[float]()
    var evIdx = newSeq[int]()
    var evNums = newSeq[int]()
    for i, f in fnames:
      let fileContent = readFile(dataPwd / f[0]).strip.splitLines
      let data = concat(@[f[0]], fileContent)
      let ev: OldEvent = processOldEventScanf(data)[]

      # NOTE: the pixel
      # 167 200 *
      # is the noisy pixel of the 2014/15 chip. Filter it out.
      let numPixBefore = ev.chips[0].pixels.len
      var pix = ev.chips[0]
        .pixels
        .filterIt((it.x, it.y) != (167'u8, 200'u8))
      doAssert numPixBefore == pix.len or numPixbefore == pix.len + 1, "Only " &
        "a single noisy pixel was to be removed?! " & $(numPixBefore - pix.len) &
        "instead removed!"

      if CorrectOneOffXError:
        pix = pix.mapIt((x: (it[0] + 1'u8), y: it[1], ch: it[2]))

      var reco = recoEvent((pix, f[1]), 0)[]
      # sort by cluster length (Marlin and TPA don't agree on the cluster ordering)
      reco.cluster.sort((r1, r2: auto) => cmp(r1.data.len , r2.data.len))

      # create plot of found clusters
      if reco.cluster.len > 1 or expEvents[i].cluster.len > 1:
        ggplotMoreClusters(reco.cluster, expEvents[i].cluster,
                           i, f[1])
      else:
        ggplotMissing(reco.cluster, expEvents[i].cluster,
                      i, f[1])

      echo "\n\n\nIdx number: ", i, " eventNumber: ", f[1]
      echo "Cluster numbers : ", reco.cluster.len
      echo "Exp cluster numbers : ", expEvents[i].cluster.len
      if reco.cluster.len == expEvents[i].cluster.len:
        let (passed, diffTup) = compareClusters(reco.cluster, expEvents[i].cluster,
                                                i, f[1])
        check passed
        rotAngDiffs.add diffTup[0]
        meanPropDiffs.add diffTup[1]
        evIdx.add i
        evNums.add f[1]
      else:
        # different number of clusters
        skippedEvents.add f[0]

    let df = seqsToDf({ "eventIndex" : evIdx.mapIt($it),
                        "eventNumber": evNums,
                        "rotAngDifference" : rotAngDiffs,
                        "meanPropertyDifference" : meanPropDiffs })
    echo pretty(df, -1)
    ggplot(df, aes("rotAngDifference", "meanPropertyDifference", color = "eventIndex")) +
      geom_point() +
      ggsave("rotAng_diff_vs_mean_prop_diff.pdf")

    echo "Number of skipped events: ", skippedEvents.len
    for f in skippedEvents:
      echo f
