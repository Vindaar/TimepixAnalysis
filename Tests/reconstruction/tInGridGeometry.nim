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
const dataPwd = pwd / "data/marlinTestEvents"
const jsonData = dataPwd / "marlinEvents.json"

import times
let plotSuffix = $getTime().toUnix & ".pdf"

#proc checkEqual(ingrid: ): bool =
  # checks whther the ingrid event is equal to our expectation
#  return true

func almostEqual(x, y: float, ep = 1e-5): bool =
  ## checks very roughly if the values match. Anything beyond
  ## 1e-5 should be of no issue for us
  result = (x - y) < ep

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
                    ("data034772_1_021331054.txt", 34772)]
    # reconstruct each file in `fnames` and compare with the `RecoEvent[Pix]` from `marlinEvent.json`
    let expEventsJ = jsonData.readFile.parseJson
    var expEvents = newSeq[RecoEvent[Pix]]()
    for exp in expEventsJ:
      var expReco = exp.to(RecoEvent[Pix])
      # sort by cluster length (Marlin and TPA don't agree on the cluster ordering)
      expReco.cluster.sort((r1, r2: auto) => cmp(r1.data.len , r2.data.len))
      expEvents.add expReco

    for i, f in fnames:
      let fileContent = readFile(dataPwd / f[0]).strip.splitLines
      let data = concat(@[f[0]], fileContent)
      let ev: OldEvent = processOldEventScanf(data)[]

      # NOTE: the pixel
      # 167 200 *
      # is the noisy pixel of the 2014/15 chip. Filter it out.
      let pix = ev.chips[0].pixels.filterIt(it.x != 167'u8 and it.y != 200'u8)

      var reco = recoEvent((pix, f[1]), 0)[]
      # sort by cluster length (Marlin and TPA don't agree on the cluster ordering)
      reco.cluster.sort((r1, r2: auto) => cmp(r1.data.len , r2.data.len))

      # create plot of found clusters
      var dfReco = DataFrame()
      dfReco.bindToDf(reco.cluster)
      var dfExp = DataFrame()
      dfExp.bindToDf(expEvents[i].cluster)
      let df = bind_rows([dfReco, dfExp], id = "origin")
      echo dfReco
      echo dfExp
      if "from" in df:
        # work around issue in ggplotnim
        ggplot(df, aes("x", "y", color = "from")) +
          geom_point() +
          facet_wrap(~origin) +
          ggtitle("Is event " & $i & " with clusters: " & $reco.cluster.len) +
          ggsave(&"recoed_cluster_{i}.pdf")
      else:
        ggplot(df, aes("x", "y")) +
          geom_point() +
          facet_wrap(~origin) +
          ggsave(&"recoed_cluster_{i}.pdf")

      echo "Cluster numbers : ", reco.cluster.len
      echo "Exp cluster numbers : ", expEvents[i].cluster.len
      # TODO: "fix" the code below. Does not quite work, simply because Marlin and TPA
      # produce different resutls on the cluster algorithm already.
      # Somehow fix that future me.
      #for j in 0 ..< reco.cluster.len:
      #  echo "Comparing cluster ", j, " of event ", i
      #  echo reco.cluster[j].data
      #  echo expEvents[i].cluster[j].data
      #  check reco.cluster[j].hits == expEvents[i].cluster[j].hits
      #  check reco.cluster[j].data.len == expEvents[i].cluster[j].data.len
      #  # sort cluster content by pixels x, y
      #  reco.cluster[j].data = reco.cluster[j].data.sortedByIt((it[0], it[1]))
      #  expEvents[i].cluster[j].data = expEvents[i].cluster[j].data.sortedByIt((it[0], it[1]))
      #
      #  # compare content by x, y pixels
      #  for k in 0 ..< reco.cluster[j].data.len:
      #    echo "Comparing ", reco.cluster[j].data[k]
      #    echo "With ", expEvents[i].cluster[j].data[k]
      #    check reco.cluster[j].data[k].x + 1 == expEvents[i].cluster[j].data[k].x # Marlin x coordinate is off by 1!
      #    check reco.cluster[j].data[k].y     == expEvents[i].cluster[j].data[k].y
          # cannot compare charge yet
