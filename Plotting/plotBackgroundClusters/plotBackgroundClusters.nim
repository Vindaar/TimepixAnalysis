import plotly, nimhdf5, docopt, arraymancer, seqmath
import strformat, sequtils, strutils, os
import helpers / utils
import ingrid / tos_helpers

const docStr = """
Usage:
  plotBackgroundClusters <likelihoodFile> [options]
  plotBackgroundClusters -h | --help
  plotBackgroundClusters --version

Options:
  -h, --help   Show this help
  --version    Show the version number

Given some Likelihood File (the output of the likelihood.nim), this script
creates a heatmap of the cluster centers that still remain visible on the plot.
In addition it highligths the gold region on the plot.
"""
const doc = withDocopt(docStr)

iterator extractClusters(h5f: var H5FileObj): (seq[float], seq[float]) =
  for num, grp in runs(h5f, likelihoodBase()):
    var mgrp = h5f[grp.grp_str]
    let centerChip = mgrp.attrs["centerChip", int]
    let chipGrp = h5f[(grp / "chip_" & $centerChip).grp_str]
    let cX = h5f[chipGrp.name / "centerX", float]
    let cY = h5f[chipGrp.name / "centerY", float]
    yield (cX, cY)

proc main =
  let args = docopt(doc)
  let h5file = $args["<likelihoodFile>"]

  var h5f = H5file(h5file, "r")

  var
    cX: seq[float]
    cY: seq[float]
  for centerX, centerY in extractClusters(h5f):
    cX.add centerX
    cY.add centerY

  # convert center positions to a 256x256 map
  var occ = zeros[int]([256, 256])
  for i in 0 ..< cX.len:
    var
      xInd = (256.0 * cX[i] / 14.0).round.int
      yInd = (256.0 * cY[i] / 14.0).round.int
    echo xInd, " for corresponding ", cX[i]
    echo yInd, " for corresponding ", cY[i]
    if xInd == 256:
      xInd = 255
    if yInd == 256:
      yInd = 255
    occ[yInd, xInd] += 1

  let perc99 = occ.toRawSeq.percentile(99)
  echo perc99
  occ = occ.clamp(0, 6)

  heatmap(occ.toSeq2D)
    .width(1600)
    .height(1600)
    .show()



when isMainModule:
  main()
