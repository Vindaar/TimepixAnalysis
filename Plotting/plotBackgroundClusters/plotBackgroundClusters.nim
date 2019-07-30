import plotly, nimhdf5, docopt, arraymancer, seqmath
import strformat, sequtils, strutils, os
import helpers / utils
import ingrid / tos_helpers
import json

import colorMaps

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

proc goldRegionOutline(maxVal: int): Tensor[int] =
  ## returns the outline of the gold region
  let min = (256.0 * 4.5 / 14.0).round.int
  let max = (256.0 * 9.5 / 14.0).round.int
  result = zeros[int]([256, 256])
  for coord, val in result:
    if (coord[0] in @[min, min + 1, min - 1] and
        coord[1] in {min .. max}) or
       (coord[0] in @[max, max + 1, max - 1] and
        coord[1] in {min .. max}) or
       (coord[1] in @[min, min + 1, min - 1] and
        coord[0] in {min .. max}) or
       (coord[1] in @[max, max + 1, max - 1] and
        coord[0] in {min .. max}):
      result[coord[0], coord[1]] = maxVal

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
    if xInd == 256:
      xInd = 255
    if yInd == 256:
      yInd = 255
    occ[yInd, xInd] += 1

  let perc99 = occ.toRawSeq.percentile(99)
  echo perc99
  occ = occ.clamp(0, 6)

  let outline = goldRegionOutline(5)
  occ = occ .+ outline

  let plt = heatmap(occ.toSeq2D)
    .width(1600)
    .height(1600)
    .toPlotJson

  template createPlot(cmap: untyped): untyped =
    plt.traces[0]["colorscale"] = cmap
    plt.traces[0]["zmax"] = % 6
    plt.traces[0]["zauto"] = % false
    plt.show(h5file.extractFilename & "_" & astToStr(cmap) & ".svg")
  createPlot(viridisZeroWhitePlotly)
  createPlot(viridisPlotly)
  createPlot(plasmaPlotly)
  createPlot(plasmaZeroWhitePlotly)

when isMainModule:
  main()
