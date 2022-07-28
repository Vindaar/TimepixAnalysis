import nimhdf5, arraymancer, seqmath, chroma, cligen
import std / [strformat, sequtils, strutils, os, json, sets]
import helpers / utils
import ingrid / tos_helpers
#import plotly
#import plotly / color
import ggplotnim
import ggplotnim / [ggplot_vegatex]

const docStr = """
Given some Likelihood File (the output of the likelihood.nim), this script
creates a heatmap of the cluster centers that still remain visible on the plot.
In addition it highligths the gold region on the plot.
"""

iterator extractClusters(h5f: var H5FileObj): (seq[float], seq[float], seq[float]) =
  for num, grp in runs(h5f, likelihoodBase()):
    var mgrp = h5f[grp.grp_str]
    var centerChip: int
    try:
      centerChip = mgrp.attrs["centerChip", int]
    except KeyError:
      echo "WARNING: could not find `centerChip` attribute. Will assume 3!"
      centerChip = 3
    var chipGrp: H5Group
    try:
      chipGrp = h5f[(grp / "chip_" & $centerChip).grp_str]
    except KeyError:
      echo "INFO: no events left in run number " & $num & " for chip " & $centerChip
      continue
    let cX = h5f[chipGrp.name / "centerX", float]
    let cY = h5f[chipGrp.name / "centerY", float]
    let cE = h5f[chipGrp.name / "energyFromCharge", float]
    yield (cX, cY, cE)

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

when false:
  const noisyPixels = [(64, 109),
                       (64, 110),
                       (67, 112),
                       (65, 108),
                       (66, 108),
                       (67, 108),
                       (65, 109),
                       (66, 109),
                       (67, 109),
                       (68, 109),
                       (65, 110),
                       (66, 110),
                       (67, 110),
                       (65, 111),
                       (66, 111),
                       (67, 111),
                       (68, 110),
                       (68, 109),
                       (68, 111),
                       (68, 108),
                       (67, 107),
                       (66, 111),
                       (69, 110)]

var noisyPixels = newSeq[(int, int)]()
for x in 150 .. 250:
  for y in 130 .. 162:
    noisyPixels.add (x, y)

for x in 125 .. 135:
  for y in 110 .. 120:
    noisyPixels.add (x, y)

proc readClusters(file: string, cTab: var CountTable[(int, int)],
                  filterNoisyPixels: bool,
                  filterEnergy: float) =
  var h5f = H5open(file, "r")
  func toPixel(s: float): int = (256.0 * s / 14.0).round.int
  for centerX, centerY, energy in extractClusters(h5f):
    for (cX, cY, cE) in zipEm(centerX, centerY, energy):
      # drop clusters with energy too large
      if filterEnergy > 0.0 and cE > filterEnergy: continue

      let (pX, pY) = (cX.toPixel, cY.toPixel)
      #echo (pX, pY) in noisyPixels, " pos ", (pX, pY)
      if filterNoisyPixels and (pX, pY) notin noisyPixels:
        cTab.inc((cX.toPixel, cY.toPixel))
      elif not filterNoisyPixels:
        cTab.inc((cX.toPixel, cY.toPixel))
  doAssert h5f.close() >= 0

proc readClusters(file: string, filterNoisyPixels: bool,
                  filterEnergy: float): CountTable[(int, int)] =
  result = initCountTable[(int, int)]()
  file.readClusters(result, filterNoisyPixels, filterEnergy)

proc readClusters(files: seq[string], filterNoisyPixels: bool,
                  filterEnergy: float): CountTable[(int, int)] =
  result = initCountTable[(int, int)]()
  for f in files:
    f.readClusters(result, filterNoisyPixels, filterEnergy)

proc toDf(cTab: CountTable[(int, int)]): DataFrame =
  var
    cX, cY, cC: seq[int]
  for (pos, count) in pairs(cTab):
    cX.add pos[0]
    cY.add pos[1]
    cC.add count
  result = toDf({"x" : cX, "y" : cY, "count" : cC})
    .arrange("count", order = SortOrder.Ascending)

proc writeNoisyClusters(cTab: CountTable[(int, int)],
                        threshold: int) =
  ## Writes a `noisyPixels` array to file with all pixels that have
  ## more counts than `threshold`.
  cTab.toDf().writeCsv("/t/data_tab.csv")
  var f = open("/t/noisy_clusters.txt", fmWrite)
  f.write("const noisyPixels = [\n")
  for tup, val in cTab:
    if val > threshold:
      f.write($tup & ",\n")
  f.write("]")
  f.close()

proc main(
  files: seq[string], suffix = "", title = "",
  names: seq[string] = @[], # can be used to create facet plot. All files with same names will be stacked
  filterNoisyPixels = false,
  filterEnergy = 0.0,
  useTikZ = false,
  writeNoisyClusters = false,
  threshold = 2) =

  let cTab = readClusters(files, filterNoisyPixels, filterEnergy)
  if writeNoisyClusters:
    cTab.writeNoisyClusters(threshold)

  var df = newDataFrame()
  if names.len == 0 or names.deduplicate.len == 1:
    echo files
    df = cTab.toDf()
  else:
    doAssert names.len == files.len, "Need same number of names as input files!"
    for i in 0 ..< names.len:
      var dfLoc = readClusters(@[files[i]], filterNoisyPixels, filterEnergy).toDf()
      dfLoc["Type"] = names[i]
      df.add dfLoc
    # now sort all again
    df = df.arrange("count", order = SortOrder.Ascending)

  createDir("plots")
  if names.len == 0 or names.deduplicate.len == 1:
    let totalEvs = df["count", int].sum()
    let plt = ggplot(df, aes("x", "y", color = "count")) +
      geom_point(size = some(1.0)) +
      xlim(0, 256) + ylim(0, 256) +
      xlab("x [Pixel]") + ylab("y [Pixel]") +
      margin(top = 1.75) +
      scale_color_continuous(scale = (low: 0.0, high: 5.0))
    if not useTikZ:
      let fname = &"plots/background_cluster_centers{suffix}.pdf"
      echo "[INFO]: Saving plot to ", fname
      plt + ggtitle(title & &". Total # cluster = {totalEvs}") +
        ggsave(fname)
    else:
      let fname = &"/home/basti/phd/Figs/backgroundClusters/background_cluster_centers{suffix}"
      echo "[INFO]: Saving plot to ", fname
      plt + ggtitle(title & r". Total \# cluster = " & $totalEvs) +
        ggvegatex(fname)
  else:
    var totalEvs = newSeq[string]()
    for tup, subDf in groups(df.group_by("Type")):
      let numEvs = subDf["count", int].sum()
      totalEvs.add &"{tup[0][1]}: {numEvs}"
    let fname = &"/home/basti/org/Figs/statusAndProgress/IAXO_TDR/background_cluster_centers{suffix}.pdf"
    echo "[INFO]: Saving plot to ", fname
    let ticks = arange(0, 260, 5).mapIt(it.float)
    ggplot(df, aes("x", "y", color = "count")) +
      facet_wrap("Type") +
      geom_point(size = some(1.0)) +
      #xlim(0, 256) + ylim(0, 256) +
      #xlab("x [Pixel]", tickMargin = 2.0) + ylab("y [Pixel]", margin = 2.0, tickMargin = -0.5) + # for TikZ
      xlab("x [Pixel]") + ylab("y [Pixel]") +
      scale_x_continuous(breaks = ticks) + scale_y_continuous(breaks = ticks) +
      margin(top = 1.75) +
      scale_color_continuous(scale = (low: 0.0, high: 15.0)) +
      ggtitle(r"Total # clusters " & $totalEvs.join(", ")) +
      ggsave(fname,
             width = 900, height = 480)
             #useTeX = true, standalone = true) # onlyTikZ = true)

  when false:
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

    var plt = heatmap(occ.toSeq2D)
      .width(1600)
      .height(1600)

    template createPlot(cmap: untyped): untyped =
      plt = plt.zmax(6)
      plt = plt.colormap(cmap)
      plt.show(h5file.extractFilename & "_" & astToStr(cmap) & ".svg")
    createPlot(Viridis)
    createPlot(ViridisZeroWhite)
    createPlot(Plasma)
    createPlot(PlasmaZeroWhite)

when isMainModule:
  dispatch main
