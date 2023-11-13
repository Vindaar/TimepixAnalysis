import nimhdf5, seqmath, chroma, cligen
import std / [strformat, sequtils, strutils, os, sets]
import helpers / utils
import ingrid / tos_helpers
import ggplotnim
import ggplotnim / [ggplot_vegatex]

from ginger import transparent

const docStr = """
Given some Likelihood File (the output of the likelihood.nim), this script
creates a heatmap of the cluster centers that still remain visible on the plot.
In addition it highligths the gold region on the plot.
"""

type
  ## `ClusterTable` maps (x, y) positions to energies of all clusters at that pixel
  ##
  ## Can be converted to a pure `CountTable[(int, int)]`.
  ClusterTable = Table[(int, int), seq[float]]

  ColorBy = enum
    count, energy

proc initClusterTable(): ClusterTable = initTable[(int, int), seq[float]]()

proc toCountTable(tab: ClusterTable): CountTable[(int, int)] =
  result = initCountTable[(int, int)]()
  for (k, v) in pairs(tab):
    for i in 0 ..< v.len:
      result.inc(k)

iterator extractClusters(h5f: var H5File, chip: int): (int, seq[float], seq[float], seq[float]) =
  const attrName = "Total number of cluster"
  for num, grp in runs(h5f, likelihoodBase()):
    var mgrp = h5f[grp.grp_str]
    var targetChip: int
    if chip >= 0: # use given chip if any (default -1)
      targetChip = chip
    else: # try to determine the center chip
      try: #
        targetChip = mgrp.attrs["centerChip", int]
      except KeyError: # else use 3
        echo "WARNING: could not find `centerChip` attribute. Will assume 3!"
        targetChip = 3
    var chipGrp: H5Group
    try:
      chipGrp = h5f[(grp / "chip_" & $targetChip).grp_str]
    except KeyError:
      echo "INFO: no events left in run number " & $num & " for chip " & $targetChip
      continue
    let cX = h5f[chipGrp.name / "centerX", float]
    let cY = h5f[chipGrp.name / "centerY", float]
    let cE = h5f[chipGrp.name / "energyFromCharge", float]
    let totalClusters = chipGrp.attrs[attrName, int]
    yield (totalClusters, cX, cY, cE)

proc goldRegionOutline(): DataFrame =
  ## returns the outline of the gold region
  let min = (256.0 * 4.5 / 14.0).round.int
  let max = (256.0 * 9.5 / 14.0).round.int

  let xs = [min, max, max, min, min]
  let ys = [min, min, max, max, min]
  result = toDf(xs, ys)

  ## This was for the very old version using `plotly` that took a 256x256 tensor
  when false:
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

## XXX: Handle the noisy pixels better!
when true:
  const noisyPixels = [
    # "Main" noise cluster in Run-2
    (64, 109),
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
    (66, 107),
    (67, 107),
    (66, 111),
    (69, 110),
    # Secondary cluster in Run-2
    (75, 197),
    (76, 197),
    (74, 196),
    (75, 196),
    (76, 196),
    (76, 195),
    (77, 196),
    # Deich bottom left
    (1 , 83),
    (2 , 83),
    (3 , 83),
    (4 , 83),
    (5 , 83),
    (6 , 83),
    # Deich middle left
    (1 , 166),
    (2 , 166),
    (3 , 166),
    (4 , 166),
    (5 , 166),
    # Deich middle right
    (251 , 166),
    (252 , 166),
    (253 , 166),
    (254 , 166),
    (255 , 166)]

when false:
  var noisyPixels = newSeq[(int, int)]()
  for x in 150 .. 250:
    for y in 130 .. 162:
      noisyPixels.add (x, y)

  for x in 125 .. 135:
    for y in 110 .. 120:
      noisyPixels.add (x, y)

proc readClusters(file: string, cTab: var ClusterTable,
                  filterNoisyPixels: bool,
                  energyMin, energyMax: float,
                  chip: int,
                  switchAxes: bool): int =
  ## Fills the `cTab` of those pixels that are hit and also returns the total number of
  ## clusters on the target chip.
  echo "reading: ", file
  var h5f = H5open(file, "r")
  func toPixel(s: float): int = min((255.0 * s / 14.0).round.int, 255)
  for totalNum, centerX, centerY, energy in extractClusters(h5f, chip):
    result += totalNum # count total number of events
    for (cX, cY, cE) in zipEm(centerX, centerY, energy):
      # drop clusters with energy too large
      if   energyMin > 0.0 and cE < energyMin: continue
      elif energyMax > 0.0 and cE > energyMax: continue

      #echo (pX, pY) in noisyPixels, " pos ", (pX, pY)
      # Need original position for noisy pixels!
      let posOriginal = (cX.toPixel, cY.toPixel)
      let pos = if not switchAxes: posOriginal
                else: (cY.toPixel, cX.toPixel)
      if pos notin cTab:
        cTab[pos] = newSeq[float]()

      if filterNoisyPixels and posOriginal notin noisyPixels:
        cTab[pos].add ce
      elif not filterNoisyPixels:
        cTab[pos].add ce
  doAssert h5f.close() >= 0

proc readClusters(file: string, filterNoisyPixels: bool,
                  energyMin, energyMax: float, chip: int,
                  switchAxes: bool): (int, ClusterTable) =
  var cTab = initClusterTable()
  let totalNum = file.readClusters(cTab, filterNoisyPixels, energyMin, energyMax, chip, switchAxes)
  result = (totalNum, cTab)

proc readClusters(files: seq[string], filterNoisyPixels: bool,
                  energyMin, energyMax: float, chip: int,
                  switchAxes: bool): (int, ClusterTable) =
  var cTab = initClusterTable()
  var totalNum = 0
  for f in files:
    totalNum += f.readClusters(cTab, filterNoisyPixels, energyMin, energyMax, chip, switchAxes)
  result = (totalNum, cTab)

proc toDf(cTabTup: (int, ClusterTable)): DataFrame =
  let (totalNum, cTab) = cTabTup
  var
    cX, cY: seq[int]
    cE: seq[float]
  for (pos, energies) in pairs(cTab):
    for E in energies:
      cX.add pos[0]
      cY.add pos[1]
      cE.add E
  ## XXX: `toDf` is broken here due to proc of same name!!!
  result = seqsToDf({"x" : cX, "y" : cY, "Energy [keV]" : cE, "TotalNumber" : totalNum})

proc toDf(cTabTup: (int, CountTable[(int, int)])): DataFrame =
  let (totalNum, cTab) = cTabTup
  var
    cX, cY, cC: seq[int]
  for (pos, count) in pairs(cTab):
    cX.add pos[0]
    cY.add pos[1]
    cC.add count
  ## XXX: `toDf` is broken here due to proc of same name!!!
  result = seqsToDf({"x" : cX, "y" : cY, "count" : cC, "TotalNumber" : totalNum})
    .arrange("count", order = SortOrder.Ascending)

proc toDf(cTab: CountTable[(int, int)]): DataFrame =
  result = (-1, cTab).toDf()

proc toDf(cTab: ClusterTable): DataFrame =
  result = (-1, cTab).toDf()

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

proc plotClusters(df: DataFrame, names: seq[string], useTikZ: bool, zMax: float,
                  colorBy: ColorBy,
                  energyText: bool,
                  energyTextRadius: float,
                  suffix, title, outpath, axionImage: string,
                  scale: float,
                  preliminary: bool,
                  showGoldRegion: bool) =
  let outpath = if outpath.len > 0: outpath else: "plots"
  var colorCol: string
  let totalEvs = df.len
  var df = df
  # Set the `color` aesthetic and potentially convert input DF to count form
  case colorBy
  of energy:
    colorCol = "Energy [keV]"
  else:
    colorCol = "count"
    # generate a version with unique positions but a `count` column
    var countGroups = 0
    for (tup, subDf) in groups(df.group_by(["x", "y"])):
      inc countGroups
    df = df.group_by(["x", "y"]).summarize(f{"count" << len(col("Energy [keV]"))})

  createDir(outpath)
  let fname = &"{outpath}/background_cluster_centers{suffix}.pdf"
  if names.len == 0 or names.deduplicate.len == 1:
    let maxCount = if colorBy == count: df[colorCol, int].max
                   else: 0
    echo df
    var plt = ggplot(df, aes("x", "y")) +
      xlim(0, 256) + ylim(0, 256) +
      xlab("x [Pixel]") + ylab("y [Pixel]")
      #margin(top = 1.75 * scale)

    if axionImage.len > 0:
      var customInferno = inferno()
      customInferno.colors[0] = 0 # transparent
      let rasterData = readCsv(axionImage)
      let zCol = if "z" in rasterData: "z" else: "photon flux"
      plt = plt +
        geom_raster(data = rasterData, aes = aes("x", "y", fill = zCol), alpha = 0.3) +
        minorGridLines() +
        scale_fill_gradient(customInferno)

    if colorBy == energy and energyText:
      var dfText = df
      if energyTextRadius > 0.0:
        dfText = df.filter(f{float -> bool: (128.0 - `x`)^2 + (128.0 - `y`)^2 < energyTextRadius^2})
      plt = plt + geom_text(data = dfText,
                            aes = aes(y = f{`y` + 3}, text = "Energy [keV]"),
                            font = font(8.0, alignKind = taLeft))

    # Add the main point geom
    if colorBy == count and  maxCount < 10:
      plt = plt + geom_point(aes = aes(color = factor(colorCol)), size = some(scale * 1.0))
    else:
      plt = plt + geom_point(aes = aes(color = colorCol), size = some(scale * 1.0)) +
        scale_color_continuous(scale = (low: 0.0, high: zMax))

    if preliminary:
      let red = color(1.0, 0.0, 0.0)
      plt = plt + annotate("Preliminary",
                           backgroundColor = transparent,
                           x = 25, y = 70, rotate = 45.0,
                           font = font(32.0, color = red))

    if showGoldRegion:
      plt = plt + geom_line(data = goldRegionOutline(), aes = aes("xs", "ys"), color = "red")

    if not useTikZ:
      echo "[INFO]: Saving plot to ", fname
      plt + theme_scale(scale, family = "serif") +
        ggtitle(title & &". Total # clusters = {totalEvs}") +
        ggsave(fname, width = 640.0 * scale, height = 480 * scale)
        #ggsave(fname, width = 640, height = 480)#width = 1200, height = 800)
    else:
      #let fname = &"/home/basti/phd/Figs/backgroundClusters/background_cluster_centers{suffix}"
      echo "[INFO]: Saving plot to ", fname
      plt + ggtitle(title & r". Total \# clusters = " & $totalEvs) +
        theme_scale(scale) +
        ggsave(fname, width = 800, height = 600, useTeX = true, standalone = true)
        #ggsave(fname, width = 600, height = 450, useTeX = true, standalone = true)
        #ggvegatex(fname)
  else:
    var totalEvs = newSeq[string]()
    for tup, subDf in groups(df.group_by("Type")):
      let numEvs = if colorBy == count: subDf[colorCol, int].sum()
                   else: subDf.len
      totalEvs.add &"{tup[0][1]}: {numEvs}"
    #let fname = &"/home/basti/org/Figs/statusAndProgress/IAXO_TDR/background_cluster_centers{suffix}.pdf"
    echo "[INFO]: Saving plot to ", fname
    let ticks = arange(0, 260, 5).mapIt(it.float)
    ggplot(df, aes("x", "y", color = colorCol)) +
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

proc plotSuppression(cTab: CountTable[(int, int)], totalNum, tiles: int,
                     outpath: string,
                     showGoldRegion: bool) =
  ## Plots a tilemap of background suppressions. It uses `tiles` elements in each axis.
  proc toTensor(cTab: CountTable[(int, int)]): Tensor[int] =
    result = zeros[int]([256, 256])
    for tup, val in cTab:
      let
        x = tup[0]
        y = tup[1]
      result[x, y] = val
  let countT = cTab.toTensor()
  let numPerTile = totalNum / (tiles * tiles) ## Note: This assumes the original background is homogeneous!
  let step = 256.0 / tiles.float
  var
    xI = newSeq[int]()
    yI = newSeq[int]()
    cI = newSeq[int]()
    sI = newSeq[float]()
  for i in 0 ..< tiles:
    let xStart = (i.float * step).floor.int
    let xStop = ((i+1).float * step).floor.int
    for j in 0 ..< tiles:
      let yStart = (j.float * step).floor.int
      let yStop = ((j+1).float * step).floor.int
      var counts = 0
      for x in xStart ..< xStop:
        for y in yStart ..< yStop:
          counts += countT[x, y]
      xI.add((i.float * step).floor.int)
      yI.add((j.float * step).floor.int)
      cI.add counts
      sI.add (numPerTile.float / counts.float)

  let dfTile = seqsToDf(xI, yI, cI, sI)
  echo dfTile

  let size = step.ceil
  let outfile = if outpath.len > 0: outpath / "background_suppression_tile_map.pdf"
                else: "plots/background_suppression_tile_map.pdf"
  ## XXX: FIX UP title for TeX backend (unicode) & make # of total cluster dependent on code, not
  ## hardcoded!
  ggplot(dfTile, aes("xI", "yI", fill = "sI", width = size, height = size)) +
    geom_tile() +
    geom_text(aes = aes(x = f{`xI` + size / 2.0}, y = f{`yI` + size / 2.0}, text = "sI"), color = "white") +
    xlim(0, 256) + ylim(0, 256) +
    xlab("x [Pixel]") + ylab("y [Pixel]") +
    margin(top = 1.75) +
    ggtitle(&"Local background suppression compared to {totalNum.float:.2g} raw clusters") +
    ggsave(outfile, useTeX = true, standalone = true)

proc main(
  files: seq[string], suffix = "", title = "",
  names: seq[string] = @[], # can be used to create facet plot. All files with same names will be stacked
  filterNoisyPixels = false,
  energyMin = 0.0, energyMax = 0.0,
  useTikZ = false,
  writeNoisyClusters = false,
  energyText = false,
  energyTextRadius = -1.0, ## Radius in which around the center to print energy as text
  colorBy: ColorBy = count,
  zMax = 5.0,
  threshold = 2,
  chip = -1, # chip can be used to overwrite center chip reading
  tiles = 7,
  outpath = "",
  axionImage = "", # "/home/basti/org/resources/axion_images/axion_image_2018_1487_93_0.989AU.csv"
  preliminary = false,
  backgroundSuppression = false,
  showGoldRegion = false,
  scale = 1.0, # Scale the output image and all texts etc. by this amount. Default is 640x480
  switchAxes = false # if true, will replace X by Y (to effectively rotate the clusters into CAST setup)
     ) =

  if energyText and colorBy == count:
    raise newException(ValueError, "`energyText` incompatible with `colorBy = count`.")

  # Read the data
  let (totalNum, cTab) = readClusters(files, filterNoisyPixels, energyMin, energyMax, chip, switchAxes)
  if writeNoisyClusters:
    cTab.toCountTable.writeNoisyClusters(threshold)

  var df = newDataFrame()
  if names.len == 0 or names.deduplicate.len == 1:
    echo files
    df = cTab.toDf()
  else:
    doAssert names.len == files.len, "Need same number of names as input files!"
    for i in 0 ..< names.len:
      var dfLoc = readClusters(@[files[i]], filterNoisyPixels, energyMin, energyMax, chip, switchAxes).toDf()
      dfLoc["Type"] = names[i]
      df.add dfLoc

  plotClusters(df, names, useTikZ, zMax, colorBy, energyText, energyTextRadius, suffix, title, outpath, axionImage, scale, preliminary, showGoldRegion)
  # `df.len` is total number clusters
  if backgroundSuppression:
    doAssert names.len == 0, "Suppression plot when handing multiple files that are not combined not supported."
    plotSuppression(cTab.toCountTable(), totalNum, tiles, outpath, showGoldRegion)

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
