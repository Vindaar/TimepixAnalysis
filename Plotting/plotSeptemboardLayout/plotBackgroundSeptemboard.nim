import std / [strutils, os, strformat]
import ingrid / [ingrid_types, tos_helpers]
import ggplotnim, nimhdf5, seqmath

## Just a test of plotting the septemboard layout

let UseRealLayout = true
const chipOutlineX = @[
  (x: 0,   yMin: 0, yMax: 255),
  (x: 255, yMin: 0, yMax: 255)
]
const chipOutlineY = @[
  (y: 0  , xMin: 0, xMax: 255),
  (y: 255, xMin: 0, xMax: 255)
]

const Pix = 900
proc readFile(file: string, onlyCount: bool, chargeCut: int): Tensor[float] =
  result = zeros[float]([Pix, Pix]) # real layout fits into (900, 900) tensor
  withH5(file, "r"):
    let fileInfo = getFileInfo(h5f)
    for run in fileInfo.runs:
      let df = getSeptemDataFrame(h5f, run, ToT = true, charge = false)
        .chpPixToRealPix(realLayout = true)
      let
        xs = df["x", int]
        ys = df["y", int]
        zs = df["ToT", int]
      if onlyCount:
        for i in 0 ..< xs.len:
          result[xs[i], ys[i]] += 1.0
      else:
        for i in 0 ..< xs.len:
          if zs[i] < chargeCut:
            result[xs[i], ys[i]] += zs[i].float
      #break

proc toLongDf(t: Tensor[float]): DataFrame =
  var
    xs = newSeq[int]()
    ys = newSeq[int]()
    zs = newSeq[float]()
  for x in 0 ..< Pix:
    for y in 0 ..< Pix:
      xs.add x
      ys.add y
      zs.add t[x, y]
  result = toDf(xs, ys, zs)

proc main(file: string,
          onlyCount = true,
          chargeCut = int.high, # cut charge values above this (highest ToT is 11810)
          quantile = 95,
          title = "",
          zCol = "",
          outfile = "/tmp/occupancy_map.pdf",
         ) =
  let zCol = if zCol.len > 0: zCol
             elif onlyCount: "Counts"
             else: "Charge"
  # 0. generate a temporary file name (to avoid rereading data when changing the plot)
  let tmpFile = "/tmp" / file.extractFilename.replace(".h5", "") & &"_onlyCount_{onlyCount}_chargeCut_{chargeCut}.csv"
  var df = newDataFrame()
  if existsFile tmpFile:
    df = readCsv(tmpFile)
  else:
    # 1. read input data (x, y, ToT) into septem layout tensor & compute occupancy
    let data = readFile(file, onlyCount, chargeCut)
    # 2. merge all into a long DF

    df = toLongDf(data)
      .rename(f{zCol <- "zs"})
    echo df
    # write temp CSV
    df.writeCsv(tmpFile)

  # 3. set up custom color scale
  var customColors = viridis()
  customColors.colors[0] = 0
  # 4. title suffix
  let cs = if onlyCount: "(by count)" else: "(by charge)"
  let title = title & " " & cs & " cut to " & $quantile & "-th percentile"
  # 5. plot into single raster plot with septemboard layout
  let zMax = df[zCol, float].percentile(quantile)
  var plt = ggplot(df, aes()) +
    geom_raster(aes = aes(xs, ys, fill = zCol)) +
    coord_fixed(1.0) +
    xlab("x [pixel]") + ylab("y [pixel]") +
    scale_x_continuous() + scale_y_continuous() +
    scale_fill_gradient(customColors, dataScale = (low: 0.0, high: zMax)) +
    ggtitle(title) +
    themeLatex(fWidth = 0.9, width = 600, baseTheme = singlePlot)
  plt.addSeptemboardOutline(useRealLayout = true)
  plt + ggsave(outfile, width = 800)

when isMainModule:
  import cligen
  dispatch main
