import sequtils, plotly, tables, seqmath, nimhdf5, strformat, chroma
import os

const GridColor = color(0.8, 0.8, 0.8, 0.8)
const LegendBg = color(1.0, 1.0, 1.0, 0.5)
const LegendBorder = color(0.6, 0.6, 0.6, 0.6)
const Color1 = color(1.0, 0.0, 102.0 / 256.0)
const Color2 = color(0.0, 153.0 / 256.0, 204 / 256.0)

## This file uses the `XrayReferenceDatasets.h5` file to compare properties
## for different targets (-> energies)

proc makePlot(h5f: var H5FileObj, grp, dset: string, filter: (float, float)): Plot[float32] =
  let data1d = h5f[grp / dset, float32]
  let
    (binsA, countsA) = data1d.reshape2D([data1d.len div 2, 2]).split(SplitSeq.Seq2Col)
    bins = binsA.filterIt(it >= filter[0] and it <= filter[1])
    counts = countsA[0 .. bins.high]
  result = barPlot(bins, counts)
    .title(&"Comparison of distributions for {dset} at different energies")
    .name(&"{grp} {dset}") # "Outer chips # pix hit")
    .xlabel(&"{dset}")
    .ylabel("Counts normalized by PDF")
    .gridColor(GridColor)
    .markerColor(@[Color2], idx = 0)
    .width(800)
    .height(500)
  result.traces[0].autoWidth = true

proc main =

  if paramCount() < 1:
    quit()

  let h5file = paramStr(1)
  var h5f = H5file(h5file, "r")
  defer: discard h5f.close()

  const groups = ["Cu-EPIC-0.9kV", "Mn-Cr-12kV"]
  const dsets = ["excentricity", "fractionwithinrmsy"]
  const filters = [(0.0, 3.0),
                   (0.0, 0.8)]

  var data = initTable[string, seq[float32]]()

  for i, dset in dsets:
    var plts = newSeq[Plot[float32]]()
    for grp in groups:
      plts.add h5f.makePlot(grp, dset, filters[i])
    let plt = plts[0]
      .addTrace(plts[1].traces[0])
      .legendLocation(0.63, 0.95)
      .legendBgColor(LegendBg)
      .legendBorderColor(LegendBorder)
      .legendBorderWidth(1)
      .markerColor(@[Color1], idx = 0)
      .markerColor(@[Color2], idx = 1)
    plt.traces[1].opacity = 0.5
    plt.layout.showLegend = true
    plt.show(&"{dset}_comparison.svg")

when isMainModule:
  main()
