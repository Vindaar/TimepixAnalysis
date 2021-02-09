import nimhdf5, ggplotnim, os, strformat, strutils, sequtils, algorithm, seqmath, arraymancer
import ingrid / tos_helpers, times

import cligen

import numericalnim except linspace

type
  HistTuple = tuple[bins: seq[float64], hist: seq[float64]]

  YearKind = enum
    yr2014 = "2014"
    yr2018 = "2018"

proc splitSeq[T, U](s: seq[seq[T]], dtype: typedesc[U]): (seq[U], seq[U]) =
  ## splits a (N, 2) nested seq into two seqs
  result[0] = newSeq[dtype](s.len)
  result[1] = newSeq[dtype](s.len)
  for i in 0..s.high:
    result[0][i] = s[i][0].U
    result[1][i] = s[i][1].U

proc readRefDsets(refFile: string,
                  dkKind: InGridDsetKind,
                  yearKind: YearKind): DataFrame =
  ## reads the reference datasets from the `refFile` and returns them.
  var h5ref = H5open(refFile, "r")
  # create a table, which stores the reference datasets from the ref file
  let dset = dkKind.toDset(fkTpa)
  var tab = initTable[string, HistTuple]()
  let xrayRef = getXrayRefTable()
  let energies = getXrayFluorescenceLines()
  let ranges = getEnergyBinning()
  for idx in 0 ..< xray_ref.len:
    let dset_name = xray_ref[idx]
    var data = h5ref[(dset_name / dset).dset_str]
    tab[dset_name] = data.readAs(float64).reshape2D(data.shape).splitSeq(float64)
    var dfDset = seqsToDf({ "Bins" : tab[dset_name].bins, "Hist" : tab[dset_name].hist })
    dfDset["Dset"] = constantColumn(dset_name, dfDset.len)

    dfDset["Energy"] = constantColumn(energies[idx], dfDset.len)
    dfDset["yMin"] = if idx == 0: constantColumn(0, dfDset.len)
                     else: constantColumn(ranges[idx-1], dfDset.len)
    dfDset["yMax"] = if idx == ranges.high: constantColumn(10.0, dfDset.len)
                     else: constantColumn(ranges[idx], dfDset.len)
    result.add dfDset
  result["runType"] = constantColumn("CDL", result.len)

proc filterKdeMorph(df: DataFrame): DataFrame =
  let xrayRef = getXrayRefTable()
  let energies = getXrayFluorescenceLines()
  let eKde = df.filter(f{`runType` == "Morph"})["Energy", float].toRawSeq
  result = df
  for idx, e in energies:
    let eEx = eKde[min(eKde.lowerBound(e), eKde.high)]
    var modDf = df.filter(f{`runType` == "Morph"}, f{float: `Energy` == eEx})
      .mutate(f{float -> float: "Hist" ~ `Hist` / sum(df["Hist"])})
    let name: string = xrayRef[idx]
    modDf["Dset"] = constantColumn(name, modDf.len)
    result.add modDf
  result = result.filter(f{`Dset` != "Morph"})

proc plotRef(df: DataFrame,
             dset: string,
             refFile: string, yearKind: YearKind,
             isKde: static bool = false) =
  let xrayRef = getXrayRefTable()
  var labelOrder = initTable[Value, int]()
  for idx, el in xrayRef:
    labelOrder[%~ el] = idx
  #labelOrder[%~ "Morph"] = labelOrder.len
  ggplot(df.filter(f{string: `runType` == "Morph"}), aes("Bins", "Energy", fill = "Hist")) +
    #ggridges("Dset", overlap = 1.0, labelOrder = labelOrder) +
    #facet_wrap("Dset", scales = "free") +
    geom_raster() +
    #geom_histogram(stat = "identity", position = "identity",
    #               alpha = some(0.5)) +
    ggtitle(&"{dset} of reference file, year: {yearKind}") +
    ggsave(&"out/{dset}_ridgeline_tile_{refFile.extractFilename}_{yearKind}.pdf",
            width = 800, height = 480)

  when isKde:
    let df = filterKdeMorph(df)
  #let df = df.filter(f{float: `Energy` in energies or `Energy` in eExctract})
  #echo df.filter(f{float: `Energy` in energies}, f{string: `runType` == "Morph"}).pretty(-1)

  ggplot(df, aes("Bins", "Hist", fill = "runType")) +
    ggridges("Dset", overlap = 1.0, labelOrder = labelOrder) +
    #facet_wrap("Dset", scales = "free") +
    geom_histogram(stat = "identity", position = "identity",
                   alpha = some(0.5)) +
    ggtitle(&"{dset} of reference file, year: {yearKind}") +
    ggsave(&"out/{dset}_ridgeline_{refFile.extractFilename}_{yearKind}.pdf",
            width = 800, height = 480)

proc normalize(df: DataFrame): DataFrame =
  result = df
  result = result.group_by("Dset")
    .mutate(f{float -> float: "Hist" ~ `Hist` / sum(df["Hist"])})

import sugar
proc morph(df: DataFrame, energy: float, offset = 1): Tensor[float] =
  ## generates a distribution for `col` at `energy`
  let lineEnergies = getXrayFluorescenceLines()
  let idx = max(lineEnergies.lowerBound(energy) - 1, 0)
  # need idx and idx+offset
  dump idx
  let xrayRef = getXrayRefTable()
  let refLow = xrayRef[idx]
  let refHigh = xrayRef[idx+offset]
  echo refLow, "  ", refHigh, " at ", energy
  let refLowT = df.filter(f{string -> bool: `Dset` == refLow and `runType` == "CDL"})["Hist", float]
  let refHighT = df.filter(f{string -> bool: `Dset` == refHigh and `runType` == "CDL"})["Hist", float]
  let bins = df.filter(f{string -> bool: `Dset` == refLow and `runType` == "CDL"})["Bins", float]
  result = zeros[float](refLowT.size.int)
  # walk over each bin and compute linear interpolation between
  let Ediff = abs(lineEnergies[idx] - lineEnergies[idx+offset])
  dump Ediff
  dump (abs(lineEnergies[idx] - energy)) / Ediff
  dump (abs(lineEnergies[idx+offset] - energy)) / Ediff
  for i in 0 ..< bins.size:
    result[i] = refLowT[i] * (1 - (abs(lineEnergies[idx] - energy)) / Ediff) +
      refHighT[i] * (1 - (abs(lineEnergies[idx+offset] - energy)) / Ediff)

proc linKernel(x, x_i, bw: float): float {.inline.} =
  result = (1 - abs(x - x_i) / bw)

proc morphKde(df: DataFrame, dset: string): DataFrame =
  ## morph the given df from the reference datasets bin by bin using a KDE
  #let bins = df.filter(f{string -> bool: `Dset` == refLow and `runType` == "CDL"})["Bins", float]
  #let df = df.filter(f{string: `Dset` == dset})
  result = df.select(@["Bins", "Hist", "Dset", "Energy", "runType"])
  for tup, subDf in groups(df.group_by("Bins")):
    var energies = subDf["Energy", float]
    var hist = subDf["Hist", float]

    let kd = kde(energies, weights = hist, bw = 0.3)#, kernel = linKernel, kernelKind = knCustom, cutoff = 3.0) #, bw = 0.3)
    let es = linspace(min(energies), max(energies), 1000)
    var dfLoc = seqsToDf({"Energy" : es, "Hist" : kd })
    let bin = subDf["Bins", float][0]
    dfLoc["Bins"] = constantColumn(bin, dfLoc.len)
    dfLoc["Dset"] = constantColumn("Morph", dfLoc.len)
    dfLoc["runType"] = constantColumn("Morph", dfLoc.len)

    when false:
      ggplot(dfLoc, aes("Energy", "Hist")) +
        geom_line() + ggsave(&"out/line_{dset}_{bin}.pdf")

    result.add dfLoc

proc morphSpline(df: DataFrame, exclude: int): DataFrame =
  ## morph the given df from the reference datasets bin by bin using a spline interpolation
  result = df.select(@["Bins", "Hist", "Dset", "Energy", "runType"])
  for tup, subDf in groups(df.group_by("Bins")):
    let subDf = subDf.arrange("Energy")

    var energies = newSeq[float]() # = subDf["Energy", float]
    var hist = newSeq[float]()# = subDf["Hist", float]
    for idx in 0 ..< subDf.len:
      if idx == exclude:
        echo "Excluding index ", idx, " which is ", subDf["Energy", float][idx]
        continue
      energies.add subDf["Energy", float][idx]
      hist.add subDf["Hist", float][idx]
    echo energies
    #echo hist
    let kd = newCubicSpline(energies, hist)
    let es = linspace(min(energies), max(energies), 1000)
    var dfLoc = seqsToDf({"Energy" : es, "Hist" : es.mapIt(kd.eval(it)) })
    let bin = subDf["Bins", float][0]
    dfLoc["Bins"] = constantColumn(bin, dfLoc.len)
    dfLoc["Dset"] = constantColumn("Morph", dfLoc.len)
    dfLoc["runType"] = constantColumn("Morph", dfLoc.len)

    when false:
      ggplot(dfLoc, aes("Energy", "Hist")) +
        geom_line() + ggsave(&"out/line_{dset}_{bin}.pdf")

    result.add dfLoc

proc filterSplineToEnergy(df: DataFrame, idx: int, E: float): DataFrame =
  let energies = df.filter(f{`runType` == "Morph"})["Energy", float].toRawSeq
  let eExclude = energies[min(energies.lowerBound(E), energies.high)]
  result = df.filter(f{`runType` == "Morph"}, f{float: `Energy` == eExclude})
    .mutate(f{float -> float: "Hist" ~ `Hist` / sum(df["Hist"])})
  let xrayRef = getXrayRefTable()
  let name = xrayRef[idx]
  result["Dset"] = constantColumn(name, result.len)

proc getMorphedDfSpline(df: DataFrame): DataFrame =
  let energies = getXrayFluorescenceLines()
  result = df
  for idx in 1 ..< energies.high: # exclude first and last
    let dfSpl = morphSpline(df, idx)
      .filterSplineToEnergy(idx, energies[idx])
    result.add dfSpl

proc morphDf(df: DataFrame, idx: int, energy: float): DataFrame =
  ## returns morphed as added thing in df
  result = df
  #const idx = 1
  let xrayRef = getXrayRefTable()
  let energies = getEnergyBinning()
  let res = morph(df, energy, offset = 2)
  let bins = df["Bins", float][0 ..< res.size]
  var dfMorph = seqsToDf({"Bins" : bins, "Hist" : res})
  dfMorph["Dset"] = constantColumn(xrayRef[idx], dfMorph.len)
  dfMorph["runType"] = constantColumn("Morph", dfMorph.len)
  dfMorph["Energy"] = constantColumn(energy, dfMorph.len)
  dfMorph["yMin"] = constantColumn(energies[idx-1], dfMorph.len)
  dfMorph["yMax"] = constantColumn(energies[idx], dfMorph.len)
  result.add dfMorph

proc getMorphedDf(df: DataFrame): DataFrame =
  let energies = getXrayFluorescenceLines()
  result = df
  for idx in 1 ..< energies.high: # exclude first and last
    result = result.morphDf(idx, energies[idx])
    echo result

proc getInterpolatedDf(df: DataFrame, num = 1000): DataFrame =
  ## returns a DF with `num` interpolated distributions using next neighbors
  ## for linear interpolation (to be used to plot as a raster)
  let energiesLines = getXrayFluorescenceLines()
  result = df.select(["Bins", "Hist", "runType", "Energy", "Dset"])
  let energies = linspace(energiesLines[0], energiesLines[^1], num)
  let xrayRef = getXrayRefTable()
  #let energies = getEnergyBinning()
  var binWidth: float
  for idx, E in energies:
    let res = morph(df, E, offset = 1)
    let bins = df["Bins", float][0 ..< res.size]
    if binWidth == 0.0:
      binWidth = bins[1] - bins[0]
    var dfMorph = seqsToDf({"Bins" : bins, "Hist" : res})
    dfMorph["runType"] = constantColumn("Morph", dfMorph.len)
    dfMorph["Energy"] = constantColumn(E, dfMorph.len)
    dfMorph["Dset"] = constantColumn("None", dfMorph.len)
    result.add dfMorph
    echo result
  result["binWidth"] = constantColumn(binWidth, result.len)
  result["binHeight"] = constantColumn(energies[1] - energies[0], result.len)

proc plotTile(df: DataFrame, dset: string) =
  let binWidth = block:
                   let t = df["Bins", float]
                   t[1] - t[0]
  var plt = ggplot(df, aes("Bins", "Energy")) +
    geom_tile(aes = aes(fill = "Hist", yMin = "yMin", yMax = "yMax",
                        xMin = f{`Bins` - binWidth / 2.0}, xMax = f{`Bins` + binWidth / 2.0})) +
    geom_line(aes = aes("Bins", "Energy", color = "Dset")) +
    scale_y_continuous() +
    ylim(0, 10) +
    #xlim(0, 15) +
    ggtitle(&"2D heatmap plot of CDL data for {dset}")
  let lines = getXrayFluorescenceLines()
  var lineDf = newDataFrame()
  for idx in 0 ..< lines.len:
    var bins = df.filter(f{string -> bool: `Dset` == "Cu-Ni-15kV" and `runType` == "CDL"})["Bins", float]
    let line = lines[idx]
    let minVal = bins.min - (bins.max - bins.min) * 0.05
    let maxVal = bins.min
    echo minVal
    echo maxval
    let lStr = &"{line:.3g}"
    plt = plt + geom_text(aes = aes(x = minVal, y = line, text = lStr),
                          color = some(parseHex("FF0000")),
                          size = some(4.0))

  plt + ggsave(&"out/cdl_as_tile_{dset}.pdf")

proc plotInterpolatedRaster(df: DataFrame, dset: string, interpolated = false) =
  var dfMorph: DataFrame
  if not interpolated:
    dfMorph = df.getInterpolatedDf.filter(f{string: `runType` == "Morph"})
  else:
    dfMorph = df.filter(f{string: `runType` == "Morph"})
  var plt = ggplot(dfMorph, aes("Bins", "Energy", fill = "Hist")) +
    geom_raster(aes = aes(width = "binWidth", height = "binHeight")) +
    scale_y_continuous() +
    ylim(0, 10) +
    #xlim(0, 15) +
    ggtitle(&"2D heatmap plot of linearly interpolated CDL data for {dset}")
  when false:
    let lines = getXrayFluorescenceLines()
    var lineDf = newDataFrame()
    for idx in 0 ..< lines.len:
      var bins = df.filter(f{string -> bool: `Dset` == "Cu-Ni-15kV" and `runType` == "CDL"})["Bins", float]
      let line = lines[idx]
      let minVal = bins.min - (bins.max - bins.min) * 0.05
      let maxVal = bins.min
      let lStr = &"{line:.3g}"
      plt = plt + geom_text(data = newDataFrame(), aes = aes(x = minVal, y = line, text = lStr),
                            color = some(parseHex("FF0000")),
                            size = some(4.0))
  plt + ggsave(&"out/cdl_as_raster_interpolated_{dset}.pdf")

proc main(files: seq[string],
          refFile: string = "/mnt/1TB/CAST/CDL_2019/XrayReferenceFile2018.h5",
          fkKind: FrameworkKind = fkTpa) =

  const dkKinds = [igEccentricity, igLengthDivRmsTrans, igFractionInTransverseRms]
  var df = newDataFrame()
  for dkKind in dkKinds:
    let dfRef = readRefDsets(refFile, dkKind, yr2018)
      .normalize
    let dset = dkKind.toDset(fkKind)
    let res = morphKde(dfRef, dset)
    #let res = getMorphedDfSpline(dfRef)
    echo res
    #var res = getMorphedDf(dfRef)
    plotRef(res,
            dset,
            refFile, yr2018, isKde = true)
    #res["Variable"] = constantColumn(dset, res.len)
    #
    #ggplot(res, aes("Bins", "Energy", color = "Hist")) +
    #  geom_point(size = some(3.0)) +
    #  scale_y_continuous() +
    #  ggtitle(&"Scatter plot of CDL data for {dset}") +
    #  ggsave(&"out/cdl_as_scatter_{dset}.pdf")
    #
    ##plotTile(res, dset)
    #plotInterpolatedRaster(dfRef, dset, interpolated = true)

  #echo df.pretty(-1)


when isMainModule:
  dispatch main
