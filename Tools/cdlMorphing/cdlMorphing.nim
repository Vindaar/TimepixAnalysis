import nimhdf5, ggplotnim, os, strformat, strutils, sequtils, algorithm, seqmath, arraymancer
import ingrid / [tos_helpers, ingrid_types], times
import helpers / utils

import cligen

import numericalnim except linspace

type
  YearKind = enum
    yr2014 = "2014"
    yr2018 = "2018"

  MorphTechnique = enum
    mtNone, mtLinear, mtKde, mtSpline

proc readRefDf(ctx: LikelihoodContext): DataFrame =
  ## Get the full reference data as a DF with added `yMin` and `yMax`
  ## containing the energy boundaries of each target / filter combination (`Dset`)
  let df = readRefDsetsDf(ctx)
  result = newDataFrame()
  let ranges = getEnergyBinning()
  let tab = getInverseXrayRefTable()
  for (tup, subDf) in groups(df.group_by("Dset")):
    var dfLoc = subDf
    let idx = tab[tup[0][1].toStr] # get energy index of this target / filter
    dfLoc["yMin"] = if idx == 0: 0.0
                    else: ranges[idx-1]
    dfLoc["yMax"] = if idx == ranges.high: 10.0
                    else: ranges[idx]
    result.add dfLoc
  result["runType"] = "CDL"

proc normalize(df: DataFrame): DataFrame =
  result = df
  result = result.group_by("Dset")
    .mutate(f{float -> float: "Hist" ~ `Hist` / sum(df["Hist"])})

proc morphDf(df: DataFrame, idx: int, energy: float, igDset: string): DataFrame =
  ## returns morphed as added thing in df
  #result = df# .select(@["Bins", "Hist", "Dset", "Energy", "runType"])
  let xrayRef = getXrayRefTable()
  let energies = getEnergyBinning()
  ## NOTE: `morph` nowadays in `likelihood_utils.nim` because we want to check the code we actually
  ## run in practice!
  let (bins, res) = morph(df, energy, offset = 2)
  result = toDf({"Bins" : bins, "Hist" : res})
  result["Dset"] = xrayRef[idx]
  result["runType"] = "Morph"
  result["Energy"] = energy
  result["yMin"] = energies[idx-1]
  result["yMax"] = energies[idx]
  result["Variable"] = igDset

proc getMorphedDf(df: DataFrame, igDset: string): DataFrame =
  let energies = getXrayFluorescenceLines()
  result = df # keep `CDL` data
  for idx in 1 ..< energies.high: # exclude first and last
    result.add df.morphDf(idx, energies[idx], igDset)
  result = result.filter(f{`Variable` == igDset})

proc linKernel(x, x_i, bw: float): float {.inline.} =
  result = (1 - abs(x - x_i) / bw)

proc morphedKde(df: DataFrame, exclude: int, igDset: string): DataFrame =
  ## morph the given df from the reference datasets bin by bin using a KDE
  ##
  ## XXX: This entirey approach is pretty much broken.
  ## When I came up with this I didn't have in mind that we would be working
  ## with the distinct energies as inputs. What we should do to actually see
  ## if a KDE approach can work is to
  ##
  ## 1. take the *raw* clusters that go into the reference DF
  ## 2. split (but not histogram!) them into bins corresponding to the final bins
  ##   of each reference histogram
  ## 3. apply the KDE to all those clusters *along their energy* (their value
  ##   in the property is irrelevant!)
  ## 4. use the resulting KDE to 'predict' any point on the energy scale
  ##
  ## I assume this approach might work reasonably well actually, because then we
  ## have statistics to use a KDE with.
  ##
  ## ```nim
  ##   var eccs = newSeq[float]()
  ##   var ldiv = newSeq[float]()
  ##   var frac = newSeq[float]()
  ##   var Es   = newSeq[float]()
  ##   #withXrayRefCuts(cdlFile, dset, year, energyDset):
  ##   withLogLFilterCuts(cdlFile, dset, year, energyDset, LogLCutDsets):
  ##     # read filtered data
  ##     eccs.add data[igEccentricity][i]
  ##     ldiv.add data[igLengthDivRmsTrans][i]
  ##     frac.add data[igFractionInTransverseRms][i]
  ##     Es.add   data[igEnergyFromCharge][i]
  ## ```
  ##
  ## to get the raw data. Then get the *bins* for each property and walk
  ## all bins, get all clusters for each bin, apply KDE, done.
  ## *However*, with this we still cannot 'exclude' an interval properly
  ## as the KDE doesn't work like that (missing clusters 'in the middle' just
  ## leads to a drop in density). And more importantly we would then
  ## break an implicit assumption we make, namely that the actual cluster energy
  ## is irrelevant, only the associated energy of the fluorescence line counts!
  # Filter out only subset for this `IngridDsetKind`
  let df = df.filter(f{string: `Variable` == igDset})
  result = df.select(@["Bins", "Hist", "Dset", "Energy", "runType", "Variable"])
  for tup, subDf in groups(df.group_by("Bins")):
    let subDf = subDf.arrange("Energy", SortOrder.Ascending)
    var energies = newSeq[float]()
    var hist = newSeq[float]()

    for idx in 0 ..< subDf.len:
      if idx == exclude:
        echo "Excluding index ", idx, " which is ", subDf["Energy", float][idx]
        continue
      energies.add subDf["Energy", float][idx]
      hist.add subDf["Hist", float][idx]

    let kd = kde(energies.toTensor, weights = hist.toTensor, bw = 0.3)#, kernel = linKernel, kernelKind = knCustom, cutoff = 3.0) #, bw = 0.3)
    let es = linspace(min(energies), max(energies), 1000)
    var dfLoc = toDf({"Energy" : es, "Hist" : kd })
    let bin = subDf["Bins", float][0]
    dfLoc["Dset"] = "None"
    dfLoc["Bins"] = bin
    dfLoc["Variable"] = igDset
    dfLoc["runType"] = "Morph"

    when false:
      ggplot(dfLoc, aes("Energy", "Hist")) +
        geom_line() + ggsave(&"out/line_{igDset}_{bin}.pdf")

    result.add dfLoc

proc filterKdeMorph(df: DataFrame, idx: int): DataFrame =
  let xrayRef = getXrayRefTable()
  let energies = getXrayFluorescenceLines()
  let energy = energies[idx]
  let eKde = df.filter(f{`runType` == "Morph"})["Energy", float].toRawSeq
  let eEx = eKde[min(eKde.lowerBound(energy), eKde.high)]
  result = df.filter(f{`runType` == "Morph"}, f{float: `Energy` == eEx})
    .mutate(f{float -> float: "Hist" ~ `Hist` / sum(df["Hist"])})
  result["Dset"] = xrayRef[idx]
  result["yMin"] = energies[idx-1]
  result["yMax"] = energies[idx]

proc getMorphedDfKde(df: DataFrame, igDset: string): DataFrame =
  result = df
  let energies = getXrayFluorescenceLines()
  for idx in 1 ..< energies.high: # leave out first and last
    let dfM = morphedKde(df, idx, igDset)
      .filterKdeMorph(idx)
    result.add dfM

proc morphSpline(df: DataFrame, exclude: int, igDset: string): DataFrame =
  ## morph the given df from the reference datasets bin by bin using a spline interpolation
  result = df.select(@["Bins", "Hist", "Dset", "Energy", "runType", "Variable"])
  for tup, subDf in groups(df.group_by("Bins")):
    let subDf = subDf.arrange("Energy", SortOrder.Ascending)

    var energies = newSeq[float]() # = subDf["Energy", float]
    var hist = newSeq[float]()# = subDf["Hist", float]
    for idx in 0 ..< subDf.len:
      if idx == exclude:
        echo "Excluding index ", idx, " which is ", subDf["Energy", float][idx]
        continue
      energies.add subDf["Energy", float][idx]
      hist.add subDf["Hist", float][idx]
    let kd = newCubicSpline(energies, hist)
    let es = linspace(min(energies), max(energies), 1000)
    var dfLoc = toDf({"Energy" : es, "Hist" : es.mapIt(kd.eval(it)) })
    let bin = subDf["Bins", float][0]
    dfLoc["Bins"] = bin
    dfLoc["Dset"] = "None" # spline interpolation in all energies. `Dset` filled in `filterSplineToEnergy`
    dfLoc["runType"] = "Morph"
    dfLoc["Variable"] = igDset
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
  result["Dset"] = name
  result["yMin"] = energies[idx-1]
  result["yMax"] = energies[idx]

proc getMorphedDfSpline(df: DataFrame, igDset: string): DataFrame =
  let energies = getXrayFluorescenceLines()
  result = df
  for idx in 1 ..< energies.high: # exclude first and last
    let dfSpl = morphSpline(df, idx, igDset)
      .filterSplineToEnergy(idx, energies[idx])
    result.add dfSpl

proc getLinearInterpolatedDf(df: DataFrame, num = 1000): DataFrame =
  ## returns a DF with `num` interpolated distributions using next neighbors
  ## for linear interpolation (to be used to plot as a raster)
  let energiesLines = getXrayFluorescenceLines()
  result = df.select(["Bins", "Hist", "runType", "Energy", "Dset"])
  let energies = linspace(energiesLines[0], energiesLines[^1], num)
  let xrayRef = getXrayRefTable()
  #let energies = getEnergyBinning()
  var binWidth: float
  var dfs = newSeq[DataFrame]()
  for idx, E in energies:
    let (bins, res) = morph(df, E, offset = 1)
    if binWidth == 0.0:
      binWidth = bins[1] - bins[0]
    var dfMorph = toDf({"Bins" : bins, "Hist" : res})
    dfMorph["runType"] = "Morph"
    dfMorph["Energy"] = E
    dfMorph["Dset"] = "None"

    dfs.add dfMorph
  result = assignStack(dfs)
  result["binWidth"] = binWidth
  result["binHeight"] = energies[1] - energies[0]

proc plotRef(df: DataFrame,
             dset: string,
             refFile: string, yearKind: YearKind,
             mtKind: MorphTechnique,
             outpath: string) =
  let xrayRef = getXrayRefTable()
  var labelOrder = initTable[Value, int]()
  for idx, el in xrayRef:
    labelOrder[%~ el] = idx
  #labelOrder[%~ "Morph"] = labelOrder.len
  #ggplot(df.filter(f{string: `runType` == "Morph"}), aes("Bins", "Energy", fill = "Hist")) +
  #  #ggridges("Dset", overlap = 1.0, labelOrder = labelOrder) +
  #  #facet_wrap("Dset", scales = "free") +
  #  geom_raster() +
  #  #geom_histogram(stat = "identity", position = "identity",
  #  #               alpha = some(0.5)) +
  #  ggtitle(&"{dset} of reference file, year: {yearKind}") +
  #  ggsave(&"out/{dset}_ridgeline_tile_{refFile.extractFilename}_{yearKind}.pdf",
  #          width = 800, height = 480)
  #let df = df.filter(f{float: `Energy` in energies or `Energy` in eExctract})
  #echo df.filter(f{float: `Energy` in energies}, f{string: `runType` == "Morph"}).pretty(-1)

  ggplot(df, aes("Bins", "Hist", fill = "runType")) +
    ggridges("Dset", overlap = 1.0, labelOrder = labelOrder) +
    #facet_wrap("Dset", scales = "free") +
    geom_histogram(stat = "identity", position = "identity",
                   alpha = some(0.5)) +
    ggtitle(&"{dset} of reference file, year: {yearKind}") +
    ggsave(&"{outpath}/{dset}_ridgeline_morph_{mtKind}_{refFile.extractFilename}_{yearKind}.pdf",
            width = 800, height = 480)

proc plotTile(df: DataFrame, dset: string, mtKind: MorphTechnique, outpath: string) =
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
    let lStr = &"{line:.3g}"
    plt = plt + geom_text(aes = aes(x = minVal, y = line, text = lStr),
                          color = "#FF0000",
                          size = 4.0)

  plt + ggsave(&"{outpath}/cdl_as_tile_morph_{mtKind}_{dset}.pdf")

proc plotScatter(df: DataFrame, igDset: string, mtKind: MorphTechnique, outpath: string) =
  ggplot(df, aes("Bins", "Energy", color = "Hist")) +
    geom_point(size = 3.0) +
    scale_y_continuous() +
    ggtitle(&"Scatter plot of CDL data for {igDset}") +
    ggsave(&"{outpath}/cdl_as_scatter_morph_{mtKind}_{igDset}.pdf")

proc plotInterpolatedRaster(df: DataFrame, dset: string, mtKind: MorphTechnique, outpath: string) =
  var dfMorph: DataFrame
  ## XXX: The morphing kind here is not really used / does not make much sense as it is.
  ## (only `mtNone` does). If the input data is already morphed, it is not morphed in the
  ## intended way:
  ## we would want to fully interpolate. Instead of excluding one energy each, interpolate
  ## between energies, which we don't do. To make that plot, need to change the `getMorphed*Df`
  ## procs to have the option to fully interpolate in addition.
  case mtKind
  of mtNone: dfMorph = df.getLinearInterpolatedDf.filter(f{string: `runType` == "Morph"})
  else: dfMorph = df.filter(f{string: `runType` == "Morph"})

  let (f, t) = (dfMorph["Bins", float].min, dfMorph["Bins", float].max)
  var dfE = newDataFrame()
  let xrayRef = getXrayRefTable()
  for idx, E in getXrayFluorescenceLines():
    dfE.add toDf({"Bins" : @[f, t], "Energy": @[E, E], "Dset" : xrayRef[idx]})
  var plt = ggplot(dfMorph, aes("Bins", "Energy")) +
    geom_raster(aes = aes(fill = "Hist")) + # aes = aes(width = "binWidth", height = "binHeight")) +
    geom_line(data = dfE, aes = aes("Bins", "Energy", color = "Dset")) +
    scale_y_continuous() +
    ylim(0, 10) +
    ylab("Energy [keV]") + xlab(dset) +
    #xlim(0, 15) +
    ggtitle(&"Linearly interpolated CDL data for {dset}")
  let lines = getXrayFluorescenceLines()
  var lineDf = newDataFrame()
  for idx in 0 ..< lines.len:
    var bins = df.filter(f{string -> bool: `Dset` == "Cu-Ni-15kV" and `runType` == "CDL"})["Bins", float]
    let line = lines[idx]
    let minVal = bins.min - (bins.max - bins.min) * 0.05
    let maxVal = bins.min
    let lStr = &"{line:.3g}"
    plt = plt + geom_text(data = newDataFrame(), aes = aes(x = minVal, y = line, text = lStr),
                          color = "#FF0000",
                          size = 4.0)
  plt + ggsave(&"{outpath}/cdl_as_raster_interpolated_morph_{mtKind}_{dset}.pdf")

proc computeMorphingSystematics(df: DataFrame): Table[string, float] =
  ## Computes statistics for the systematics introduced via morphing.
  ##
  ## We compute the sumo of squared differences between the real line & the morphed
  ## result for the input. The input contains a column `runType` that designates whether
  ## it's real data (CDL) or morphed (morph).
  ## For each target/filter (aside the outer two) we compute the difference in each bin.

  let dfCdl = df.filter(f{`runType` == "CDL"})
  let dfMorph = df.filter(f{`runType` == "Morph"})

  let xrayRef = getXrayRefTable()
  var diffs = initTable[string, float]()
  for idx in 0 ..< xrayRef.len:
    let dset_name = xrayRef[idx]
    if dsetName notin ["Cu-Ni-15kV", "C-EPIC-0.6kV"]:
      # compute differences
      let dfC = dfCdl.filter(f{`Dset` == dsetName})
      let dfM = dfMorph.filter(f{`Dset` == dsetName})

      let cH = dfC["Hist", float]
      let mH = dfM["Hist", float]
      doAssert cH.size == mH.size
      diffs[dsetName] = abs(cH -. mH).sum() / cH.size.float # normalize to number of bins
  echo "Diffs: ", diffs
  result = diffs

proc main(#files: seq[string],
          cdlFile: string = "/mnt/1TB/CAST/CDL_2019/calibration-cdl-2018.h5",
          fkKind: FrameworkKind = fkTpa,
          outpath = "out") =
  const dkKinds = [igEccentricity, igLengthDivRmsTrans, igFractionInTransverseRms]

  # Read a likelihood context for convenience to access the reference data via the
  # `calibration-cdl` file instead of the `XrayReferenceFile`. The former is saner,
  # because we don't rely on some weird state in the reference file.
  # We hand `morphKind = mkLinear` for no real reason. We don't use the morphed fields
  # of the Context here in this code. But we need it to call `readRefDsetsDf`, which
  # returns a DF of the reference data (replacing our previous `readRefDsets` proc.
  let ctx = initLikelihoodContext(cdlFile,
                                  year = yr2018,
                                  energyDset = igEnergyFromCharge,
                                  region = crSilver,
                                  timepix = Timepix1,
                                  morphKind = mkLineart,
                                  useLnLCut = true)

  var df = newDataFrame()
  var diffTab = initTable[string, Table[string, float]]()
  let dfRefAll = readRefDf(ctx)
  for dkKind in dkKinds:
    let igDset = dkKind.toDset(fkKind)
    let dfRef = dfRefAll.filter(f{`Variable` == igDset})
      .normalize

    ## "Interpolated" raster without interpolation
    proc allPlots(df: DataFrame, igDset: string, mt: MorphTechnique) =
      plotRef(df, igDset, cdlFile.extractFilename, yr2018, mt, outpath)
      plotScatter(df, igDset, mt, outpath)
      plotTile(df, igDset, mt, outpath)
      if mt == mtNone:
        plotInterpolatedRaster(dfRef, igDset, mt,outpath)

    allPlots(dfRef, igDset, mtNone)

    let dfLinear = getMorphedDf(dfRef, igDset)
    let dfSpline = getMorphedDfSpline(dfRef, igDset)
    ## NOTE: KDE does not make sense, see docstring above
    #let dfKde    = getMorphedDfKde(dfRef, igDset)

    # for systematics we only care about what we actually use:
    diffTab[$dkKind] = computeMorphingSystematics(dfLinear)

    allPlots(dfLinear, igDset, mtLinear)
    allPlots(dfSpline, igDset, mtSpline)
    ## NOTE: KDE does not make sense, see docstring above

  # print info about morphing systematic differences
  let xrayRef = getXrayRefTable()
  var meanDiff = 0.0
  for idx in 0 ..< xrayRef.len:
    let dsetName = xrayRef[idx]
    if dsetName in ["Cu-Ni-15kV", "C-EPIC-0.6kV"]: continue
    var diff = 0.0
    for dkKind in dkKinds:
      diff += diffTab[$dkKind][dsetName]
    diff /= 3.0
    echo "Target/Filter: ", dsetName, " = ", diff
    meanDiff += diff
  echo "Mean difference ", meanDiff / 6.0

when isMainModule:
  dispatch main
