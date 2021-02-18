import nimhdf5, ggplotnim, os, strformat, strutils, sequtils, algorithm, seqmath, times
import ingrid / [ingrid_types, tos_helpers]
import helpers / utils

import cligen

type
  HistTuple = tuple[bins: seq[float64], hist: seq[float64]]

  SupportedRead = SomeFloat | SomeInteger | string | bool | Value

const CommonNames = ["totalCharge", "eventNumber", "hits", "energyFromCharge"]

proc readData(h5f: H5FileObj,
              evNames: seq[string]): DataFrame =
  const allNames = ["timestamp", "eventNumber"]

  template readIt(names, baseName, df: untyped): untyped =
    for name in names:
      let dsetName = baseName / name
      let dsetH5 = h5f[dsetName.dset_str]
      withDset(dsetH5):
        when type(dset) is seq[SupportedRead]:
          df[name] = dset

  for run, grp in runs(h5f, recoBase()):
    var dfEv = newDataFrame()
    var dfAll = newDataFrame()
    let group = h5f[grp.grp_str]
    readIt(evNames, grp / "chip_3", dfEv)
    readIt(allNames, grp, dfAll)
    var dfJoined = inner_join(dfEv, dfAll, "eventNumber")
    dfJoined["runNumber"] = constantColumn(run, dfJoined.len)
    result.add dfJoined

proc classifyClusters(df: var DataFrame) =
  let energies = df["energyFromCharge"].toTensor(float).toRawSeq
  var dsets = newSeq[string](energies.len)
  for i, e in energies:
    dsets[i] = energies[i].toRefDset()
  df["Dset"] = dsets

proc readFiles(files: seq[string], runType: string,
               dsets: seq[string]): DataFrame =
  for file in files:
    let h5f = H5open(file, "r")
    result.add readData(h5f, dsets)
    discard h5f.close()
  result = result.arrange("timestamp")
  result["runType"] = constantColumn(runType, result.len)
  result.classifyClusters()

proc binIt(dfLoc: DataFrame, dset: string): (seq[int], seq[float]) =
  let (numBins, minVal, maxVal) = cdlToXrayBinning2018(dset)
  result = histogram(dfLoc[dset].toTensor(float).clone.toRawSeq,
                     numBins,
                     range = (minVal, maxVal),
                     upperRangeBinRight = false)

proc binDf(df: DataFrame, dset: string): DataFrame =
  for tup, subDf in groups(df.group_by("Dset")):
    let (hist, bins) = binIt(subDf, dset)
    var dfToAdd = seqsToDf({ "Bins" : bins[0 ..< ^1],
                             "Hist" : hist })
    dfToAdd["runType"] = constantColumn(df["runType", 0].toStr, dfToAdd.len)
    dfToAdd["Dset"] = constantColumn(tup[0][1].toStr, dfToAdd.len)
    result.add dfToAdd

proc readRefDsets(refFile: string,
                  dkKind: InGridDsetKind,
                  yearKind: YearKind): DataFrame =
  ## reads the reference datasets from the `refFile` and returns them.
  var h5ref = H5open(refFile, "r")
  # create a table, which stores the reference datasets from the ref file
  let dset = dkKind.toDset(fkTpa)
  var tab = initTable[string, HistTuple]()
  let xrayRef = getXrayRefTable()
  for dset_name in values(xray_ref):
    var data = h5ref[(dset_name / dset).dset_str]
    tab[dset_name] = data.readAs(float64).reshape2D(data.shape).splitSeq(float64)
    var dfDset = seqsToDf({ "Bins" : tab[dset_name].bins, "Hist" : tab[dset_name].hist })
    dfDset["Dset"] = constantColumn(dset_name, dfDset.len)
    result.add dfDset
  result["runType"] = constantColumn("CDL", result.len)

proc plotRef(df: DataFrame,
             dset: string,
             refFile: string, yearKind: YearKind) =
  #block RefPlots:
  #  ggplot(df, aes("Eccentricity", "Ecc #", fill = "Dset")) +
  #    geom_histogram(stat = "identity", position = "identity", alpha = some(0.5)) +
  #    ggtitle(&"Eccentricity of reference file, year: {yearKind}") +
  #    ggsave(&"out/eccentricity_{refFile.extractFilename}_{yearKind}.pdf",
  #            width = 800, height = 480)
  #  ggplot(df, aes("L / RMS_trans", "L / RMS_trans #", fill = "Dset")) +
  #    geom_histogram(stat = "identity", position = "identity", alpha = some(0.5)) +
  #    ggtitle(&"L / RMS_trans of reference file, year: {yearKind}") +
  #    ggsave(&"out/lengthDivRmsTrans_{refFile.extractFilename}_{yearKind}.pdf",
  #            width = 800, height = 480)
  #  ggplot(data = df.filter(f{Value: isNull(df["fracRmsTrans"][idx]) == (%~ false)}),
  #         aes("fracRmsTrans", "fracRmsTrans #", fill = "Dset")) +
  #    geom_histogram(stat = "identity", position = "identity", alpha = some(0.5)) +
  #    ggtitle(&"fracRmsTrans of reference file, year: {yearKind}") +
  #    ggsave(&"out/fracRmsTrans_{refFile.extractFilename}_{yearKind}.pdf",
  #            width = 800, height = 480)

  let xrayRef = getXrayRefTable()
  var labelOrder = initTable[Value, int]()
  for idx, el in xrayRef:
    labelOrder[%~ el] = idx

  ggplot(df, aes("Bins", "Hist", fill = "runType")) +
    ggridges("Dset", overlap = 1.0, labelOrder = labelOrder) +
    #facet_wrap("Dset", scales = "free") +
    geom_histogram(stat = "identity", position = "identity",
                   alpha = some(0.5)) +
    ggtitle(&"{dset} of reference file, year: {yearKind}") +
    ggsave(&"out/{dset}_ridgeline_{refFile.extractFilename}_{yearKind}.pdf",
            width = 800, height = 480)

proc normalize(df: DataFrame, grp = "Dset", values = "Hist"): DataFrame =
  result = df
  result = result.group_by("Dset")
    .mutate(f{float -> float: "Hist" ~ `Hist` / sum(df["Hist"])})

proc plotCdlFile(cdlFile, refFile: string) =
  let xrayTab = getXrayRefTable()
  var df = newDataFrame()
  for idx, key in xrayTab:
    let (logL, energies) = buildLogLHist(cdlFile, refFile, key,
                                         year = yr2018,
                                         region = crGold)
    var dfLoc = seqsToDf({"logL" : logL, "Energy" : energies})
    dfLoc["Dset"] = constantColumn(key, dfLoc.len)
    df.add dfLoc

  let xrayRef = getXrayRefTable()
  var labelOrder = initTable[Value, int]()
  for idx, el in xrayRef:
    labelOrder[%~ el] = idx

  let breaks = linspace(0, 30.0, 201)
  ggplot(df, aes("logL", fill = "Dset")) +
    ggridges("Dset", overlap = 1.5, labelOrder = labelOrder) +
    geom_histogram(breaks = breaks, position = "identity",
                   hdKind = hdOutline,
                   density = true) +
    ggsave(&"out/logL_ridgeline.pdf")

  ggplot(df, aes("logL", color = "Dset")) +
    geom_histogram(binWidth = 0.15, position = "identity", alpha = some(0.0),
                   lineWidth = some(1.5),
                   hdKind = hdOutline) +
    ggtitle("LogL distributions from CDL data") +
    ggsave(&"out/logL_outline.pdf")

proc main(files: seq[string] = @[],
          refFile: string = "/home/basti/CastData/data/CDL_2019/XrayReferenceFile2018.h5",
          cdlFile: string = "/home/basti/CastData/data/CDL_2019/calibration-cdl-2018.h5",
          fkKind: FrameworkKind = fkTpa) =

  if files.len == 0:
    plotCdlFile(cdlFile, refFile)
  else:
    let dfBack = readFiles(files,
                           "background",
                           concat(@CommonNames, @["eccentricity",
                                                  "fractionInTransverseRms",
                                                  "lengthDivRmsTrans"]))
    const dkKinds = [igEccentricity, igLengthDivRmsTrans, igFractionInTransverseRms]
    for dkKind in dkKinds:
      let dfRef = readRefDsets(refFile, dkKind, yr2018)
      #echo dfRef.pretty(-1)

      let dset = dkKind.toDset(fkKind)
      let dfB = dfBack.binDf(dset)
      var df: DataFrame
      df.add dfRef.normalize
      df.add dfB.normalize
      plotRef(df,
              dset,
              refFile, yr2018)

when isMainModule:
  dispatch main
