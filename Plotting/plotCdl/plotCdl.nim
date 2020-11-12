import nimhdf5, ggplotnim, os, strformat, strutils, sequtils, algorithm, seqmath
import ingrid / tos_helpers, times

import cligen

type
  HistTuple = tuple[bins: seq[float64], hist: seq[float64]]

  SupportedRead = SomeFloat | SomeInteger | string | bool | Value

  YearKind = enum
    yr2014 = "2014"
    yr2018 = "2018"

const CommonNames = ["totalCharge", "eventNumber", "hits", "energyFromCharge"]

proc splitSeq[T, U](s: seq[seq[T]], dtype: typedesc[U]): (seq[U], seq[U]) =
  ## splits a (N, 2) nested seq into two seqs
  result[0] = newSeq[dtype](s.len)
  result[1] = newSeq[dtype](s.len)
  for i in 0..s.high:
    result[0][i] = s[i][0].U
    result[1][i] = s[i][1].U

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
    dfJoined["runNumber"] = constantColumn(run.parseInt, dfJoined.len)
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
    let h5f = H5File(file, "r")
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
  var h5ref = H5file(refFile, "r")
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

proc normalize(df: DataFrame): DataFrame =
  #let sum = result["Hist"].toTensor(float).sum
  result = df
  result = result.mutate(f{float -> float: "Hist" ~ `Hist` / max(df["Hist"])})

proc main(files: seq[string],
          refFile: string = "/mnt/1TB/CAST/CDL_2019/XrayReferenceFile2018.h5",
          fkKind: FrameworkKind = fkTpa) =
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



  #echo df.pretty(-1)


when isMainModule:
  dispatch main
