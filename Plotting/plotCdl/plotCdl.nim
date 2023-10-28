import nimhdf5, ggplotnim, os, strformat, strutils, sequtils, algorithm, seqmath, times
import ingrid / [ingrid_types, tos_helpers]
import helpers / utils

import cligen

import ginger

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
    dfJoined["runNumber"] = run
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
  result["runType"] = runType
  result.classifyClusters()

proc readRefDf(ctx: LikelihoodContext): DataFrame =
  ## Use existing logic in `likelihood_utils.nim`
  result = readRefDsetsDf(ctx)
  result["runType"] = "CDL"

proc binIt(dfLoc: DataFrame, dset: string): (seq[int], seq[float]) =
  let (numBins, minVal, maxVal) = cdlToXrayBinning2018(dset)
  result = histogram(dfLoc[dset].toTensor(float).clone.toRawSeq,
                     numBins,
                     range = (minVal, maxVal),
                     upperRangeBinRight = false)

proc binDf(df: DataFrame, dset: string): DataFrame =
  for tup, subDf in groups(df.group_by("Dset")):
    let (hist, bins) = binIt(subDf, dset)
    var dfToAdd = toDf({ "Bins" : bins[0 ..< ^1],
                         "Hist" : hist })
    dfToAdd["runType"] = df["runType", 0].toStr
    dfToAdd["Dset"] = tup[0][1].toStr
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
    var dfDset = toDf({ "Bins" : tab[dset_name].bins, "Hist" : tab[dset_name].hist })
    dfDset["Dset"] = dset_name
    result.add dfDset
  result["runType"] = "CDL"

proc plotRef(df: DataFrame,
             dset: string,
             refFile: string, yearKind: YearKind,
             outpath: string) =
  #block RefPlots:
  #  ggplot(df, aes("Eccentricity", "Ecc #", fill = "Dset")) +
  #    geom_histogram(stat = "identity", position = "identity", alpha = some(0.5)) +
  #    ggtitle(&"Eccentricity of reference file, year: {yearKind}") +
  #    ggsave(&"{outpath}/eccentricity_{refFile.extractFilename}_{yearKind}.pdf",
  #            width = 800, height = 480)
  #  ggplot(df, aes("L / RMS_trans", "L / RMS_trans #", fill = "Dset")) +
  #    geom_histogram(stat = "identity", position = "identity", alpha = some(0.5)) +
  #    ggtitle(&"L / RMS_trans of reference file, year: {yearKind}") +
  #    ggsave(&"{outpath}/lengthDivRmsTrans_{refFile.extractFilename}_{yearKind}.pdf",
  #            width = 800, height = 480)
  #  ggplot(data = df.filter(f{Value: isNull(df["fracRmsTrans"][idx]) == (%~ false)}),
  #         aes("fracRmsTrans", "fracRmsTrans #", fill = "Dset")) +
  #    geom_histogram(stat = "identity", position = "identity", alpha = some(0.5)) +
  #    ggtitle(&"fracRmsTrans of reference file, year: {yearKind}") +
  #    ggsave(&"{outpath}/fracRmsTrans_{refFile.extractFilename}_{yearKind}.pdf",
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
    ggsave(&"{outpath}/{dset}_ridgeline_{refFile.extractFilename}_{yearKind}.pdf",
            width = 800, height = 480)

proc normalize(df: DataFrame, grp = "Dset", values = "Hist"): DataFrame =
  result = df
  result = result.group_by("Dset")
    .mutate(f{float -> float: "Hist" ~ `Hist` / sum(df["Hist"])})

proc plotCdlFile(ctx: LikelihoodContext, outpath: string) =
  var df = newDataFrame()
  let xrayTab = getXrayRefTable()

  var refDf = newDataFrame()
  for idx, key in xrayTab:
    let (logL, energies) = buildLogLHist(key, ctx)
    var dfLoc = toDf({"logL" : logL, "Energy" : energies})
    dfLoc["Dset"] = key
    df.add dfLoc

    ## read the data from the CDL file and generate the reference data using cuts
    let (eccs, ldiv, frac) = key.genRefData(ctx)

    proc dsetToDf(dset: HistTuple, property: string): DataFrame =
      result = toDf({ "Bins" : dset[0],
                      "Counts" : dset[1],
                      "Dset" : property })
    var propDf = newDataFrame()
    propDf.add dsetToDf(eccs, "eccentricity")
    propDf.add dsetToDf(ldiv, "lengthDivRmsTrans")
    propDf.add dsetToDf(frac, "fractionInTransverseRms")

    propDf["tfKind"] = key
    refDf.add propDf

  let xrayRef = getXrayRefTable()
  var labelOrder = initTable[Value, int]()
  for idx, el in xrayRef:
    labelOrder[%~ el] = idx

  var img = initViewport(wImg = 1200, hImg = 600, backend = bkCairo)
  img.layout(3, rows = 1)
  var idx = 0
  let title = "Overview of all reference distributions for each of the three cluster properties"
  for (tup, subDf) in refDf.group_by("Dset").groups:
    let dset = tup[0][1].toStr
    proc filterDset(df: DataFrame, dset: string): DataFrame =
      case dset
      of "eccentricity":            result = df.filter(f{`Bins` <= 4.0 and `Bins` >= 0.0})
      of "lengthDivRmsTrans":       result = df.filter(f{`Bins` <= 15.0 and `Bins` >= 0.0})
      of "fractionInTransverseRms": result = df#.filter(f{`Bins` <= 15.0 and `Bins` >= 0.0})
      else: doAssert false, "Dset: " & $dset
    let dfF = filterDset(subDf, dset)
    var plt = ggplot(dfF, aes("Bins", "Counts", fill = "tfKind"), backend = bkCairo) +
      ggridges("tfKind", overlap = 1.5, labelOrder = labelOrder) +
      geom_histogram(position = "identity", stat = "identity", hdKind = hdOutline, color = "black", lineWidth = 1.0) +
      xlab(dset, tickMargin = 1.5) +
      margin(left = 4.0, right = 1, top = 1.5) +
      hideLegend()
    img.embedAsRelative(idx, ggcreate(plt).view)
    inc idx

    ggplot(dfF, aes("Bins", "Counts", fill = "tfKind")) +
      facet_wrap("tfKind", scales = "free") +
      geom_histogram(stat = "identity", position = "identity") +
      ggtitle(&"{dset} of reference file") +
      xlab(dset) +
      ggsave(&"{outpath}/{dset}_facet_{ctx.cdlFile.extractFilename}.pdf",
              width = 1200, height = 1000)

  var area = img.addViewport()
  let text = area.initText(c(0.5, 0.05, ukRelative), title, goText, taCenter, font = some(font(16.0)))
  area.addObj text
  img.children.add area
  img.draw(&"{outpath}/ridgeline_all_properties_side_by_side.pdf")

  ggplot(refDf.filter(f{`Bins` <= 10.0 and `Bins` >= 0.0}), aes("Bins", "Counts", fill = "Dset")) +
    ggridges("tfKind", overlap = 1.5, labelOrder = labelOrder) +
    geom_histogram(position = "identity", stat = "identity", hdKind = hdOutline, color = "black", lineWidth = 1.0) +
    xlim(0, 10) +
    ggsave(&"{outpath}/ridgeline_all_properties_same_ridge.pdf", width = 1000, height = 600)


  let breaks = linspace(0, 30.0, 201)
  ggplot(df, aes("logL", fill = "Dset")) +
    ggridges("Dset", overlap = 1.5, labelOrder = labelOrder) +
    geom_histogram(breaks = breaks, position = "identity",
                   hdKind = hdOutline,
                   color = "black", lineWidth = 1.0,
                   density = true) +
    xlab("-ln L") +
    ggtitle("-ln L distributions for each target/filter combination") +
    ggsave(&"{outpath}/logL_ridgeline.pdf")

  ggplot(df, aes("logL", color = "Dset")) +
    geom_histogram(binWidth = 0.15, position = "identity", alpha = some(0.0),
                   lineWidth = some(1.5),
                   hdKind = hdOutline) +
    ggtitle("LogL distributions from CDL data") +
    ggsave(&"{outpath}/logL_outline.pdf")

  ## XXX: replace this by the version that we generate in `cdl_spectrum_creation`? We'd need to
  ## write the energies we calculate there to the output file though! Anyhow this plot is also
  ## useful to have.
  ggplot(df.filter(f{`Energy` < 20.0}), aes("Energy", fill = "Dset")) +
    geom_histogram(position = "identity", alpha = 0.5, bins = 300, hdKind = hdOutline) +
    xlab("Energy [keV]") + ylab("#") +
    ggtitle("Energy spectra of all GridPix Feb 2019 X-ray tube data") +
    ggsave(&"{outpath}/cdl_energies.pdf")

  # XXX: Well, it doesn't work in the way I thought, because obviously we cannot just compute
  # the likelihood distributions directly! We can compute *a* likelihood value for a cluster
  # at an arbitrary energy, but to get a likelihood distribution we'd need actual X-ray data
  # at all energies! The closest we have to that is the general CDL data (not just the main
  # peak & using correct energies for each cluster)
  # XXX 2: Not quite true. We *DO* compute the interpolated likelihood *DISTRIBUTIONS* already,
  # namely to define the cut values!
  # ## And now a heatmap of the interpolated *likelihood* distributions including the cut values
  # let cutVs = ctx.calcCutValueTab()

proc main(files: seq[string] = @[],
          # refFile: string = "/home/basti/CastData/data/CDL_2019/XrayReferenceFile2018.h5", ## <-- deprecated
          cdlFile: string = "/home/basti/CastData/data/CDL_2019/calibration-cdl-2018.h5",
          fkKind: FrameworkKind = fkTpa,
          outpath = "out") =
  let ctx = initLikelihoodContext(cdlFile,
                                  year = yr2018,
                                  energyDset = igEnergyFromCharge,
                                  region = crGold,
                                  timepix = Timepix1,
                                  morphKind = mkLinear) # morphing to plot interpolation
  if files.len == 0:
    ctx.plotCdlFile(outpath)
  else:
    # If other input files given we produce a comparison ridge line plot of each property
    # between the input and the reference data.
    let dfBack = readFiles(files,
                           "background",
                           concat(@CommonNames, @["eccentricity",
                                                  "fractionInTransverseRms",
                                                  "lengthDivRmsTrans"]))
    const dkKinds = [igEccentricity, igLengthDivRmsTrans, igFractionInTransverseRms]
    let dfRefAll = ctx.readRefDf()
    for dkKind in dkKinds:
      let dset = dkKind.toDset(fkKind)
      let dfRef = dfRefAll.filter(f{`Variable` == dset})
      let dfB = dfBack.binDf(dset)
      let df = assignStack(@[dfRef.normalize, dfB.normalize])
      plotRef(df, dset, ctx.cdlFile, yr2018, outpath)

when isMainModule:
  dispatch main
