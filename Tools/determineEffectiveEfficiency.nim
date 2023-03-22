import std / [os, strutils, random, sequtils, stats, strformat]
import nimhdf5, unchained
import numericalnim except linspace


import ingrid / private / [likelihood_utils, hdf5_utils, ggplot_utils, geometry, cdl_cuts]
import ingrid / calibration
import ingrid / calibration / [fit_functions]
import ingrid / [ingrid_types, tos_helpers]
import ingridDatabase / [databaseRead, databaseDefinitions, databaseUtils]

import arraymancer / stats / kde

# cut performed regardless of logL value on the data, since transverse
# rms > 1.5 cannot be a physical photon, due to diffusion in 3cm drift
# distance
const RmsCleaningCut = 1.5

#let CdlFile = "/mnt/1TB/CAST/CDL_2019/calibration-cdl-2018.h5" #
let CdlFile = "/home/basti/CastData/data/CDL_2019/calibration-cdl-2018.h5"

proc drawNewEvent(rms, energy: seq[float]): int =
  let num = rms.len - 1
  var idx = rand(num)
  while rms[idx] >= RmsCleaningCut or
        (energy[idx] <= 4.5 or energy[idx] >= 7.5):
    idx = rand(num)
  result = idx

proc computeEnergy(h5f: H5File, pix: seq[Pix], group: string, centerChip: int,
                   a, b, c, t, bL, mL: float,
                   capacitance: FemtoFarad): float =
  let totalCharge = pix.mapIt(calibrateCharge(it.ch.float, capacitance, a, b, c, t)).sum
  # compute mean of all gas gain slices in this run (most sensible)
  let gain = h5f[group / &"chip_{centerChip}/gasGainSlices", GasGainIntervalResult].mapIt(it.G).mean
  let calibFactor = linearFunc(@[bL, mL], gain) * 1e-6
  # now calculate energy for all hits
  result = totalCharge * calibFactor

proc generateFakeData(ctx: LikelihoodContext, h5f: H5File, nFake: int, energy = 3.0): DataFrame =
  ## For each run generate `nFake` fake events
  let refSetTuple = ctx.readRefDsets()
  result = newDataFrame()
  let fileInfo = h5f.getFileInfo()
  let tpx = fileInfo.timepix
  let centerChip = fileInfo.centerChip
  let capacitance = tpx.getCapacitance()
  for (num, group) in runs(h5f):
    # first read all x / y / tot data
    let xs = h5f[group / &"chip_{centerChip}/x", special_type(uint8), uint8]
    let ys = h5f[group / &"chip_{centerChip}/y", special_type(uint8), uint8]
    let ts = h5f[group / &"chip_{centerChip}/ToT", special_type(uint16), uint16]
    let rms = h5f[group / &"chip_{centerChip}/rmsTransverse", float]
    let cX = h5f[group / &"chip_{centerChip}/centerX", float]
    let cY = h5f[group / &"chip_{centerChip}/centerY", float]
    let energyInput = h5f[group / &"chip_{centerChip}/energyFromCharge", float]
    let chipGrp = h5f[(group / &"chip_{centerChip}").grp_str]
    let chipName = chipGrp.attrs["chipName", string]
    # get factors for charge calibration
    let (a, b, c, t) = getTotCalibParameters(chipName, num)
    # get factors for charge / gas gain fit
    let (bL, mL) = getCalibVsGasGainFactors(chipName, num, suffix = $gcIndividualFits)
    var count = 0
    var evIdx = 0

    when false:
      for i in 0 ..< xs.len:
        if xs[i].len < 150 and energyInput[i] > 5.5:
          # recompute from data
          let pp = toSeq(0 ..< xs[i].len).mapIt((x: xs[i][it], y: ys[i][it], ch: ts[i][it]))
          let newEnergy = h5f.computeEnergy(pp, group, centerChip, a, b, c, t, bL, mL, capacitance)
          echo "Length ", xs[i].len , " w/ energy ", energyInput[i], " recomp ", newEnergy
          let df = toDf({"x" : pp.mapIt(it.x.int), "y" : pp.mapIt(it.y.int), "ch" : pp.mapIt(it.ch.int)})
          ggplot(df, aes("x", "y", color = "ch")) +
            geom_point() +
            ggtitle("funny its real") +
            ggsave("/tmp/fake_event_" & $i & ".pdf")
          sleep(200)
      if true: quit()

    # to store fake data
    var energies = newSeqOfCap[float](nFake)
    var logLs = newSeqOfCap[float](nFake)
    var rmss = newSeqOfCap[float](nFake)
    var eccs = newSeqOfCap[float](nFake)
    var ldivs = newSeqOfCap[float](nFake)
    var frins = newSeqOfCap[float](nFake)
    var cxxs = newSeqOfCap[float](nFake)
    var cyys = newSeqOfCap[float](nFake)
    var lengths = newSeqOfCap[float](nFake)
    while count < nFake:
      # draw index from to generate a fake event
      evIdx = drawNewEvent(rms, energyInput)
      # draw number of fake pixels
      # compute ref # pixels for this event taking into account possible double counting etc.
      let basePixels = (energy / energyInput[evIdx] * xs[evIdx].len.float)
      let nPix = round(basePixels + gauss(sigma = 10.0)).int  # ~115 pix as reference in 3 keV (26 eV), draw normal w/10 around
      if nPix < 4:
        echo "Less than 4 pixels: ", nPix, " skipping"
        continue
      var pix = newSeq[Pix](nPix)
      var seenPix: set[uint16] = {}
      let evNumPix = xs[evIdx].len

      if nPix >= evNumPix:
        echo "More pixels to draw than available! ", nPix, " vs ", evNumPix, ", skipping!"
        continue
      if not inRegion(cX[evIdx], cY[evIdx], crSilver):
        echo "Not in silver region. Not a good basis"
        continue

      var pIdx = rand(evNumPix - 1)
      for j in 0 ..< nPix:
        # draw pix index
        while pIdx.uint16 in seenPix:
          pIdx = rand(evNumPix - 1)
        seenPix.incl pIdx.uint16
        pix[j] = (x: xs[evIdx][pIdx], y: ys[evIdx][pIdx], ch: ts[evIdx][pIdx])
      # now draw
      when false:
        let df = toDf({"x" : pix.mapIt(it.x.int), "y" : pix.mapIt(it.y.int), "ch" : pix.mapIt(it.ch.int)})
        ggplot(df, aes("x", "y", color = "ch")) +
          geom_point() +
          ggsave("/tmp/fake_event.pdf")
        sleep(200)

      # reconstruct event
      let inp = (pixels: pix, eventNumber: 0, toa: newSeq[uint16](), toaCombined: newSeq[uint64]())
      let recoEv = recoEvent(inp, -1,
                             num, searchRadius = 50,
                             dbscanEpsilon = 65,
                             clusterAlgo = caDefault)
      if recoEv.cluster.len > 1 or recoEv.cluster.len == 0:
        echo "Found more than 1 or 0 cluster! Skipping"
        continue
      # compute charge
      let energy = h5f.computeEnergy(pix, group, centerChip, a, b, c, t, bL, mL, capacitance)

      # puhhh, now the likelihood...
      let ecc = recoEv.cluster[0].geometry.eccentricity
      let ldiv = recoEv.cluster[0].geometry.lengthDivRmsTrans
      let frin = recoEv.cluster[0].geometry.fractionInTransverseRms
      let logL = calcLikelihoodForEvent(energy,
                                        ecc,
                                        ldiv,
                                        frin,
                                        refSetTuple)
      # finally done
      energies.add energy
      logLs.add logL
      rmss.add recoEv.cluster[0].geometry.rmsTransverse
      eccs.add ecc
      ldivs.add ldiv
      frins.add frin
      cxxs.add recoEv.cluster[0].centerX
      cyys.add recoEv.cluster[0].centerY
      lengths.add recoEv.cluster[0].geometry.length
      inc count
    let df = toDf({ "energyFromCharge" : energies,
                    "likelihood" : logLs,
                    "runNumber" : num,
                    "rmsTransverse" : rmss,
                    "eccentricity" : eccs,
                    "lengthDivRmsTrans" : ldivs,
                    "centerX" : cxxs,
                    "centerY" : cyys,
                    "length" : lengths,
                    "fractionInTransverseRms" : frins })
    result.add df

proc applyLogLCut(df: DataFrame, cutTab: CutValueInterpolator): DataFrame =
  result = df.mutate(f{float -> bool: "passLogL?" ~ (
    block:
      #echo "Cut value: ", cutTab[idx(igEnergyFromCharge.toDset())], " at dset ", toRefDset(idx(igEnergyFromCharge.toDset())), " at energy ", idx(igEnergyFromCharge.toDset())
      idx(igLikelihood.toDset()) < cutTab[idx(igEnergyFromCharge.toDset())])})

proc readRunData(h5f: H5File, chip: int): DataFrame =
  result = h5f.readDsets(chipDsets =
    some((chip: chip,
          dsets: @[igEnergyFromCharge.toDset(),
                   igRmsTransverse.toDset(),
                   igLengthDivRmsTrans.toDset(),
                   igFractionInTransverseRms.toDset(),
                   igEccentricity.toDset(),
                   igCenterX.toDset(),
                   igCenterY.toDset(),
                   igLength.toDset(),
                   igLikelihood.toDset()]))
  )

proc filterEvents(df: DataFrame, energy: float = Inf,
                  eccFilters = none[tuple[e, p: float]]()): DataFrame =
  let xrayCutsTab {.global.} = getXrayCleaningCuts()
  template applyFilters(dfI: typed, eccCutOverride = Inf): untyped {.dirty.} =
    let minRms = xrayCuts.minRms
    let maxRms = xrayCuts.maxRms
    let maxLen = xrayCuts.maxLength
    let maxEcc = if classify(eccCutOverride) != fcInf: eccCutOverride else: xrayCuts.maxEccentricity
    dfI.filter(f{float -> bool: idx(igRmsTransverse.toDset()) < RmsCleaningCut and
      inRegion(idx("centerX"), idx("centerY"), crSilver) and
      idx("rmsTransverse") >= minRms and
      idx("rmsTransverse") <= maxRms and
      idx("length") <= maxLen and
      idx("eccentricity") <= maxEcc
    })
  if "Peak" in df:
    doAssert classify(energy) == fcInf
    result = newDataFrame()
    for (tup, subDf) in groups(df.group_by("Peak")):
      case tup[0][1].toStr
      of "Escapepeak":
        let dset = 2.9.toRefDset()
        let xrayCuts = xrayCutsTab[dset]
        let dfF = if eccFilters.isSome:
                    applyFilters(subDf, eccFilters.get.e)
                  else:
                    applyFilters(subDf)
        result.add dfF
      of "Photopeak":
        let dset = 5.9.toRefDset()
        let xrayCuts = xrayCutsTab[dset]
        let dfF = if eccFilters.isSome:
                    applyFilters(subDf, eccFilters.get.p)
                  else:
                    applyFilters(subDf)
        result.add dfF
      else: doAssert false, "Invalid name"
  else:
    doAssert classify(energy) != fcInf
    let dset = energy.toRefDset()
    let xrayCuts = xrayCutsTab[dset]
    result = applyFilters(df)

proc splitPeaks(df: DataFrame): DataFrame =
  let eD = igEnergyFromCharge.toDset()
  echo df
  result = df.mutate(f{float -> string: "Peak" ~ (
    if idx(eD) < 3.5 and idx(eD) > 2.5:
      "Escapepeak"
    elif idx(eD) > 4.5 and idx(eD) < 7.5:
      "Photopeak"
    else:
      "None")})
    .filter(f{`Peak` != "None"})

proc handleFile(fname: string, cutTab: CutValueInterpolator,
                eccFilters = none[tuple[e, p: float]]()): DataFrame =
  ## Given a single input file, performs application of the likelihood cut for all
  ## runs in it, split by photo & escape peak. Returns a DF with column indicating
  ## the peak, energy of each event & a column whether it passed the likelihood cut.
  ## Only events that are pass the input cuts are stored.
  let h5f = H5open(fname, "r")
  let fileInfo = h5f.getFileInfo()
  randomize(423)
  result = newDataFrame()
  let data = h5f.readRunData(fileInfo.centerChip)
    .splitPeaks()
    .filterEvents(eccFilters = eccFilters)
    .applyLogLCut(cutTab)
  result.add data
  when false:
    ggplot(result, aes("energyFromCharge")) +
      geom_histogram(bins = 200) +
      ggsave("/tmp/ugl.pdf")
  discard h5f.close()

proc handleFakeData(ctx: LikelihoodContext, fname: string, energy: float, cutTab: CutValueInterpolator): DataFrame =
  let h5f = H5open(fname, "r")
  var data = generateFakeData(ctx, h5f, 5000, energy = energy)
    .filterEvents(energy)
    .applyLogLCut(cutTab)
  result = data
  discard h5f.close()

proc getIndices(dset: string): seq[int] =
  ## XXX: This is not safe and does not work if `fitByRun` is used!!
  result = newSeq[int]()
  withLogLFilterCuts(CdlFile, dset, yr2018, igEnergyFromCharge, LogLCutDsets):
  #withXrayRefCuts(CdlFile, dset, yr2018, igEnergyFromCharge):
    result.add i

proc readCdlData(df: DataFrame, energy: float, cutTab: CutValueInterpolator): DataFrame =
  # map input fake energy to reference dataset
  let grp = energy.toRefDset()

  let passedInds = getIndices(grp)
  let h5f = H5open(CdlFile, "r")
  const xray_ref = getXrayRefTable()

  var dfR = newDataFrame()
  for dset in IngridDsetKind:
    try:
      let d = dset.toDset()
      if d notin df: continue # skip things not in input
      ## first read data from CDL file (exists for sure)
      ## extract all CDL data that passes the cuts used to generate the logL histograms
      var cdlFiltered = newSeq[float](passedInds.len)
      let cdlRaw = h5f[cdlGroupName(grp, "2019", d), float]
      for i, idx in passedInds:
        cdlFiltered[i] = cdlRaw[idx]
      echo "Total number of elements ", cdlRaw.len, " filtered to ", passedInds.len
      dfR[d] = cdlFiltered
      when false:
        ## XXX: if we want to keep the following, we need to use the `withXrayRefCuts`
        ## template to generate each dataset on the fly. 'Problem' is only that the
        ## template needs to open the file itself, but the file is already open
        ## in this scope. So need to restructure the code for that.

        ## now read histograms from RefFile, if they exist (not all datasets do)
        if grp / d in h5f:
          let dsetH5 = h5f[(grp / d).dset_str]
          let (bins, data) = dsetH5[float].reshape2D(dsetH5.shape).split(Seq2Col)
          let fname = &"/tmp/{grp}_{d}_energy_{energy:.1f}.pdf"
          echo "Storing histogram in : ", fname
          # now add fake data
          let dataSum = simpson(data, bins)
          let refDf = toDf({"bins" : bins, "data" : data})
            .mutate(f{"data" ~ `data` / dataSum})
          let df = df.filter(f{float: idx(d) <= bins[^1]})
          ggplot(refDf, aes("bins", "data")) +
            geom_histogram(stat = "identity", hdKind = hdOutline, alpha = 0.5) +
            geom_histogram(data = df, aes = aes(d), bins = 200, alpha = 0.5,
                           fillColor = "orange", density = true, hdKind = hdOutline) +
            ggtitle(&"{d}. Orange: fake data from 'reducing' 5.9 keV data @ {energy:.1f}. Black: CDL ref {grp}") +
            ggsave(fname, width = 1000, height = 600)
    except AssertionError:
      continue

  discard h5f.close()

  result = dfR.applyLogLCut(cutTab)

proc plotRefHistos(df: DataFrame, energy: float, cutTab: CutValueInterpolator,
                   dfAdditions: seq[tuple[name: string, df: DataFrame]] = @[]) =
  let grp = energy.toRefDset()
  # get effect of logL cut on CDL data
  let dfR = df.readCdlData(energy, cutTab)
  var dfs = @[("Fake", df), ("Real", dfR)]

  if dfAdditions.len > 0:
    dfs = concat(dfs, dfAdditions)
  var dfPlot = bind_rows(dfs, "Type")
  echo "Rough filter removes: ", dfPlot.len
  dfPlot = dfPlot.filter(f{`lengthDivRmsTrans` <= 50.0 and `eccentricity` <= 5.0})
  echo "To ", dfPlot.len, " elements"
  ggplot(dfPlot, aes("lengthDivRmsTrans", "fractionInTransverseRms", color = "eccentricity")) +
    facet_wrap("Type") +
    geom_point(size = 1.0, alpha = 0.5) +
    ggtitle(&"Fake energy: {energy:.2f}, CDL dataset: {grp}") +
    ggsave(&"/tmp/scatter_colored_fake_energy_{energy:.2f}.png", width = 1200, height = 800)

  # plot likelihood histos
  ggplot(dfPlot, aes("likelihood", fill = "Type")) +
    geom_histogram(bins = 200, alpha = 0.5, hdKind = hdOutline) +
    ggtitle(&"Fake energy: {energy:.2f}, CDL dataset: {grp}") +
    ggsave(&"/tmp/histogram_fake_energy_{energy:.2f}.pdf", width = 800, height = 600)

  #echo "DATASET : ", grp, "--------------------------------------------------------------------------------"
  echo "Efficiency of logL cut on filtered CDL data (should be 80%!) = ", dfR.filter(f{idx("passLogL?") == true}).len.float / dfR.len.float
  echo "Elements passing using `passLogL?` ", dfR.filter(f{idx("passLogL?") == true}).len, " vs total ", dfR.len
  let (hist, bins) = histogram(dfR["likelihood", float].toRawSeq, 200, (0.0, 30.0))
  ggplot(toDf({"Bins" : bins, "Hist" : hist}), aes("Bins", "Hist")) +
    geom_histogram(stat = "identity") +
    ggsave("/tmp/usage_histo_" & $grp & ".pdf")
  let cutval = determineCutValue(hist, eff = 0.8)
  echo "Effficiency from `determineCutValue? ", bins[cutVal]

proc plotDset(dfIn: DataFrame, peakDf, dfR: DataFrame, dset, typ: string,
              energy: float, suffix = "") =
  ## XXX: remove dfIn
  let dM = peakDf[dset, float].mean()
  let cM = dfR[dset, float].mean()
  echo "Mean of ", dset, " = ", dM, " vs CDL = ", cM
  let yM = 1.5 #(dfIn.len div 100).float
  let dfLine = toDf({ "x" : [dM, dM, cM, cM], "y" : [0.0, yM, 0.0, yM],
                      "Type" : ["55Fe", "55Fe", "CDL", "CDL"]})
  let fname = "/t/" & dset & "_" & $typ & $suffix & "_compare.pdf"
  echo "Saving plot: ", fname
  var dfDens = newDataFrame()
  for t, s in groups(dfIn.group_by("Type")):
    let x = s[dset, float]
    let xs = linspace(x.min, x.max, 1000)
    let ys = kde(x, normalize = true, adjust = 2.5)
    dfDens.add toDf({"x" : xs, "y" : ys, "Type" : t[0][1]})
  ggplot(dfIn, aes(dset, fill = "Type")) +
    geom_histogram(bins = 40, hdKind = hdOutline, alpha = 0.5, position = "identity", density = true) +
    geom_line(data = dfDens, aes = aes(x = "x", y = "y", color = "Type"), fillColor = color(0,0,0,0)) +
    geom_line(data = dfLine, aes = aes(x = "x", y = "y")) +
    ggtitle("Property: " & dset & ", Energy: " & $typ) +
    ggsave(fname)

proc mutateCdl(fnTab: var Table[string, FormulaNode], df, peakDf: DataFrame, dset, typ: string,
               dsetMin, dsetMean, dsetMax, cdlMin, cdlMean, cdlMax: float,
               energy: float, suffix: string, genPlots: bool): DataFrame =
  if true: # dset == "eccentricity":
     #let fn = f{float: dset ~ (idx(dset) - col(dset).min) / (col(dset).max - col(dset).min) * (Max - Min) + Min}
    let fn = f{float: dset ~ (idx(dset) - cdlMin) / (cdlMax - cdlMin) * (dsetMax - dsetMin) + dsetMin}
    #let fn = f{float: dset ~ idx(dset)}
    fnTab[dset] = fn
    result = df.clone().mutate(fn)   #idx(dset) + shift})
    # echo "80 PERCENTILE FOR ECC ", result[dset, float].toSeq1D.percentile(80), " and data ", peakDf[dset, float].toSeq1D.percentile(80)
  else:
    let shift = dsetMean - cdlMean
    let fn = f{float: dset ~ idx(dset) + shift}
    #let fn = f{float: dset ~ idx(dset)}
    fnTab[dset] = fn
    result = df.clone().mutate(fn)

  ## XXX: move the following out here
  if genPlots:
    let dfL = bind_rows([("55Fe", peakDf), ("CDL", result)], "Type")
    dfL.plotDset(peakDf, result, dset, typ, energy, suffix)
    #result = dfS.clone() ## XXX: without clone we run into some ref semantics issue!!!

proc filterNaN(s: seq[float]): seq[float] =
  result = s.filterIt(classify(it) notin {fcInf, fcNan, fcNegInf})

proc computeLogL(fnTab: Table[string, FormulaNode], df: DataFrame, dset: string): seq[float] =
  let E = df["energyFromCharge", float]
  let e = df["eccentricity", float]
  let l = df["lengthDivRmsTrans", float]
  let f = df["fractionInTransverseRms", float]
  var fn = newSeq[FormulaNode]()
  for val in values(fnTab):
    fn.add val

  let (eccs, ldiv, frac) = genRefData(CdlFile, dset, yr2018, igEnergyFromCharge, fn)
  let (h, b) = histogram(e.toSeq1D.filterNaN(),
                         bins = 100)
  let dfC = toDf({"b" : eccs[0], "h" : eccs[1]})
  echo dfC

  ggplot(toDf({"h" : h, "b" : b}), aes("b", "h")) +
    geom_histogram(alpha = 0.5, hdKind = hdOutline, stat = "identity", fillColor = "blue") +
    geom_histogram(data = dfC, aes = aes("b", "h"),
                   stat = "identity", alpha = 0.5, hdKind = hdOutline, fillColor = "red") +
    xlim(0.0, 2.0) +
    ggsave("/t/compare_ecc_data_ref.pdf")
  result = newSeqOfCap[float](e.size)
  var num = 0
  for i in 0 ..< e.size:
    let logL = calcLogL(e[i], l[i], f[i], eccs, ldiv, frac)
    if classify(logL) notin {fcInf, fcNan, fcNegInf}:
      result.add logL
      inc num
  echo "COMPARISON ", num, " vs ", e.size

proc applyMutation(fnTab: var Table[string, FormulaNode],
                   cdlDf, peakDf: DataFrame,
                   peak: string, energy: float,
                   genPlots = true): DataFrame =
  result = cdlDf.clone()

  for dset in keys(fnTab):
    let data = peakDf[dset, float]
    let cdl = cdlDf[dset, float]
    let
      dataMin = data.percentile(1) #.min
      dataMean = data.mean
      dataMax = data.percentile(99) # max
      cdlMin = cdl.percentile(1) #min
      cdlMean = cdl.mean
      cdlMax = cdl.percentile(99) # max
    result = fnTab.mutateCdl(
      result, peakDf,
      dset, peak,
      dataMin, dataMean, dataMax, cdlMin, cdlMean, cdlMax,
      energy, "_shifted",
      genPlots
    )

proc calcEff(dataLogL: seq[float], passed: Tensor[bool], cdlCutVal, dataCutVal: float,
             peak: string): float =
  var numLeft = 0
  var numLeftIt = 0
  var numPassed = 0
  for i, l in dataLogL:
    if l <= cdlCutVal:
      inc numLeft
    if l <= dataCutVal:
      inc numLeftIt
    if passed[i]:
      inc numPassed
  echo "==================== Effective efficiency with shifted logL ===================="
  echo "Cut value at = ", cdlCutVal
  echo "Data cut value = ", dataCutVal
  echo peak, " = ", numLeft.float / dataLogL.len.float
  echo peak, " based on data = ", numLeftIt.float / dataLogL.len.float
  echo peak, " passed = ", numPassed.float / dataLogL.len.float
  result = numLeft.float / dataLogL.len.float

proc calcCut(fnTab: Table[string, FormulaNode], df: DataFrame, energy: float): (seq[float], float) =
  let logL = fnTab.computeLogL(df, energy.toRefDset())
  let logLSorted = logL.sorted
  let cutIdx = determineCutValueData(logLSorted, 0.8)
  let cutVal = logLSorted[cutIdx]
  result = (logLSorted, cutVal)

proc computeMeans(df: DataFrame, cutTab: CutValueInterpolator) =
  #var shifts = newDataFrame()
  echo df
  for (tup, peakDf) in groups(df.group_by("Peak")):
    var dfR: DataFrame
    var energy: float
    let peak = tup[0][1].toStr
    case peak
    of "Escapepeak":
      energy = 2.9
      dfR = peakDf.readCdlData(2.9, cutTab)
        .mutate(f{"Peak" <- "Escapepeak"})
    of "Photopeak":
      energy = 5.9
      dfR = peakDf.readCdlData(5.9, cutTab)
        .mutate(f{"Peak" <- "Photopeak"})
    else: doAssert false
    echo "------------------------------ Means of ", peak, " ------------------------------"
    ggplot(peakDf, aes("eccentricity")) +
      geom_histogram(bins = 100, hdKind = hdOutline, position = "identity") + #, density = true
      ggtitle("Eccentricity of 55Fe data alone") +
      ggsave("/t/eccentricity_" & $tup[0][1] & "_55fe_alone.pdf")
    let dfL = bind_rows([("55Fe", peakDf), ("CDL", dfR)], "Type")

    dfL.plotDset(peakDf, dfR, "eccentricity", peak, energy)
    dfL.plotDset(peakDf, dfR, "fractionInTransverseRms", peak, energy)
    dfL.plotDset(peakDf, dfR, "lengthDivRmsTrans", peak, energy)

    # 1. compute means, min, max
    # 2. apply mutation of data based on those
    # 3. compute logL

    ## now plot each dset again, but shifted
    echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Shifted >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    var fnTab = { "eccentricity" : f{""}, "lengthDivRmsTrans" : f{""},
                  "fractionInTransverseRms" : f{""} }.toTable

    let dfShifted = fnTab.applyMutation(dfR, peakDf, peak, energy)
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Shifted <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

    let (cdlLogL, cdlCutVal) = fnTab.calcCut(dfShifted, energy)
    let (dataLogL, dataCutVal) = fnTab.calcCut(peakDf, energy)

    let eff = calcEff(dataLogL, peakDf["passLogL?", bool], cdlCutVal, dataCutVal, peak)

    echo "EFF ", eff, "s"
    proc toDfH(s: seq[float]): DataFrame =
      let (h, b) = histogram(s.filterNaN(), bins = 50)
      result = toDf({"bins" : b, "hist" : h})

    ggplot(cdlLogL.toDfH, aes("bins", "hist")) +
           #aes("bins", "hist")) +
      geom_histogram(stat = "identity", alpha = 0.5, hdKind = hdOutline, fillColor = "blue") +
      #geom_histogram(alpha = 0.5, hdKind = hdOutline, fillColor = "blue") +
      geom_histogram(
        data = toDfH(dataLogL), stat = "identity", alpha = 0.5, fillColor = "red",
        #data = toDf(dataLogL), aes = aes("dataLogL"), alpha = 0.5, fillColor = "red",
        hdKind = hdOutline
      ) +
      geom_linerange(aes = aes(x = cdlCutVal, yMin = 0.0, yMax = 1e4), color = "blue") +
      geom_linerange(aes = aes(x = dataCutVal, yMin = 0.0, yMax = 1e4), color = "red") +
      scale_y_continuous() +
      ggtitle("Typ " & $peak) +
      ggsave("/t/logl_" & $peak & ".pdf")
    echo "trying to write, plot done"
    dfShifted.writeCsv("/t/data_" & $peak & ".csv")
    plotRefHistos(peakDf.filter(f{`eccentricity` < 1.35}), energy,
                  cutTab, dfAdditions = @[("ShiftedCDL", dfShifted.filter(f{`eccentricity` < 1.35}))])
    # now compute logL for each input event & compute effective efficiency

  when false:
    ggplot(shifts, aes("energy", "shift", color = "Type")) +
      facet_wrap("Dset", scales = "free")  +
      facet_margin(0.5) +
      geom_point() +
      ggtitle("Difference between 55Fe and CDL for each dataset at Photo-&Escapepeak") +
      ggsave("/t/shift_cdl_55fe.pdf")

proc toEdf(x: seq[float], bins: seq[float]): seq[float] =
  ## Computes the EDF of the input data `x` given some `bins`.
  ##
  ## The bins are used as boundaries to count elements below each bin edge. The
  ## result contains `bins.len` elements, all in [0, 1]
  let xS = x.sorted
  var i = 0
  result = newSeq[float](bins.len)
  for j, b in bins:
    while i < xS.len and xS[i] <= b:
      inc i
    result[j] = i.float / x.len.float
  doAssert result.allIt(it <= 1.0)

proc kolmogorovSmirnov(x1, x2: seq[float]): float =
  ## Compute the Kolmogorov-Smirnov test for two datasets.
  ##
  ## The data is binned first to min & max of the combined data range and based on the
  ## associated EDF the KS test is performed.
  ##
  ## ``KS(x) = sup_x | EDF₁(x) - EDF₂(x) |``
  ##
  ## or in ~code
  ##
  ## ``KS(x) = min(abs(EDF₁ -. EDF₂(x)))``
  let range = (min: min(x1.min, x2.min), max: max(x1.max, x2.max))
  let bins = linspace(range[0], range[1], 200)
  let x1Edf = x1.toEdf(bins)
  let x2Edf = x2.toEdf(bins)
  result = 0.0
  for i, b in bins:
    result = max( result, abs(x1Edf[i] - x2Edf[i]) )

proc plotEccs(ecVal, ksVal: seq[float], peak: string) =
  let ecDf = toDf({"ecc" : ecVal, "b" : toSeq(0 ..< ecVal.len), "ks" : ksVal})
  ggplot(ecDf, aes("ecc", "ks", color = "b")) +
    geom_point() + geom_line() +
    ggsave("/t/eccs_tested_" & $peak & ".pdf")

proc determineEccentricityCutoff(data, dataCdl: DataFrame): tuple[ecc_esc, ecc_pho: float] =
  for (tup, peakDf) in groups(data.group_by("Peak")):
    let peak = tup[0][1].toStr
    let energy = if peak == "Escapepeak": 2.9 else: 5.9
    var ecc = 1.9 ## Turn this into a parameter
    var lastKs = Inf
    var ks = Inf
    var lastEcc = Inf
    var fnTab = { "eccentricity" : f{""}, "lengthDivRmsTrans" : f{""},
                  "fractionInTransverseRms" : f{""} }.toTable

    var ecVal = newSeq[float]()
    var ksVal = newSeq[float]()
    const Δε = 0.1
    const α = 0.1
    const absKS = 0.005
    const ΔKS = 1e-3

    var n = 0

    var bestKs = Inf
    var bestEcc = Inf
    var sign = 1.0
    while ks > absKS:
      let peakDf = peakDf.filter(f{`eccentricity` < ecc})
      echo "MAX PEAK ", peakDf["eccentricity", float].max
      var dfCdl = dataCdl.filter(f{`Peak` == peak})
      dfCdl = fnTab.applyMutation(dfCdl, peakDf, peak, energy, genPlots = false)
      ks = kolmogorovSmirnov(peakDf["eccentricity", float].toSeq1D,
                             dfCdl["eccentricity", float].toSeq1D)
      let dfL = bind_rows([("55Fe", peakDf), ("CDL", dfCdl)], "Type")
      dfL.plotDset(peakDf, dfCdl, "eccentricity", peak, energy, "_1")

      ecVal.add ecc
      ksVal.add ks
      if ks < absKS or # stop early and not adjust `ecc`
        abs(ks - lastKs) < ΔKS: break

      if ks > lastKs: #ks > bestKs and ks > lastKs: # only possibly change sign if we're worse than before!
        sign = -sign
      let adj = sign * Δε * exp(- n * α)
      lastEcc = ecc
      ecc = ecc - adj

      if ks < bestKs:
        bestEcc = ecc
        bestKs = ks
      lastKs = ks
      inc n

      plotEccs(ecVal, ksVal, peak)

    if peak == "Escapepeak":
      result.ecc_esc = ecc
    else:
      result.ecc_pho = ecc

    plotEccs(ecVal, ksVal, peak)

proc matchDistributions(files: seq[string], cutTab: CutValueInterpolator) =
  # 1. read raw data for CDL & 55Fe (up to some larger number, e.g. 4 in eccentricity)
  # 2. compute effective eff
  # 3. binary search with new eccentricity filter until close to expected
  var data = newDataFrame()
  var eccFilters = (e: 5.0, p: 5.0)
  for f in files:
    data.add handleFile(f, cutTab, eccFilters = some(eccFilters))
  var dataCdl = newDataFrame()
  dataCdl.add(data.readCdlData(2.9, cutTab).mutate(f{"Peak" <- "Escapepeak"}))
  dataCdl.add(data.readCdlData(5.9, cutTab).mutate(f{"Peak" <- "Photopeak"}))

  let (ecc_esc, ecc_pho) = determineEccentricityCutoff(data, dataCdl)

  for (tup, peakDf) in groups(data.group_by("Peak")):
    ## now use final `ecc` to compute  efficiency
    let peak = tup[0][1].toStr
    let ecc = if peak == "Escapepeak": ecc_esc else: ecc_pho
    let energy = if peak == "Escapepeak": 2.9 else: 5.9
    let peakDf = peakDf.filter(f{`eccentricity` < ecc})
    echo "MAX PEAK ", peakDf["eccentricity", float].max
    var dfCdl = dataCdl.filter(f{`Peak` == peak})
    var fnTab = { "eccentricity" : f{""}, "lengthDivRmsTrans" : f{""},
                  "fractionInTransverseRms" : f{""} }.toTable

    dfCdl = fnTab.applyMutation(dfCdl, peakDf, peak, energy)
    let (cdlLogL, cdlCutVal) = fnTab.calcCut(dfCdl, energy)
    let (dataLogL, dataCutVal) = fnTab.calcCut(peakDf, energy)
    let eff = calcEff(dataLogL, peakDf["passLogL?", bool], cdlCutVal, dataCutVal, peak)
    echo "RESULTING EFFICIENCY : ", eff

  # look at resulting efficiencies and compare with cuts
  let xrayCuts = getXrayCleaningCuts()
  let energies = getXrayFluorescenceLines()
  let eccCuts = toSeq(values(xrayCuts)).mapIt(it.maxEccentricity)
  echo eccCuts
  var df = toDf({"energies" : energies, "cuts" : eccCuts, "Type" : "CDL"})
  df.add toDf({"energies" : [2.9, 5.9], "cuts" : [ecc_esc, ecc_pho], "Type" : "55Fe"})
  let m = (ecc_pho / 1.3 - ecc_esc / 1.4) / (5.9 - 2.9)
  let b = ecc_pho / 1.3  - m * 5.9
  var extraCuts = newSeq[float]()
  for i, cut in eccCuts:
    extraCuts.add (energies[i] * m + b) * cut
  df.add toDf({"energies" : energies, "cuts" : extraCuts, "Type" : "Extra"})
  df = df.filter(f{`cuts` != Inf})
  ggplot(df, aes("energies", "cuts", color = "Type")) +
    geom_point() +
    ggsave("/t/cut_values.pdf")

proc main(files: seq[string], fake = false, real = false, refPlots = false,
          computeMeans = false,
          matchDistributions = false,
          energies: seq[float] = @[]) =
  ## given the input files of calibration runs, walks all files to determine the
  ## 'real' software efficiency for them & generates a plot
  let ctx = initLikelihoodContext(CdlFile, yr2018, crSilver, igEnergyFromCharge,
                                  Timepix1, mkLinear)
  let cutTab = ctx.calcCutValueTab()
  var df = newDataFrame()
  if real and not fake:
    for f in files:
      df.add handleFile(f, cutTab)
    var effEsc = newSeq[float]()
    var effPho = newSeq[float]()
    var nums = newSeq[int]()
    for (tup, subDf) in groups(df.group_by(@["runNumber", "Peak"])):
      echo "------------------"
      echo tup
      #echo subDf
      let eff = subDf.filter(f{idx("passLogL?") == true}).len.float / subDf.len.float
      echo "Software efficiency: ", eff
      if tup[1][1].toStr == "Escapepeak":
        effEsc.add eff
      elif tup[1][1].toStr == "Photopeak":
        effPho.add eff
        # only add in one branch
        nums.add tup[0][1].toInt
      echo "------------------"
    let dfEff = toDf({"Escapepeak" : effEsc, "Photopeak" : effPho, "RunNumber" : nums})
    echo dfEff.pretty(-1)
    let stdEsc = effEsc.standardDeviationS
    let stdPho = effPho.standardDeviationS
    let meanEsc = effEsc.mean
    let meanPho = effPho.mean
    echo "Std Escape = ", stdEsc
    echo "Std Photo = ", stdPho
    echo "Mean Escape = ", meanEsc
    echo "Mean Photo = ", meanPho
    ggplot(dfEff.gather(["Escapepeak", "Photopeak"], "Type", "Value"), aes("Value", fill = "Type")) +
      geom_histogram(bins = 20, hdKind = hdOutline, alpha = 0.5) +
      ggtitle(&"σ_escape = {stdEsc:.4f}, μ_escape = {meanEsc:.4f}, σ_photo = {stdPho:.4f}, μ_photo = {meanPho:.4f}") +
      ggsave("/tmp/software_efficiencies_cast_escape_photo.pdf", width = 800, height = 600)

    for (tup, subDf) in groups(df.group_by("Peak")):
      case tup[0][1].toStr
      of "Escapepeak": plotRefHistos(subDf, 2.9, cutTab)
      of "Photopeak": plotRefHistos(subDf, 5.9, cutTab)
      else: doAssert false, "Invalid data: " & $tup[0][1].toStr
  if fake and not real:
    var effs = newSeq[float]()
    for e in energies:
      if e > 5.9:
        echo "Warning: energy above 5.9 keV not allowed!"
        return
      df = newDataFrame()
      for f in files:
        df.add ctx.handleFakeData(f, e, cutTab)
      plotRefHistos(df, e, cutTab)

      echo "Done generating for energy ", e
      effs.add(df.filter(f{idx("passLogL?") == true}).len.float / df.len.float)
    let dfL = toDf({"Energy" : energies, "Efficiency" : effs})
    echo dfL
    ggplot(dfL, aes("Energy", "Efficiency")) +
      geom_point() +
      ggtitle("Software efficiency from 'fake' events") +
      ggsave("/tmp/fake_software_effs.pdf")
  if fake and real:
    doAssert files.len == 1, "Not more than 1 file supported!"
    let f = files[0]
    let dfCast = handleFile(f, cutTab)

    for (tup, subDf) in groups(dfCast.group_by("Peak")):
      case tup[0][1].toStr
      of "Escapepeak":
        plotRefHistos(ctx.handleFakeData(f, 2.9, cutTab), 2.9, cutTab, @[("CAST", subDf)])
      of "Photopeak":
        plotRefHistos(ctx.handleFakeData(f, 5.9, cutTab), 5.9, cutTab, @[("CAST", subDf)])
      else: doAssert false, "Invalid data: " & $tup[0][1].toStr

  if computeMeans:
    for f in files:
      df.add handleFile(f, cutTab, eccFilters = some((e: 1.53, p: 1.43)))
    computeMeans(df, cutTab)

  if matchDistributions:
    matchDistributions(files, cutTab)
  #if refPlots:
  #  plotRefHistos()

when isMainModule:
  import cligen
  dispatch main
