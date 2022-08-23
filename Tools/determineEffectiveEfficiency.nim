import std / [os, strutils, random, sequtils, stats, strformat]
import nimhdf5, unchained
import numericalnim except linspace


import ingrid / private / [likelihood_utils, hdf5_utils, ggplot_utils, geometry, cdl_cuts]
import ingrid / calibration
import ingrid / calibration / [fit_functions]
import ingrid / ingrid_types
import ingridDatabase / [databaseRead, databaseDefinitions, databaseUtils]

# cut performed regardless of logL value on the data, since transverse
# rms > 1.5 cannot be a physical photon, due to diffusion in 3cm drift
# distance
const RmsCleaningCut = 1.5

let CdlFile = "/mnt/1TB/CAST/CDL_2019/calibration-cdl-2018.h5" # "/home/basti/CastData/data/CDL_2019/calibration-cdl-2018.h5"
let RefFile = "/mnt/1TB/CAST/CDL_2019/XrayReferenceFile2018.h5" # "/home/basti/CastData/data/CDL_2019/XrayReferenceFile2018.h5"

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

proc generateFakeData(h5f: H5File, nFake: int, energy = 3.0): DataFrame =
  ## For each run generate `nFake` fake events
  let refSetTuple = readRefDsets(RefFile, yr2018)
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
  result = df.mutate(f{float: "passLogL?" ~ (block:
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

proc filterEvents(df: DataFrame, energy: float = Inf): DataFrame =
  let xrayCutsTab {.global.} = getXrayCleaningCuts()
  template applyFilters(dfI: untyped): untyped {.dirty.} =
    let minRms = xrayCuts.minRms
    let maxRms = xrayCuts.maxRms
    let maxLen = xrayCuts.maxLength
    let maxEcc = xrayCuts.maxEccentricity
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
        let dset = 5.9.toRefDset()
        let xrayCuts = xrayCutsTab[dset]
        result.add applyFilters(df)
      of "Photopeak":
        let dset = 2.9.toRefDset()
        let xrayCuts = xrayCutsTab[dset]
        result.add applyFilters(df)
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

proc handleFile(fname: string, cutTab: CutValueInterpolator): DataFrame =
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
    .filterEvents()
    .applyLogLCut(cutTab)
  result.add data
  when false:
    ggplot(result, aes("energyFromCharge")) +
      geom_histogram(bins = 200) +
      ggsave("/tmp/ugl.pdf")
  discard h5f.close()

proc handleFakeData(fname: string, energy: float, cutTab: CutValueInterpolator): DataFrame =
  let h5f = H5open(fname, "r")
  var data = generateFakeData(h5f, 5000, energy = energy)
    .filterEvents(energy)
    .applyLogLCut(cutTab)
  result = data
  discard h5f.close()

proc getIndices(dset: string): seq[int] =
  result = newSeq[int]()
  applyLogLFilterCuts(CdlFile, dset, yr2018, igEnergyFromCharge):
    result.add i

proc readCdlData(df: DataFrame, energy: float, cutTab: CutValueInterpolator): DataFrame =
  # map input fake energy to reference dataset
  let grp = energy.toRefDset()

  let passedInds = getIndices(grp)
  let h5f = H5open(RefFile, "r")
  let h5fC = H5open(CdlFile, "r")
  const xray_ref = getXrayRefTable()
  #for (i, grp) in pairs(xray_ref):

  var dfR = newDataFrame()
  for dset in IngridDsetKind:
    try:
      let d = dset.toDset()
      if d notin df: continue # skip things not in input
      ## first read data from CDL file (exists for sure)
      ## extract all CDL data that passes the cuts used to generate the logL histograms
      var cdlFiltered = newSeq[float](passedInds.len)
      let cdlRaw = h5fC[cdlGroupName(grp, "2019", d), float]
      for i, idx in passedInds:
        cdlFiltered[i] = cdlRaw[idx]
      echo "Total number of elements ", cdlRaw.len, " filtered to ", passedInds.len
      dfR[d] = cdlFiltered
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
  discard h5fC.close()

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

proc computeMeans(df: DataFrame, cutTab: CutValueInterpolator) =
  for (tup, subDf) in groups(df.group_by("Peak")):
    var dfR: DataFrame
    case tup[0][1].toStr
    of "Escapepeak":
      dfR = df.readCdlData(2.9, cutTab)
    of "Photopeak":
      dfR = df.readCdlData(5.9, cutTab)
    else: doAssert false
    echo "------------------------------ Means of ", tup[0][1], " ------------------------------"
    echo "Mean of eccentricity = ", subDf["eccentricity", float].mean, " vs CDL = ", dfR["eccentricity", float].mean()
    echo "Mean of lengthDivRmsTrans = ", subDf["lengthDivRmsTrans", float].mean, " vs CDL = ", dfR["lengthDivRmsTrans", float].mean()
    echo "Mean of fractionInTransverseRms = ", subDf["fractionInTransverseRms", float].mean, " vs CDL = ", dfR["fractionInTransverseRms", float].mean()

proc main(files: seq[string], fake = false, real = false, refPlots = false,
          computeMeans = false,
          energies: seq[float] = @[]) =
  ## given the input files of calibration runs, walks all files to determine the
  ## 'real' software efficiency for them & generates a plot
  let cutTab = calcCutValueTab(CdlFile, yr2018, igEnergyFromCharge)
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
      of "Escapepeak": plotRefHistos(df, 2.9, cutTab)
      of "Photopeak": plotRefHistos(df, 5.9, cutTab)
      else: doAssert false, "Invalid data: " & $tup[0][1].toStr
  if fake and not real:
    var effs = newSeq[float]()
    for e in energies:
      if e > 5.9:
        echo "Warning: energy above 5.9 keV not allowed!"
        return
      df = newDataFrame()
      for f in files:
        df.add handleFakeData(f, e, cutTab)
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
        plotRefHistos(handleFakeData(f, 2.9, cutTab), 2.9, cutTab, @[("CAST", subDf)])
      of "Photopeak":
        plotRefHistos(handleFakeData(f, 5.9, cutTab), 5.9, cutTab, @[("CAST", subDf)])
      else: doAssert false, "Invalid data: " & $tup[0][1].toStr

  if computeMeans:
    for f in files:
      df.add handleFile(f, cutTab)
    computeMeans(df, cutTab)
  #if refPlots:
  #  plotRefHistos()

when isMainModule:
  import cligen
  dispatch main
