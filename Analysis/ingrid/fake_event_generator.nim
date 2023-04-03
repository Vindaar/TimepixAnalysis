import std / [strformat, os, random, sequtils, stats, math, options, strutils]
import pkg / [nimhdf5, ggplotnim, unchained, seqmath]

from ingrid / calibration / fit_functions import linearFunc
import ingridDatabase / [databaseRead, databaseDefinitions, databaseUtils]
import ingrid / [calibration, ingrid_types, tos_helpers, gas_physics]

import helpers / sampling_helper

type
  FakeDataKind* = enum
    fkRemovePixels,     ## target a lower energy than input & produce events by randomly selecting subsection of pixels
    fkHighThreshold,    ## target same energy events of a higher threshold, by removing pixels below certain thresholds
    fkDiffusionFixedLoc ## target a different diffusion, by 'pulling' electrons towards the cluster center or away to fake
                        ## a higher/lower level of diffusion. This version is for a *fixed* 'conversion position' (or in
                        ## other words a distance left to drift)
    fkDiffusion         ## target a different diffusion, by 'pulling' electrons towards the cluster center or away to fake
                        ## a higher/lower level of diffusion. This version draws from an exponential distribution given
                        ## a specific absorption length


  FakeEvent = object
    valid: bool
    runNumber: int
    cluster: ClusterObject[Pix] # contains all geometry & energy
    totalCharge: float # the total charge of the cluster
    logL: float # the likelihood value

  # Parameters needed to reconstruct
  CalibInfo = object
    a, b, c, t: float
    mL, bL: float
    capacitance: FemtoFarad

  FakeDesc* = object
    case kind*: FakeDataKind
    of fkRemovePixels: discard
    of fkDiffusionFixedLoc:
      zDrift*: float ## conversion position as distance from the *readout* (not the cathode!), i.e. this is
                     ## the distance still to drift!
    of fkDiffusion:
      λ*: cm ## absorption length to use. Will use exponential distribution to draw a position for each event
             ## for this absorption length
             ## NOTE: If `λ` is not set (i.e. 0) we will compute it for CAST conditions automatically using
             ## the desired target energy!
      σT*: float # the transverse diffusion to use
      sampler*: Sampler
    of fkHighThreshold:
      threshold*: float

const RmsCleaningCut = 1.5

proc drawNewEvent(rms, energy, cX, cY: seq[float], energyTarget: float): int =
  let num = rms.len - 1
  var idx = rand(num)
  while rms[idx] >= RmsCleaningCut or
        energy[idx] < energyTarget or #(energy[idx] <= 4.5 or energy[idx] >= 7.5) or
        not inRegion(cX[idx], cY[idx], crSilver):
    idx = rand(num)
  result = idx

proc computeEnergy(h5f: H5File, pix: seq[Pix], group: string, centerChip: int,
                   calibInfo: CalibInfo): tuple[totalCharge, energy: float] =
  let totalCharge = pix.mapIt(
    calibrateCharge(
      it.ch.float,
      calibInfo.capacitance,
      calibInfo.a,
      calibInfo.b,
      calibInfo.c,
      calibInfo.t)
  ).sum
  # compute mean of all gas gain slices in this run (most sensible)
  let gain = h5f[group / &"chip_{centerChip}/gasGainSlices", GasGainIntervalResult].mapIt(it.G).mean
  let calibFactor = linearFunc(@[calibInfo.bL, calibInfo.mL], gain) * 1e-6
  # now calculate energy for all hits
  result = (totalCharge: totalCharge, energy: totalCharge * calibFactor)

proc reconstructFakeEvent(h5f: H5File, group: string, pix: seq[Pix],
                          runNumber: int,
                          centerChip: int,
                          calibInfo: CalibInfo,
                          ctx: LikelihoodContext): FakeEvent =
  let inp = (pixels: pix, eventNumber: 0, toa: newSeq[uint16](), toaCombined: newSeq[uint64]())
  let recoEv = recoEvent(inp, -1,
                         runNumber, searchRadius = 50,
                         dbscanEpsilon = 65,
                         clusterAlgo = caDefault)
  if recoEv.cluster.len > 1 or recoEv.cluster.len == 0:
    echo "Found more than 1 or 0 cluster! Skipping"
    return

  # get the actual cluster
  var cluster = recoEv.cluster[0]
  # compute charge
  let (totalCharge, energy) = h5f.computeEnergy(pix, group, centerChip, calibInfo)
  cluster.energy = energy # assign the energy

  # puhhh, now the likelihood...
  let ecc = cluster.geometry.eccentricity
  let ldiv = cluster.geometry.lengthDivRmsTrans
  let frin = cluster.geometry.fractionInTransverseRms
  let logL = ctx.calcLikelihoodForEvent(energy,
                                        ecc,
                                        ldiv,
                                        frin)

  result = FakeEvent(valid: true, runNumber: runNumber, cluster: cluster, logL: logL, totalCharge: totalCharge)

proc toDf(fakeEvs: seq[FakeEvent]): DataFrame =
  ## Serializes all information of the fake events into a data frame
  ##
  ## DF contains all InGridDsetKinds that are actually computed in TPA.
  result = newDataFrame()
  for dset in TPADatasets - {igNumClusters}:
    result[dset.toDset] = newColumn(colFloat, fakeEvs.len)
  for i, ev in fakeEvs:
    let cl = ev.cluster
    result[igHits.toDset][i]             = cl.hits
    result[igCenterX.toDset][i]          = cl.centerX
    result[igCenterY.toDset][i]          = cl.centerY
    result[igEnergyFromCharge.toDset][i] = cl.energy
    result[igTotalCharge.toDset][i]      = ev.totalCharge
    result[igLikelihood.toDset][i]       = ev.logL
    for field, val in fieldPairs(cl.geometry):
      result[field][i] = val
  result["runNumber"] = fakeEvs.mapIt(it.runNumber)

randomize()
proc getDiffusion(fakeDesc: FakeDesc): float =
  result = rand(510.0 .. 700.0)
  #if fakeDesc.σT > 0.0:
  #  result = fakeDesc.σT
  #else:
  #  ## Random diffusion!!!!
  #  result = rand(510.0 .. 700.0)

    #let df = readCsv("/home/basti/org/resources/ar_iso_97_7_2_3_septemboard_cast_different_temps.csv")
    #result = df.mutate(f{"Diff" ~ abs(idx("T [K]") - 303.0)})
    #  .arrange("Diff", SortOrder.Ascending)["σT [μm/√cm]", float][0]
  #echo "DIFFUSION: ", result
  #result = result * result
  #result = 911.0 #4696.96n
  #result = 400.0
  #result = 560.0 # <- something like 560 seems to match the CDL data at ~1 keV best :/

proc genRemovePixelsEvent(xs, ys: seq[uint8], ts: seq[uint16], energyInput: float,
                          energy: float): seq[Pix] =
  # draw number of fake pixels
  # compute ref # pixels for this event taking into account possible double counting etc.
  let basePixels = (energy / energyInput * xs.len.float)

  ## XXX: only needed for fkRemovePixels!
  let nPix = round(basePixels + gauss(sigma = basePixels * 0.1)).int  # ~115 pix as reference in 3 keV (26 eV), draw normal w/10 around
  echo "BASE PIXELS: ", basePixels, " TARGET PIX : ", nPix, " from: ", energyInput, " and nHits ", xs.len
  if nPix < 4:
    echo "Less than 4 pixels: ", nPix, " skipping"
    return @[]
  result = newSeq[Pix](nPix)
  var seenPix: set[uint16] = {}
  let evNumPix = xs.len

  if nPix >= evNumPix:
    echo "More pixels to draw than available! ", nPix, " vs ", evNumPix, ", skipping!"
    return @[]
  var pIdx = rand(evNumPix - 1)
  for j in 0 ..< nPix:
    # draw pix index
    while pIdx.uint16 in seenPix:
      pIdx = rand(evNumPix - 1)
    seenPix.incl pIdx.uint16
    result[j] = (x: xs[pIdx], y: ys[pIdx], ch: ts[pIdx])
  # now draw
  when false:
    let df = toDf({"x" : pix.mapIt(it.x.int), "y" : pix.mapIt(it.y.int), "ch" : pix.mapIt(it.ch.int)})
    ggplot(df, aes("x", "y", color = "ch")) +
      geom_point() +
      ggsave("/tmp/fake_event.pdf")
    sleep(200)

proc expFn(x: float, λ: float): float =
  ## Exponential distribution for the absorption length `λ`
  result = 1.0 / λ * exp(- x / λ)

proc genDiffusionEvent(xs, ys: seq[uint8], ts: seq[uint16], energyInput: float,
                       energy: float, fakeDesc: FakeDesc): seq[Pix] =
  var # mutable copies for lower energy event if needed
    xs = xs
    ys = ys
    ts = ts
  if energy < energyInput * 0.85: ## allow 15% range to smaller values!
    echo "GENERATING LOWER ENERGY EVENT: ", energy, " compare ", energyInput
    # call `genRemovePixelsEvent` to generate target lower energy event
    let initNumPix = xs.len
    let pix = genRemovePixelsEvent(xs, ys, ts, energyInput, energy)
    # copy over data to `xs`, `ys`, `ts`
    doAssert xs.len >= pix.len
    for i in 0 ..< pix.len:
      xs[i] = pix[i].x
      ys[i] = pix[i].y
      ts[i] = pix[i].ch
    xs.setLen(pix.len)
    ys.setLen(pix.len)
    ts.setLen(pix.len)

  let nPix = xs.len
  result = newSeq[Pix](nPix)
  ## Change the effective diffusion of the drawn event.
  ## 1. get the diffusion (transverse) coefficient
  #let σT = getDiffusion(fakeDesc)
  #echo "DIFFUSION COEFF : ", σT
  ## 2. using the existing energy and target energies compute the required distance we 'move'
  ##    the cluster from and to.
  ##    -> For now the start and stop positions (mean pos!) will be an argument, but in the future
  ##    we'll also draw from the exponential distribution that describes the likely absorption position!
  ##
  ## 3. define the a gaussian with mean of resulting diffusion around that distance
  ##    (what sigma does it have?)
  ##    -> Or: define gaussian of the transverse diffusion coefficient and simply multiply!
  # -> this code corresponds to a _conversion location_ of course! Need to adjust by plugging
  # into exponential distr to get _a_ conversion point given a conversion length
  let P1 =
    if fakeDesc.kind == fkDiffusionFixedLoc:
      (proc(): float =
        let x = gauss(mu = 0.0, sigma = fakeDesc.σT * sqrt(fakeDesc.zDrift)) #sqrt(6.0 * σT / 2.262))#σT * sqrt(zDrift))
        let y = gauss(mu = 0.0, sigma = fakeDesc.σT * sqrt(fakeDesc.zDrift)) #sqrt(6.0 * σT / 2.262))#σT * sqrt(zDrift))
        result = sqrt(x*x + y*y)
      )
    else:
      # sample a *single* conversion point for this event
      let conversionPoint = fakeDesc.sampler.sample() # point 'behind' cathode at which conversion takes place
      (proc(): float =
        let zDrift = 3.0 - conversionPoint # use conversion point to define a closure to sample for each electron
        let x = gauss(mu = 0.0, sigma = fakeDesc.σT * sqrt(zDrift)) #sqrt(6.0 * σT / 2.262))#σT * sqrt(zDrift))
        let y = gauss(mu = 0.0, sigma = fakeDesc.σT * sqrt(zDrift)) #sqrt(6.0 * σT / 2.262))#σT * sqrt(zDrift))
        result = sqrt(x*x + y*y)
      )
  ## 4. for each pixel sample it and move each pixel the resulting distance towards the
  ##    center / away from the center (depending)
  let xc = xs.mapIt(it.float).mean
  let yc = ys.mapIt(it.float).mean
  result = newSeq[Pix](nPix)
  for j in 0 ..< nPix:
    let x0 = xs[j].float
    let y0 = ys[j].float
    # draw a new distance from the center for this pixel
    let pixelPos = P1().μm
    let dist = clamp(pixelPos.to(mm) / 14.1.mm * 256, -255.0, 255.0)
    # transport (x0, y0) to (x1, y1) along the line [(x0, y0), (xc, yc)] where
    # `xc`, `yc` are the cluster center coordinates.
    # First compute distance between pixel and cluster center
    let x0c = x0 - xc; let y0c = y0 - yc # ; let m = y0c / x0c; let b = y0 / (m * x0)
    # let r = sqrt( x0c*x0c + y0c*y0c )
    let φ = arctan2(y0c, x0c)
    ## XXX: think about if we want this!
    #let φ = rand(0.0 .. 2*PI)
    let xp = dist * cos(φ) + xc; let yp = dist * sin(φ) + yc
    # compute point `dist` away from `(xc, yc)`
    result[j] = (x: xp.round.uint8, y: yp.round.uint8, ch: ts[j])

  when false:
    ## XXX: check that these events look correct!
    var count {.global.} = 0
    block:
      let xs = xs.mapIt(it.float); let ys = ys.mapIt(it.float)
      let dfOrig = seqstoDf({ "x" : xs, "y" : ys })
      let dfNew =  seqstoDf({ "x" : result.mapIt(it.x.float),
                              "y" : result.mapIt(it.y.float) })
      let df = bind_rows([("Orig", dfOrig), ("Fake", dfNew)],
                         "Type")
      ggplot(df, aes("x", "y", color = "Type")) +
        geom_point() +
        xlim(0, 256) + ylim(0, 256) +
        ggsave("/tmp/events/fake_event_" & $count & ".pdf")
      inc count

proc sanityCheckSampler(fakeDesc: FakeDesc) =
  var samples = newSeq[float]()
  for i in 0 ..< 10_000_000:
    samples.add fakeDesc.sampler.sample()
  ggplot(toDf(samples), aes("samples")) +
    geom_histogram(bins = 1000, hdKind = hdOutline) +
    ggsave("/tmp/samples.pdf", width = 1200, height = 800)

  ggplot(seqstoDf({"x" : fakeDesc.sampler.xs, "y" : fakeDesc.sampler.ys}), aes("x", "y")) +
    geom_line() +
    ggsave("/tmp/sampplerrrr.pdf")
  if true: quit()

import ../../Tools/determineDiffusion/determineDiffusion

proc generateFakeData*(ctx: LikelihoodContext, h5f: H5File, nFake: int, energy = 3.0,
                       fakeDesc: FakeDesc,
                       tfKind = none[TargetFilterKind](),
                       energyDset = igEnergyFromCharge,
                       run = -1): DataFrame =

  ## For each run generate `nFake` fake events
  ##
  ## XXX:
  ## Insert the other kinds of fake data:
  ## - stretched / compressed data (different absorption length)
  ## - same energy but higher threshold (dropping pixels below a threshold)
  result = newDataFrame()
  let fileInfo = h5f.getFileInfo()
  let tpx = fileInfo.timepix
  let centerChip = fileInfo.centerChip
  let capacitance = tpx.getCapacitance()

  var fakeDesc = fakeDesc
  if fakeDesc.kind == fkDiffusion:
    ## Generate the exponential distribution to sample from based on the
    ## given absorption length
    let λ = if fakeDesc.λ > 0.0.cm: fakeDesc.λ else: absorptionLengthCAST(energy.keV)
    let fnSample = (proc(x: float): float =
                      result = expFn(x, λ.float)
    )
    # we want to be able to sample between 0 and 3 cm
    fakeDesc.sampler = sampler(fnSample, 0.0, 3.0, num = 5000)

    when false:
      sanityCheckSampler(fakeDesc)

  for (num, group) in runs(h5f):
    if run > 0 and run != num: continue # skip if not the desired run
    if tfKind.isSome:
      let tf = tfKind.get
      let runGrp = h5f[group.grp_str]
      if "tfKind" in runGrp.attrs:
        let tfGet = parseEnum[TargetFilterKind](runGrp.attrs["tfKind", string])
        if tf != tfGet:
          echo "Skipping run : ", num, " as it does not contain tfKind ", tf, " data."
          continue

    # first read all x / y / tot data
    let xs = h5f[group / &"chip_{centerChip}/x", special_type(uint8), uint8]
    let ys = h5f[group / &"chip_{centerChip}/y", special_type(uint8), uint8]
    let ts = h5f[group / &"chip_{centerChip}/ToT", special_type(uint16), uint16]
    let rms = h5f[group / &"chip_{centerChip}/rmsTransverse", float]
    let cX = h5f[group / &"chip_{centerChip}/centerX", float]
    let cY = h5f[group / &"chip_{centerChip}/centerY", float]
    let energyInput = h5f[group / &"chip_{centerChip}/{energyDset.toDset()}", float]
    let chipGrp = h5f[(group / &"chip_{centerChip}").grp_str]
    let chipName = chipGrp.attrs["chipName", string]
    # get factors for charge calibration
    let (a, b, c, t) = getTotCalibParameters(chipName, num)
    # get factors for charge / gas gain fit
    let (bL, mL) = getCalibVsGasGainFactors(chipName, num, suffix = $gcIndividualFits)
    let calibInfo = CalibInfo(a: a, b: b, c: c, t: t, mL: mL, bL: bL, capacitance: capacitance)
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
    var fakeEvents = newSeqOfCap[FakeEvent](nFake)

    # get one diffusion value for this run
    fakeDesc.σT = getDiffusion(rms) #getDiffusion(fakeDesc)
    echo "THIS DIFFUSION: ", fakeDesc.σT
    while count < nFake:
      # draw index from to generate a fake event
      evIdx = drawNewEvent(rms, energyInput, cX, cY, energy)
      var pix: seq[Pix]
      case fakeDesc.kind
      of fkRemovePixels:
        pix = genRemovePixelsEvent(xs[evIdx], ys[evIdx], ts[evIdx],
                                   energyInput[evIdx],
                                   energy)
      of fkHighThreshold:
        doAssert false
      of fkDiffusion, fkDiffusionFixedLoc:
        pix = genDiffusionEvent(xs[evIdx], ys[evIdx], ts[evIdx],
                                energyInput[evIdx],
                                energy,
                                fakeDesc)
      if pix.len == 0:
        echo "[INFO] Skipping invalid event with 0 entries. ", count, " events done"
        continue # invalid event

      # reconstruct event
      let fakeEv = h5f.reconstructFakeEvent(group,
                                            pix, num,
                                            centerChip,
                                            calibInfo,
                                            ctx)
      if not fakeEv.valid:
        continue
      # finally done
      ## replace by better logic without so much duplicity?
      fakeEvents.add fakeEv
      inc count

    result.add toDf(fakeEvents)
