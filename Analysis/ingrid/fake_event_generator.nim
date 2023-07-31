import std / [strformat, os, random, sequtils, stats, math, strutils]
import pkg / [nimhdf5, ggplotnim, unchained, seqmath, xrayAttenuation]

from ingrid / calibration / fit_functions import linearFunc, polyaImpl
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
    fkDiffusionFromData ## target a different diffusion, by 'pulling' electrons towards the cluster center or away to fake
                        ## a higher/lower level of diffusion. This version draws from an exponential distribution given
                        ## a specific absorption length, but keeps the center positions of *real data*
    fkGainDiffusion     ## Simulates fake data purely based on the gas gain of a run and its diffusion determined by
                        ## the fit to the transverse RMS of the data.

  FakeEvent* = object
    valid*: bool
    runNumber*: int
    cluster*: ClusterObject[Pix] # contains all geometry & energy
    totalCharge*: float # the total charge of the cluster
    σT*: float
    gasGain*: float
    logL*: float # the likelihood value

  ## XXX: add target energy to be separate from target? Doesn't really make sense at least
  ## for `fkDiffusion`
  FakeDesc* = object
    nFake*: int ## Number of fake events to generate
    tfKind*: TargetFilterKind ## the desired target for which to generate fake events
    energyDset*: InGridDsetKind = igEnergyFromCharge ## dataset to use for energy
    case kind*: FakeDataKind
    of fkRemovePixels: discard
    of fkDiffusionFixedLoc:
      zDrift*: float ## conversion position as distance from the *readout* (not the cathode!), i.e. this is
                     ## the distance still to drift!
    of fkDiffusionFromData, fkGainDiffusion:
      λ*: cm ## absorption length to use. Will use exponential distribution to draw a position for each event
             ## for this absorption length
             ## NOTE: If `λ` is not set (i.e. 0) we will compute it for CAST conditions automatically using
             ## the desired target energy!
      σT*: float # the transverse diffusion to use
      gasMixture*: GasMixture
      sampler*: Sampler
      sampleTab: Table[string, Sampler] # samplers for each fluorescence line
    of fkHighThreshold:
      threshold*: float

  GainInfo* = object # a helper object that only stores information to define the
                    # gas gain information
    N*: float = 100000.0 # amplitude of Pólya to draw from (irrelevant)
    G*: float = 3500.0# the absolute gas gain
    theta*: float = 2.3 # the width of the distribution

const RmsCleaningCut = 1.5

proc drawNewEvent(rnd: var Rand, energy: seq[float], energyTarget: keV): int =
  let num = energy.len - 1
  var idx = rnd.rand(num)
  ## XXX: replace this. Instead pre select all data that is allowed to be used
  ## for new event selection
  #while #rms[idx] >= RmsCleaningCut or
  #      energy[idx] < energyTarget: # or #(energy[idx] <= 4.5 or energy[idx] >= 7.5) or
  #      #not inRegion(cX[idx], cY[idx], crSilver):
  #  idx = rnd.rand(num)
  result = idx

proc computeEnergy(pix: Pixels, calibInfo: CalibInfo, gain: float): tuple[totalCharge, energy: float] =
  let totalCharge = pix.mapIt(
    calibrateCharge(
      it.ch.float,
      calibInfo.capacitance,
      calibInfo.a,
      calibInfo.b,
      calibInfo.c,
      calibInfo.t)
  ).sum
  # use the gain to calculate the calibration factor
  let calibFactor = linearFunc(@[calibInfo.bL, calibInfo.mL], gain) * 1e-6
  # now calculate energy for all hits
  result = (totalCharge: totalCharge, energy: totalCharge * calibFactor)

type
  Missing = object
  MaybeContext = LikelihoodContext | Missing
func missing(): Missing = discard

proc reconstructFakeEvent[C: MaybeContext](
  pix: Pixels,
  σT: float,
  gasGain: float,
  runNumber: int,
  calibInfo: CalibInfo,
  ctx: C = missing()
                                         ): FakeEvent =
  let inp = (pixels: pix, eventNumber: 0, toa: newSeq[uint16](), toaCombined: newSeq[uint64]())
  let recoEv = recoEvent(inp, -1,
                         runNumber, searchRadius = 50,
                         dbscanEpsilon = 65,
                         clusterAlgo = caDefault)
  if recoEv.cluster.len > 1 or recoEv.cluster.len == 0:
    echo "Found more than 1 or 0 cluster! Skipping. Number of clusters: ", recoEv.cluster.len
    return

  # get the actual cluster
  var cluster = recoEv.cluster[0]
  # compute charge
  let (totalCharge, energy) = computeEnergy(pix, calibInfo, gasGain)
  cluster.energy = energy # assign the energy

  # puhhh, now the likelihood...
  let ecc = cluster.geometry.eccentricity
  let ldiv = cluster.geometry.lengthDivRmsTrans
  let frin = cluster.geometry.fractionInTransverseRms
  when C isnot Missing:
    let logL = ctx.calcLikelihoodForEvent(energy,
                                          ecc,
                                          ldiv,
                                          frin)
  else:
    let logL = NaN

  result = FakeEvent(valid: true, runNumber: runNumber, cluster: cluster, σT: σT, gasGain: gasGain, logL: logL, totalCharge: totalCharge)

proc fakeToDf*(fakeEvs: seq[FakeEvent]): DataFrame =
  ## Serializes all information of the fake events into a data frame
  ##
  ## DF contains all InGridDsetKinds that are actually computed in TPA.
  result = newDataFrame()
  for dset in TPADatasets - {igNumClusters}:
    result[dset.toDset] = newColumn(colFloat, fakeEvs.len)
  result["σT"] = newColumn(colFloat, fakeEvs.len)
  result["gasGain"] = newColumn(colFloat, fakeEvs.len)
  for i, ev in fakeEvs:
    let cl = ev.cluster
    result[igHits.toDset][i]             = cl.hits
    result[igCenterX.toDset][i]          = cl.centerX
    result[igCenterY.toDset][i]          = cl.centerY
    result[igEnergyFromCharge.toDset][i] = cl.energy
    result[igTotalCharge.toDset][i]      = ev.totalCharge
    result[igLikelihood.toDset][i]       = ev.logL
    result["σT"][i]                      = ev.σT / 1000.0 # from μm/√cm to mm/√cm
    result["gasGain"][i]                 = ev.gasGain / 1000.0 # to 1k electrons
    for field, val in fieldPairs(cl.geometry):
      result[field][i] = val
  result["runNumber"] = fakeEvs.mapIt(it.runNumber)

proc getDiffusion(rnd: var Rand, fakeDesc: FakeDesc): float =
  result = rnd.rand(510.0 .. 700.0)
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

proc genRemovePixelsEvent(rnd: var Rand, xs, ys: seq[uint8], ts: seq[uint16], energyInput: keV,
                          energy: keV): Pixels =
  # draw number of fake pixels
  # compute ref # pixels for this event taking into account possible double counting etc.
  let basePixels = (energy / energyInput * xs.len.float)

  ## XXX: only needed for fkRemovePixels!
  let nPix = round(basePixels + rnd.gauss(sigma = basePixels * 0.1)).int  # ~115 pix as reference in 3 keV (26 eV), draw normal w/10 around
  #echo "BASE PIXELS: ", basePixels, " TARGET PIX : ", nPix, " from: ", energyInput, " and nHits ", xs.len
  if nPix < 4:
    echo "Less than 4 pixels: ", nPix, " skipping"
    return @[]
  result = newSeq[Pix](nPix)
  var seenPix: set[uint16] = {}
  let evNumPix = xs.len

  if nPix >= evNumPix:
    #echo "More pixels to draw than available! ", nPix, " vs ", evNumPix, ", skipping!"
    discard
    #return @[]
  #var pIdx = rand(evNumPix - 1)
  var pIdx: int # = rand(nPix - 1)
  var toDraw = toSeq(0 ..< evNumPix) # to draw is all valid indices of the original event
  for j in 0 ..< nPix:
    # draw pix index from `toDraw`
    when false:
      while pIdx.uint16 in seenPix:
        pIdx = rnd.rand(evNumPix - 1)
      seenPix.incl pIdx.uint16
    else:
      let dIdx = rnd.rand(0 .. toDraw.high) # the draw index
      pIdx = toDraw[dIdx]
      toDraw.del(dIdx) # remove this index
    # insert pix
    result[j] = (x: xs[pIdx], y: ys[pIdx], ch: ts[pIdx])
    if toDraw.len == 0: # drawn all pixels, but need more. Allow every pixel again
                        # it's important that we only do this when using diffusion events!
      toDraw = toSeq(0 ..< evNumPix)
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

proc expShift(x: float, u, l, p: float): float =
  ## Shifted exponential distribution that satisfies:
  ## f(l) = p, f(u) = 1.0
  let λ = (u - l) / ln(p)
  result = exp(- (x - u) / λ)

proc genDiffusionEvent(rnd: var Rand, xs, ys: seq[uint8], ts: seq[uint16], energyInput: keV,
                       energy: keV, fakeDesc: FakeDesc): Pixels =
  # first compute center of the *original* cluster!
  let xc = xs.mapIt(it.float).mean
  let yc = ys.mapIt(it.float).mean
  var # mutable copies for lower energy event if needed
    xs = xs
    ys = ys
    ts = ts
  if false: # true: #energy < energyInput * 0.85: ## allow 15% range to smaller values!
    #echo "GENERATING LOWER ENERGY EVENT: ", energy, " compare ", energyInput
    # call `genRemovePixelsEvent` to generate target lower energy event
    let initNumPix = xs.len
    let pix = rnd.genRemovePixelsEvent(xs, ys, ts, energyInput, energy)
    # copy over data to `xs`, `ys`, `ts`
    xs.setLen(pix.len)
    ys.setLen(pix.len)
    ts.setLen(pix.len)
    for i in 0 ..< pix.len:
      xs[i] = pix[i].x
      ys[i] = pix[i].y
      ts[i] = pix[i].ch

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
      (proc(rnd: var Rand): float =
        let x = rnd.gauss(mu = 0.0, sigma = fakeDesc.σT * sqrt(fakeDesc.zDrift)) #sqrt(6.0 * σT / 2.262))#σT * sqrt(zDrift))
        let y = rnd.gauss(mu = 0.0, sigma = fakeDesc.σT * sqrt(fakeDesc.zDrift)) #sqrt(6.0 * σT / 2.262))#σT * sqrt(zDrift))
        result = sqrt(x*x + y*y)
      )
    else:
      # sample a *single* conversion point for this event
      let conversionPoint = rnd.sample(fakeDesc.sampler) # point 'behind' cathode at which conversion takes place
      (proc(rnd: var Rand): float =
        let zDrift = 3.0 - conversionPoint # use conversion point to define a closure to sample for each electron
        let x = rnd.gauss(mu = 0.0, sigma = fakeDesc.σT * sqrt(zDrift)) #sqrt(6.0 * σT / 2.262))#σT * sqrt(zDrift))
        let y = rnd.gauss(mu = 0.0, sigma = fakeDesc.σT * sqrt(zDrift)) #sqrt(6.0 * σT / 2.262))#σT * sqrt(zDrift))
        result = sqrt(x*x + y*y)
      )
  ## 4. for each pixel sample it and move each pixel the resulting distance towards the
  ##    center / away from the center (depending)
  #let xc = xs.mapIt(it.float).mean
  #let yc = ys.mapIt(it.float).mean
  ## XXX: Maybe our data with cuts applied (readCdlCutsIdx thing) has a bias that is not
  ## present in the distribution we plot the fake data against?
  result = newSeq[Pix](nPix)
  for j in 0 ..< nPix:
    let x0 = xs[j].float
    let y0 = ys[j].float
    # draw a new distance from the center for this pixel
    let pixelPos = P1(rnd).μm
    let dist = clamp(pixelPos.to(mm) / 14.1.mm * 256, -255.0, 255.0)
    # transport (x0, y0) to (x1, y1) along the line [(x0, y0), (xc, yc)] where
    # `xc`, `yc` are the cluster center coordinates.
    # First compute distance between pixel and cluster center
    let x0c = x0 - xc; let y0c = y0 - yc # ; let m = y0c / x0c; let b = y0 / (m * x0)
    # let r = sqrt( x0c*x0c + y0c*y0c )
    #let φ = arctan2(y0c, x0c)
    ## XXX: think about if we want this!
    let φ = rnd.rand(0.0 .. 2*PI)
    let xp = dist * cos(φ) + xc; let yp = dist * sin(φ) + yc
    # compute point `dist` away from `(xc, yc)`
    result[j] = (x: xp.round.uint8, y: yp.round.uint8, ch: ts[j])

func toIdx[L: Length | float](arg: L): int =
  when L is float:
    let argMm = arg.mm
  else:
    let argMm = arg.to(mm)
  result = (argMm / 14.1 * 256.0).float.round.int.clamp(0, 255)

proc sampleInSilver(rnd: var Rand): (int, int) =
  ## Returns a coordinate in the silver region
  let reg = getRegionCut(crSilver)
  var radius = Inf
  var x, y: float
  while radius > reg.radius:
    x = rnd.rand(-4.5 .. 4.5)
    y = rnd.rand(-4.5 .. 4.5)
    radius = sqrt(x*x + y*y)
  result = (toIdx (x + TimepixSize / 2.0), toIdx (y + TimepixSize / 2.0))

proc initPolyaSampler(params: seq[float],
                      frm = 0.0, to = 20000.0): Sampler =
  ## Helper to sample from Pólya distribution with the given parameters
  ## `(N, G, θ)`
  let fnSample = (
    proc(x: float): float =
      result = polyaImpl(params, x)
  )
  # we want to be able to sample between 0 and 20k
  result = sampler(fnSample, frm, to, num = 1000)

proc sampleTargetCharge(rnd: var Rand, targetEnergyCharge: float): float =
  ## Samples a target charge for this event based on the given charge
  ## `targetEnergyCharge` which is the equivalent energy that corresponds
  ## to the desired taget energy in charge based on the gas gain of the
  ## current data slice.
  ##
  ## We sample a normal distribution around the target with a sigma that
  ## corresponds to 10%.
  result = rnd.gauss(mu = targetEnergyCharge, sigma = targetEnergyCharge * 0.075)

proc genGainDiffusionEvent(rnd: var Rand, gain: GainInfo,
                           calibInfo: CalibInfo,
                           energy: keV, fakeDesc: FakeDesc): Pixels =
  ## 0. compute the equivalent charge of the target energy using the "gas gain vs energy calibration" relation
  ##   using `E = totalCharge · calibrationFactor`
  ##   with  `calibrationFactor = m · gain + b` with `m, b` determined by the gas gain vs energy calib fit.
  let gainToUse = gain.G * 0.9

  let calibFactor = linearFunc(@[calibInfo.bL, calibInfo.mL], gain.G) * 1e-6
  let targetEnergyCharge = energy.float / calibFactor
  #echo "TARGET CHARGE: ", targetEnergyCharge
  #echo "Calibration afctor: ", calibFactor, " for gain ", gain.G, " and energy ", energy, " is ", targetEnergyCharge
  ## 1. Sample a conversion point and generate a sampler for diffused pixels for that point
  let P1 = block:
    # sample a *single* conversion point for this event
    let conversionPoint = rnd.sample(fakeDesc.sampler) # point 'behind' cathode at which conversion takes place
    (proc(rnd: var Rand): float =
      let zDrift = 3.0 - conversionPoint # use conversion point to define a closure to sample for each electron
      let x = rnd.gauss(mu = 0.0, sigma = fakeDesc.σT * sqrt(zDrift)) #sqrt(6.0 * σT / 2.262))#σT * sqrt(zDrift))
      let y = rnd.gauss(mu = 0.0, sigma = fakeDesc.σT * sqrt(zDrift)) #sqrt(6.0 * σT / 2.262))#σT * sqrt(zDrift))
      result = sqrt(x*x + y*y)
    )

  template invert(x, cInfo: untyped): untyped =
    invertCharge(x, cInfo.capacitance,
                 cInfo.a, cInfo.b, cInfo.c, cInfo.t).float
  template calib(x, cInfo: untyped): untyped =
    calibrateCharge(x, cInfo.capacitance,
                    cInfo.a, cInfo.b, cInfo.c, cInfo.t).float

  ## 2. Sample a target charge from the equivalent charge that corresponds to our
  ##   target energy
  let targetCharge = rnd.sampleTargetCharge(targetEnergyCharge)
  ## 3. Define a sampler to sample from the Pólya distribution required for the
  ##   current gain parameters
  #let gInv = invert(gain.G * 0.6, calibInfo)
  #let gInv = invert(gain.G, calibInfo)
  #let params = @[gain.N, gInv, gain.theta]
  let params = @[gain.N, gainToUse, gain.theta / 3.0]
  let psampler = initPolyaSampler(params, frm = 0.0, to = 20000.0) #invert(20000.0, calibInfo))
  ## 4. Sample a target center x and center y position within crSilver
  let (xc, yc) = rnd.sampleInSilver()
  ## 5. Generate pixels until desired target charge accumulated
  var totalCharge = 0.0
  # store pixels in a table for now to detect and deal with duplicate hits!
  var pixTab = initTable[(uint8, uint8), float]()

  ## XXX: make dependent on Run-2 or Run-3 data!
  ## XXX: These cutoffs are not super important because they correspond to a
  ## ToT value of ~0...
  const cutoffs = @[1027.450870326596, # Run-2
                    893.4812944899318] # Run-3
  let cutoff = cutoffs[1] * 1.15

  ## XXX: for the time being we simply always accept any ToT value above cutoff!
  let actSampler = (proc(rnd: var Rand, x: float): bool =
                      let activateThreshold = expShift(x, gainToUse, cutoff, 0.3)
                      result = x > cutoff # and rnd.rand(1.0) < activateThreshold
  )
  let neighSampler = (proc(rnd: var Rand, x: float): int =
    ## Returns the number of neighbors to activate
    ## Slope and offset set: f(1000) = 0.0, f(10000) = 1.0 (outside none and always)
    let m = 1.0 / 9000.0
    let activateThreshold = m * x - 1000 * m
    let val = rnd.rand(1.0)
    if val * 4.0 < activateThreshold:
      result = 4
    elif val * 3.0 < activateThreshold:
      result = 3
    elif val * 2.0 < activateThreshold:
      result = 2
    elif val       < activateThreshold:
      result = 1
    else:
      result = 0
  )

  ## 6. Begin the generation of individual pixels in the cluster
  while totalCharge < targetCharge:
    # sample a charge from Pólya first (to know if to continue)
    #let ToT = rnd.sample(psampler)
    #let charge = calib(ToT, calibInfo)
    let charge = rnd.sample(psampler)
    let ToT = invert(charge, calibInfo)

    if abs(targetCharge - (totalCharge + charge)) > abs(targetCharge - totalCharge):
      break # stop early, do not use this charge
    # cutoff: pixel threshold!
    #echo "charge ", charge, " from ToT ", ToT, " for 1 is ", calib(1, calibInfo), " for cutoff: ", cutoff, " inverted cutoff: ", invert(cutoff, calibInfo)
    let activatePixel = actSampler(rnd, charge)
    if not activatePixel: #charge < cutoff:
      #totalCharge += charge # increase total charge anyway! We simply didn't record it, but was there
      #                      # *BUT* the total charge we target *ALSO* did not record these!
      continue
    # draw a new distance from the center for this pixel
    let pixelPos = abs(P1(rnd).μm) # take `abs` as we only care about distance from center
    let dist = (toIdx pixelPos).float # convert to pixel index, but use as floa
    # sample a random angle
    let φ = rnd.rand(0.0 .. 2*PI)
    # compute the new x, y positions based on distance of this pixel from center
    let xp = (dist * cos(φ) + xc).round.uint8
    let yp = (dist * sin(φ) + yc).round.uint8

    proc insert(x, y: uint8, tot: float) =
      if tot > 0.0: # do not insert if ToT is 0
        if (x, y) notin pixTab:
          pixTab[(x, y)] = 0.0
        pixTab[(x, y)] += tot

    ## Insert the pixel
    insert(xp, yp, ToT)
    totalCharge += charge
    ## Now possibly generate some neighbors
    let numNeighbors = neighSampler(rnd, charge)
    if numNeighbors > 0: # charge > 3500.0: # whatever
      # possibly activate a neighbor pixel!
      #let activateNeighbor = rnd.rand(1.0) < 0.5
      var count = 0
      type Neighbor = enum Right, Left, Up, Down
      var seen: array[Neighbor, bool]
      while count < numNeighbors:
        let chargeNeighbor = rnd.sample(psampler) / 3.0 # reduce amount
        #let activateNeighbor = actSampler(rnd, chargeNeighbor)
        let neighbor = block:
          var num = Neighbor(rnd.rand(3))
          while seen[num]:
            num = Neighbor(rnd.rand(3)) # [right, left, up, down]
          seen[num] = true
          num
        #let chargeNeighbor = rnd.sample(psampler) / 2.0 # reduce amount
        let totNeighbor = invert(chargeNeighbor, calibInfo)
        case neighbor
        of Right: insert(xp + 1, yp,     totNeighbor)
        of Left:  insert(xp - 1, yp,     totNeighbor)
        of Up:    insert(xp,     yp + 1, totNeighbor)
        of Down:  insert(xp,     yp - 1, totNeighbor)
        if totNeighbor > 0.0:
          totalCharge += chargeNeighbor
        inc count
  ## Finally: convert the table to our pixel sequence and convert the charge back
  ## into ToT counts!
  result = newSeq[Pix](pixTab.len)
  var i = 0
  for pix, tot in pairs(pixTab):
    ## XXX: we need to think about whether we should sample in charge space or instead in
    ## ToT space. In the former case we need to convert the charge to ToT or have a
    ## way to reconstruct events containing `(uint8, uint8, float)` (which should be fine
    ## in principle, but requires care in the types allowed)
    ## For now we try back conversion no ToT.
    #doAssert charge < uint16.high.float, " charge was: " & $charge & " and the pixtab " & $pixTab & " of len " & $pixTab.len & " at diffusion " & $fakeDesc.σT
    #result[i] = (x: pix[0], y: pix[1], ch: charge.uint16) #tot)
    #let tot = invert(charge, calibInfo)
    #doAssert tot < uint16.high.float, " value is ? " & $tot # & " from charge " & $charge
    if tot < uint16.high.float:
      result[i] = (x: pix[0], y: pix[1], ch: tot.uint16)
    else:
      result[i] = (x: pix[0], y: pix[1], ch: 11810.uint16)
    inc i

proc sanityCheckSampler(rnd: var Rand, fakeDesc: FakeDesc, run: int) =
  var samples = newSeq[float]()
  for i in 0 ..< 10_000_000:
    samples.add rnd.sample(fakeDesc.sampler)
  ggplot(toDf(samples), aes("samples")) +
    geom_histogram(bins = 1000, hdKind = hdOutline) +
    ggsave(&"/tmp/samples_{run}.pdf", width = 1200, height = 800)

  ggplot(seqstoDf({"x" : fakeDesc.sampler.xs, "y" : fakeDesc.sampler.ys}), aes("x", "y")) +
    geom_line() +
    ggsave(&"/tmp/sampplerrrr_{run}.pdf")
  #if true: quit()

import ../../Tools/determineDiffusion/determineDiffusion

proc getFluorescenceLines(tfKind: TargetFilterKind): seq[FluorescenceLine] =
  case tfKind
  of tfCEpic0_6:  result = getFluorescenceLines(Carbon.init())
  of tfCuEpic0_9: result = getFluorescenceLines(Oxygen.init())
  of tfCuEpic2:   result = getFluorescenceLines(Copper.init()).filterIt(it.energy < 1.0.keV)
  of tfAlAl4:     result = getFluorescenceLines(Aluminium.init())
  of tfAgAg6:     result = getFluorescenceLines(Silver.init()).filterIt(it.energy < 4.0.keV)
  of tfTiTi9:     result = getFluorescenceLines(Titanium.init()).filterIt(it.energy > 1.0.keV)
  of tfMnCr12:    result = getFluorescenceLines(Manganese.init()).filterIt(it.energy > 1.0.keV)
  of tfCuNi15:    result = getFluorescenceLines(Copper.init()).filterIt(it.energy > 1.0.keV)

proc toIntensityScaled(lines: seq[FluorescenceLine]): seq[FluorescenceLine] =
  result = newSeqOfCap[FluorescenceLine](200)
  for l in lines:
    for j in 0 ..< l.intensity.round.int: # intensities are essentially all ints anyway
      result.add l

proc sampleFluorescenceLine(rnd: var Rand, lines: seq[FluorescenceLine]): FluorescenceLine =
  ## We'll do the sampling based on duplicating the informatino in the input
  ## seq `intensity` times for each line, then uniform sampling from the result
  ##
  ## The argument should already be 'intensity scaled' like that
  doAssert lines.len == 1 or lines.len > 20, "Input is not 'intensity scaled'. Call `toIntensityScaled`"
  # draw uniform sample of all indices
  let idx = rnd.rand(0 ..< lines.len)
  result = lines[idx]

proc generateSampler(fakeDesc: FakeDesc, targetEnergy: keV): Sampler =
  ## Generate the exponential distribution to sample from based on the
  ## given absorption length
  let λ = if fakeDesc.λ > 0.0.cm: fakeDesc.λ
          else: absorptionLengthCAST(fakeDesc.gasMixture, targetEnergy) #  / 2.0
  let fnSample = (proc(x: float): float =
                    result = expFn(x, λ.float)
  )
  # we want to be able to sample between 0 and 3 cm
  ## IMPORTANT NOTE: The fact that we only sample between 0 and 3 cm means we do not
  ## correctly recover the X-ray energy distributio one might expect for a given
  ## energy assuming the absorption probability. For that we would have to sample above
  ## 3 cm, too and simply treat it as an invalid cluster in case the sampled value is
  ## above 3 cm (outside chamber). But we don't care about that here!
  echo "USING λ == ", λ, " to sample! Target energy ", targetEnergy, " and gasMix ", fakeDesc.gasMixture, " gives ", absorptionLengthCAST(fakeDesc.gasMixture, targetEnergy)

  result = sampler(fnSample, 0.0, 3.0, num = 1000)

proc assignSampler(rnd: var Rand, fakeDesc: var FakeDesc, targetEnergy: keV, lines: seq[FluorescenceLine]): keV =
  if fakeDesc.kind in {fkDiffusionFromData, fkGainDiffusion}:
    ## Draw a fluorescence line from which to sample from
    let line = rnd.sampleFluorescenceLine(lines)
    ## Assign the target energy
    result = line.energy
    ## Generate the sampler sampling the conversion point
    if line.name notin fakeDesc.sampleTab:
      if fakeDesc.gasMixture.gases.len == 0:
        fakeDesc.gasMixture = initCASTGasMixture()
      fakeDesc.sampleTab[line.name] = fakeDesc.generateSampler(result)
    ## Assign the (possibly new) sampler as the sampler for this event!
    fakeDesc.sampler = fakeDesc.sampleTab[line.name]
    #echo "Sampling from line: ", line
    when false: # true:
      rnd.sanityCheckSampler(fakeDesc, runNumber)

proc generateRemoveOrDiff(rnd: var Rand, fakeDesc: var FakeDesc, lines: seq[FluorescenceLine],
                          xs, ys: seq[seq[uint8]], ts: seq[seq[uint16]],
                          targetEnergy: keV, energyInput: seq[float]
                          ): Pixels =
  # assign the correct target sampler based on fluorescence lines of the target
  # (modifying `targetEnergy` as well)
  let targetEnergy = rnd.assignSampler(fakeDesc, targetEnergy, lines)
  # draw index from to generate a fake event
  let evIdx = rnd.drawNewEvent(energyInput, targetEnergy)
  var pix: Pixels
  case fakeDesc.kind
  of fkRemovePixels:
    result = rnd.genRemovePixelsEvent(xs[evIdx], ys[evIdx], ts[evIdx],
                                      energyInput[evIdx].keV,
                                      targetEnergy)
  of fkDiffusionFromData, fkDiffusionFixedLoc:
    result = rnd.genDiffusionEvent(xs[evIdx], ys[evIdx], ts[evIdx],
                                   energyInput[evIdx].keV,
                                   targetEnergy,
                                   fakeDesc)
  of fkHighThreshold, fkGainDiffusion:
    doAssert false

type
  ReturnType = DataFrame | seq[Pixels]

proc generateFakeEvent[T: ReturnType, C: MaybeContext](
  rnd: var Rand,
  fakeDesc: FakeDesc,
  runNumber: int,
  xs, ys: seq[seq[uint8]], ts: seq[seq[uint16]],
  targetEnergy: keV, energyInput: seq[float],
  calibInfo: CalibInfo,
  _: typedesc[T],
  ctx: C = missing()
                                              ): T =
  ## Generates fake events based on the input data
  var count = 0
  # to store fake data
  var fakeEvents = newSeqOfCap[FakeEvent](fakeDesc.nFake)

  let lines = getFluorescenceLines(fakeDesc.tfKind)
    .toIntensityScaled()

  # mutable copy for sampling
  var fakeDesc = fakeDesc
  when T is Pixels:
    result = newSeqOfCap[Pix](fakeDesc.nFake)
  while count < fakeDesc.nFake:
    let pix = generateRemoveOrDiff(rnd, fakeDesc, lines, xs, ys, ts, targetEnergy, energyInput)
    if pix.len == 0:
      echo "[INFO] Skipping invalid event with 0 entries. ", count, " events done"
      continue # invalid event
    when T is Pixels:
      result.add pix
    else:
      # reconstruct event
      doAssert false, "THISCODENEEDS TO BE FIXED"
      let fakeEv = reconstructFakeEvent(pix,
                                        fakeDesc.σT,
                                        0.0, ## XXX: FIX ME FROM REAL DATA
                                        runNumber,
                                        calibInfo,
                                        ctx)
      if not fakeEv.valid:
        continue
      when false:
        # plotFakeEvent()
        discard
      fakeEvents.add fakeEv
      inc count
  when T is DataFrame:
    result = fakeToDf(fakeEvents)

proc generateFakeFromData[T: ReturnType, C: MaybeContext](
  rnd: var Rand, h5f: H5File,
  run, chipNumber: int,
  calibInfo: CalibInfo,
  fakeDesc: FakeDesc,
  runType: RunTypeKind,
  _: typedesc[T],
  ctx: C = missing()
                                       ): T =
  let group = (recoBase() & $run)
  let targetEnergy = fakeDesc.tfKind.toXrayLineEnergy().keV
  # filter to clusters that pass our preselection cuts (applies the Xray cleaning cuts)
  let passIdx =
    if true:  # isCdlFile:
      ## XXX: this must be the energy of the main fluorescence line!
      let bins = concat(@[0.0], getEnergyBinning())
      let eIdx = bins.lowerBound(targetEnergy.float)
      let eMin = bins[eIdx-1]
      let eMax = bins[eIdx]
      echo "MIN ENERGY: ", eMin, " to ", eMax, "\n\n\n"
      getCdlCutIdxs(h5f, run, chipNumber, fakeDesc.tfKind,
                    eMin, eMax, fakeDesc.energyDset)
    else:
      getCdlCutIdxs(h5f, run, chipNumber, fakeDesc.tfKind)
  # now need `rms` to calculate diffusion
  let rms = h5f[group / &"chip_{chipNumber}/rmsTransverse", passIdx, float]

  # get one diffusion value for this run
  var fakeDesc = fakeDesc
  let isBackground = runType == rtBackground
  var loss: float
  (fakeDesc.σT, loss) = getDiffusion(rms, isBackground, run) #getDiffusion(fakeDesc)
  echo "THIS DIFFUSION: ", fakeDesc.σT, " has # events : ", passIdx.len

  # and energy of each cluster
  let energyInput = h5f[group / &"chip_{chipNumber}/{fakeDesc.energyDset.toDset()}", passIdx, float]

  # now read x, y, charge data
  let xs = h5f[group / &"chip_{chipNumber}/x", passIdx, special_type(uint8), uint8]
  let ys = h5f[group / &"chip_{chipNumber}/y", passIdx, special_type(uint8), uint8]
  let ts = h5f[group / &"chip_{chipNumber}/ToT", passIdx, special_type(uint16), uint16]
  # use existing data to generate fake data
  result = rnd.generateFakeEvent(fakeDesc, run, xs, ys, ts, targetEnergy, energyInput, calibInfo, T, ctx)

  when false:
    let cy = h5f[group / &"chip_{chipNumber}/centerY", passIdx, float]
    ggplot(toDf(cy), aes("cy")) +
      geom_histogram(bins = 100) +
      ggsave("/t/read_data_in_fake_gen.pdf")
    let cysFromRaw = ys.mapIt(it.mapIt(it.float).mean)
    ggplot(toDf(cysFromRaw), aes("cysFromRaw")) +
      geom_histogram(bins = 100) +
      ggsave("/t/read_data_in_fake_gen_from_raw.pdf")

var rms = newSeq[float]()
proc generateGainDiffusion*(rnd: var Rand,
                            fakeDesc: var FakeDesc,
                            lines: seq[FluorescenceLine],
                            gainInfo: GainInfo,
                            calibInfo: CalibInfo,
                            targetEnergy: keV,
                            runNumber = -1): Pixels =
  # assign the correct target sampler based on fluorescence lines of the target
  # (modifyig `targetEnergy` as well)
  let targetEnergy = rnd.assignSampler(fakeDesc, targetEnergy, lines)
  # using correct sampler & target energy to generate a new diffusion event
  result = rnd.genGainDiffusionEvent(gainInfo,
                                     calibInfo,
                                     targetEnergy,
                                     fakeDesc)
  let xps = result.mapIt(applyPitchConversion(it.x.float, 0.0, 256)[0])
  rms.add(xps.standardDeviation())

proc generateAndReconstruct*[C: MaybeContext](rnd: var Rand,
                             fakeDesc: var FakeDesc,
                             lines: seq[FluorescenceLine],
                             gainInfo: GainInfo,
                             calibInfo: CalibInfo,
                             targetEnergy: keV,
                             ctx: C = missing(),
                             runNumber = -1): FakeEvent =
  let pix = generateGainDiffusion(rnd, fakeDesc, lines, gainInfo, calibInfo, targetEnergy, runNumber)
  if pix.len == 0:
    echo "[INFO] Skipping invalid event with 0 entries. "
    return FakeEvent(valid: false) # invalid event
  # reconstruct event
  result = reconstructFakeEvent(pix,
                                fakeDesc.σT,
                                gainInfo.G * 0.9, # * 0.75, #  * 0.85, ## XXX: FIX ME
                                runNumber,
                                calibInfo,
                                ctx)

proc generateAndReconstruct*[C: MaybeContext](rnd: var Rand,
                             fakeDesc: var FakeDesc,
                             nFake: int,
                             lines: seq[FluorescenceLine],
                             gain: GainInfo,
                             calibInfo: CalibInfo,
                             targetEnergy: keV,
                             ctx: C = missing(),
                             runNumber = -1): seq[FakeEvent] =
  ## XXX: ADD RUN NUMBER HERE!
  result = newSeqOfCap[FakeEvent](nFake)
  var count = 0
  while count < nFake:
    let fakeEv = rnd.generateAndReconstruct(fakeDesc, lines, gain, calibInfo, targetEnergy, ctx, runNumber)
    if not fakeEv.valid:
      continue
    result.add fakeEv
    inc count
  let df = toDf(rms)
  ggplot(df, aes("rms")) +
    geom_histogram(bins = 100, hdKind = hdOutline) +
    ggsave("/home/basti/Sync/histo_rms_from_fake.pdf")

proc generateRawGainDiff*(rnd: var Rand,
                          fakeDesc: var FakeDesc,
                          nFake: int,
                          lines: seq[FluorescenceLine],
                          gain: GainInfo,
                          calibInfo: CalibInfo,
                          targetEnergy: keV,
                          runNumber = -1): seq[Pixels] =
  ## XXX: ADD RUN NUMBER HERE!
  result = newSeqOfCap[Pixels](nFake)
  var count = 0
  while count < nFake:
    let pix = rnd.generateGainDiffusion(fakeDesc, lines, gain, calibInfo, targetEnergy, runNumber)
    if pix.len == 0: continue
    result.add pix
    inc count

proc generateFullFakeEvents[T: ReturnType, C: MaybeContext](
  rnd: var Rand,
  fakeDesc: FakeDesc,
  runNumber: int,
  targetEnergy: keV,
  gains: seq[GasGainIntervalResult],
  calibInfo: CalibInfo,
  _: typedesc[T],
  ctx: C = missing()
                                         ): T =
  ## This proc generates fake events using the gas gain and diffusion determined by
  ## the `rmsTransverse` data as a starting point. The events
  let lines = getFluorescenceLines(fakeDesc.tfKind)
    .toIntensityScaled()

  # mutable copy for sampling
  var fakeDesc = fakeDesc
  var targetEnergy = targetEnergy # mutable to use energies of all fluorescence lines

  ## Need to sample nFake div slices events!
  ## Q: More important to sample exactly nFake or same number for each slice?
  doAssert fakeDesc.kind == fkGainDiffusion
  let numSlices = gains.len
  when T is DataFrame:
    # to store fake data
    var fakeEvents = newSeqOfCap[FakeEvent](fakeDesc.nFake)
  else:
    result = newSeqOfCap[Pixels](fakeDesc.nFake)
  for gain in gains:
    let gainInfo = GainInfo(N: gain.N, G: gain.G, theta: gain.theta)
    let nFakeThisGain = fakeDesc.nFake div numSlices
    when T is DataFrame:
      fakeEvents.add rnd.generateAndReconstruct(fakeDesc, nFakeThisGain, lines, gainInfo, calibInfo, targetEnergy, ctx, runNumber)
    else:
      result.add rnd.generateRawGainDiff(fakeDesc, nFakeThisGain, lines, gainInfo, calibInfo, targetEnergy, runNumber)
  when T is DataFrame:
    result = fakeToDf(fakeEvents)

proc generateFakeFromGain[T: ReturnType, C: MaybeContext](
  rnd: var Rand, h5f: H5File,
  run, chipNumber: int,
  calibInfo: CalibInfo,
  fakeDesc: FakeDesc,
  runType: RunTypeKind,
  _: typedesc[T],
  ctx: C = missing(),
  useCache = true
                          ): T =
  ## Generates fake events for the desired CDL target using the diffusion and gas gains
  ## as determined by the input data.
  let group = (recoBase() & $run)
  let targetEnergy = fakeDesc.tfKind.toXrayLineEnergy().keV
  # need rmsTransverse to determine diffusion
  let grp = h5f[(group / "chip_" & $chipNumber).grp_str]
  let passIdx = if run > 306:
      ## XXX: this must be the energy of the main fluorescence line!
      let bins = concat(@[0.0], getEnergyBinning())
      let eIdx = bins.lowerBound(targetEnergy.float)
      let eMin = bins[eIdx-1]
      let eMax = bins[eIdx]
      echo "MIN ENERGY: ", eMin, " to ", eMax, "\n\n\n"
      getCdlCutIdxs(h5f, run, chipNumber, fakeDesc.tfKind,
                    eMin, eMax, fakeDesc.energyDset)
    else: # very weak cuts that should also apply well for background data that remove
                  # obvious events that are noisy, and not well defined photons or tracks!
      h5f.cutOnProperties(grp,
                          crSilver,
                          ("width", 0.0, 6.0),
                          ("rmsTransverse", 0.0, 1.3),
                          ("hits", 20.0, Inf)) # to exclude noise clusters in Run-2 data
  # now need `rms` to calculate diffusion
  let rms = h5f[group / &"chip_{chipNumber}/rmsTransverse", passIdx, float]

  #let rms = h5f[group / &"chip_{chipNumber}/rmsTransverse", float]
  # get one diffusion value for this run
  var fakeDesc = fakeDesc
  let isBackground = runType == rtBackground
  ## XXX: WARNING the `targetEnergy` here is not strictly the correct energy to use.
  ## THe energy mainly depends on the energy that best produces the `rmsT` data in
  ## the real data! We only care about the enegy in the context of `getDiffusion` to
  ## simulate similar data to get σT. Nothing else.
  var loss: float
  (fakeDesc.σT, loss) = getDiffusion(rms, isBackground, run, targetEnergy, useCache)
  echo "Using a diffusion value of ", fakeDesc.σT
  # now extract gas gain values and polya fit parameters
  let gains = h5f[group / &"chip_{chipNumber}/gasGainSlices", GasGainIntervalResult]
  # use existing data to generate fake data
  result = rnd.generateFullFakeEvents(fakeDesc, run, targetEnergy, gains, calibInfo, T, ctx)

proc generateRunFakeData*[T: ReturnType, C: MaybeContext](
  rnd: var Rand, h5f: H5File,
  run, chipNumber: int,
  chipName: string,
  capacitance: FemtoFarad,
  fakeDesc: FakeDesc,
  runType: RunTypeKind,
  _: typedesc[T],
  ctx: C = missing(),
  useCache = true
                                       ): T =
  ## Generates fake data for the given `run` in the input file `h5f` using the
  ## given `fakeDesc` and the desired target `tfKind`.
  ##
  let group = (recoBase() & $run)
  let runGrp = h5f[group.grp_str]
  var isCdlFile = false
  if "tfKind" in runGrp.attrs:
    isCdlFile = true
    let tfGet = parseEnum[TargetFilterKind](runGrp.attrs["tfKind", string])
    if fakeDesc.tfKind != tfGet:
      echo "Skipping run : ", run, " as it does not contain tfKind ", fakeDesc.tfKind, " data, but ", tfGet, " instead."
      when T is DataFrame:
        return newDataFrame()
      else:
        return @[]

  let calibInfo = h5f.initCalibInfo(runNumber = run,
                                    chipName = chipName,
                                    chipNumber = chipNumber,
                                    capacitance = capacitance)
  if fakeDesc.kind in {fkRemovePixels, fkHighThreshold, fkDiffusionFixedLoc, fkDiffusionFromData}:
    result = generateFakeFromData(rnd, h5f, run, chipNumber, calibInfo, fakeDesc, runType, T, ctx)
  else:
    result = generateFakeFromGain(rnd, h5f, run, chipNumber, calibInfo, fakeDesc, runType, T, ctx, useCache)

proc generateFakeData*[C: MaybeContext](
  rnd: var Rand, h5f: H5File,
  fakeDesc: FakeDesc,
  run = -1,
  useCache = true,
  ctx: C = missing()
                                      ): DataFrame =
  ## For each run generate `nFake` fake events
  ##
  ## XXX:
  ## Insert the other kinds of fake data:
  ## - same energy but higher threshold (dropping pixels below a threshold)
  result = newDataFrame()
  let fileInfo = h5f.getFileInfo()
  let tpx = fileInfo.timepix
  let centerChip = fileInfo.centerChip
  let chipName = fileInfo.centerChipName
  let capacitance = tpx.getCapacitance()

  for (num, group) in runs(h5f):
    if run > 0 and run != num: continue # skip if not the desired run
    let dfRun = generateRunFakeData(rnd, h5f, num,
                                    centerChip, chipName,
                                    capacitance,
                                    fakeDesc,
                                    fileInfo.runType,
                                    DataFrame,
                                    ctx,
                                    useCache)
    if dfRun.len > 0:
      result.add dfRun

proc generateFakeData*(h5f: H5File, fakeDesc: FakeDesc,
                       run = -1,
                       seed = 42): DataFrame =
  ## Convenience version to generate fake data based on an input file.
  const CdlFile = "/home/basti/CastData/data/CDL_2019/calibration-cdl-2018.h5"
  let ctx = initLikelihoodContext(CdlFile,
                                  year = yr2018,
                                  energyDset = igEnergyFromCharge,
                                  region = crSilver,
                                  timepix = Timepix1,
                                  morphKind = mkLinear) # morphing to plot interpolation
  var rnd = initRand(seed)
  result = generateFakeData(rnd, h5f, fakeDesc, run, ctx = ctx)

proc generateRawFakeData*(h5f: H5File, fakeDesc: FakeDesc,
                          run = -1,
                          seed = 42,
                          useCache = true): Table[int, seq[Pixels]] =
  ## Convenience version to generate raw fake data based on an input file, i.e. not
  ## reconstructed but raw `Pix` data.
  ##
  ## Data is returned as a table with run number of the original gain / diffusion parameter
  ## based on and `Pixels` data for each.
  var rnd = initRand(seed)
  let fileInfo = h5f.getFileInfo()
  let tpx = fileInfo.timepix
  let centerChip = fileInfo.centerChip
  let chipName = fileInfo.centerChipName
  let capacitance = tpx.getCapacitance()

  result = initTable[int, seq[Pixels]]()
  for (num, group) in runs(h5f):
    if run > 0 and run != num: continue # skip if not the desired run
    result[num] = generateRunFakeData(rnd, h5f, num,
                                      centerChip, chipName,
                                      capacitance,
                                      fakeDesc,
                                      fileInfo.runType,
                                      seq[Pixels],
                                      useCache = useCache)

proc user(energy: seq[float], nmc: int,
          gains: seq[float], diffusions: seq[float],
          run = -1) = #tfKind = XXX
  ## Either generate data like an existing run
  ## or generate data according to gas gain, diffusion and energies
  # 1. generate fake events according to input
  discard

proc like(path: string, ## Input file
          run: int,
          outpath: string,
          outRun: int,
          tfKind: TargetFilterKind,
          nmc: int,
          chip = 3) =
  ## Generates fake data like `run` in `path` for chip `chip`
  let fakeDesc = FakeDesc(nFake: nmc,
                          tfKind: tfKind,
                          kind: fkGainDiffusion)
  let filter = H5Filter(kind: fkZlib, zlibLevel: 4)
  var rnd = initRand(42)
  withH5(path, "r"):
    let dfFake = generateFakeData(rnd, h5f, fakeDesc, run)
    #let dfFake = h5f.generateFakeData(fakeDesc, run = run)
    var h5fout = H5open(outpath, "rw")

    ggplot(dfFake, aes("energyFromCharge")) +
      geom_histogram(bins = 100) +
      ggsave("/tmp/energy_dist.pdf")
    for c in getKeys(dfFake):
      var dset = h5fout.create_dataset(
        recoBase() & $outRun / "chip_" & $chip / c,
        dfFake.len,
        dtype = float,
        chunksize = @[10000],
        maxshape = @[int.high],
        filter = filter)
      dset.unsafeWrite(cast[ptr float](dfFake[c, float].unsafe_raw_offset()), dfFake.len)

when isMainModule:
  import cligen
  dispatchMulti([like], [user])
