import std / [random, strformat, os]
import unchained, ggplotnim
import nimhdf5 except copyFlat, fromFlat, calcSize, newBuffer
from std / sequtils import toSeq, mapIt

import ingrid / [fake_event_generator, gas_physics, ingrid_types, tos_helpers]
from pkg / xrayAttenuation import FluorescenceLine
import helpers / sampling_helper

import ./io_helpers

proc generateSingleEvent(rnd: var Rand,
                         fakeDesc: var FakeDesc,
                         calibInfo: CalibInfo,
                         energySampler: sampling_helper.Sampler,
                         gains: seq[float],
                         diffusion: seq[float]): FakeEvent =
  # 1. sample an energy from the given sampler
  let energy = rnd.sample(energySampler).keV
  let lineName = &"{energy.float:.2f}" ## <-- This should yield 1000 different absorption length samplers (if going to 10 keV)!
  let lines = @[FluorescenceLine(name: lineName, energy: energy, intensity: 1.0)]
  # 2. sample a gas gain
  var G = Inf
  while G > 4500:
    G = rnd.gauss(mu = (gains[1] + gains[0]) / 2.0, sigma = (gains[1] - gains[0]) / 3.0)
  #let gain = GainInfo(N: 100_000.0, G: G, theta: rnd.rand(0.4 .. 2.4))
  ## Maybe better matching gain?
  let gain = GainInfo(N: 900_000.0, G: G, theta: rnd.rand(1.4 .. 2.4))
  # 3. sample a diffusion
  # let σT = rnd.gauss(mu = 660.0, sigma = (diffusion[1] - diffusion[0]) / 4.0) ## <- had a parenthesis bug!
  let σT = rnd.rand(diffusion[0] .. diffusion[1])
  # update σT field of FakeDesc
  fakeDesc.σT = σT
  # now generate and add
  while not result.valid:
    result = rnd.generateAndReconstruct(fakeDesc, lines, gain, calibInfo, energy)

type
  SimContext = object
    rnd: Rand
    nFake: int
    fakeDesc: FakeDesc
    calibInfo: CalibInfo
    energySampler: Sampler
    gains: seq[float]
    diffusion: seq[float]
    outfile: string ## path where we write `events` as strings `/dev/shm/...`
    events: seq[FakeEvent]

import flatBuffers

#var dataChan: Channel[string]
#dataChan.open()
#proc generateParallelEvent(ctx: ptr SimContext) {.gcsafe.} =
#  {.cast(gcsafe).}:
#    ctx.events = newSeq[FakeEvent](ctx.nFake)
#    var rnd = ctx.rnd
#    echo "Generating: ", ctx.nFake, " events"
#    for i in 0 ..< ctx.nFake:
#      if i mod 5000 == 0:
#        echo "Generated ", i, " / ", ctx.nFake, " events."
#      ctx.events[i] = rnd.generateSingleEvent(ctx.fakeDesc, ctx.calibInfo,
#                                              ctx.energySampler, ctx.gains,
#                                              ctx.diffusion)
#    echo "Sending data"
#    dataChan.send(asFlat(ctx.events).toString())

proc initSimContexts(num: int, nFake: int,
                     calibInfo: CalibInfo,
                     energySampler: sampling_helper.Sampler,
                     gains: seq[float],
                     diffusion: seq[float]): seq[SimContext] =
  const baseSeed = 0x12312
  let nFakeEach = nFake div num
  for i in 0 ..< num:
    let fakeDesc = FakeDesc(kind: fkGainDiffusion, gasMixture: initCASTGasMixture())
    let rnd = initRand(baseSeed + num)
    let each = if i != num - 1: nFakeEach else: nFakeEach + nFake mod num
    result.add SimContext(rnd: rnd,
                          nFake: each,
                          fakeDesc: fakeDesc,
                          calibInfo: calibInfo,
                          energySampler: energySampler,
                          outfile: "/dev/shm/data_buffer_" & $i & ".dat",
                          gains: gains, diffusion: diffusion)

#import pkg / malebolgia
import cligen / [procpool, mslice, osUt]

from std / strutils import parseBool, parseInt

proc generateFakeEvents(rnd: var Rand,
                        calibInfo: CalibInfo,
                        energySampler: sampling_helper.Sampler,
                        gains: seq[float],
                        diffusion: seq[float],
                        nFake = 100_000): seq[FakeEvent] =
  let parallel = getEnv("MT", "false").parseBool
  if not parallel:
  #when true:
    var count = 0
    result = newSeqOfCap[FakeEvent](nFake)
    # Note: tfKind and nFake irrelevant for us!
    var fakeDesc = FakeDesc(kind: fkGainDiffusion, gasMixture: initCASTGasMixture())
    while count < nFake:
      if count mod 5000 == 0:
        echo "Generated ", count, " events."
      result.add rnd.generateSingleEvent(fakeDesc, calibInfo, energySampler, gains, diffusion)
      inc count
  #else:
  #when false:
  else:
    let nJobs = getEnv("JOBS", "28").parseInt
    var ctxs = initSimContexts(nJobs, nFake, calibInfo, energySampler,
                               gains, diffusion)
    var pp = initProcPool((
      proc(r, w: cint) =
        const baseSeed = 0xAFFE1337
        let i = open(r)
        var val: int
        while i.uRd(val):
          echo "Generating fake events with seed: ", val
          var ctx = ctxs[val]
          ctx.events = newSeq[FakeEvent](ctx.nFake)
          for j in 0 ..< ctx.nFake:
            if j mod 5000 == 0:
              echo "Generated ", j, " / ", ctx.nFake, " events."
            ctx.events[j] = ctx.rnd.generateSingleEvent(ctx.fakeDesc, ctx.calibInfo,
                                                        ctx.energySampler, ctx.gains,
                                                        ctx.diffusion)
          let buf = asFlat(ctx.events)
          buf.writeBuffer(ctx.outfile)
          ## NOTE: writing the buffer directly breaks for reasons I don't yet understand
          discard w.wrLenBuf("wrote buffer to " & $ctx.outfile & " of length: " & $buf.size)
    ),
                        framesLenPfx,
                        nJobs)
    let nums = toSeq(0 ..< nJobs)
    var readRes = proc(s: MSlice) = echo $s
    pp.evalOb nums, readRes
    for ctx in ctxs:
      let dat = readFile(ctx.outfile)
      echo "Read data from ", ctx.outfile, " of len ", dat.len
      let evs = flatTo[seq[FakeEvent]](fromString(dat))
      result.add evs

  when false:
    let nJobs = ThreadPoolSize
    var ctxs = initSimContexts(nJobs, nFake, calibInfo, energySampler,
                               gains, diffusion)
    var m = createMaster()
    m.awaitAll:
      for i in 0 ..< nJobs:
        m.spawn generateParallelEvent(ctxs[i].addr)
    #for ctx in ctxs:
    #  result.add ctx.events
    var (av, dat) = tryRecv(dataChan)
    while av:
      result.add(flatTo[seq[FakeEvent]](fromString(dat)))
      (av, dat) = tryRecv(dataChan)


proc initEnergySampler(energyMin, energyMax: float,
                       yEnergyMin, yEnergyMax: float): sampling_helper.Sampler =
  let fnSample = (
    proc(x: float): float =
      let m = (yEnergyMin - yEnergyMax) / (energyMin - energyMax)
      let b = yEnergyMax - m * energyMax
      result = m * x + b
  )
  result = sampler(fnSample, energyMin, energyMax, num = 1000)

proc generateFakeEventsDf*(rnd: var Rand,
                           calibInfo: CalibInfo,
                           gains: seq[float],
                           diffusion: seq[float],
                           nFake = 100_000): DataFrame =
  let eS = initEnergySampler(0.1, 10.0, 1.0, 1.0) # uniform sampler
  result = fakeToDf( rnd.generateFakeEvents(calibInfo, eS, gains, diffusion, nFake) )
  result["Type"] = $dtSignal

import ingridDatabase / databaseRead
proc simulateEvents(calibFiles: seq[string], readRaw: bool, rnd: var Rand,
                    nFake: int,
                    energyMin, energyMax: float,
                    yEnergyMin, yEnergyMax: float,
                    outfile: string,
                    note: string
             ): seq[FakeEvent] =
  ## Generates synthetic data to train. The calibration input files are only used to determine the
  ## boundaries of the gas gain and diffusion in which to generate events in and the calibration
  ## function for energy calibration.
  ##
  ## The output HDF5 file `outfile` will contain `calibFiles.len * nFake` events.
  # for now we just define hardcoded bounds
  let gains = @[2400.0, 4000.0]
  let diffusion = @[550.0, 700.0] # μm/√cm, will be converted to `mm/√cm` when converted to DF
  # the theta parameter describing the Pólya also needs to be varied somehow

  # define the energy sampler as a linear function between (eMin,
  let eS = initEnergySampler(energyMin, energyMax, yEnergyMin, yEnergyMax)
  for c in calibFiles:
    withH5(c, "r"):
      let calibInfo = h5f.initCalibInfo()
      result.add rnd.generateFakeEvents(calibInfo, eS, gains, diffusion, nFake)
  echo "result: ", result.len

proc simulate(calib: seq[string],
              nFake: int,
              energyMin, energyMax: float,
              yEnergyMin, yEnergyMax: float,
              outfile: string,
              note: string,
              plotEnergy = true,
              chip = 3,
              run = 0
             ) =
  ## Simulates data to train a network on after.
  ## Energies for X-rays will be within energyMin to energyMax
  ## with linear change from yEnergyMin to yEnergyMax.
  ##
  ## The data will be written to `outfile` and `note` is an additional
  ## note added to the root group.
  var rnd = initRand(1349)
  let data = simulateEvents(calib, false, rnd, nFake, energyMin, energyMax,
                            yEnergyMin, yEnergyMax, outfile, note)
  # now write to output file
  var dfFake = fakeToDf(data)
  dfFake["eventNumber"] = toSeq(0 ..< dfFake.len)
  var h5f = H5open(outfile, "rw")
  let filter = H5Filter(kind: fkZlib, zlibLevel: 4)
  if plotEnergy:
    ggplot(dfFake, aes("energyFromCharge")) +
      geom_histogram(bins = 100) +
      ggsave("/tmp/energy_dist.pdf")
  for c in getKeys(dfFake):
    var dset = h5f.create_dataset(
      recoBase() & $run / "chip_" & $chip / c,
      dfFake.len,
      dtype = float,
      chunksize = @[10000],
      maxshape = @[int.high],
      filter = filter)
    dset.unsafeWrite(cast[ptr float](dfFake[c, float].unsafe_raw_offset()), dfFake.len)

  template fields(field, fieldName: untyped): untyped =
    var dset = h5f.create_dataset(
      recoBase() & $run / "chip_" & $chip / fieldName,
      data.len,
      dtype = float,
      chunksize = @[10000],
      maxshape = @[int.high],
      filter = filter)
    dset[dset.all] = data.mapIt(it.cluster.field)

  fields(centerX, "centerX")
  fields(centerY, "centerY")

  # And now write the raw x, y, ToT, charge data
  # First turn into 3 1D seqs
  var
    xs = newSeq[seq[uint8]](data.len)
    ys = newSeq[seq[uint8]](data.len)
    ch = newSeq[seq[uint16]](data.len)
  for i in 0 ..< data.len:
    let cl = data[i].cluster
    xs[i] = newSeq[uint8](cl.data.len)
    ys[i] = newSeq[uint8](cl.data.len)
    ch[i] = newSeq[uint16](cl.data.len)
    for j in 0 ..< cl.data.len:
      xs[i][j] = cl.data[j].x
      ys[i][j] = cl.data[j].y
      ch[i][j] = cl.data[j].ch
  # now writep
  template createWrite(name: string, data, spType: untyped): untyped =
    var dset = h5f.create_dataset(
      recoBase() & $run / "chip_" & $chip / name,
      data.len,
      dtype = spType,
      chunksize = @[10000],
      maxshape = @[int.high],
      filter = filter)
    dset[dset.all] = data
  createWrite("x", xs, special_type(uint8))
  createWrite("y", ys, special_type(uint8))
  createWrite("ToT", ch, special_type(uint16))

  var root = h5f["/".grp_str]
  root.attrs["Note"] = note
  ## Write some default attributes we might read via `getFileInfo`
  proc writeAttrs(h5f: H5File, name: string) =
    let grp = h5f[name.grp_str]
    grp.attrs["centerChip"] = 3
    grp.attrs["centerChipName"] = "H 10 W69"
    grp.attrs["numChips"] = 7 # we pretend
    grp.attrs["runType"] = $rtCalibration
    grp.attrs["TimepixVersion"] = "Timepix1"
    grp.attrs["runNumber"] = run
  h5f.writeAttrs("/reconstruction")
  h5f.writeAttrs("/reconstruction/run_" & $run)

when isMainModule:
  import cligen
  dispatch simulate
