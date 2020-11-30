import mclimit
import ggplotnim, seqmath, sequtils, tables, os, options
import ingrid / [tos_helpers]
from arraymancer import tensor

import nlopt, nimhdf5, numericalnim
import shell, strutils

import cligen

type
  FitObject = object
    back: Histogram
    cand: Histogram
    axModel: DataFrame
    gae: ptr float
    rnd: ptr Random
    obsCLs: ptr float
    obsCLb: ptr float
    obsCLsb: ptr float
    opKind: OptimizeByKind

  OptimizeByKind = enum
    opCLs = "CLs"
    opCLb = "CLb" ## NOTE: optimizing by this probably does not make sense in most use cases!
    opCLsb = "CLsb"

func toTensor[T](t: Tensor[T]): Tensor[T] = t
template toHisto(arg, binsArg: typed): untyped =
  var counts = arg.toTensor.asType(float)
  # why cannot directly map_inline with sqrt :(
  var err = newSeq[float](counts.len)
  for i in 0 ..< counts.len:
    err[i] = sqrt(counts[i])
  Histogram(ndim: 1,
            bins: binsArg.toTensor.asType(float),
            counts: concat([counts, zeros[float](1)], axis = 0),
            err: err.toTensor)

proc readDsets(h5f: H5FileObj, names: varargs[string]): DataFrame =
  ## reads all likelihood data in the given `h5f` file as well as the
  ## corresponding energies. Flattened to a 1D seq.
  ## This proc is for TPA generated H5 files! (i.e. containing run_* groups, ...)
  # iterate over all groups, read all likelihood and energy dsets
  result = newDataFrame()
  for name in names:
    var data = newSeq[float]()
    for run, grp in runs(h5f, likelihoodBase()):
      let group = h5f[grp.grp_str]
      let centerChip = if "centerChip" in group.attrs: "chip_" & $group.attrs["centerChip", int]
                       else: "chip_3"
      if grp / centerChip in h5f:
        doAssert grp / centerChip / name in h5f[(group.name / centerChip).grp_str]
        data.add h5f[grp / centerChip / name, float64]
      else:
        echo &"INFO: Run {run} does not have any candidates"
    result[name] = toColumn data

proc rescale(flux: var Tensor[float], gae_new, gae_current: float) =
  echo "gae current ", gae_current
  flux.apply_inline(x * pow(gae_new / gae_current, 2.0))

proc drawLimitPlot(flux, energy: Tensor[float], param: float,
                   backHist, candHist: Histogram,
                   fIsSB: bool) =
  # TODO: get rid of hardcoded value here!
  const g_agamma = 1e-12
  let fluxPlot = if fIsSB: flux +. backHist.counts else: flux.clone
  let axLab = if fIsSB: "ax. sig+back" else: "axion signal"
  var df = seqsToDf({ axLab : fluxPlot,
                      "Energy" : energy,
                      "background" : backHist.counts,
                      "exp. cand." : candHist.counts })
  let suff = if fIsSB: "_sb" else: ""
  template cc(h, df, col, op): untyped =
    let x = df[col].toTensor(float)
    let err = h.err
    var res = zeros[float](df.len)
    for i in 0 ..< df.len:
      res[i] = op(x[i], err[i])
    res
  df.write_csv(&"/tmp/current_data{suff}.csv")

  var yMin = zeros[float](df.len * 3)
  yMin[0 ..< df.len] = cc(backHist, df, "background", `-`)
  yMin[df.len ..< 2 * df.len] = cc(candHist, df, "exp. cand.", `-`)
  var yMax = zeros[float](df.len * 3)
  yMax[0 ..< df.len] = cc(backHist, df, "background", `+`)
  yMax[df.len ..< 2 * df.len] = cc(candHist, df, "exp. cand.", `+`)
  df = df.gather(["background", "exp. cand.", axLab], "Type", "y")
  df["yMin"] = yMin.map_inline(
    if x < 0.0: 0.0
    else: x
  )
  df["yMax"] = yMax
  # echo df.pretty(-1)
  ggplot(df, aes(Energy, "y", fill = "Type", color = "Type")) +
    geom_histogram(stat = "identity", position = "identity", alpha = some(0.5)) +
    geom_point(binPosition = "center") +
    geom_errorbar(data = df.filter(f{`Type` == "background"}),
                  aes = aes(yMin = "yMin", yMax = "yMax"), binPosition = "center") +
    ylab("#") +
    xlab("Energy [keV]") +
    # ggtitle(&"Expected g_ae·g_aγ = {param * g_agamma:.2e}, CLs = {obsCLs:.4f}, CLb = {obsCLb:.4f}") +
    ggtitle(&"Expected g_ae·g_aγ = {param * g_agamma:.2e} GeV⁻¹ at 95% CLs") +
    ggsave(&"/tmp/current_flux{suff}.pdf", width = 800, height = 480)


proc runLimitCalc(p: float, data: FitObject) =
  const nmc = 100_000
  var (backHist, candHist, rnd) = (data.back, data.cand, data.rnd)

  let axModel = data.axModel
  # NOTE: since the flux tensor is a reference it's enough to change it. We don't have
  # to reassign it as a Column to the DF
  var (flux, energy) = (axModel["Flux"].toTensor(float), axModel["Energy"].toTensor(float))

  let param = if classify(p) != fcNan: p else: 1e-22
  echo "Scaling to ", param, " from ", data.gae[]
  # now rescale the flux
  flux.rescale(param, data.gae[])
  data.gae[] = param

  var
    obsCLs: float
    obsCLb: float
    obsCLsb: float
  when not defined(useRoot):
    let sigHist = toHisto(flux, energy)
    let ch = mclimit.Channel(sig: sigHist, back: backHist, cand: candHist,
                     systErr: { "Software" : SystematicError(cand: 0.05, back: 0.05),
                                "Stat" :  SystematicError(cand: 0.3, back: 0.1),
                                "Tel" : SystematicError(cand: 0.05, back: 0.05),
                                "Window" : SystematicError(cand: 0.10, back: 0.10)
                              }.toOrderedTable)
    var rand = wrap(initMersenneTwister(49))
    let limit = computeLimit(@[ch], rand, stat = false, nmc = nmc,
                             verbose = false)
    obsCLs = limit.CLs()
    obsCLb = limit.CLb()
    obsCLsb = limit.CLsb()
    echo "<CLb>  = ", limit.getExpectedCLb_b()
    echo "<CLsb> = ", limit.getExpectedCLsb_b()
    echo "<CLs>  = ", limit.getExpectedCLs_b()
  else:
    let res = shellVerbose:
      "../../../mclimit/tools/calcLimit /tmp/current_data.csv true"
    obsCLs = res[0].splitLines[^2].parseFloat
    obsCLb = res[0].splitLines[^1].parseFloat
  data.obsCLs[] = obsCLs
  data.obsCLb[] = obsCLb
  data.obsCLsb[] = obsCLsb
  block Plot:
    # plot current model
    drawLimitPlot(flux, energy, data.gae[], backHist, candHist, true)
    drawLimitPlot(flux, energy, data.gae[], backHist, candHist, false)

proc calcCL95(p: seq[float], data: FitObject): float =
  runLimitCalc(p[0], data)
  case data.opKind
  of opCLs:
    result = data.obsCLs[]
  of opCLb:
    result = data.obsCLb[]
  of opCLsb:
    result = data.obsCLsb[]

proc constrainCL95(p: seq[float], data: FitObject): float =
  ## TODO: instead of calling runLimitCalc again here, we should add
  ## the results to the FitObject and then just read them here
  var
    obsCLs = data.obsCLs[]
    obsCLb = data.obsCLb[]
    obsCLsb = data.obsCLsb[]
  case data.opKind
  of opCLs:
    result = abs(obsCLs - 0.05 - 1e-3) + 1e-3
  of opCLb:
    result = abs(obsCLb - 0.05 - 1e-3) + 1e-3
  of opCLsb:
    result = abs(obsCLsb - 0.05 - 1e-3) + 1e-3
  echo result, " at a param ", p
  echo "CLb    = ", obsCLb
  echo "CLs    = ", obsCLs
  echo "CLsb   = ", obsCLsb

proc readAxModel(f: string, scale: float, limit2013 = false): DataFrame =
  ## scale is the scaling required from a purely weight based flux
  ## to one corresponding to the tracking time. The flux by itself is
  ## essentially corresponding to a fixed, very short time, based on the
  ## number of axions being simulated. For a fixed coupling constant we get
  ## a fixed number of axions per year and cm²
  let
    gaeDf = toDf(readCsv(f))
    energy = gaeDf["Axion energy [keV]"].toTensor(float).toRawSeq
    weights = gaeDf["Flux after experiment"].toTensor(float).toRawSeq
  var
    val, bins: seq[float]
  if limit2013:
    (val, bins) = histogram(energy, weights = weights, range = (0.8, 6.2285), bins = 20)
  else:
    (val, bins) = histogram(energy, weights = weights, range = (0.0, 10.0), bins = 50)
    #(val, bins) = histogram(energy, range = (0.0, 10.2), bins = 51)
  val.add 0.0
  var flux = val.toTensor.asType(float)
  # scale to tracking time
  flux.apply_inline(x * scale)
  result = seqsToDf({ "Energy" : bins[0 .. ^1],
                      "Flux" : flux })
  ggplot(result, aes("Energy", "Flux")) +
    geom_histogram(stat = "identity") +
    ggsave("/tmp/axionModel.pdf")

proc drawExpCand(h: Histogram): Histogram =
  ## given a histogram as input, draws a new histogram using Poisson
  ## statistics
  var pois: Poisson
  var rnd = wrap(initMersenneTwister(0x1337))
  result = h.clone()
  for i in 0 ..< h.counts.len:
    let cnt = h.counts[i]
    pois = poisson(cnt)
    let cntDraw = rnd.sample(pois)
    result.counts[i] = cntDraw
    result.err[i] = sqrt(cntDraw)

proc main(backFiles, candFiles: seq[string], axionModel: string,
          optimizeBy: string = $opCLs,
          limit2013: bool = false) =
  var
    h5Backs = backFiles.mapIt(H5open(it, "r"))
    h5Cands = candFiles.mapIt(H5open(it, "r"))
  let searchStr = "flux_after_exp_N_"
  let idxN = find(axionModel, searchStr) + searchStr.len
  let N_sim = parseInt(axionModel[idxN ..< ^4])

  # FIX ME! has to be cleaned up
  let resPath = "../../../AxionElectronLimit"
  let diffFluxDf = toDf(readCsv(resPath / "axion_diff_flux_gae_1e-13_gagamma_1e-12.csv"))
  let totalFluxPerYear = simpson(diffFluxDf["Flux / keV⁻¹ m⁻² yr⁻¹"].toTensor(float).toRawSeq,
                                 diffFluxDf["Energy / eV"].toTensor(float).map_inline(x * 1e-3).toRawSeq)

  #let gaeRawDf = toDf(readCsv(axionModel)).rename(f{"Energy" <- "Axion energy [keV]"})
  proc flatten(dfs: seq[DataFrame]): DataFrame =
    ## flatten a seq of DFs, which are identical by stacking them
    for df in dfs:
      result.add df.clone
  let backEnergy = h5Backs.mapIt(
    it.readDsets("energyFromCharge")
      .rename(f{"Energy" <- "energyFromCharge"})).flatten
  let candEnergy = h5Cands.mapIt(
    it.readDsets("energyFromCharge")
      .rename(f{"Energy" <- "energyFromCharge"})).flatten

  let df2013 = toDf(readCsv("../../resources/background_rate_cast_gae_2013.csv"))

  let ratePlot = bind_rows([("back", backEnergy), ("cand", candEnergy)], "Type")
  ggplot(ratePlot, aes(Energy, fill = "Type")) +
    geom_histogram(bin_width = 0.2, position = "identity", alpha = some(0.5)) +
    xlim(0, 11) +
    ggsave("/tmp/back_histo.pdf")

  var rnd = wrap(initMersenneTwister(49))
  let back = backEnergy["Energy"].toTensor(float).toRawSeq
  let cand = candEnergy["Energy"].toTensor(float).toRawSeq
  let (backHI, binsB) = histogram(back, range = (0.0, 10.0), bins = 50)
  let (candHI, binsC) = histogram(cand, range = (0.0, 10.0), bins = 50)

  proc readDuration(h5f: H5FileObj): float =
    let lhGrp = h5f["/likelihood".grp_str]
    result = lhGrp.attrs["totalDuration", float]

  var
    backHist: Histogram
    candHist: Histogram
    backTime: float
    trackToBackRatio: float
  if limit2013:
    backHist = toHisto(df2013["Background"].toTensor(float).toRawSeq,
                           df2013["Energy"].toTensor(float).toRawSeq)
    backHist.err = backHist.err.map_inline(x.float / 4.0)

    backTime = 1890 * 3600 # background time of 2013 paper
    trackToBackRatio = backTime / (197 * 3600)
    candHist = toHisto(df2013["Candidates"].toTensor(float).toRawSeq,
                       df2013["Energy"].toTensor(float).toRawSeq)
  else:
    # var backH = backHI
    # var candH = candHI
    backHist = toHisto(backHI, binsB)
    trackToBackRatio = 19.56 # this is the ratio of background to
                                   # tracking time in the combined Run 2 and 3
                                   # only required for poisson sampled candidates
    backHist.counts = backHist.counts.map_inline(x.float / trackToBackRatio)
    backHist.err = backHist.err.map_inline(x.float / trackToBackRatio)

    backTime = h5Backs.mapIt(it.readDuration).sum
    #let candHist = toHisto(candH, binsC)
    candHist = backHist.drawExpCand()

  # compute signal
  let trackingTime = backTime / trackToBackRatio
  echo "Total background time ", backTime / 3600.0, " h"
  echo "Total tracking time ", trackingTime / 3600.0, " h"
  let secondsOfSim = N_sim.float / totalFluxPerYear * 86400 * 365
  echo &"secondsOfSim = {secondsOfSim}"
  let areaBore = PI * pow(2.15, 2.0) # area of bore in cm²
  echo &"areaBore = {areaBore} cm²"
  # - calculate how much more time is in tracking than simulation
  # - convert from m² to cm²
  # - multiply by area of bore
  #let scale = totalFluxPerYear / N_sim.float * 5.0 / (100 * 100) * areaBore * (trackingTime / (86400 * 365))
  let scale = trackingTime.float / secondsOfSim / (100 * 100) * areaBore

  echo &"Scale = {scale}"
  let gaeDf = readAxModel(axionModel, scale, limit2013)

  # 1e-13 is the value, which was used to calculate the currently used flux
  var gae = 1e-13
  var obsCLs: float
  var obsCLb: float
  var obsCLsb: float
  let fitObj = FitObject(back: backHist, cand: candHist,
                         axModel: gaeDf,
                         gae: gae.addr,
                         rnd: rnd.addr,
                         obsCLs: obsCLs.addr,
                         obsCLb: obsCLb.addr,
                         obsCLsb: obsCLsb.addr,
                         opKind: parseEnum[OptimizeByKind](optimizeBy))
  var opt = newNloptOpt[FitObject](LN_COBYLA, 1, @[(l: 1e-14, u: 1e-8)])
  let varStruct = newVarStruct(calcCL95, fitObj)
  opt.setFunction(varStruct)
  var constrainVarStruct = newVarStruct(constrainCL95, fitObj)
  opt.addInequalityConstraint(constrainVarStruct)
  opt.maxtime = 30.0
  opt.initialStep = 1e-10
  let optRes = opt.optimize(@[2.1e-10])
  echo opt.status

  echo optRes
  destroy(opt)

  ## TODO: closing produces attribute closing errors, investigate!
  #for h5f in concat(h5Backs, h5Cands):
  #  discard h5f.close()

  # finally run root as comp:
  let res = shellVerbose:
    "../../../mclimit/tools/calcLimit /tmp/current_data.csv"


when isMainModule:
  dispatch(main, echoResult = false, noAutoEcho = true)
