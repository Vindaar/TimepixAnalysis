import mclimit
import ggplotnim, seqmath, sequtils, tables, os, options
import ingrid / [tos_helpers]
from arraymancer import tensor

import nlopt, nimhdf5
import shell, strutils

import cligen

type
  FitObject = object
    back: Histogram
    cand: Histogram
    axModel: DataFrame
    gae: ptr float
    rnd: ptr Random

func toTensor[T](t: Tensor[T]): Tensor[T] = t
template toHisto(arg, binsArg: typed): untyped =
  let counts = arg.toTensor.asType(float)
  # why cannot directly map_inline with sqrt :(
  var err = newSeq[float](counts.len)
  for i in 0 ..< counts.len:
    err[i] = sqrt(counts[i])
  Histogram(ndim: 1,
            bins: binsArg.toTensor.asType(float),
            counts: counts,
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
  let fluxPlot = if fIsSB: flux +. backHist.counts else: flux
  let axLab = if fIsSB: "ax. sig+back" else: "axion signal"
  var df = seqsToDf({ axLab : fluxPlot,
                      "Energy" : energy,
                      "background" : backHist.counts,
                      "exp. cand." : candHist.counts })
  template cc(h, df, col, op): untyped =
    let x = df[col].toTensor(float)
    let err = h.err
    var res = zeros[float](df.len)
    for i in 0 ..< df.len:
      res[i] = op(x[i], err[i])
    res
  df.write_csv("/tmp/current_data.csv")

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
  let suff = if fIsSB: "_sb" else: ""
  ggplot(df, aes(Energy, "y", fill = "Type", color = "Type")) +
    geom_histogram(stat = "identity", position = "identity", alpha = some(0.5)) +
    geom_point(binPosition = "center") +
    geom_errorbar(data = df.filter(f{`Type` == "background"}),
                  aes = aes(yMin = "yMin", yMax = "yMax"), binPosition = "center") +
    ylab("#") +
    xlab("Energy [keV]") +
    # ggtitle(&"Expected g_ae·g_aγ = {param * g_agamma:.2e}, CLs = {obsCLs:.4f}, CLb = {obsCLb:.4f}") +
    ggtitle(&"Expected g_ae·g_aγ = {param * g_agamma:.2e} GeV⁻¹ at 95% CLs") +
    ggsave(&"/tmp/current_flux{suff}.pdf")


proc runLimitCalc(p: float, data: FitObject) =
  const nmc = 100_000
  var (backHist, candHist, rnd) = (data.back, data.cand, data.rnd)

  let axModel = data.axModel
  var (flux, energy) = (axModel["Flux"].toTensor(float), axModel["Energy"].toTensor(float))

  let param = if classify(p) != fcNan: p else: 1e-22
  echo "Scaling to ", param
  # now rescale the flux
  flux.rescale(param, data.gae[])
  data.gae[] = param

  #var obsCLs: float
  #var obsCLb: float
  when not defined(useRoot):
    let sigHist = toHisto(flux, energy)
    let ch = mclimit.Channel(sig: sigHist, back: backHist, cand: candHist,
                     systErr: { "Tel" : SystematicError(cand: 0.1, back: 0.0),
                                "Rad" : SystematicError(cand: 0.2, back: 0.0) }.toOrderedTable)
    var rand = wrap(initMersenneTwister(49))
    let limit = computeLimit(@[ch], rand, stat = false, nmc = nmc,
                             verbose = false)
    obsCLs = limit.CLs()
    obsCLb = limit.CLb()
  else:
    let res = shellVerbose:
      "../../../mclimit/tools/calcLimit /tmp/current_data.csv true"
    obsCLs = res[0].splitLines[^2].parseFloat
    obsCLb = res[0].splitLines[^1].parseFloat

  block Plot:
    # plot current model
    drawLimitPlot(flux, energy, param, backHist, candHist, true)
    drawLimitPlot(flux, energy, param, backHist, candHist, false)

proc calcCL95(p: seq[float], data: FitObject): float =
  var
    obsCLs: float
    obsCLb: float
  runLimitCalc(p[0], data, obsCLs, obsCLb)
  result = obsCLs

proc constrainCL95(p: seq[float], data: FitObject): float =
  var
    obsCLs: float
    obsCLb: float
  runLimitCalc(p[0], data, obsCLs, obsCLb)
  result = abs(obsCLs - 0.5 - 1e-3) - 1e-3
  echo result, " at a CLs: ", obsCLs, "  and CLb  ", obsCLb, " for param ", p

proc readAxModel(f: string): DataFrame =
  let
    gaeDf = toDf(readCsv(f))
    energy = gaeDf["Axion energy [keV]"].toTensor(float).toRawSeq
    weights = gaeDf["Flux after experiment"].toTensor(float).toRawSeq
    (val, bins) = histogram(energy, weights = weights, range = (0.0, 10.2), bins = 51)
  var flux = val.toTensor
  flux.rescale(1e6, 1.0)
  result = seqsToDf({ "Energy" : bins[0 .. ^2],
                      "Flux" : flux })

proc main(backFile, candFile, axionModel: string) =

  echo "file: ", backFile, " exists? ", existsFile(backFile)
  let h5Back = H5File(backFile, "r")
  let h5Cand = H5File(candFile, "r")
  let gaeDf = readAxModel(axionModel)
  echo gaeDf

  #let gaeRawDf = toDf(readCsv(axionModel)).rename(f{"Energy" <- "Axion energy [keV]"})
  let backEnergy = h5Back.readDsets("energyFromCharge").rename(f{"Energy" <- "energyFromCharge"})
  let candEnergy = h5Cand.readDsets("energyFromCharge").rename(f{"Energy" <- "energyFromCharge"})
  let ratePlot = bind_rows([("back", backEnergy), ("cand", candEnergy)], "Type")
  ggplot(ratePlot, aes(Energy, fill = "Type")) +
    geom_histogram(bin_width = 0.2, position = "identity", alpha = some(0.5)) +
    xlim(0, 11) +
    ggsave("/tmp/back_histo.pdf")

  var rnd = wrap(initMersenneTwister(49))
  let back = backEnergy["Energy"].toTensor(float).toRawSeq
  let cand = candEnergy["Energy"].toTensor(float).toRawSeq
  let (backHI, binsB) = histogram(back, range = (0.0, 10.2), bins = 51)
  let (candHI, binsC) = histogram(cand, range = (0.0, 10.2), bins = 51)

  proc scaleDset(h5f: H5FileObj, data: seq[int]): seq[float] =
    let lhGrp = h5f["/likelihood".grp_str]
    let time_back = lhGrp.attrs["totalDuration", float]
    let area = pow(0.95 - 0.45, 2)
    const bin_width = 0.392
    const shutter_open = 1.0
    const factor = 1.0
    let scale = factor / (time_back * shutter_open * area * bin_width) * 86400 #* 1e5
    result = data.mapIt(it.float * scale)

  # NOTE: scaling ``before`` running the calculation does not make sense!
  # We need the counts for the poisson statistics / smearing.
  let backH = backHI#scaleDset(h5Back, backHI)
  let candH = candHI#scaleDset(h5Cand, candHI)

  let backHist = toHisto(backH, binsB[0 .. ^2])
  let candHist = toHisto(candH, binsC[0 .. ^2])

  # 1e-13 is the value, which was used to calculate the currently used flux
  var gae = 1e-13
  let fitObj = FitObject(back: backHist, cand: candHist,
                         axModel: gaeDf,
                         gae: gae.addr,
                         rnd: rnd.addr)
  var opt = newNloptOpt[FitObject](LN_COBYLA, 1, @[(l: 1e-12, u: 1e-8)])
  let varStruct = newVarStruct(calcCL95, fitObj)
  opt.setFunction(varStruct)
  var constrainVarStruct = newVarStruct(constrainCL95, fitObj)
  opt.addInequalityConstraint(constrainVarStruct)
  #opt.xtol_rel = 1e-10
  #opt.ftol_rel = 1e-10
  opt.maxtime = 600.0
  opt.initialStep *= 1e-10
  let optRes = opt.optimize(@[2.1e-9])
  echo opt.status

  echo optRes
  destroy(opt)

when isMainModule:
  dispatch(main, echoResult = false, noAutoEcho = true)
