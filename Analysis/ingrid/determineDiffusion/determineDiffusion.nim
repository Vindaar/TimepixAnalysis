import std / [os, sequtils, tables, stats, strformat, strutils, options]
import ingrid / [tos_helpers, ingrid_types]
import pkg / [nimhdf5, datamancer, ggplotnim, seqmath]
import unchained


## IMPORTANT NOTICE:
## The code here (and in related files touching on the diffusion constant)
## uses the variable `σT` to talk about the diffusion constant `D_T`:
## `D_T = σ_T / √( drift distance )`
## or written as the definition of `σ_T`:
## `σ_T = D_T · √( drift distance )`
##
## This is of course *INACCURATE*, but done for "historic" reasons. I started with
## an implementation that directly computed σ_T, but then switched to extracting
## the D_T value from it. We could update the code, but it would break the existing
## data files that contain D_T data, but under the term "σ_T".
##
## The reason for this is that as part of the fake data simulation here in
## `simulateRmsTrans`, we *first* sample a distance a cluster will diffuse
## and then generate from a normal distribution with `σ = D_T * √x`:
##
##   let x = rnd.gauss(mu = 0.0, sigma = σT * sqrt(zDrift))
##   let y = rnd.gauss(mu = 0.0, sigma = σT * sqrt(zDrift))
##
## while this should clearly say "D_T". Just be aware of it!

## XXX: a bit annoying that this is here...
const CdlFile = "/home/basti/CastData/data/CDL_2019/calibration-cdl-2018.h5"
const RmsCleaningCut = 1.5


proc linear(p: Tensor[float], x: float): float =
  result = p[0] * x + p[1]

import pkg / numericalnim except linspace

proc resultStr[T: Tensor | seq](pRes, pErr: T, χ²dof: float): string =
  result.add &"χ²/dof = {χ²dof:.4f}\n"
  for i in 0 ..< pRes.len.int:
    result.add &"p[{i}] = {pRes[i]:.4e} ± {pErr[i]:.4e}\n"
  result.strip()

proc echoResult(pRes, pErr: Tensor[float], χ²dof: float) =
  echo resultStr(pRes, pErr, χ²dof)

proc high[T](x: Tensor[T]): int = x.size.int - 1

const CutVal = "CutVal" # the actual cut value required to get `ε`!
const Err = 0.025

proc gauss[T: Tensor[float] | seq[float]](p: T, x: float): float =
  let N = p[0]
  let mean = p[1]
  let sigma = p[2]
  let
    arg = (x - mean).float / sigma.float
    res = exp(-0.5 * arg * arg)
  result = N / sqrt(2*PI) * res

proc gaussSeq(p: seq[float], x: float): float =
  result = gauss(p.toTensor, x)

proc toFitDf[T: Tensor | seq](params: T, xmin, xmax: float): DataFrame =
  let xFit = seqmath.linspace(xmin, xmax, 1000)
  let yFit = xFit.mapIt(gauss(params, it))
  result = toDf({ "xFit" : xFit,
                  "yFit" : yFit })

import mpfit
proc fitRmsTransverse(data: Tensor[float], bins: int): (DataFrame, seq[float], string) =
  var (hist, bins) = histogram(data.toSeq1D, bins = bins)
  let binWidth = bins[1] - bins[0]
  bins = bins[0 ..< ^1].mapIt(it + binWidth / 2.0) # remove last bin edge, move half bin right to have centers
  let hArgMax = hist.argmax
  let binArgMax = bins[hist.argmax]
  let histMax = hist[hArgMax].float

  let binLowerIdx = max(0, (hArgMax - (0.03 * bins.len.float)).round.int)
  let binUpperIdx = bins.high #(hArgMax + (0.055 * bins.len.float)).round.int
  let binUpper = bins[^1] # bins.lowerBound(data.percentile(98))]

  let params = @[histMax, binArgMax, abs(binArgMax - binUpper)]
  let xT = bins[binLowerIdx .. binUpperIdx]
  let yT = hist[binLowerIdx .. binUpperIdx].mapIt(it.float)
  let eT = yT.mapIt(1.0) #sqrt(x))

  let (pRes, res) = fit(gaussSeq, params, xT, yT, eT)
  let xFit = linspace(binArgMax, binUpper, 200)
  let yFit = xFit.mapIt(gauss(pRes, it))
  let resStr = resultStr(pRes, res.error, res.chiSq)

  echo resStr
  result = (toFitDf(pRes, xT.min, xT.max), pRes, resStr)

proc determineFromData(x: seq[float]): float =
  var (hist, bins) = histogram(x, bins = 100)
  let binWidth = bins[1] - bins[0]
  bins = bins[0 ..< ^1].mapIt(it + binWidth / 2.0) # remove last bin edge, move half bin right to have centers
  let hArgMax = hist.argmax
  let binArgMax = bins[hist.argmax]
  let histMax = hist[hArgMax].float
  var idx = hArgMax
  while hist[idx] > histMax * 0.1:
    inc idx
  result = bins[idx]

import ingrid / gas_physics
import helpers / [sampling_helper, stats_helpers]
import xrayAttenuation
import std / random

import os
import arraymancer
var HistoPath = "/home/basti/Sync/"
let HistoFile = "test_histo.pdf"
proc simulateRmsTrans(rnd: var Rand, rmsReal: seq[float], gasMixture: GasMixture, σT: float, energy: keV,
                      isBackground: bool, run: int): float =
  ## 2. get absorption length for this energy
  #rnd = initRand(42) # same RND for every run!
  let λ = if not isBackground: absorptionLengthCAST(gasMixture, energy)
          else: 0.cm
  ## 3. generate a sampler for the conversion point using λ
  proc expFn(x: float, λ: float): float =
    if not isBackground:
      result = 1.0 / λ * exp(- x / λ)
    else:
      result = 1.0 # uniform for backgound. Mostly muons parallel to readout. At any distance from cathode
  let fnSample = (proc(x: float): float =
    result = expFn(x, λ.float)
  )
  let distSampler = sampler(fnSample, 0.0, 3.0, num = 1000)

  var energySampler = if isBackground:
                        sampler(
                          (proc(x: float): float =
                            let λ = -10.0 / (ln(0.4))
                            result = exp(-x / λ)),
                          0.1, 10.0, num = 1000
                        )
                      else:
                        Sampler() # not needed

  ## 4. now sample N events with each `num` electrons, evaluate `rmsTransverse` for each
  const N = 5000
  var rms = newSeq[float](N)
  for i in 0 ..< N:
    ## 5. define number of electrons expected for this energy
    ## (or also receive `hits` dataset? Problematic due to background requirement!
    ##  instead use heuristic?)
    let num = if isBackground: rnd.sample(energySampler) / 0.026
              else: energy / 26.eV
    ## 6. sample a conversion position for this event, generate sampler for this events
    ##  possible diffusion
    let P =
      block:
        let conversionPoint = rnd.sample(distSampler) # point 'behind' cathode at which conversion takes place
        (proc(rnd: var Rand): float =
             let zDrift = 3.0 - conversionPoint # use conversion point to define a closure to sample for each electron
             let x = rnd.gauss(mu = 0.0, sigma = σT * sqrt(zDrift))
             let y = rnd.gauss(mu = 0.0, sigma = σT * sqrt(zDrift))
             result = sqrt(x*x + y*y)
        )
    let energySigma = if isBackground: 0.4 # larger range of energies allowed for background (but already uniform)
                      else: 0.1
    let numThis = max(rnd.gauss(mu = num, sigma = num * energySigma).floor.int, 15)
    # 2D Tensor of the pixel positions to later compute covariance matrix
    var posT = zeros[float]([numThis, 2])
    for j in 0 ..< numThis:
      let dist = P(rnd).μm.to(mm).float
      let φ = rnd.rand(0.0 .. 2*PI)
      let xp = (dist * cos(φ) + 7.0)
      let yp = (dist * sin(φ) + 7.0)
      posT[j, _] = [xp, yp].toTensor.unsqueeze(0)
    ## Compute rotated coordinate system using eigenvectors of covariance matrix
    let cov = covariance_matrix(posT, posT)
    let (eigenval, eigenvec) = symeig(cov, return_eigenvectors = true)
    let R = eigenvec
    # rotate the existing positions into the rotated coordinates using the eigenvectors
    let posTrot = posT * R
    # 0 corresponds to the shorter axis
    let rmsT = posTrot[_, 0].std()
    rms[i] = rmsT
  ## 6. now compute the loss
  result = cramerVonMises(rmsReal, rms) #kolmogorovSmirnov(rmsReal, rms)
  let df = bind_rows([("real", toDf({"rmsT" : rmsReal})), ("sim", toDf({"rmsT" : rms}))], "typ")
  ggplot(df, aes("rmsT", fill = "typ")) +
    geom_histogram(bins = 100, hdKind = hdOutline, alpha = 0.5, position = "identity", density = true) +
    ggtitle(&"Run: {run} with Cramér-von Mises = {result}") +
    ggsave(HistoPath / HistoFile)

proc determineMonteCarlo(rmsT: seq[float], energy: keV, isBackground = false, run = -1): (float, float) =
  ## It determines the best diffusion constant, `D_T`. That is:
  ##
  ## `σ_T = D_T · √( drift distance )`
  ##
  ## and not the diffusion coefficient `σ_T`!
  ##
  ## Determines the best diffusion matching the data using a MonteCarlo approach
  ## combined with a simple non linear optimization strategy that checks the
  ## improvement via Kolmogorov-Smirnov. We aim to reproduce the best possible
  ## match of the rmsTransverse data using the known data and the target energy.
  ## 1. define gas mixture
  let gm = initCASTGasMixture()
  ## 2. start diffusion value
  var diffusion = 660.0
  ## 3. define dual number
  let rmsReal = rmsT
  var diffD = diffusion
  var rnd = initRand(1212)
  #var ks = rnd.simulateRmsTrans(rmsReal, gm, diffD, energy)

  var ks = rnd.simulateRmsTrans(rmsReal, gm, diffD, energy, isBackground, run)

  # define a function that computes the derivative of ks with respect to diffD
  # this can be done either analytically or numerically
  # here we use a simple finite difference approximation
  proc dks_ddiffD(diffD: float): float =
    let h = 1.0 # big enough step size to overcome randomness of MC process!
    result = (rnd.simulateRmsTrans(rmsReal, gm, diffD + h, energy, isBackground, run) -
              rnd.simulateRmsTrans(rmsReal, gm, diffD - h, energy, isBackground, run)) / (2 * h)

  # define a stopping criterion
  let max_iter = 100 # maximum number of iterations
  let tol = 1e-6 # tolerance for change in diffD

  # iterate Newton's method until convergence or maximum iterations reached
  var iter = 0 # iteration counter
  var bestEstimate = diffD
  var bestLoss = ks

  proc lr(iter: int): float =
    let a = 10.0
    let c = 1.0 / 100.0
    let λ = -ln(c) * c
    result = a * exp(-λ * iter.float)
  #var lr = 5.0

  var numBad = 0 # keeps track of number of iters without improvements. Alternative: adjust gradient?
  const MaxConsBad = 10 # maximum consecutive bad samples
  while iter < max_iter and numBad < MaxConsBad:
    let delDiffD = dks_ddiffD(diffD)
    let diffD_new = diffD - lr(iter) * delDiffD # update diffD using Newton's formula
    echo "iter = ", iter, ", ∂σT = ", delDiffD, " and new diff ", diffD_new, " KS diff? ", ks, " and best: ", bestLoss

    if abs(diffD_new - diffD) < tol: # check if change is small enough
      echo "Stopping iteratoin due to difference = ", abs(diffD_new - diffD)
      break # stop the loop
    diffD = clamp(abs(diffD_new), 400.0, 800.0) # assign new value to diffD
    ks = rnd.simulateRmsTrans(rmsReal, gm, diffD, energy, isBackground, run) # update ks using new diffD
    iter += 1 # increment iteration counter

    if ks < bestLoss:
      bestLoss = ks
      bestEstimate = diffD
      numBad = 0

      # copy the last plot (and potential CSVs) to one based on this run number
      copyFile(HistoPath / HistoFile, HistoPath / &"histo_sim_vs_real_run_{run}_loss_{bestLoss}_sigma_{bestEstimate}.pdf")
      for f in walkFiles(HistoPath / HistoFile & "*csv"):
        let (path, file, ext) = f.splitFile() # e.g. `/foo/bar/baz/run_89_rmsTransverseFit.pdf_geom_1_line.csv`
        let (_, filePdf, extPdf) = file.splitFile() # e.g. `run_89_rmsTransverseFit.pdf_geom_1_line`
        copyFile(f, HistoPath / &"histo_sim_vs_real_run_{run}_loss_{bestLoss}_sigma_{bestEstimate}{extPdf}.csv")
    else:
      inc numBad
      if numBad > MaxConsBad div 2:
        # reset estimate back to best estimate, i.e. try again with smaller 'learning rate'
        echo "Resetting current estimate: ", diffD, " back to previous best: ", bestEstimate
        diffD = bestEstimate

  result = (bestEstimate, bestLoss) # diffD # return the final value of diffD as the optimal diffusion parameter

  # now apply correction from 1D to 2D?
  #result = sqrt(2.0) * result

when false:
  while abs(0.1 * diffD) > 1e-2:
    echo "Starting with diffusion : ", diffD, " and ks result ", ks
    diffusion = diffusion - ks.d * 1e19
    diffD = D(diffusion, 1.0)
    ks = rnd.simulateRmsTrans(rmsReal, gm, diffD, energy)
    echo "which is now: ", diffD

  # 1. MC simulate rms transverse based on energy & diffusion
  # 2. compute KS test of real and simulated
  # 3. update diffusion using autograd

## XXX: MAYBE adjust the scaling parameter based on the whole range? Or
## generally on the 'absolute' scale of sigma? using 2 sigma if sigma
## is huge is dumb!
#const Scale = 2.0
import nimhdf5 / serialize_tables
const Scale = 1.65
const MaxZ = 3.0
const CacheTabFile = "/dev/shm/cacheTab_diffusion_runs.h5"
type
  TabKey = int
  #        ^-- run number
  TabVal = (float, float)
  #         ^-- diffusion value
  #                ^-- loss of gradient descent via Cramér-von Mises
  CacheTabTyp = Table[TabKey, TabVal]
## Global helper table that stores diffusion values determined from
## the RMS data of this run
var CacheTab =
  if fileExists(CacheTabFile):
    tryDeserializeH5[CacheTabTyp](CacheTabFile)
  else:
    initTable[TabKey, TabVal]()
proc runAvailable(run: int): bool =
  if run in CacheTab:
    result = true
  else:
    if fileExists(CacheTabFile):
      let tab = tryDeserializeH5[CacheTabTyp](CacheTabFile)
      # merge `tab` and `CacheTab`
      for k, v in tab:
        CacheTab[k] = v # overwrite possible existing keys in table
      # write merged file
      CacheTab.tryToH5(CacheTabFile)
    result = run in CacheTab # still not in: not available

const BackgroundDiffusionCorrection = 40.0
proc getDiffusionForRun*(run: int, isBackground: bool): float =
  ## Attempts to return the diffusion determined for run `run`. If the CacheTab
  ## does not contain the run yet, it raises an exception.
  if run in CacheTab:
    result = CacheTab[run][0]
  else:
    raise newException(ValueError, "Diffusion for run " & $run & " not determined yet. Please " &
      "generate the cache table of diffusion values using `determineDiffusion`.")
  if isBackground:
    result = result - BackgroundDiffusionCorrection

proc getDiffusion*(rmsT: seq[float],
                   isBackground: bool,
                   run = -1,
                   energy = 0.0.keV,
                   useCache = true): (float, float) =
  ## DataFrame must contain `rmsTransverse` column!
  ##
  ## It returns the diffusion constant, `D_T`. That is:
  ##
  ## `σ_T = D_T · √( drift distance )`
  ##
  ## and not the diffusion coefficient `σ_T`!
  ##
  ## XXX: The `energy` argument here is a bit dangerous. We shouldn't use the
  ## energy of the target energy we want to simulate later for an NN cut value,
  ## because the energy that is needed is the one that actually reproduces the input
  ## data best to determine its diffusion only!
  var (df, pRes, _) = fitRmsTransverse(rmsT.toTensor, 100)
  let scaleShift = Scale * abs(pRes[2])
  #let scaleShift = 0.25
  var arg = pRes[1] + scaleShift # 0.15 #scaleShift
  #echo "Diffusion to use from pRes + 0.15 : ", arg / sqrt(MaxZ) * 1000.0
  #echo "Diffusion to use from pRes + scaleShift : ", arg / sqrt(MaxZ) * 1000.0
  #arg = determineFromData(rmsT)
  if useCache and run > 0 and runAvailable(run):
    result = CacheTab[run] # potentially re-read, but now available
  else:
    # try to reload first, maybe available then
    result = determineMonteCarlo(rmsT, energy, isBackground = isBackground, run = run)
    if useCache and run > 0:
      CacheTab[run] = result
      # serialize file
      CacheTab.tryToH5(CacheTabFile)
  #result = arg / sqrt(MaxZ) * 1000.0
  if isBackground:
    ## All our background determinations have a bias to about ~30 larger than their closest 5.9 keV
    ## calibration run. Hence for now just subtract that to get the ~likely correct~ value.
    ## Run this file as a standalone on all data files and see the generated plot.
    result[0] = result[0] - BackgroundDiffusionCorrection

  echo "USED RESULT: ", result
  if run > 0:
    let dfD = toDf(rmsT)
    echo "GENERATING PLOT:::::::::::\n\n\n\n"
    let yMax = df["yFit", float].max
    ggplot(dfD, aes("rmsT")) +
      geom_histogram(bins = 100, hdKind = hdOutline) +
      geom_line(data = df, aes = aes("xFit", "yFit"), color = "purple") +
      geom_linerange(aes = aes(x = arg, y = 0.0, yMin = 0.0, yMax = yMax), color = "green") +
      scale_y_continuous() +
      ylab("rmsT") +
      theme_font_scale(1.0, family = "serif") +
      ggsave(&"{HistoPath}/run_{run}_rmsTransverseFit.pdf")
    echo "RESULT COMES OUT TO ", result, " from fit parameters? ? ", pRes, " FOR RUN: ", run
  #if true: quit()

## XXX: THIS IS OUTDATED and dangerous! It gives the wrong numbers to use for the
## DF containing the real data read via `readValidDsets`, i.e. give the wrong inputs
## to the NN when evaluating the effective efficiency (or anything working with network
## prediction of real data!)
#proc getDiffusion*(df: DataFrame): float =
#  ## DataFrame must contain `rmsTransverse` column!
#  let (df, pRes, _) = fitRmsTransverse(df["rmsTransverse", float], 100)
#  result = (pRes[1] + Scale * abs(pRes[2])) / sqrt(MaxZ) * 1000.0

proc determineDiffusion(df: DataFrame, outpath: string, runInput: int, useCache, useTeX: bool) =
  var dfAll = newDataFrame()
  for (tup, subDf) in groups(df.group_by("run")):
    let run = tup[0][1].toInt
    if runInput > 0 and run != runInput: continue
    let runType = parseEnum[RunTypeKind](subDf["typ"].unique.item(string))
    let isBackground = runType == rtBackground
    let energy = if not isBackground: subDf["energy"].unique.item(float).keV
                 else: 0.keV

    # or use cache ? At least if available :)
    let (σT, loss) = getDiffusion(subDf["rms", float].toSeq1D,
                                  isBackground = isBackground,
                                  run = run,
                                  energy = energy,
                                  useCache = useCache)
    dfAll.add toDf( {"run" : run, "σT" : σT, "loss" : loss, "isBackground" : isBackground })

  let ylabel = if useTeX: r"$D_T$ [$\si{μm.cm^{-1/2}}$]" else: "D_T [μm/√cm]"
  ggplot(dfAll, aes("run", "σT", color = "loss", shape = "isBackground")) +
    geom_point() +
    xlab("Run number") + ylab(ylabel) +
    ggtitle("Determined diffusion constant of all CAST runs") +
    themeLatex(fWidth = 0.9, width = 600, baseTheme = singlePlot) +
    ggsave(&"{outpath}/σT_per_run.pdf", useTeX = useTeX, standalone = useTeX)

when isMainModule:
  import ../nn/io_helpers
  proc readData(fname: string): DataFrame =
    withH5(fname, "r"):
      let fileInfo = h5f.getFileInfo()
      let chipNumber = fileInfo.centerChip
      result = newDataFrame()
      for run in fileInfo.runs:
        let group = (recoBase() & $run)
        # need rmsTransverse to determine diffusion
        let grp = h5f[(group / "chip_" & $chipNumber).grp_str]
        # very weak cuts that should also apply well for background data that remove
        # obvious events that are noisy, and not well defined photons or tracks!
        ## NOTE: These do not apply to the CDL datasets!
        # try to read attribute and if not given `tfKind` skip this run
        var tfKind = tfMnCr12 ## Default we assume
        let runGrp = h5f[group.grp_str]
        if "tfKind" in runGrp.attrs:
          tfKind = parseEnum[TargetFilterKind](runGrp.attrs["tfKind", string])
        let passIdx = h5f.cutOnProperties(
          grp,
          crSilver,
          ("width", 0.0, 6.0),
          ("rmsTransverse", 0.0, 1.3),
          ("hits", 20.0, Inf))
        # now need `rms` to calculate diffusion
        let rms = h5f[group / &"chip_{chipNumber}/rmsTransverse", passIdx, float]
        let df = toDf({ "rms" : rms,
                        "run" : run,
                        "typ" : $fileInfo.runType,
                        "energy" : toXrayLineEnergy(tfKind)
        })
        result.add df

  proc main(fnames: seq[string],
            plotPath: string = "/home/basti/Sync",
            histoPlotPath: string = "/home/basti/Sync",
            run = -1,
            useCache = true,
            useTeX = false) =
    ## `plotPath` is where the D_T for all runs plot will be placed.
    ## `histoPlotPath` is the path for all the histograms of the distributions
    ## created during optimization.

    # Adjust the plotting path for the histograms
    HistoPath = histoPlotPath

    var df = newDataFrame()
    for f in fnames:
      df.add readData(f)
    determineDiffusion(df, plotPath, run, useCache, useTeX)

  #when isMainModule:
  import cligen
  dispatch main
