import std / [os, sequtils, tables, stats, strformat, strutils, options]
import ingrid / [tos_helpers, ingrid_types]
from ingrid / calibration import iterGainSlicesFromDset
import pkg / [nimhdf5, datamancer, ggplotnim, seqmath]
import unchained


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

import ../NN_playground/io_helpers

proc readRealCalib(f: string, typ: string, eLow, eHigh: float,
                   tfKind: TargetFilterKind,
                   validDsets: set[InGridDsetKind] = {}): DataFrame =
  if validDsets.card == 0:
    result = readCalibData(f, typ, eLow, eHigh, tfKind = some(tfKind))
  else:
    result = readCalibData(f, typ, eLow, eHigh, tfKind = some(tfKind), validDsets = validDsets)
  result = result.rename(f{"DataType" <- "CalibType"})
  result["isFake?"] = false
  ## Overwrite the `Target` field in the read data, because that is based on calling `toRefDset`
  ## for each cluster. Some may just slightly cross a boundary causing the rest of our code
  ## that uses the target as a factor to break
  result["Target"] = $tfKind

proc readCdlData(cdlFile: string,
                 tfKind: TargetFilterKind,
                 eLow, eHigh: float,
                 typeSuffix = "_CDL"): DataFrame =
  ## Reads the data fro the given tfKind filtering to the desired energy range
  const dsets = ValidReadDsets + {igEnergyFromCdlFit} - {igLikelihood}
  proc filterAfter(df: DataFrame, frm, to: float): DataFrame =
    let E = igEnergyFromCdlFit.toDset()
    df.filter(f{float -> bool: idx(E) >= frm and idx(E) <= to})
      .drop([E])
  let energy = toXrayLineEnergy(tfKind)
  result = filterAfter(
    readRealCalib(cdlFile, $energy & $typeSuffix,
                  0.0, Inf,
                  tfKind = tfKind,
                  validDsets = dsets),
    eLow, eHigh
  )

import mpfit
proc fitRmsTransverse(data: Tensor[float], bins: int): (DataFrame, seq[float], string) =
  var (hist, bins) = histogram(data.toSeq1D, bins = bins)
  let binWidth = bins[1] - bins[0]
  bins = bins[0 ..< ^1].mapIt(it + binWidth / 2.0) # remove last bin edge, move half bin right to have centers
  let hArgMax = hist.argmax
  let binArgMax = bins[hist.argmax]
  let histMax = hist[hArgMax].float

  let binLowerIdx = (hArgMax - (0.05 * bins.len.float)).round.int
  let binUpper = bins[^1] # bins.lowerBound(data.percentile(98))]

  let params = @[histMax, binArgMax, abs(binArgMax - binUpper)]
  let xT = bins[binLowerIdx .. ^1]
  let yT = hist[binLowerIdx .. ^1].mapIt(it.float)
  let eT = yT.mapIt(1.0) #sqrt(x))

  let (pRes, res) = fit(gaussSeq, params, xT, yT, eT)
  let xFit = linspace(binArgMax, binUpper, 200)
  let yFit = xFit.mapIt(gauss(pRes, it))
  let resStr = resultStr(pRes, res.error, res.chiSq)

  echo resStr
  result = (toFitDf(pRes, xT.min, xT.max), pRes, resStr)

proc getDiffusion*(rmsT: seq[float]): float =
  ## DataFrame must contain `rmsTransverse` column!
  let (df, pRes, _) = fitRmsTransverse(rmsT.toTensor, 100)
  result = (pRes[1] + 1.5 * pRes[2]) / sqrt(3.0) * 1000.0

proc getDiffusion*(df: DataFrame): float =
  ## DataFrame must contain `rmsTransverse` column!
  let (df, pRes, _) = fitRmsTransverse(df["rmsTransverse", float], 100)
  result = (pRes[1] + pRes[2]) / sqrt(3.0) * 1000.0

proc determineDiffusion(df: DataFrame) =
  var dfAll = newDataFrame()
  for (tup, subDf) in groups(df.group_by("runNumber")):
    let run = tup[0][1].toInt
    let dfL = subDf # .filter(f{`rmsTransverse` < percentile(`rmsTransverse`, 99)})
    const bins = 100
    let (dfFit, pRes, text) = dfL["rmsTransverse", float].fitRmsTransverse(bins)

    ggplot(dfL, aes("rmsTransverse")) +
      geom_histogram(bins = bins, hdKind = hdOutline) +
      geom_line(data = dfFit, aes = aes("xFit", "yFit"), color = "purple") +
      ggsave(&"/tmp/run_{run}_rmsTransverseFit.pdf")

    dfAll.add toDf({ "rmsT" : pRes[1] + pRes[2],
                     "run" : run })

  dfAll = dfAll.mutate(f{"σT" ~ `rmsT` / sqrt(3.0)})
  echo dfAll.pretty(-1)

  ggplot(dfAll, aes("run", "rmsT")) +
    geom_point() +
    ggsave("/tmp/rmsT_per_run.pdf")

  ggplot(dfAll, aes("run", "σT")) +
    geom_point() +
    ggsave("/tmp/sigmaT_per_run.pdf")

proc main(fnames: seq[string],
          cdlFile: string = "",
          plotPath: string = "") =
  ## The `cdlFile` currently must be the `CDL_2019_Reco.h5` file!
  let ctx = initLikelihoodContext(CdlFile,
                                  year = yr2018,
                                  energyDset = igEnergyFromCharge,
                                  region = crSilver,
                                  timepix = Timepix1,
                                  morphKind = mkLinear) # morphing to plot interpolation
  var df = newDataFrame()
  for c in fnames:
    # 3 keV
    #df.add readRealCalib(c, "3.0", 2.5, 3.5, tfAgAg6)
    # 5.9 keV
    if "Calibration" in c:
      df.add readRealCalib(c, "5.9", 5.0, 7.0, tfMnCr12)
    else:
      df.add prepareAllBackground(c, false)

    determineDiffusion(df)

  #if cdlFile.len > 0:
  #  df.add readCdlData(cdlFile, tfAgAg6, 2.5, 3.5)
  #  df.add readCdlData(cdlFile, tfMnCr12, 5.0, 7.0)
  #  df.add readCdlData(cdlFile, tfCEpic0_6, 0.0, 0.4)

when isMainModule:
  import cligen
  dispatch main
