import os, nimhdf5, datamancer, strutils, sugar, sequtils, stats, strformat, measuremancer
import ingrid / ingrid_types

#[
This is a companion tool to `TimepixAnalysis/Analysis/runLimits.nim`.

Given the path to where the former outputs the log and H5 files to, this script reads
all files from the path and generates an Org table of the different setups and their
respective expected limits.
]#

type
  LimitData = object
    nmc: int # Number of toy limits
    expectedLimit: float
    expLimitVariance: float
    expLimitStd: float
    l5: float
    l16: float
    l25: float
    l75: float
    l84: float
    l95: float
    limitNoSignal: float
    vetoes: set[LogLFlagKind]
    eff: Efficiency

  Efficiency = object
    totalEff: float # total efficiency multiplier based on signal efficiency of lnL cut, FADC & veto random coinc rate
    signalEff: float # the lnL cut signal efficiency used in the inputs
    nnSignalEff: float # target signal efficiency of MLP
    nnEffectiveEff: float # effective efficiency based on
    nnEffectiveEffStd: float
    eccLineVetoCut: float # the eccentricity cutoff for the line veto (affects random coinc.)
    vetoPercentile: float # if FADC veto used, the percentile used to generate the cuts
    septemVetoRandomCoinc: float # random coincidence rate of septem veto
    lineVetoRandomCoinc: float # random coincidence rate of line veto
    septemLineVetoRandomCoinc: float # random coincidence rate of septem + line veto

  CouplingKind = enum
    ck_g_ae²       ## We vary the `g_ae²` and leave `g_aγ²` fully fixed
    ck_g_aγ⁴       ## We vary the `g_aγ⁴` and leave `g_ae²` fully fixed (and effectively 'disabled'); for axion-photon searches
    ck_g_ae²·g_aγ² ## We vary the *product* of `g_ae²·g_aγ²`, i.e. direct `g⁴` proportional search.
                   ## Note that this is equivalent in terms of the limit!
    ck_g_ae·g_aγ   ## We vary the *product* of `g_ae·g_aγ`, but in `g²` form. This is *NOT* equivalent and gives the ``WRONG``
                   ## limit. Only for illustration!
    ck_β⁴          ## We vary `β`, the chameleon coupling. For chameleon searches.


proc convertLimit(limit: float, coupling: CouplingKind): float =
  case coupling
  of ck_g_ae²: result = sqrt(limit) * 1e-12
  of ck_g_ae²·g_aγ²: result = sqrt(limit)
  of ck_g_aγ⁴: result = pow(limit, 0.25)
  of ck_β⁴: result = pow(limit, 0.25)
  else: doAssert false

proc expLimit(limits: seq[float], coupling: CouplingKind): float =
  result = limits.median.convertLimit(coupling)

import random

template withBootstrap(rnd: var Rand, samples: seq[float], num: int, body: untyped): untyped =
  let N = samples.len
  for i in 0 ..< num:
    # resample
    var newSamples {.inject.} = newSeq[float](N)
    for j in 0 ..< N:
      newSamples[j] = samples[rnd.rand(0 ..< N)] # get an index and take its value
    # compute our statistics
    body

proc expLimitVarStd(limits: seq[float], coupling: CouplingKind): (float, float) =
  var rnd = initRand(12312)
  let limits = limits.mapIt(convertLimit(it, coupling)) # rescale limits
  const num = 1000
  var medians = newSeqOfCap[float](num)
  withBootstrap(rnd, limits, num):
    medians.add median(newSamples)
  #echo "Medians? ", medians
  result = (variance(medians), standardDeviation(medians))

proc readVetoes(h5f: H5File): set[LogLFlagKind] =
  let flags = h5f["/ctx/logLFlags", string]
  for f in flags:
    result.incl parseEnum[LogLFlagKind](f)

import std / strscans
proc tryParseEccLine(s: string): float =
  let (success, _, val) = scanTuple(s, "$*_eccCutoff_$f")
  if success:
    result = val

proc readEfficiencies(h5f: H5File): Efficiency =
  let eff = h5f["/ctx/eff".grp_str]
  if "eccLineVetoCut" in eff.attrs:
    result = h5f.deserializeH5[:Efficiency](eff.name)
  else:
    result = h5f.deserializeH5[:Efficiency](eff.name, exclude = @["eccLineVetoCut"])
    result.eccLineVetoCut = tryParseEccLine(h5f.name)

proc readLimit(fname: string, coupling: CouplingKind): LimitData =
  var h5f = H5open(fname, "r")
  let nmc = h5f["/".grp_str].attrs["nmc", int]
  let limits = h5f["/limits", float]
  let noCands = h5f.attrs["limitNoSignal", float]
  let vetoes = readVetoes(h5f)
  let effs = readEfficiencies(h5f)
  let (variance, std) = expLimitVarStd(limits, coupling)
  echo "Standard deviation of existing limits: ", limits.standardDeviation
  let
    l5  = convertLimit(limits.percentile(5), coupling)
    l16 = convertLimit(limits.percentile(16), coupling)
    l25 = convertLimit(limits.percentile(25), coupling)
    l75 = convertLimit(limits.percentile(75), coupling)
    l84 = convertLimit(limits.percentile(84), coupling)
    l95 = convertLimit(limits.percentile(95), coupling)
  result = LimitData(nmc: nmc,
                     expectedLimit: expLimit(limits, coupling),
                     expLimitVariance: variance,
                     expLimitStd: std,
                     limitNoSignal: convertLimit(noCands, coupling),
                     l5: l5, l16: l16, l25: l25, l75: l75, l84: l84, l95: l95,
                     vetoes: vetoes,
                     eff: effs)

proc getSuffix(coupling: CouplingKind, power = "¹"): string =
  case coupling
  of ck_g_ae²:       result = &" [GeV⁻{power}]"
  of ck_g_aγ⁴:       result = &" [GeV⁻{power}]"
  of ck_g_ae²·g_aγ²: result = &" [GeV⁻{power}]"
  of ck_β⁴: discard # nothing to add!
  else: doAssert false

proc expLimitCol(coupling: CouplingKind): string =
  result = "Expected limit"
  result.add getSuffix(coupling)

proc limitNoSigCol(coupling: CouplingKind): string =
  result = "Limit no signal"
  result.add getSuffix(coupling)

proc expLimitVarCol(coupling: CouplingKind): string =
  result = "Exp. limit variance"
  result.add getSuffix(coupling, power = "²")

proc expLimitStdCol(coupling: CouplingKind): string =
  result = "Exp. limit σ"
  result.add getSuffix(coupling)

proc asDf(limit: LimitData, coupling: CouplingKind): DataFrame =
  ## Calling it `toDf` causes issues...
  let typ = if fkMLP in limit.vetoes: "MLP"
             else: "LnL"
  let eff = if fkMLP in limit.vetoes: limit.eff.nnEffectiveEff
            else: limit.eff.signalEff
  let septem = fkSeptem in limit.vetoes
  let line = fkLineVeto in limit.vetoes
  let fadc = fkFadc in limit.vetoes
  result = toDf({ "ε_eff" : eff,
                  "nmc" : limit.nmc,
                  "Type" : typ,
                  "Scinti" : fkScinti in limit.vetoes,
                  "FADC" : fadc,
                  "ε_FADC" : 1.0 - (1.0 - limit.eff.vetoPercentile) * 2.0,
                  "Septem" : septem,
                  "Line" : line,
                  "eccLineCut" : limit.eff.eccLineVetoCut,
                  "ε_Septem" : if septem and not line: limit.eff.septemVetoRandomCoinc else: 1.0,
                  "ε_Line" : if line and not septem: limit.eff.lineVetoRandomCoinc else: 1.0,
                  "ε_SeptemLine" : if septem and line: limit.eff.septemLineVetoRandomCoinc else: 1.0,
                  "ε_total" : limit.eff.totalEff,
                  limitNoSigCol(coupling) : limit.limitNoSignal,
                  expLimitCol(coupling) : limit.expectedLimit,
                  expLimitVarCol(coupling) : limit.expLimitVariance,
                  expLimitStdCol(coupling) : limit.expLimitStd,
                  "P_5": limit.l5, "P_16": limit.l16, "P_25": limit.l25, "P_75": limit.l75, "P_84": limit.l84, "P_95": limit.l95,
  })

proc toTightLayout(df: DataFrame, coupling: CouplingKind, precision: int): DataFrame =
  ## Turns the table into a tight layout version. This means dropping some columns
  ## and merging others.
  # 1. merge septem & line veto to one
  result = df.mutate(f{bool -> string: "Veto" ~ (
    if `Septem` and `Line`: "SL"
    elif `Septem`: "S"
    elif `Line`: "L"
    else: "-"
  )})
  # 2. merge Type with veto
  result = result.mutate(f{string -> string: "Type" ~ (
    `Type` & " " & `Veto`
  )})
  # 3. merge exp limit & exp limit σ into one
  proc mergeConvert(a, b: float): string =
    pretty( a ± b, precision = precision, merge = true )
  result = result.mutate(f{float -> string: "Expected" ~ (
    mergeConvert(idx(expLimitCol(coupling)), idx(expLimitStdCol(coupling)))
  )})
  result = result.drop([expLimitCol(coupling), expLimitStdCol(coupling), "Septem", "Line", "Scinti", "FADC", "eccLineCut",
                        "ε_FADC", "ε_Line", "ε_SeptemLine", "ε_Septem",
                        expLimitVarCol(coupling), limitNoSigCol(coupling),
                        "Veto"])

proc main(path: seq[string] = @[],
          prefix: seq[string] = @[],
          coupling = ck_g_ae²,
          precision = 4) =
  var df = newDataFrame()
  doAssert path.len == prefix.len, "Need one prefix for each path!"
  for i, p in path:
    let pref = prefix[i]
    for f in walkDirRec(p):
      let fname = extractFilename(f)
      if fname.startsWith(pref) and fname.endsWith(".h5"):
        echo "File: ", fname
        let limit = readLimit(f, coupling)
        df.add asDf(limit, coupling)

  let dfA = df.arrange(expLimitCol(coupling))
  echo dfA.toOrgTable(precision = precision)
  echo dfA.toTightLayout(coupling, precision).toOrgTable(precision = precision)

when isMainModule:
  import cligen
  dispatch main
