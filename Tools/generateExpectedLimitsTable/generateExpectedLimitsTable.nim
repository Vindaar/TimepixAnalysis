
import os, nimhdf5, datamancer, strutils, sugar
import ingrid / ingrid_types

#[
This is a companion tool to `TimepixAnalysis/Analysis/runLimits.nim`.

Given the path to where the former outputs the log and H5 files to, this script reads
all files from the path and generates an Org table of the different setups and their
respective expected limits.
]#

type
  LimitData = object
    expectedLimit: float
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

proc expLimit(limits: seq[float]): float =
  result = sqrt(limits.percentile(50)) * 1e-12

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

proc readLimit(fname: string): LimitData =
  var h5f = H5open(fname, "r")
  let limits = h5f["/limits", float]
  let noCands = h5f.attrs["limitNoSignal", float]
  let vetoes = readVetoes(h5f)
  let effs = readEfficiencies(h5f)
  result = LimitData(expectedLimit: expLimit(limits),
                     limitNoSignal: noCands,
                     vetoes: vetoes,
                     eff: effs)

proc asDf(limit: LimitData): DataFrame =
  ## Calling it `toDf` causes issues...
  let typ = if fkMLP in limit.vetoes: "MLP"
             else: "LnL"
  let eff = if fkMLP in limit.vetoes: limit.eff.nnEffectiveEff
            else: limit.eff.signalEff
  result = toDf({ "ε" : eff,
                  "Type" : typ,
                  "Scinti" : fkScinti in limit.vetoes,
                  "FADC" : fkFadc in limit.vetoes,
                  "ε_FADC" : 1.0 - (1.0 - limit.eff.vetoPercentile) * 2.0,
                  "Septem" : fkSeptem in limit.vetoes,
                  "Line" : fkLineVeto in limit.vetoes,
                  "eccLineCut" : limit.eff.eccLineVetoCut,
                  "ε_Septem" : limit.eff.septemVetoRandomCoinc,
                  "ε_Line" : limit.eff.lineVetoRandomCoinc,
                  "ε_SeptemLine" : limit.eff.septemLineVetoRandomCoinc,
                  "Total eff." : limit.eff.totalEff,
                  "Limit no signal" : limit.limitNoSignal,
                  "Expected Limit" : limit.expectedLimit })

proc main(path = @["/t/lhood_outputs_adaptive_fadc_limits/"],
          prefix = @["mc_limit_lkMCMC_skInterpBackground_nmc_1000"]) =
  var df = newDataFrame()
  doAssert path.len == prefix.len, "Need one prefix for each path!"
  for i, p in path:
    let pref = prefix[i]
    for f in walkFiles(p / pref & "*.h5"):
      echo "File: ", f
      let limit = readLimit(f)
      df.add asDf(limit)
  echo df.arrange("Expected Limit").toOrgTable()

when isMainModule:
  import cligen
  dispatch main
