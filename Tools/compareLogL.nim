import std / [os, strutils, sequtils, strformat]
import nimhdf5
import seqmath
import helpers / utils

import ingrid / private / [likelihood_utils, hdf5_utils, ggplot_utils, geometry, cdl_cuts]
import ingrid / ingrid_types

when true:
  const CdlFile = "/mnt/1TB/CAST/CDL_2019/calibration-cdl-2018.h5"
  const RefFile = "/mnt/1TB/CAST/CDL_2019/XrayReferenceFile2018.h5"
  const Suffix = "Local"
else:
  const CdlFile = "/t/calibration-cdl-2018.h5"
  const RefFile = "/t/XrayReferenceFile2018.h5"
  const Suffix = "Void"

proc readRefEx(): tuple[eccs, ldiv, frac: histTuple] =
  var h5f = H5open(RefFile, "r")
  let d = toRefDset(5.9)
  proc r(x: string): histTuple =
    let dset = h5f[(d / x).dset_str]
    let data = dset.readAs(float).reshape2D(dset.shape).splitSeq(float)
    result = (bins: data[0], hist: data[1])
  result = (eccs: r("eccentricity"), ldiv: r("lengthDivRmsTrans"), frac: r("fractionInTransverseRms"))
  discard h5f.close()

proc readNew(): (histTuple, histTuple, histTuple) =
  result = genRefData(CdlFile, toRefDset(5.9), yr2018, igEnergyFromCharge)

block ReferenceComparison:
  let (eccs, ldiv, frac) = readRefEx()
  let (ecc2, ldi2, fra2) = readNew()
  proc p(s: string, o, n: histTuple) =
    doAssert o.bins.len == n.bins.len
    let df = toDf({"old" : o.hist, "new" : n.hist, "bins" : o.bins})
      .gather(["old", "new"], "type", "val")
    echo df.pretty(-1)
    ggplot(df, aes("bins", "val", fill = "type")) +
      geom_histogram(bins = 100, position = "identity", stat = "identity", alpha = 0.5, hdKind = hdOutline) +
      ggsave(&"/t/{s}_old_new_{Suffix}.pdf")
  p("ecc", eccs, ecc2)
  p("ldiv", ldiv, ldi2)
  p("frac", frac, fra2)

proc readExisting(): seq[float] =
  var h5f = H5open(CdlFile, "r")
  result = h5f[cdlPath(toRefDset(5.9)) / "likelihood", float]
  discard h5f.close()

proc getNew(): seq[float] =
  let (logL, E) = buildLogLHist(CdlFile, toRefDset(5.9), yr2018, igEnergyFromCharge)
  result = logL

proc calcLikeOld(): seq[float] =
  # get old reference
  let (eccs, ldiv, frac) = readRefEx()
  when false:
    withLogLFilterCuts(CdlFile, toRefDset(5.9), yr2018, igEnergyFromCharge, LogLCutDsets):
      result.add calcLogL(data[igEccentricity][i],
                          data[igLengthDivRmsTrans][i],
                          data[igFractionInTransverseRms][i],
                          eccs, ldiv, frac)
  else:
    ## XXX: CAREFUL: The following applies the charge cut around the main peak for the data
    ##. Think about if this is desired here!
    withCdlData(CdlFile, toRefDset(5.9), yr2018, igEnergyFromCharge, LogLCutDsets):
      for i {.inject.} in 0 ..< energy.len:
        result.add calcLogL(data[igEccentricity][i],
                            data[igLengthDivRmsTrans][i],
                            data[igFractionInTransverseRms][i],
                            eccs, ldiv, frac)

block LogLComparisons:
  let logLEx = readExisting()
  let logLNew = getNew()

  let df = toDf({"old" : logLEx, "new" : logLNew, "likeOld" : calcLikeOld()}).gather(["old", "new", "likeOld"], "type", "val")
    .dropNull("val")
    .filter(f{`val` < 30.0})
  echo df

  ggplot(df, aes("val", fill = "type")) +
    geom_histogram(bins = 100, position = "identity", alpha = 0.5, hdKind = hdOutline) +
    ggsave(&"/t/compare_old_logl_new_logl_{Suffix}.pdf")
