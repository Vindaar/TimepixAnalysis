from random import shuffle
import std / [sequtils, os, strutils]

import ingrid / [tos_helpers, ingrid_types]
import pkg / [nimhdf5, datamancer]

const ValidReadDSets* = XrayReferenceDsets - { igNumClusters,
                                               igFractionInHalfRadius,
                                               igRadiusDivRmsTrans,
                                               igRadius, igBalance,
                                               igLengthDivRadius } + {
                                                 igCenterX, igCenterY }

const ValidDsets* = ValidReadDSets - { igLikelihood, igCenterX, igCenterY, igHits, igEnergyFromCharge, igTotalCharge }

proc shuff*(x: seq[int]): seq[int] =
  ## Not in place shuffle
  result = x
  result.shuffle()

proc readDsets*(h5f: H5File, cdlPath, dset: string): DataFrame =
  result = newDataFrame()
  let dsets = concat(toSeq(ValidReadDsets - { igLikelihood } ), @[igEventNumber])
  var passData: array[InGridDsetKind, seq[float]]
  for s in mitems(passData):
    s = newSeqOfCap[float](50_000)
  withLogLFilterCuts(cdlPath, dset, yr2018, igEnergyFromCharge, dsets):
    for d in dsets:
      passData[d].add data[d][i]
  for d in dsets:
    result[d.toDset(fkTpa)] = passData[d]

proc readRaw*(h5f: H5File, grpName: string, idxs: seq[int] = @[]): DataFrame =
  ## XXX: Need to implement filtering to suitable data for CDL data!
  result = newDataFrame()
  let
    xs = h5f[grpName / "x", special_type(uint8), uint8]
    ys = h5f[grpName / "y", special_type(uint8), uint8]
    ev = h5f.readAs(grp_name / "eventNumber", int)
  doAssert xs.len == ev.len
  var xsAll = newSeqOfCap[int](xs.len * 100)
  var ysAll = newSeqOfCap[int](xs.len * 100)
  var evAll = newSeqOfCap[int](xs.len * 100)
  for i in idxs:
    for j in 0 ..< xs[i].len:
      xsAll.add xs[i][j].int
      ysAll.add ys[i][j].int
      evAll.add ev[i]
  result = toDf({"x" : xsAll, "y" : ysAll, "eventNumber" : evAll})

proc prepareCDL*(readRaw: bool,
                 cdlPath = "/home/basti/CastData/data/CDL_2019/calibration-cdl-2018.h5"
               ): DataFrame =
  var h5f = H5file(cdlPath, "r")
  let tb = getXrayRefTable()
  var df = newDataFrame()
  for k, bin in tb:
    var dfLoc = newDataFrame()
    if not readRaw:
      # read dsets only returns those indices that pass
      dfLoc = h5f.readDsets(cdlPath, bin)
      dfLoc["pass?"] = true
    else:
      doAssert false, "Raw data reading currently not supported!"
      when false:
        ## XXX: We need a way to read the indices _for each run_ of the CDL data such that
        ## we can read the raw x,y,ToT data (cannot be handled by our `array[InGridDsetKind, seq[float]]`
        ## data type we use now.
        ## Maybe an iterator similar which yields the
        ## - (tfKind, runNumber)
        ## - event indices that pass the cuts in each?
        let pass = h5f.buildLogLHist(bin)
        var idxs = newSeqOfCap[int](pass.len)
        for i, p in pass:
          if p:
            idxs.add i
        dfLoc = h5f.readRaw(cdlPrefix($yr2018) & bin, idxs)
        dfLoc["pass?"] = true
    dfLoc["Target"] = bin
    # remove all that don't pass
    ## XXX: this should now be unnecessary
    dfLoc = dfLoc.filter(f{bool: idx("pass?") == true})
    df.add dfLoc
  discard h5f.close()
  df["Type"] = "signal"
  result = df
  result["Idx"] = toSeq(0 ..< result.len).shuff()
  result = result.arrange("Idx")
  #result = result[0 ..< 55000]
  result.drop("Idx")
  #echo result

proc prepareBackground*(h5f: H5File, run: int, readRaw: bool): DataFrame =
  let path = "/reconstruction/run_$#/chip_3" % $run
  let dsets = toSeq(ValidReadDsets - { igLikelihood }).mapIt(it.toDset(fkTpa))
  let evNumDset = "eventNumber"
  let grp = h5f[path.grp_str]
  if not readRaw:
    for dset in dsets:
      result[dset] = h5f.readAs(grp.name / dset, float)
    result["eventNumber"] = h5f.readAs(grp.name / "eventNumber", int)
  else:
    result = h5f.readRaw(grp.name)
  result["pass?"] = true # does not matter!
  #result["runNumber
  result["Type"] = "back"

  let energyBins = getEnergyBinning()
  let targetTab = getXrayRefTable()
  ## TODO: the following is broken? CHECK!
  # result = result.mutate(f{"Target" ~ targetTab[energyBins.lowerBound(idx("energyFromCharge"))]})
  if not readRaw:
    result = result.mutate(f{"Target" ~ `energyFromCharge`.toRefDset}) # ["energyFromCharge", float].toSeq1D.mapIt(targetTab[energyBins.lowerBound(it)])

proc prepareBackground*(fname: string, run: int, readRaw: bool): DataFrame =
  var h5f = H5open(fname, "r")
  result = h5f.prepareBackground(run, readRaw)
  discard h5f.close()

proc prepareAllBackground*(fname: string, readRaw: bool): DataFrame =
  var h5f = H5open(fname, "r")
  for run, grp in runs(h5f):
    var df = h5f.prepareBackground(run, readRaw)
    df["runNumber"] = run
    result.add df

  # filter gold region
  result = result.filter(f{float -> bool: inRegion(`centerX`, `centerY`, crGold)})
  discard h5f.close()

{.experimental: "views".}
import flambeau / [flambeau_raw, tensors]
proc toInputTensor*(df: DataFrame): (RawTensor, RawTensor) {.noInit.} =
  ## Converts an appropriate data frame to a tuple of input / target tensors
  let cols = ValidDsets.card
  var input = rawtensors.zeros(df.len * cols).reshape(sizes = [df.len.int64, cols].asTorchView())
  for i, c in ValidDsets.toSeq.mapIt(it.toDset(fkTpa)).sorted:
    let xp = fromBlob[float](cast[pointer](df[c, float].unsafe_raw_offset()), df.len).convertRawTensor()
    input[_, i] = xp
  var target = rawtensors.zeros(df.len * 2).reshape([df.len.int64, 2].asTorchView())
  let typ = df["Type", string]
  for i in 0 ..< typ.size:
    if typ[i] == "signal":
      target[i, _] = [1, 0].toRawTensorFromScalar #.toTensor.convertRawTensor()
    else:
      target[i, _] = [0, 1].toRawTensorFromScalar #toTensor.convertRawTensor()
  result = (input, target)
