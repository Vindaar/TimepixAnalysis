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

type
  DataType* = enum
    dtSignal = "signal"
    dtBack = "back"

proc shuff*(x: seq[int]): seq[int] =
  ## Not in place shuffle
  result = x
  result.shuffle()

proc readCdlDset*(h5f: H5File, cdlPath, dset: string): DataFrame =
  ## Reads the given dataset (target/filter kind) `dset` according
  ## to the required cuts for X-rays for this target/filter.
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

proc shuffle*(df: DataFrame): DataFrame =
  ## Shufles the input data frame
  result = df.shallowCopy()
  result["Idx"] = toSeq(0 .. df.high).shuff()
  result = result.arrange("Idx")
  result.drop("Idx")

proc randomHead*(df: DataFrame, head: int): DataFrame =
  ## Returns the `head` elements of the input data frame shuffled
  #if head > df.len:
  #  echo "Dataframe: ", df
  result = df.shuffle().head(min(head, df.len))

proc readValidDsets*(h5f: H5File, path: string, readRaw = false,
                     typ = dtBack,
                     subsetPerRun = 0): DataFrame =
  ## Reads all data for the given run `path` (must be a chip path)
  ##
  ## `subsetPerRun` is an integer which if given only returns this many entries (random)
  ## from each run.
  let dsets = toSeq(ValidReadDsets - { igLikelihood }).mapIt(it.toDset(fkTpa))
  let evNumDset = igEventNumber.toDset(fkTpa)
  let grp = h5f[path.grp_str]
  if not readRaw:
    for dset in dsets:
      result[dset] = h5f.readAs(grp.name / dset, float)
    result[evNumDset] = h5f.readAs(grp.name / evNumDset, int)
  else:
    result = h5f.readRaw(grp.name)
  result["pass?"] = true # does not matter!
  ## XXX: Add run number?
  #result["runNumber
  result["Type"] = $typ
  if not readRaw:
    result = result.mutate(f{"Target" ~ `energyFromCharge`.toRefDset})
  if subsetPerRun > 0:
    result = result.randomHead(subsetPerRun)

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
      dfLoc = h5f.readCdlDset(cdlPath, bin)
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
  df["Type"] = $dtSignal
  result = df.shuffle()

proc prepareData*(h5f: H5File, run: int, readRaw: bool, typ = dtBack, subsetPerRun = 0): DataFrame =
  let path = "/reconstruction/run_$#/chip_3" % $run
  result = h5f.readValidDsets(path, readRaw, typ, subsetPerRun)

proc prepareBackground*(fname: string, run: int, readRaw: bool, subsetPerRun = 0): DataFrame =
  var h5f = H5open(fname, "r")
  result = h5f.prepareData(run, readRaw, subsetPerRun = subsetPerRun)
  discard h5f.close()

proc prepareAllBackground*(fname: string, readRaw: bool, subsetPerRun = 0): DataFrame =
  var h5f = H5open(fname, "r")
  for run, grp in runs(h5f):
    var df = h5f.prepareData(run, readRaw, subsetPerRun = subsetPerRun)
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
    if typ[i] == $dtSignal:
      target[i, _] = [1, 0].toRawTensorFromScalar #.toTensor.convertRawTensor()
    elif typ[i] == $dtBack:
      target[i, _] = [0, 1].toRawTensorFromScalar #toTensor.convertRawTensor()
    else:
      doAssert false, "Invalid type field " & $typ[i]
  result = (input, target)

proc toNimSeq*[T](t: RawTensor): seq[T] =
  doAssert t.sizes().len == 1
  result = newSeq[T](t.size(0))
  for i in 0 ..< result.len:
    result[i] = t[i].item(T)
