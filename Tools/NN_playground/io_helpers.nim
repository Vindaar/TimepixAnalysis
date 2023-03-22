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

# include total charge
#const ValidDsets* = { igEccentricity, igLengthDivRmsTrans, igFractionInTransverseRms }
const ValidDsets* = ValidReadDSets - { igLikelihood, igCenterX, igCenterY, igHits, igEnergyFromCharge } # , igTotalCharge }

## `CurrentDsets` is the variable in use to read the correct data / extract the correct data
## from a given DF and turn it into a Torch Tensor.
## For now a global variable, that might likely become a field of `MLPDesc`.
var CurrentDsets* = ValidDsets.toSeq

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
  ## XXX: Add run number?
  #result["runNumber
  result["Type"] = $typ
  if not readRaw:
    result = result.mutate(f{"Target" ~ `energyFromCharge`.toRefDset})
  if subsetPerRun > 0:
    result = result.randomHead(subsetPerRun)

proc prepareData*(h5f: H5File, run: int, readRaw: bool, typ = dtBack, subsetPerRun = 0): DataFrame
proc readCalibData*(fname, calibType: string, eLow, eHigh: float,
                    subsetPerRun = 0): DataFrame =
  ## `subsetPerRun` is an integer which if given only returns this many entries (random)
  ## from each run.
  var h5f = H5open(fname, "r")
  let fileInfo = h5f.getFileInfo()

  var peakPos = newSeq[float]()
  result = newDataFrame()
  for run in fileInfo.runs:
    let xrayRefCuts = getXrayCleaningCuts()
    let runGrp = h5f[(recoBase() & $run).grp_str]
    let tfKind = if calibType == "photo": tfMnCr12 # same as 5.9 keV peak
                 else: tfMnCr12 #tfAgAg6 # ~3 keV line similar to escape
    let cut = xrayRefCuts[$tfKind]
    let grp = h5f[(recoBase() & $run / "chip_3").grp_str]
    let passIdx = cutOnProperties(
      h5f,
      grp,
      crSilver, # try cutting to silver
      (toDset(igRmsTransverse), cut.minRms, cut.maxRms),
      (toDset(igEccentricity), 0.0, cut.maxEccentricity),
      (toDset(igLength), 0.0, cut.maxLength),
      (toDset(igHits), cut.minPix, Inf),
      (toDset(igEnergyFromCharge), eLow, eHigh)
    )
    let dfChip = h5f.readRunDsets(run, chipDsets = some((chip: 3, dsets: @["eventNumber"])))
    let allEvNums = dfChip["eventNumber", int]
    let evNums = passIdx.mapIt(allEvNums[it]).toSet
    # filter to allowed events & remove any noisy events
    # note: do not use `subsetPerRun` here because we still need to extract valid event numbers!
    var df = prepareData(h5f, run, readRaw = false, typ = dtSignal)
    df = df.filter(f{int: `eventNumber` in evNums})
    df["runNumber"] = run
    if subsetPerRun > 0:
      df = df.randomHead(subsetPerRun)
    result.add df
  result["CalibType"] = calibType # Photo or escape peak

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
    dfLoc["Target"] = bin
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
  let cols = CurrentDsets.len
  var input = rawtensors.zeros(df.len * cols).reshape(sizes = [df.len.int64, cols].asTorchView())
  for i, c in CurrentDsets.mapIt(it.toDset(fkTpa)).sorted:
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

proc toTorchTensor*(df: DataFrame): RawTensor {.noInit.} =
  ## Converts an appropriate data frame to a 2D RawTensor of the data to be fed to the network
  let df = df.mutate(f{"totalCharge" ~ `totalCharge` / 1e7})
  let cols = CurrentDsets.len
  var input = rawtensors.zeros(df.len * cols).reshape(sizes = [df.len.int64, cols].asTorchView())
  for i, c in CurrentDsets.mapIt(it.toDset(fkTpa)).sorted:
    let xp = fromBlob[float](cast[pointer](df[c, float].unsafe_raw_offset()), df.len).convertRawTensor()
    input[_, i] = xp
  result = input

proc toNimSeq*[T](t: RawTensor): seq[T] =
  doAssert t.sizes().len == 1
  result = newSeq[T](t.size(0))
  for i in 0 ..< result.len:
    result[i] = t[i].item(T)
