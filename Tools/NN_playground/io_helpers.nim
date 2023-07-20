from std / random import shuffle
import std / [sequtils, os, strutils, options]

import ingrid / [tos_helpers, ingrid_types]
import pkg / [nimhdf5, datamancer]

const ValidReadDSets* = XrayReferenceDsets - { igNumClusters,
                                               igFractionInHalfRadius,
                                               igRadiusDivRmsTrans,
                                               igRadius, igBalance,
                                               igLengthDivRadius } + {
                                                 igCenterX, igCenterY, igGasGain }

# include total charge
#const ValidDsets* = { igEccentricity, igLengthDivRmsTrans, igFractionInTransverseRms }
const ValidDsets* = ValidReadDSets - { igLikelihood, igCenterX, igCenterY, igHits, igEnergyFromCharge } # , igTotalCharge }

## `CurrentDsets` is the variable in use to read the correct data / extract the correct data
## from a given DF and turn it into a Torch Tensor.
## For now a global variable, that might likely become a field of `MLPDesc`.
var CurrentDsets* = ValidDsets.toSeq.mapIt(it.toDset())

type
  DataType* = enum
    dtSignal = "signal"
    dtBack = "back"

proc shuff*(x: seq[int]): seq[int] =
  ## Not in place shuffle
  result = x
  result.shuffle()

import ../../Tools/determineDiffusion/determineDiffusion
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
  doAssert false, "This is currently not sane to use! Diffusion missing and mixing of runs!"
  #result["σT"] = getDiffusion(result) / 1000.0 # convert from μm/√cm to mm/√cm

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
                     subsetPerRun = 0,
                     validDsets: set[InGridDsetKind] = ValidReadDsets - { igLikelihood },
                     filterNan = true): DataFrame =
  ## Reads all data for the given run `path` (must be a chip path)
  ##
  ## `subsetPerRun` is an integer which if given only returns this many entries (random)
  ## from each run.
  # get run number from parent group. This group is a chip group
  let runNumber = h5f[path.parentDir.grp_str].attrs["runNumber", int]
  let dsets = toSeq(validDsets)
  let evNumDset = igEventNumber.toDset(fkTpa)
  let grp = h5f[path.grp_str]
  result = newDataFrame()
  if not readRaw:
    let data = readInGridDsetKind(h5f, grp.name, dsets)
    for d in dsets: # these are all now filled
      result[d.toDset()] = data[d]
    result[evNumDset] = h5f.readAs(grp.name / evNumDset, int)
  else:
    result = h5f.readRaw(grp.name)
  if filterNan:
    for d in dsets:
      let dset = d.toDset()
      result = result.filter(f{float: classify(idx(dset)) in {fcNormal, fcSubnormal, fcZero, fcNegZero}})
  ## XXX: Add run number?
  #result["runNumber
  result["Type"] = $typ
  if not readRaw:
    result = result.mutate(f{"Target" ~ `energyFromCharge`.toRefDset})
  if subsetPerRun > 0:
    result = result.randomHead(subsetPerRun)
  #result["Idx"] = toSeq(0 ..< result.len)
  # here we only use the cached diffusion values!
  result["σT"] = getDiffusionForRun(run = runNumber, isBackground = (typ == dtBack)) / 1000.0

proc prepareData*(h5f: H5File, run: int, readRaw: bool,
                  typ = dtBack,
                  subsetPerRun = 0,
                  validDsets: set[InGridDsetKind] = ValidReadDsets - { igLikelihood } ): DataFrame
proc readCalibData*(fname, calibType: string, eLow, eHigh: float,
                    subsetPerRun = 0,
                    tfKind = none[TargetFilterKind](),
                    validDsets: set[InGridDsetKind] = ValidReadDsets - { igLikelihood } ): DataFrame =
  ## `subsetPerRun` is an integer which if given only returns this many entries (random)
  ## from each run.
  var h5f = H5open(fname, "r")
  let fileInfo = h5f.getFileInfo()

  var peakPos = newSeq[float]()
  result = newDataFrame()
  for run in fileInfo.runs:
    let xrayRefCuts = getXrayCleaningCuts()
    let runGrp = h5f[(recoBase() & $run).grp_str]
    if tfKind.isSome:
      # try to read attribute and if not given `tfKind` skip this run
      let tf = tfKind.get
      if "tfKind" in runGrp.attrs:
        let tfKindFile = parseEnum[TargetFilterKind](runGrp.attrs["tfKind", string])
        if tf != tfKindFile: continue
      else:
        echo "[WARNING]: You asked for the target filter kind: ", tf, " but the input file does not have " &
          "this attribute in ", runGrp.name

    ## XXX: verify we really want `tfMnCr12` for escape peak and not AgAg6kV!
    let cutTfKind =
      if tfKind.isSome: tfKind.get
      elif calibType in ["photo", "5.9"]: tfMnCr12 # same as 5.9 keV peak
      else: tfAgAg6 #tfMnCr12 #tfAgAg6 # ~3 keV line similar to escape
    when false:
      let cut = xrayRefCuts[$cutTfKind]
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
    var df = prepareData(h5f, run, readRaw = false, typ = dtSignal,
                         validDsets = validDsets)
      .cutXrayCleaning(cutTfKind,
                       eLow, eHigh)
    #echo "DF IS LONG: ", df.len, " and ALLEVNUMS LONG: ", allEvNums.len
    # df = df.filter(f{int: `eventNumber` in evNums})
    df["runNumber"] = run
    if subsetPerRun > 0:
      df = df.randomHead(subsetPerRun)
    if df.filter(f{`energyFromCharge` < eLow or `energyFromCharge` > eHigh}).len > 0:
      let dfFF = df.filter(f{`energyFromCharge` < eLow or `energyFromCharge` > eHigh})
      #echo
      echo "eLow : ", eLow, " and ", eHigh
      echo "What the fuck?"

      echo df["Idx"].len
      #echo passIdx.len

      #echo "And pass idx ? ", passIdx
      echo dfFF["Idx"]



      quit()

    result.add df
  result["CalibType"] = calibType # Photo or escape peak

proc prepareCdl*(readRaw: bool,
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

proc prepareData*(h5f: H5File, run: int, readRaw: bool, typ = dtBack, subsetPerRun = 0,
                  validDsets: set[InGridDsetKind] = ValidReadDsets - { igLikelihood }
                 ): DataFrame =
  let path = "/reconstruction/run_$#/chip_3" % $run
  result = h5f.readValidDsets(path, readRaw, typ, subsetPerRun, validDsets)

proc prepareBackground*(fname: string, run: int, readRaw: bool, subsetPerRun = 0,
                        validDsets: set[InGridDsetKind] = ValidReadDsets - { igLikelihood }): DataFrame =
  var h5f = H5open(fname, "r")
  result = h5f.prepareData(run, readRaw, subsetPerRun = subsetPerRun,
                           validDsets = validDsets)
  discard h5f.close()

proc prepareAllBackground*(fname: string, readRaw: bool, subsetPerRun = 0,
                           validDsets: set[InGridDsetKind] = ValidReadDsets - { igLikelihood },
                           region = crAll
                          ): DataFrame =
  var h5f = H5open(fname, "r")
  for run, grp in runs(h5f):
    var df = h5f.prepareData(run, readRaw, subsetPerRun = subsetPerRun,
                             validDsets = validDsets)
    df["runNumber"] = run
    result.add df

  # filter to desired region
  if region != crAll:
    result = result.filter(f{float -> bool: inRegion(`centerX`, `centerY`, region)})
  discard h5f.close()

when defined(cpp):
  {.experimental: "views".}
  import flambeau / [flambeau_raw]
  proc toInputTensor*(df: DataFrame): (RawTensor, RawTensor) {.noInit.} =
    ## Converts an appropriate data frame to a tuple of input / target tensors
    let df = df.mutate(f{"totalCharge" ~ `totalCharge` / 1e7})
    let cols = CurrentDsets.len
    var input = rawtensors.zeros(df.len * cols).reshape(sizes = [df.len.int64, cols].asTorchView())
    for i, c in CurrentDsets.sorted:
      let xp = fromBlob(cast[pointer](df[c, float].unsafe_raw_offset()), df.len, float)
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
    for i, c in CurrentDsets.sorted:
      let xp = fromBlob(cast[pointer](df[c, float].unsafe_raw_offset()), df.len, float)
      input[_, i] = xp
    result = input

  proc toNimSeq*[T](t: RawTensor): seq[T] =
    doAssert t.sizes().len == 1
    result = newSeq[T](t.size(0))
    for i in 0 ..< result.len:
      result[i] = t[i].item(T)
