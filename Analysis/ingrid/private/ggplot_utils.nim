import macros, tables, strutils, os, sequtils, strformat, options

import ggplotnim
export ggplotnim

import hdf5_utils
import nimhdf5
import arraymancer

import ingrid / ingrid_types
import geometry, arraymancer_utils

macro echoType(x: typed): untyped =
  echo x.treeRepr

type SupportedRead = SomeFloat | SomeInteger | string | bool | Value


proc getDf*(h5f: H5File, path: string, keys: varargs[string]): DataFrame =
  ## read the datasets form `path` in the `h5f` file and combine them to a single
  ## DataFrame
  result = newDataFrame()
  var size = 0
  for key in keys:
    static:
      echoType(key)
    var dsetH5 = h5f[(path / key).dset_str]
    if size == 0:
      size = dsetH5.shape[0]
    withDset(dsetH5):
      when type(dset) is seq[SupportedRead]:
        result[key] = dset

proc readDsets*(h5f: H5FileObj, df: var DataFrame, names: seq[string], baseName: string) =
  ## reads the given datasets `names` from `baseName` into the existing `df`
  for name in names:
    let dsetName = baseName / name
    if dsetName in h5f:
      let dsetH5 = h5f[dsetName.dset_str]
      withDset(dsetH5):
        when type(dset) is seq[SupportedRead]:
          df[name] = dset
        else:
          doAssert false, "Invalid datatype for DataFrame! Dtype is " & $(type(dset))
    else:
      echo &"INFO: Run {baseName} does not have any data for dataset {name}"

proc readDsets*(h5f: H5FileObj, path = recoBase(),
                chipDsets = none[tuple[chip: int, dsets: seq[string]]](),
                commonDsets: openArray[string] = @[]): DataFrame =
  ## reads all desired datasets `chipDsets, commonDsets` in the given `h5f`
  ## file of `chip` under the given `path`. The result is returned as a
  ## `DataFrame`.
  ##
  ## `chipDsets` are datasets from the chip groups, whereas `commonDsets` are
  ## those from the run group (timestamp, FADC datasets etc)
  ## If input for both is given they are read as individual dataframes, which
  ## are then joined using the eventNumber dataset (which thus will always be
  ## read).
  let readChip = chipDsets.isSome
  var
    chipDsetNames: seq[string]
    chip: int
    commonDsets = @commonDsets
    evNumDset = "eventNumber"
  if readChip:
    chipDsetNames = chipDsets.get.dsets
    if evNumDset notin chipDsetNames and commonDsets.len > 0:
      chipDsetNames.add evNumDset
    if commonDsets.len > 0 and evNumDset notin commonDsets:
      commonDsets.add evNumDset
    chip = chipDsets.get.chip
  for run, grp in runs(h5f, path):
    var dfChip = newDataFrame()
    var dfAll = newDataFrame()
    if readChip:
      h5f.readDsets(dfChip, chipDsetNames, grp / "chip_" & $chip)
    h5f.readDsets(dfAll, commonDsets, grp)
    var df = if dfChip.len > 0 and dfAll.len > 0:
               innerJoin(dfChip, dfAll, evNumDset)
             elif dfChip.len > 0:
               dfChip
             else: #elif dfAll.len > 0:
               dfAll
    df["runNumber"] = constantColumn(run, df.len)
    if df.len > 0:
      result.add df

iterator getDataframes*(h5f: H5File): DataFrame =
  for num, group in runs(h5f):
    for grp in items(h5f, group):
      if "fadc" notin grp.name:
        let chipNum = grp.attrs["chipNumber", int]
        if chipNum == 3:
          let df = h5f.getDf(grp.name,
                             concat(@["eventNumber"],
                                    @(getFloatGeometryNames()),
                                    @(getIntClusterNames())))
          yield df

proc getSeptemDataFrame*(h5f: H5File, runNumber: int, allowedChips: seq[int] = @[]): DataFrame =
  ## Returns a subset data frame of the given `runNumber` and `chipNumber`, which
  ## contains only the zero suppressed event data
  var
    xs, ys: seq[uint8]
    chs: seq[float]
    evs: seq[int]
    chips = newColumn(colNone)
  let group = recoBase() & $runNumber
  for run, chip, groupName in chipGroups(h5f):
    if run == runNumber and (allowedChips.len == 0 or chip in allowedChips):
      let grp = h5f[groupName.grp_str]
      let chipNum = grp.attrs["chipNumber", int]
      echo "Reading chip ", chipNum, " of run ", runNumber
      let eventNumbersSingle = h5f[grp.name / "eventNumber", int64]
      var eventNumbers = newSeq[int]()
      let vlenXY = special_type(uint8)
      let vlenCh = special_type(float64)
      let
        x = h5f[grp.name / "x", vlenXY, uint8]
        y = h5f[grp.name / "y", vlenXY, uint8]
        ch = h5f[grp.name / "charge", vlenCh, float64]

      let
        xAll = flatten(x)
        yAll = flatten(y)
        chAll = flatten(ch)
        chipNumCol = constantColumn(chipNum, xAll.len)
      for i in 0 ..< x.len:
        # convert each event into a dataframe
        let event = eventNumbersSingle[i]
        let eventCol = toSeq(0 ..< x[i].len).mapIt(event.int)
        eventNumbers.add eventCol
      xs.add(xAll)
      ys.add(yAll)
      chs.add(chAll)
      evs.add(eventNumbers)
      chips = add(chips, chipNumCol)
  result = seqsToDf({"eventNumber" : evs, "x" : xs, "y" : ys, "charge" : chs,
                     "chipNumber" : chips})

proc getSeptemEventDF*(h5f: H5File, runNumber: int): DataFrame =
  ## Returns a DataFrame for the given `runNumber`, which contains two columns
  ## the event number and the chip number. This way we can easily extract
  ## which chips had activity on the same event.
  var
    evs, evIdx: seq[int]
    chips: Column
  let group = recoBase() & $runNumber
  for run, chip, groupName in chipGroups(h5f):
    if run == runNumber:
      echo "Reading chip ", chip, " of run ", runNumber
      let
        eventNumbers = h5f[groupName / "eventNumber", int]
        chipNumCol = constantColumn(chip, eventNumbers.len)
        evIndex = toSeq(0 ..< eventNumbers.len)
      evs.add(eventNumbers)
      chips = add(chips, chipNumCol)
      evIdx.add(evIndex)

  result = seqsToDf({"eventIndex" : evIdx, "eventNumber" : evs, "chipNumber" : chips})

iterator getSeptemDataFrame*(h5f: H5File): DataFrame =
  for num, group in runs(h5f):
    let df = h5f.getSeptemDataFrame(num)
    yield df

proc getChipOutline*(maxVal: SomeNumber): DataFrame =
  ## returns a data frame with only the outline of a Timepix chip as active pixels
  let zeroes = toSeq(0 ..< 256).mapIt(0)
  let maxvals = toSeq(0 ..< 256).mapIt(256)
  let incVals = toSeq(0 ..< 256)
  let xs = concat(zeroes, incVals, maxVals, incVals)
  let ys = concat(incVals, zeroes, incVals, maxvals)
  let ch = constantColumn(maxVal.float, xs.len)
  result = seqsToDf({"x" : xs, "y" : ys, "charge" : ch})

proc getSeptemOutlines*(maxVal: SomeNumber): Tensor[float] =
  ## returns the outline of the chips of the SeptemBoard in a
  ## full septem frame as a Tensor
  result = initSeptemFrame()
  for j in 0 ..< 7:
    let outlineDf = getChipOutline(maxVal)
    result.addChipToSeptemEvent(outlineDf, j)
  result.apply_inline:
      if x > 0.0:
        maxVal / 10.0
      else:
        x

proc getFullFrame*(maxVal: SomeNumber): DataFrame =
  ## returns a data frame with an event similar to a full timepix event, i.e. the
  ## pixels along the pad all full to 4096 pixels (after that cut off)
  let
    xData = toSeq(0 ..< 256)
    yData = toSeq(0 ..< 20)
    comb = product(@[xData, yData])
    ch = constantColumn(maxVal.float, 256 * 20)
  doAssert comb.len == ch.len
  let xy = comb.transpose
  result = seqsToDf({"x" : xy[0], "y": xy[1], "charge" : ch})

proc addChipToSeptemEvent*(occ: var Tensor[float], df: DataFrame, chipNumber: range[0 .. 6],
                           zDset = "charge") =
  doAssert zDset in df, "DataFrame has no key " & $zDset
  # now add values to correct places in tensor
  withSeptemXY(chipNumber):
    let xDf = df["x"].toTensor(int)
    let yDf = df["y"].toTensor(int)
    let zDf = df[zDset].toTensor(float)
    for i in 0 ..< df.len:
      var xIdx, yIdx: int64
      case chipNumber
      of 0, 1, 2, 3, 4:
        xIdx = x0 + xDf[i]
        yIdx = y0 + yDf[i]
      of 5, 6:
        xIdx = x0 - xDf[i]
        yIdx = y0 - yDf[i]
      occ[yIdx.int, xIdx.int] += zDf[i]
    #for i in 0 ..< chip.len:
    # instead of using the data frame, create fake data for now to test arrangment

proc dfToSeptemEvent*(df: DataFrame, zDset = "charge"): DataFrame =
  doAssert zDset in df, "DataFrame has no key " & $zDset
  # now add values to correct places in tensor
  let chipNum = df["chipNumber"].toTensor(int)
  var xDf = df["x"].toTensor(int)
  var yDf = df["y"].toTensor(int)
  let zDf = df[zDset].toTensor(float)
  for i in 0 ..< df.len:
    withSeptemXY(chipNum[i]):
      case chipNum[i]
      of 0, 1, 2, 3, 4:
        xDf[i] = x0 + xDf[i]
        yDf[i] = y0 + yDf[i]
      of 5, 6:
        xDf[i] = x0 - xDf[i]
        yDf[i] = y0 - yDf[i]
      else: doAssert false, "Invalid chip number!"
  result = seqsToDf({"x" : xDf, "y" : yDf, "charge" : zDf})
