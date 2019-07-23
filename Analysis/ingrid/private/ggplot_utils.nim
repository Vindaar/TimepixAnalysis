import macros, tables, strutils, os, sequtils

import ggplotnim
export ggplotnim

import hdf5_utils
import nimhdf5
import arraymancer

import ingrid / ingrid_types

macro echoType(x: typed): untyped =
  echo x.treeRepr

proc getDf*(h5f: var H5FileObj, path: string, keys: varargs[string]): DataFrame =
  ## read the datasets form `path` in the `h5f` file and combine them to a single
  ## DataFrame
  var data: OrderedTable[string, seq[Value]]
  var size = 0
  for key in keys:
    static:
      echoType(key)
    var dsetH5 = h5f[(path / key).dset_str]
    if size == 0:
      size = dsetH5.shape[0]
    withDset(dsetH5):
      data[key] = %~ dset
  result = toDf(data)
  result.len = size

iterator getDataframes*(h5f: var H5FileObj): DataFrame =
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

proc getSeptemDataFrame*(h5f: var H5FileObj, runNumber: int): DataFrame =
  ## Returns a subset data frame of the given `runNumber` and `chipNumber`, which
  ## contains only the zero suppressed event data
  var
    xs, ys, chs, evs, chips: seq[Value]
  let group = recoBase() & $runNumber
  for grp in items(h5f, group):
    if "fadc" notin grp.name:
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
        chipNumCol = toSeq(0 ..< xAll.len).mapIt(chipNum)
      for i in 0 ..< x.len:
        # convert each event into a dataframe
        let event = eventNumbersSingle[i]
        let eventCol = toSeq(0 ..< x[i].len).mapIt(event.int)
        eventNumbers.add eventCol

      xs.add(%~ xAll)
      ys.add(%~ yAll)
      chs.add(%~ chAll)
      evs.add(%~ eventNumbers)
      chips.add(%~ chipNumCol)

  result = seqsToDf({"eventNumber" : evs, "x" : xs, "y" : ys, "charge" : chs,
                     "chipNumber" : chips})

proc getSeptemEventDF*(h5f: var H5FileObj, runNumber: int): DataFrame =
  ## Returns a DataFrame for the given `runNumber`, which contains two columns
  ## the event number and the chip number. This way we can easily extract
  ## which chips had activity on the same event.
  var
    evs, chips, evIdx: seq[Value]
  let group = recoBase() & $runNumber
  for grp in items(h5f, group):
    if "fadc" notin grp.name:
      let chipNum = grp.attrs["chipNumber", int]
      echo "Reading chip ", chipNum, " of run ", runNumber
      let
        eventNumbers = h5f[grp.name / "eventNumber", int64]
        chipNumCol = toSeq(0 ..< eventNumbers.len).mapIt(chipNum)
        evIndex = toSeq(0 ..< eventNumbers.len)
      evs.add(%~ eventNumbers)
      chips.add(%~ chipNumCol)
      evIdx.add(%~ evIndex)

  result = seqsToDf({"eventIndex" : evIdx, "eventNumber" : evs, "chipNumber" : chips})

iterator getSeptemDataFrame*(h5f: var H5FileObj): DataFrame =
  for num, group in runs(h5f):
    let df = h5f.getSeptemDataFrame(num.parseInt)
    yield df

proc getChipOutline*(maxVal: SomeNumber): DataFrame =
  ## returns a data frame with only the outline of a Timepix chip as active pixels
  let zeroes = toSeq(0 ..< 256).mapIt(0)
  let maxvals = toSeq(0 ..< 256).mapIt(256)
  let incVals = toSeq(0 ..< 256)
  let xs = concat(zeroes, incVals, maxVals, incVals)
  let ys = concat(incVals, zeroes, incVals, maxvals)
  let ch = toSeq(0 ..< xs.len).mapIt(maxVal.float)
  result = seqsToDf({"x" : xs, "y" : ys, "charge" : ch})

proc initSeptemFrame*(): Tensor[float] {.noinit.} =
  result = zeros[float]([3 * 256, 3 * 256])

proc getFullFrame*(maxVal: SomeNumber): DataFrame =
  ## returns a data frame with an event similar to a full timepix event, i.e. the
  ## pixels along the pad all full to 4096 pixels (after that cut off)
  let
    xData = toSeq(0 ..< 256)
    yData = toSeq(0 ..< 20)
    comb = product(@[xData, yData])
    ch = toSeq(0 ..< 256 * 20).mapIt(maxVal.float)
  doAssert comb.len == ch.len
  let xy = comb.transpose
  result = seqsToDf({"x" : xy[0], "y": xy[1], "charge" : ch})

template withSeptemXY*(chipNumber: int, actions: untyped): untyped =
  ## injects the x0, y0 coordinates of the given chip number embedded into
  ## the septem frame
  var
    x0 {.inject.}: int
    y0 {.inject.}: int
  case chipNumber
  of 0:
    # chip bottom left of board
    y0 = 0 # top end of bottom row
    x0 = 128 # shifted by half a chip to the right
  of 1:
    # chip bottom right
    y0 = 0
    x0 = 128 + 256
  of 2:
    # middle left
    y0 = 256
    x0 = 0
  of 3:
    # middle middle
    y0 = 256
    x0 = 256
  of 4:
    # middle right
    y0 = 256
    x0 = 2 * 256
  of 5:
    # top right chip
    y0 = 3 * 256
    x0 = 2 * 256 + 128
  of 6:
    # top left chip
    y0 = 3 * 256
    x0 = 128 + 256
  actions

func determineChip*[T:  SomePix](p: T): int =
  ## determines which chip the given septem pixel coordinate corresponds to
  if p.y in 0 .. 255 and p.x in 128 .. 128 + 255:
    # bottom left
    result = 0
  elif p.y in 0 .. 255:
    # bottom right
    result = 1
  elif p.y in 256 .. 511 and p.x in 0 .. 255:
    # middle left
    result = 2
  elif p.y in 256 .. 511 and p.x in 256 .. 511:
    # center
    result = 3
  elif p.y in 256 .. 511:
    # middle right
    result = 4
  elif p.x in 128 + 256 .. 512 + 127:
    # top right
    result = 5
  elif p.x in 128 .. 128 + 255:
    # top left
    result = 6
  else:
    raise newException(Exception, "This chip should not exist! " & $p)

func septemPixToChpPix*[T: SomePix](p: T, chipNumber: range[0 .. 6]): T =
  ## inverse of chpPixToSeptemPix
  result = p
  case chipNumber
  of 0:
    # chip bottom left of board
    result.x = result.x - 128 # shifted by half a chip to the right
  of 1:
    # chip bottom right
    result.x = result.x - (128 + 256)
  of 2:
    # middle left
    result.y = result.y - 256
  of 3:
    # middle middle
    result.y = result.y - 256
    result.x = result.x - 256
  of 4:
    # middle right
    result.y = result.y - 256
    result.x = result.x - 512
  of 5:
    # top right chip
    result.y = -(result.y - 3 * 256)
    result.x = -(result.x - (2 * 256 + 128))
  of 6:
    # top left chip
    result.y = -(result.y - 3 * 256)
    result.x = -(result.x - (128 + 256))

func chpPixToSeptemPix*(p: Pix, chipNumber: range[0 .. 6]): PixInt =
  ## converts the given local chip pixel to the full septem frame coordinate system
  withSeptemXY(chipNumber):
    var xIdx, yIdx: int
    case chipNumber
    of 0, 1, 2, 3, 4:
      xIdx = x0 + p.x.int
      yIdx = y0 + p.y.int
    of 5, 6:
      xIdx = x0 - p.x.int
      yIdx = y0 - p.y.int
    result = (x: xIdx, y: yIdx, ch: p.ch.int)

func chpPixToSeptemPix*(pix: Pixels, chipNumber: range[0 .. 6]): PixelsInt =
  ## converts the given local chip pixels to the full septem frame coordinate system
  result.setLen(pix.len)
  for i, p in pix:
    let pp = chpPixToSeptemPix(p, chipNumber)
    result[i] = pp

proc addChipToSeptemEvent*(occ: var Tensor[float], df: DataFrame, chipNumber: range[0 .. 6],
                           zDset = "charge") =
  doAssert zDset in df, "DataFrame has no key " & $zDset
  # now add values to correct places in tensor
  withSeptemXY(chipNumber):
    for i in 0 ..< df.len:
      var xIdx, yIdx: int64
      case chipNumber
      of 0, 1, 2, 3, 4:
        xIdx = x0 + df["x"][i].toInt
        yIdx = y0 + df["y"][i].toInt
      of 5, 6:
        xIdx = x0 - df["x"][i].toInt
        yIdx = y0 - df["y"][i].toInt
      occ[yIdx.int, xIdx.int] += df[zDset][i].toFloat
    #for i in 0 ..< chip.len:
    # instead of using the data frame, create fake data for now to test arrangment
