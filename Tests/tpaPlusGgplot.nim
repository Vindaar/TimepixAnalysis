## Just a playground to play around with `ggplotnim` in combination with TimepixAnalysis.
## Might be an important basis for features required in `ggplotnim` for my work.

import ingrid / tos_helpers
import sequtils, seqmath, nimhdf5, strutils, tables, persvector
import os, sugar
import ggplotnim except `%`

import arraymancer, plotly
import helpers / utils

import macros

macro echoType(x: typed): untyped =
  echo x.treeRepr

proc getDf(h5f: var H5FileObj, path: string, keys: varargs[string]): DataFrame =
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
      data[key] = ggplotnim.`%`(dset)
  result = toDf(data)
  result.len = size

iterator getDataframes(h5f: var H5FileObj): DataFrame =
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

iterator getSeptemDataFrame(h5f: var H5FileObj): DataFrame =
  for num, group in runs(h5f):
    var df: DataFrame
    var
      xs, ys, chs, evs, chips: seq[Value]
    for grp in items(h5f, group):
      if "fadc" notin grp.name:
        let chipNum = grp.attrs["chipNumber", int]
        echo "Reading chip ", chipNum, " of run ", num
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

        xs.add(ggplotnim.`%` xAll)
        ys.add(ggplotnim.`%` yAll)
        chs.add(ggplotnim.`%` chAll)
        evs.add(ggplotnim.`%` eventNumbers)
        chips.add(ggplotnim.`%` chipNumCol)

    df = seqsToDf({"eventNumber" : evs, "x" : xs, "y" : ys, "charge" : chs,
                    "chipNumber" : chips})
    yield df

proc getChipOutline(maxVal: SomeNumber): DataFrame =
  ## returns a data frame with only the outline of a Timepix chip as active pixels
  let zeroes = toSeq(0 ..< 256).mapIt(0)
  let maxvals = toSeq(0 ..< 256).mapIt(256)
  let incVals = toSeq(0 ..< 256)
  let xs = concat(zeroes, incVals, maxVals, incVals)
  let ys = concat(incVals, zeroes, incVals, maxvals)
  let ch = toSeq(0 ..< xs.len).mapIt(maxVal.float)
  result = seqsToDf({"x" : xs, "y" : ys, "charge" : ch})

proc getFullFrame(maxVal: SomeNumber): DataFrame =
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

proc addChipToSeptemOccupancy(occ: var Tensor[float], df: DataFrame, chipNumber: range[0 .. 6]) =
  var
    x0, y0: int
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
  # now add values to correct places in tensor
  for i in 0 ..< df.len:
    var xIdx, yIdx: int64
    case chipNumber
    of 0, 1, 2, 3, 4:
      xIdx = x0 + df["x"][i].toInt
      yIdx = y0 + df["y"][i].toInt
    of 5, 6:
      xIdx = x0 - df["x"][i].toInt
      yIdx = y0 - df["y"][i].toInt
    occ[yIdx.int, xIdx.int] += df["charge"][i].toFloat
  #for i in 0 ..< chip.len:
  # instead of using the data frame, create fake data for now to test arrangment


template `[]=`(t: Tensor, stmt: untyped, val: untyped): untyped =
  t.apply(it => (if stmt: val else: it))

proc buildSeptemOccupancy(df: DataFrame) =
  ## builds an occupancy map of a single run. Data must be stored in `df`
  ## according to the following columns:
  ## chipNumber  eventNumber  X  Y  charge
  ## which gives us also a test for a *very* large dataframe.
  ## It creates a tensor of (3 * 256, 3 * 256) pixels, and maps each chip
  ## to the respective position when taking away the space between the chips.
  let dfPerEvent = df.group_by("eventNumber")
  for (pair, subgroup) in groups(dfPerEvent):
    echo "Pair ", pair
    echo "chips in subgroup ", subgroup["chipNumber"].toSeq.toSet
    echo "Number of total pixels ", subgroup.len
    # each subgroup corresponds to a full septemboard event
    # Construct a tensor covering the whole board from them
    var occ = zeros[float]([3 * 256, 3 * 256])
    # now create another subgroup of the current by the chip number
    let chips = subgroup.group_by("chipNumber")
    #for j in 0 ..< 7:
    echo subgroup
    for (chpPair, chip) in groups(chips):
      echo "In subgroup for chip ", chpPair
      let chipNumber: range[0 .. 6] = chpPair[0][1].toInt
      #let fullFrame = getFullFrame(1.0)
      occ.addChipToSeptemOccupancy(chip, chipNumber)
      echo "Added to tensor"

    # after having built the occupancy, plot it
    # create lines around the tensor
    echo "Creating outline"
    var chipOutlineOcc = zeros[float]([3 * 256, 3 * 256])
    let maxVal = max(occ)
    for j in 0 ..< 7:
      let outlineDf = getChipOutline(maxVal)
      chipOutlineOcc.addChipToSeptemOccupancy(outlineDf, j)
    echo "Clamping outline"
    #chipOutlineOcc[it > 0.0] = maxVal
    chipOutlineOcc.apply_inline:
        if x > 0.0:
          maxVal / 5.0
        else:
          x
    echo "Adding both tensors"
    occ = occ .+ chipOutlineOcc
    echo "Creating plot"
    heatmap(occ.toSeq2D).show()

proc main(fname: string) =
  var h5f = H5file(fname, "r")
  defer: discard h5f.close()

  for df in getSeptemDataFrame(h5f):
    buildSeptemOccupancy(df)
    if true:
      quit()

  for df in getDataframes(h5f):
    echo df
    let fname = "figs/" & $12
    echo "Saving as ", fname

    # small sizes
    let dfSmall = df.filter(f{"length" < 1.0})
    let dfOne = df.mutate(f{"smaller" ~ "length" < 1.0})
      #.group_by("smaller")
    echo dfOne
    ggplot(df, aes("length", "eccentricity")) +
      geom_point() +
      ggsave(fname & "le.svg")

    ggplot(dfOne, aes("length", "hits", color = "smaller")) +
      geom_point() +
      ggsave(fname & "length.svg")

    echo (dfOne.filter(f{"hits" < 300.0}).arrange("hits", order = SortOrder.Descending))
    echo f{"eccCut<1.5" ~ "eccentricity" < 1.5 and "length" < 6.0}
    echo dfOne.mutate(f{"eccCut<1.5" ~ "eccentricity" < 1.5 and "length" < 6.0})

    echo dfOne.mutate(f{"between50and100" ~ "hits" < 100.0 and "hits" > 50.0})
    echo dfOne.filter(f{"hits" < 100.0 xor "hits" > 50.0})

    dfOne.mutate(f{"eccCut<1.5" ~ "eccentricity" < 1.5 and "length" < 6.0})
                 #f{"length<6" ~ "length" < 6.0})
      .filter(f{"hits" < 300.0})
      .ggplot(aes("hits", fill = "eccCut<1.5")) +
      geom_histogram() +
      ggsave(fname & "hitsSmall.svg")


when isMainModule:
  if paramCount() > 0:
    main(paramStr(1))
  else:
    echo "Hand a filename to be read from!"
