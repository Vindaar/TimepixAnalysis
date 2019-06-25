## Just a playground to play around with `ggplotnim` in combination with TimepixAnalysis.
## Might be an important basis for features required in `ggplotnim` for my work.

import ingrid / tos_helpers
import sequtils, seqmath, nimhdf5, strutils, tables, persvector
import os, sugar

import arraymancer, plotly
import helpers / utils

import macros

template `[]=`(t: Tensor, stmt: untyped, val: untyped): untyped =
  t.apply(it => (if stmt: val else: it))

proc buildSeptemOccupancy(df: DataFrame) =
  ## builds an occupancy map of a single run. Data must be stored in `df`
  ## according to the following columns:
  ## chipNumber  eventNumber  X  Y  charge
  ## which gives us also a test for a *very* large dataframe.
  ## It creates a tensor of (3 * 256, 3 * 256) pixels, and maps each chip
  ## to the respective position when taking away the space between the chips.
  var xyz = 0
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
      echo "Chip subgroup "
      echo chip
      let chipNumber: range[0 .. 6] = chpPair[0][1].toInt
      #let fullFrame = getFullFrame(1.0)
      occ.addChipToSeptemEvent(chip, chipNumber)
      echo "Added to tensor"

    # after having built the occupancy, plot it
    # create lines around the tensor
    echo "Creating outline"
    var chipOutlineOcc = zeros[float]([3 * 256, 3 * 256])
    let maxVal = max(occ)
    for j in 0 ..< 7:
      let outlineDf = getChipOutline(maxVal)
      chipOutlineOcc.addChipToSeptemEvent(outlineDf, j)
    echo "Clamping outline"
    chipOutlineOcc.apply_inline:
        if x > 0.0:
          maxVal / 10.0
        else:
          x
    echo "Adding both tensors"
    occ = occ .+ chipOutlineOcc
    occ = occ.clamp(0.0, 100.0)#occ.toRawSeq.filterIt(it > 0.0).percentile(70))

    echo "Creating plot"
    heatmap(occ.toSeq2D)
      .title($pair)
      .width(1600)
      .height(1600)
      .show("event_" & $pair[0][1] & ".png")
    inc xyz
    if xyz > 10:
      quit()

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
