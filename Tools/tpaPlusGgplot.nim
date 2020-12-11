## Just a playground to play around with `ggplotnim` in combination with TimepixAnalysis.
## Might be an important basis for features required in `ggplotnim` for my work.

import ingrid / tos_helpers
import sequtils, seqmath, nimhdf5, strutils, tables, strformat, times
#from json import JsonNode, `[]`, `[]=`, `%`
import os, sugar

import arraymancer, plotly
import plotly / color
import helpers / utils

import macros

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
    echo "chips in subgroup ", subgroup["chipNumber"].toTensor(int).toRawSeq.toSet
    echo "Number of total pixels ", subgroup.len
    # each subgroup corresponds to a full septemboard event
    # Construct a tensor covering the whole board from them
    var occ = initSeptemFrame()
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
    var chipOutlineOcc = initSeptemFrame()
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
    when false:
      var plt = heatmap(occ.toSeq2D)
        .title($pair)
        .width(1600)
        .height(1600)
      let fname = "event_" & $pair[0][1]
      template createPlot(cmap: untyped): untyped =
        plt = plt.zmax(6)
          .colormap(cmap)
        plt.show(fname & "_" & astToStr(cmap) & ".svg")
      createPlot(Viridis)
      createPlot(Plasma)
      createPlot(WhiteToBlack)
    inc xyz
    if xyz > 10:
      quit()

proc main(fname: string) =
  var h5f = H5open(fname, "r")
  defer: discard h5f.close()

  #for df in getSeptemDataFrame(h5f):
  #  buildSeptemOccupancy(df)
  #  if true:
  #    quit()

  let df = getSeptemDataFrame(h5f, 107)
  echo df

  let outDir = "/tmp" / h5f.attrs[PlotDirPrefixAttr, string]
  createDir(outDir)
  echo "Outdir is ", outDir
  let t0 = epochTime()
  var count = 0
  #df.write_csv("/tmp/df_run_106.csv")
  for tup, dfEv in groups(df.group_by("eventNumber")):
    let evNum = tup[0][1].toInt
    let dfSeptem = dfToSeptemEvent(dfEv)
    if dfSeptem.len > 3:
      ggplot(dfSeptem, aes("x", "y", fill = "charge")) +
        geom_point() +
        ggtitle("Event number " & $evNum) +
        xlim(0, 3*256) + ylim(0, 3*256) +
        ggsave(outDir / "septem_events_run107_" & $evNum & ".pdf")
    if dfEv.len > 3:
      ggplot(dfEv, aes("x", "y", fill = "charge")) +
        facet_wrap("chipNumber") +
        geom_point() +
        ggtitle("Event number " & $evNum) +
        xlim(0, 256) + ylim(0, 256) +
        ggsave(outDir / "events_run107_" & $evNum & ".pdf")
      inc count
  let t1 = epochTime()
  echo "Created ", count, " plots in ", t1 - t0, " s. ", count.float / (t1 - t0).float, " fps"

  #for df in getDataframes(h5f):
  #  echo df
  #  let fname = "figs/"
  #  echo "Saving as ", fname
  #
  #  # small sizes
  #  let dfSmall = df.filter(f{float -> bool: `length` < 1.0})
  #  let dfOne = df.mutate(f{float -> bool: "smaller" ~ `length` < 1.0})
  #   #.group_by("smaller")
  #  echo dfOne
  #  ggplot(df, aes("length", "eccentricity")) +
  #    geom_point() +
  #    ggsave(fname & "le.svg")
  #
  #  ggplot(dfOne, aes("length", "hits", color = "smaller")) +
  #    geom_point() +
  #    ggsave(fname & "length.svg")
  #
  #  echo (dfOne.filter(f{float -> bool: `hits` < 300.0}).arrange("hits", order = SortOrder.Descending))
  #  echo f{"eccCut<1.5" ~ `eccentricity` < 1.5 and `length` < 6.0}
  #  echo dfOne.mutate(f{"eccCut<1.5" ~ `eccentricity` < 1.5 and `length` < 6.0})
  #
  #  echo dfOne.mutate(f{"between50and100" ~ `hits` < 100.0 and `hits` > 50.0})
  #  echo dfOne.filter(f{`hits` < 100.0 xor `hits` > 50.0})
  #
  #  dfOne.mutate(f{"eccCut<1.5" ~ `eccentricity` < 1.5 and `length` < 6.0})
  #               #f{"length<6" ~ "length" < 6.0})
  #    .filter(f{float -> bool: `hits` < 300.0})
  #    .ggplot(aes("hits", fill = "eccCut<1.5")) +
  #    geom_histogram() +
  #    ggsave(fname & "hitsSmall.svg")


when isMainModule:
  if paramCount() > 0:
    main(paramStr(1))
  else:
    echo "Hand a filename to be read from!"
