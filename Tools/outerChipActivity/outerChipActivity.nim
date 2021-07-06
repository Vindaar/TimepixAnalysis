import ingrid / tos_helpers
import sequtils, nimhdf5, strutils, os, cligen, sets
import ggplotnim, seqmath, chroma, ginger

proc contains[T](t: Tensor[T], x: T): bool =
  for i in 0 ..< t.size:
    if x == t[i]:
      return true

proc read(h5f: H5File, chip: int): DataFrame =
  let cDsets = some((chip: chip,
                     dsets: @["energyFromCharge"]))
  let commonDsets = @["eventDuration"]
  result = h5f.readDsets(chipDsets = cDsets,
                         commonDsets = commonDsets)
  result["chip"] = chip

#proc melt(df: DataFrame, key: string)

proc main(fname: string) =
  var h5f = H5open(fname, "r")
  defer: discard h5f.close()
  var
    totDuration = 0.0
    numEvents = 0
    numMax = 0

  let dfNoChips = h5f.readDsets(commonDsets = @["eventDuration", "eventNumber"])
  for (tup, subdf) in groups(dfNoChips.group_by("runNumber")):
    inc numEvents, subDf["eventNumber", int].max

  var df = newDataFrame()
  for chip in 0 ..< 7:
    df.add h5f.read(chip)

  var
    numNoCenter = 0
    numOnlyCenter = 0
    numCenterAnd = 0
    numAny = 0
  for (tup, subDf) in groups(df.group_by(["eventNumber", "runNumber"])):
    let chips = subDf["chip"].unique.toTensor(int)
    if 3 notin chips:
      inc numNoCenter
    if [3].toTensor == chips:
      inc numOnlyCenter
    if 3 in chips and chips.len > 0:
      inc numCenterAnd
    inc numAny
  echo df

  echo "Number of total events: ".alignLeft(50), numEvents
  echo "Number of events without center: ".alignLeft(50), numNoCenter, " | ", (numNoCenter.float / numEvents.float) * 100.0, "%"
  echo "Number of events only center: ".alignLeft(50), numOnlyCenter, " | ", (numOnlyCenter.float / numEvents.float) * 100.0, "%"
  echo "Number of events with center activity and outer: ".alignLeft(50), numCenterAnd, " | ", (numCenterAnd.float / numEvents.float) * 100.0, "%"
  echo "Number of events any hit events: ".alignLeft(50), numAny, " | ", (numAny.float / numEvents.float) * 100.0, "%"
  echo "Mean of event durations: ".alignLeft(50), dfNoChips["eventDuration", float].mean

  ggplot(df, aes("eventDuration")) +
    geom_histogram(bins = 100) +
    ggsave("/tmp/event_duration_data_full.pdf")

  ggplot(df.filter(f{`eventDuration` < 2.25 and `eventDuration` > 0.0}), aes("eventDuration")) +
    geom_histogram(bins = 100) +
    ggsave("/tmp/event_duration_data_larger.pdf")

  ggplot(df.filter(f{`eventDuration` < 1e-8}), aes("eventDuration")) +
    geom_histogram(bins = 100) +
    ggsave("/tmp/event_duration_data_smaller.pdf")

  let df0 = df.filter(f{`eventDuration` < 1e-5 and `energyFromCharge` < 10.0})
  echo df0["runNumber"].unique
  ggplot(df0, aes("energyFromCharge", color = "chip", fill = "chip")) +
    geom_histogram(position = "identity", alpha = some(0.5), hdKind = hdOutline, bins = 100) +
    ggsave("/tmp/histo_energy_no_duration.pdf")

  ggplot(df0, aes("energyFromCharge", color = "chip", fill = "chip")) +
    facet_wrap("runNumber") +
    geom_histogram(position = "identity", alpha = some(0.5), hdKind = hdOutline, bins = 50) +
    ggsave("/tmp/histo_energy_no_duration_run_facet.pdf", width = 3000, height = 2000)


when isMainModule:
  dispatch main
