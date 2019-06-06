import ingrid / tos_helpers
import sequtils, seqmath, nimhdf5, strutils, tables, persvector
import os
import ggplotnim

proc main(fname: string) =
  var h5f = H5file(fname, "r")
  defer: discard h5f.close()
  var
    totDuration = 0.0
    numEvents = 0
    numMax = 0
  for num, group in runs(h5f):
    let allEvs = h5f[group / "eventNumber", int64]
    let duration = h5f[group / "eventDuration", float64]
    #discard toTab([allEvs, duration])
    let dfGrp = seqsToDf({"eventNumber" : allEvs, "duration" : duration})
    #echo dfGrp["duration"][0 ..< 100].mapIt(it.fnum)
    for grp in items(h5f, group):
      if "fadc" notin grp.name:
        let chipNum = grp.attrs["chipNumber", int]
        if chipNum == 3:
          let evs = h5f[grp.name / "eventNumber", int64]
          let hits = h5f[grp.name / "hits", int64]
          let dfChip = seqsToDf({"eventNumber" : evs, "hits" : hits})
          let dfJoined = innerJoin(dfGrp, dfChip, by = "eventNumber")
          let maxHits = dfJoined.filter(f{"hits" > 4095})
          totDuration += duration.sum
          numEvents += allEvs.len
          numMax += maxHits.len
  echo "Total duration in seconds = ", totDuration
  echo "Number of total frames = ", numEvents
  echo "Number of full frames = ", numMax
  echo "Ratio of full frames = ", numMax / numEvents


when isMainModule:
  if paramCount() > 0:
    main(paramStr(1))
  else:
    echo "Hand a filename to be read from!"
