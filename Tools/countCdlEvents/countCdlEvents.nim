import std / [strutils, os]
import ingrid / [tos_helpers, ingrid_types]
import nimhdf5, orgtables, datamancer

proc main(fname: string) =
  var h5f = H5open(fname, "r")
  defer: discard h5f.close()
  # (runNumber, tfKind, number of events raw, number of events filtered)
  var data = newSeq[tuple[runNumber: int, tfKind: string, numRaw: int, numFiltered: int]]()

  for num, group in runs(h5f):
    let grp = h5f[group.grp_str]
    let tfKind = grp.attrs["tfKind", string]
    # we get the number of raw events from the length of the `eventNumber` dataset
    # and the number of used events (after filtering) from the `CdlSpectrumIndices` length
    let numRaw  = h5f[(group / "chip_3" / "eventNumber").dset_str].shape[0] # 1D
    let numFilt = h5f[(group / "chip_3" / "CdlSpectrumIndices").dset_str].shape[0] # 1D
    data.add (runNumber: num, tfKind: tfKind, numRaw: numRaw, numFiltered: numFilt)

  echo "Table of raw runs:"
  echo data.toOrgTable
  var df = newDataFrame()
  for d in data: # slow, but convenient here
    df.add d
  var dfR = newDataFrame()
  for (tup, subDf) in groups(df.group_by("tfKind")):
    let numRaw = subDf["numRaw", int].sum
    let numFilt = subDf["numFiltered", int].sum
    dfR.add (tfKind: tup[0][1].toStr, numRaw: numRaw, numFiltered: numFilt)

  echo dfR.toOrgTable()

when isMainModule:
  import cligen
  dispatch main
