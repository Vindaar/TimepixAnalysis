import nimhdf5
import ingrid / [tos_helpers, ingrid_types]

proc main(f: string) =
  ## Updates the run number `runNumber` attribute in the /runs or /reconstruction
  ## group attributes to match the run number of the group name.
  withH5(f, "rw"):
    let fileInfo = h5f.getFileInfo()
    for r in fileInfo.runs:
      var grp = h5f[(dataBase(fileInfo) & $r).grp_str]
      echo "[INFO] Updating run number: ", grp.attrs["runNumber", int], " to ", r
      grp.attrs["runNumber"] = r

when isMainModule:
  import cligen
  dispatch main
