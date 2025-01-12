import nimhdf5
import ingrid / [tos_helpers, ingrid_types]

proc main(fname: string, runType: RunTypeKind) =
  ## This is a convenience script if one has forgotten to set the run type
  ## in a `raw_data_manipulation` call to update the `runType` to the given
  ## target run type. This is possible for raw and reco HDF5 files.
  var h5f = H5open(fname, "rw")
  let fileInfo = h5f.getFileInfo()

  # 1. update the run type in the base group
  var baseGrp = h5f[($fileInfo.tpaFileKind).grp_str]
  baseGrp.attrs["runType"] = $runType
  for run in fileInfo.runs:
    # 2. update individual run
    var grp = h5f[(fileInfo.dataBase() & $run).grp_str]
    grp.attrs["runType"] = $runType

when isMainModule:
  import cligen
  dispatch main
