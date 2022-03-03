import nimhdf5, cligen, strformat, os

proc main(fname: string, chip: int, run: int, dset: string,
          reco = false, likelihood = false, head = 0) =
  let h5f = H5open(fname, "r")
  var path = "/runs"
  if reco:
    path = "/reconstruction"
  if likelihood:
    path = "/likelihood"
  let grp = h5f[(&"{path}/run_{run}/chip_{chip}/").grp_str]
  let h5dset = h5f[(grp.name / dset).dset_str]
  echo "Dataset: ", dset
  withDset(h5dset):
    if head > 0:
      echo dset[0 ..< head]
    else:
      echo dset

when isMainModule:
  dispatch main
