import nimhdf5, cligen, strformat, os, seqmath, stats, sequtils

proc main(fname: string, chip: int, run: int, dset: string,
          reco = false, likelihood = false, head = 0, verbose = false) =
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
    var dsetL = dset
    if head > 0:
      dsetL = dsetL[0 ..< head]
    if verbose:
      echo dsetL
    when typeof(dsetL[0]) is SomeNumber:
      let dsetF = dsetL.mapIt(it.float)
      echo "Min: ", dsetF.min
      echo "Max: ", dsetF.max
      echo "Mean: ", dsetF.mean
      echo "Sum: ", dsetF.sum
    elif typeof(dsetL[0]) is seq:
      when typeof(dsetL[0][0]) is SomeNumber:
        let dsetF = dsetL.flatten.mapIt(it.float)
        echo "Min: ", dsetF.min
        echo "Max: ", dsetF.max
        echo "Mean: ", dsetF.mean
        echo "Sum: ", dsetF.sum
    else:
      echo "No aggregates for type: ", typeof(dset)



when isMainModule:
  dispatch main
