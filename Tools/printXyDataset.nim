import nimhdf5, cligen, strformat, os

proc main(fname: string, chip: int, run: int, dset: string,
          reco = false, head = 0) =
  let h5f = H5open(fname, "r")
  var path = "/runs"
  if reco:
    path = "/reconstruction"
  let grp = h5f[(&"{path}/run_{run}/chip_{chip}/").grp_str]
  let h5dset = h5f[(grp.name / dset).dset_str]
  template readit(typ: untyped): untyped =
    if head > 0:
      echo h5Dset.readVlen(typ)[0 ..< head]
    else:
      echo h5Dset.readVlen(typ)
  case h5dset.dtypeBaseKind
  of dkUint8: readit(uint8)
  of dkUint16: readit(uint16)
  of dkUint64: readit(uint64)
  of dkInt: readit(int)
  of dkFloat: readit(float)
  else:
    echo "Unsupported data type kind ", h5dset.dtype

when isMainModule:
  dispatch main
