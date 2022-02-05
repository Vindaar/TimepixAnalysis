import nimhdf5, cligen, strformat, os

proc main(fname: string, chip: int, run: int) =
  let h5f = H5open(fname, "r")

  let grp = h5f[(&"/runs/run_{run}/chip_{chip}/").grp_str]
  echo h5f[grp.name / "raw_x", seq[uint8]]

when isMainModule:
  dispatch main
