import nimhdf5

proc main(f: string, group: string) =
  withH5(f, "rw"):
    if group in h5f:
      echo "[INFO] Deleting group: ", group
      let res = h5f.delete(group)
      echo "[INFO] Deleting of group successful? ", res
    else:
      echo "[INFO] No group of name: ", group, " in file ", f

when isMainModule:
  import cligen
  dispatch main
