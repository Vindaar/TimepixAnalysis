import nimhdf5
import std / strutils

proc writeAttr(h5o: H5Dataset | H5Group, key, val: string, dtype: string) =
  template wr(fn): untyped {.dirty.} =
    let dval = fn val
    h5o.attrs[key] = dval
  case dtype
  of "int", "i", "d", "int64": wr(parseInt)
  of "float", "f", "float64":  wr(parseFloat)
  of "string", "str":          h5o.attrs[key] = val
  of "bool", "b":              wr(parseBool)
  else:
    doAssert false, "Unsupported data type : " & $dtype

proc main(f, obj, key, val: string, dtype: string) =
  ## Writes the attribute `key` with `val` to the HDF5 object `obj`
  ## in the file `f`
  withH5(f, "rw"):
    if obj in h5f:
      if h5f.isGroup(obj):
        var grp = h5f[obj.grp_str]
        echo "[INFO] Writing attribute to group: ", grp.name
        grp.writeAttr(key, val, dtype)
      else: # must be a dataset
        doAssert h5f.isDataset(obj)
        var dset = h5f[obj.dset_str]
        echo "[INFO] Writing attribute to dataset: ", dset.name
        dset.writeAttr(key, val, dtype)
    else:
      echo "[INFO] No group or dataset of name: ", obj, " in file ", f

when isMainModule:
  import cligen
  dispatch main
