import nimpy

proc toNimSeq*(xOb: PyObject, dtype: typedesc): seq[dtype] =
  ## converts the given object stored in `x` (this should be either a Python list
  ## or a Numpy array) to a Nim sequence
  ## Recursively calls itself until the given `dtype` is `SomeNumber`
  result = newSeqOfCap[dtype](100)
  for x in xOb:
    when dtype is SomeNumber | string:
      let xNim = x.to(dtype)
      result.add xNim
    elif dtype is seq:
      # use `getInnerType` macro to extract subtype of this sequence
      # defined in `nimhdf5` as well as `utils`
      type innerDtype = utils.getInnerType(dtype)
      result.add toNimSeq(x, innerDtype)
