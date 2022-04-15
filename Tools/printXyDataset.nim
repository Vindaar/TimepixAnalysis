import nimhdf5, cligen, seqmath
import std / [strformat, os, stats, sequtils, sets, algorithm, strutils]

type
  ## A generic cut on input data using dset & low / high values
  GenericCut* = tuple[dset: string, lower, upper: float]

proc readIndices(h5f: H5File, path, dset: string, cuts: seq[GenericCut]): seq[int] =
  ## get the correct indices we need to read for the data `dset` such that it matches
  ## the given cuts. Cuts are only allowed from datasets at the same level as `path`
  ## and must be readable as `float`.
  # 1. get `dset`
  let dset = h5f[(path / dset).dset_str]
  # 2. get number of elements from shape of dset
  let num = dset.shape[0]
  # 3. create full set of indices from shape (i.e. all data)
  var idxs = toSeq(0 ..< num).toSet
  # 4. apply filters for each cut
  for cut in cuts:
    # 4.1. read data for this cut
    let data = h5f[(path / cut.dset), float] # read as float
    for i in toSeq(idxs):
      let x = data[i]
      if x < cut.lower or x > cut.upper:
        idxs.excl i
  result = toSeq(idxs).sorted

proc main(fname: string, chip: int, run: int, dset: string,
          reco = false, likelihood = false, head = 0, verbose = false,
          cuts: seq[GenericCut] = @[]) =
  let h5f = H5open(fname, "r")
  var path = "/runs"
  if reco:
    path = "/reconstruction"
  if likelihood:
    path = "/likelihood"
  let grpPath = &"{path}/run_{run}/chip_{chip}/"
  let h5dset = h5f[(grpPath / dset).dset_str]
  let idxs = readIndices(h5f, grpPath, dset, cuts)
  echo "Dataset: ", dset
  withDset(h5dset):
    var dsetL = idxs.mapIt(dset[it])
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
  import cligen/argcvt
  proc argParse(dst: var GenericCut, dfl: GenericCut,
                a: var ArgcvtParams): bool =
    echo "Parsing ", a.val
    var vals = a.val.strip(chars = {'(', ')'}).split(',')
    if vals.len != 3: return false
    try:
      echo "vals: ", vals
      dst = (dset: vals[0].strip(chars = {'"'}),
             lower: parseFloat(vals[1].strip),
             upper: parseFloat(vals[2].strip))
      echo "Parsed it to ", dst
      result = true
    except:
      result = false

  proc argHelp*(dfl: GenericCut; a: var ArgcvtParams): seq[string] =
    result = @[ a.argKeys, "(dset: string, lower, upper: float)", $dfl ]

  dispatch(main, help = {
    "fname" : "The input file name to print data from",
    "chip" : "Chip to read data from in the file",
    "run" : "Run to read data from in the file",
    "dset" : "The dataset to read",
    "reco" : "If true will read from `/reconstruction/` path",
    "likelihood" : "If true will read from `/likelihood/` path",
    "head" : "Only print `head` first elements",
    "cuts" : "Allows to cut data used for event display. Only shows events of data passing these cuts.",
    "verbose" : "Print the whole dataset (not just summaries)"
    })
