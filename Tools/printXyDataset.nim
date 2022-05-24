import nimhdf5, cligen, seqmath
import std / [strformat, os, stats, sequtils, sets, algorithm, strutils]
import ingrid / tos_helpers

type
  ## A generic cut on input data using dset & low / high values
  GenericCut* = tuple[dset: string, lower, upper: float]

  Stats = object
    min: float
    max: float
    mean: float
    sum: float

proc readIndices(h5f: H5File, path: string, cuts: seq[GenericCut]): seq[int] =
  ## get the correct indices we need to read for the data `dset` such that it matches
  ## the given cuts. Cuts are only allowed from datasets at the same level as `path`
  ## and must be readable as `float`.
  # 1. get `dset`
  let dset = h5f[path.dset_str]
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

proc singleRun(h5f: H5File, path: string,
               head: int, verbose: bool, cuts: seq[GenericCut] = @[]): Stats =
  let h5dset = h5f[path.dset_str]
  let idxs = readIndices(h5f, path, cuts)
  withDset(h5dset):
    var dsetL = idxs.mapIt(dset[it])
    if head > 0:
      dsetL = dsetL[0 ..< head]
    if verbose:
      echo dsetL
    when typeof(dsetL[0]) is SomeNumber:
      let dsetF = dsetL.mapIt(it.float)
      result = Stats(min: dsetF.min, max: dsetF.max,
                     mean: dsetF.mean, sum: dsetF.sum)
      echo "Min: ", result.min
      echo "Max: ", result.max
      echo "Mean: ", result.mean
      echo "Sum: ", result.sum
    elif typeof(dsetL[0]) is seq:
      when typeof(dsetL[0][0]) is SomeNumber:
        let dsetF = dsetL.flatten.mapIt(it.float)
        result = Stats(min: dsetF.min, max: dsetF.max,
                       mean: dsetF.mean, sum: dsetF.sum)
        echo "Min: ", result.min
        echo "Max: ", result.max
        echo "Mean: ", result.mean
        echo "Sum: ", result.sum
    else:
      echo "No aggregates for type: ", typeof(dset)

proc getPath(base: string, run: int, chip: int, dset: string, global: bool): string =
  result = ""
  if not global:
    if chip < 0:
      quit("For chip related datasets, must hand a chip number!")
    result = &"{base}{run}/chip_{chip}/{dset}"
  else:
    result = &"{base}{run}/{dset}"

proc main(fname: string, dset:string, global = false,
          run = -1, chip = -1,
          reco = false, likelihood = false, head = 0, verbose = false,
          cuts: seq[GenericCut] = @[]) =
  let h5f = H5open(fname, "r")
  var path = rawDataBase()
  if reco:
    path = recoBase()
  if likelihood:
    path = likelihoodBase()

  echo "Dataset: ", dset
  if run >= 0:
    let dsetPath = getPath(path, run, chip, dset, global)
    discard singleRun(h5f, dsetPath, head, verbose, cuts)
  else:
    var stats = newSeq[Stats]()
    for (num, grp) in runs(h5f, path):
      echo "Run: ", num
      echo "--------------------------------------------------"
      let dsetPath = getPath(path, num, chip, dset, global)
      stats.add singleRun(h5f, dsetPath, head, verbose, cuts)
    echo "=================================================="
    echo "Combined stats:"
    echo "=================================================="
    echo "Min: ", stats.mapIt(it.min).min
    echo "Max: ", stats.mapIt(it.max).max
    echo "Mean: ", stats.mapIt(it.mean).mean
    echo "Sum: ", stats.mapIt(it.sum).sum

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
    "run" : "Run to read data from in the file. If not given all runs (and summary) will be printed.",
    "dset" : "The dataset to read",
    "global" : "If set will read a 'global' (i.e. in root of run group) dataset.",
    "reco" : "If true will read from `/reconstruction/` path",
    "likelihood" : "If true will read from `/likelihood/` path",
    "head" : "Only print `head` first elements",
    "cuts" : "Allows to cut data used for event display. Only shows events of data passing these cuts.",
    "verbose" : "Print the whole dataset (not just summaries)"
    })
