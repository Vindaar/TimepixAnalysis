import ggplotnim
import nimhdf5
import os, strformat, sequtils
import ingrid / tos_helpers

proc readDset[T](h5f: H5File, dset: string, run, chip: int, eventIdxs: seq[int], _: typedesc[T]): seq[T] =
  let dset = h5f[(recoPath(run, chip).string / dset).dset_str]
  let data = dset[eventIdxs, special_type(T), T]
  doAssert data.len == 1
  result = data[0]

proc main(fname: string, run: int, eventNumber: int, chip: int,
          outpath = "out") =
  ## A *very* simple script to simply plot a single event from a single run of a single chip.
  ##
  ## This plots *all* clusters from a single event number!
  var h5f = H5open(fname, "r")
  let evNums = h5f[recoPath(run, chip).string / "eventNumber", int]
  let evNumIdx = zip(toSeq(0 ..< evNums.len), evNums).filterIt(it[1] == eventNumber).mapIt(it[0])
  let xs = h5f.readDset("x", run, chip, evNumIdx, uint8).mapIt(it.int)
  let ys = h5f.readDset("y", run, chip, evNumIdx, uint8).mapIt(it.int)
  let cs = h5f.readDset("charge", run, chip, evNumIdx, float)
  let df = toDf(xs, ys, cs)
  ggplot(df, aes("xs", "ys", color = "cs")) +
    geom_point() +
    ggtitle(&"Run {run} event number {eventNumber} for chip {chip}") +
    ggsave(outpath / &"run_{run}_eventNum_{eventNumber}.pdf")

  discard h5f.close()

when isMainModule:
  import cligen
  dispatch main
