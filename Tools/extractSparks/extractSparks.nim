import ingrid / tos_helpers
import plotly, sequtils, strutils, os
import nimhdf5, loopfusion, arraymancer, seqmath, docopt

when defined(linux):
  const commitHash = staticExec("git rev-parse --short HEAD")
else:
  const commitHash = ""
# get date using `CompileDate` magic
const currentDate = CompileDate & " at " & CompileTime

const docTmpl = """
Version: $# built on: $#
Counts number of sparks in (StartX - EndX) / (StartY / EndY), removes them
and creates an occupancy map of the rest.

Usage:
  extractSparks <H5file> [options]

Options:
  -h, --help             Show this help
  --version              Show the version number
"""
const doc = docTmpl % [commitHash, currentDate]

# NOTE: upper bounds are inclusive!
const StartX = 0
const EndX = 8
const StartY = 150
const EndY = 158
const MinPix = 10
const NPix = 256

template addPixelsToOccupancy[T](a1: var Tensor[T], a2: Tensor[T]) =
  ## template to stack two frames for an occupancy
  a1 = a1 .+ a2

proc buildEvent[T, U](x, y: seq[T], ch: seq[U]): Tensor[float] =
  ## takes all events via their *indices* (not real event numbers) and returns
  ## a (256, 256) tensor for each event
  result = newTensor[float]([NPix, NPix])
  forZip ix in x, iy in y, ich in ch:
    result[iy.int, ix.int] = ich.float

proc nonZero(t: Tensor[float]): bool =
  var count = 0
  for ix in StartX .. EndX:
    for iy in StartY .. EndY:
      if t[iy, ix] > 0:
        inc count
  result = if count > MinPix: true else: false

proc main =

  let args = docopt(doc)
  let runFile = $args["<H5file>"]
  var h5f = H5file(runFile, "r")
  let fInfo = getFileInfo(h5f)

  var count = 0
  var occupancy = newTensor[float]([NPix, NPix])

  for r in fInfo.runs:
    let chip = fInfo.centerChip
    let
      vlen = special_type(uint8)
      vlen16 = special_type(uint16)
      x = h5f[recoDataChipBase(r) & $chip / "x", vlen, uint8]
      y = h5f[recoDataChipBase(r) & $chip / "y", vlen, uint8]
      ch = h5f[recoDataChipBase(r) & $chip / "ToT", vlen, uint16]
    for i in 0 ..< x.high:
      if i mod 1000 == 0:
        echo i, " events done, current count: ", count
      let ev = buildEvent(x[i], y[i], ch[i])
      # given event check if nonzero in area
      if nonZero(ev):
        echo "Found one!"
        inc count
      else:
        occupancy.addPixelsToOccupancy(ev)

  echo "Count was ", count
  heatmap(occupancy.toRawSeq.reshape2D(@[256, 256])).show()



when isMainModule:
  main()
