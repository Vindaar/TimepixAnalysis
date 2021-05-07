import std / [sequtils, strformat]
import nimhdf5
import ggplotnim
import cligen

type
  Tpx3MetaData = object
    index_start: uint32
    index_stop: uint32
    data_length: uint32
    timestamp_start: float
    timestamp_stop: float
    scan_param_id: uint32
    discard_error: uint32
    decode_error: uint32
    trigger: float

template hasTimestamp(x: uint32): untyped =
  (x and 0xF0000000'u32) shr 28 == 0b0101

proc rawDataToDut(data: seq[uint32]) =
  echo data.len

  if data.len < 10: return # TODO: arbitrary
  var idx = 0
  var res = newSeqOfCap[uint64](data.len div 2 + 1)
  var indices = newSeqOfCap[int](data.len div 2 + 1)
  var data = newSeq[uint64](data.len div 2 + 1)

  # need notion of even, odd to know if started an element or already has entry
  # essentially like `j` for timestamps, but general for all links.

  var j = 0
  for i, el in data:
    echo j
    let linkIdx = el and 0xfe000000'u32 shr 25
    # TODO: check if linkIdx in {0 .. 7}
    if el.hasTimestamp:
      let k = el and 0xFFFFFF'u32
      if j mod 2 == 0:
        echo "no?"
        res.add k.uint64
        indices.add j
      else:
        echo "/{}"
        var x = res[idx] shl 24 + k.uint64
        x = x shl 24 + k
        res[idx] = x or 0b0101 shl 48
        inc idx
      inc j
    if i mod 2 == 0:

    data[i]
  echo res.len
  echo res[0 .. 100]

proc main(fname: string) =
  let h5f = H5file(fname, "r")
  let path = "raw_data"
  let data = h5f[path, uint32]
  let meta = h5f["meta_data", Tpx3MetaData]

  echo meta[0 .. 100]
  ## - iterate over meta data chunks, one line describes one chunk in `data`
  ## - extract timestamps and associate data in order?

  ## start new timestamp "chunk" when
  # `timestamp_combined_filter = (data_combined and 0xF000000000000) shr 48 == 0b0101`
  # in our code: if condition: append

  ## - for now: skip chunks with errors
  echo data[0 .. 100]
  #for el in data:
  for slice in meta:
    rawDataToDut(data[slice.index_start .. slice.index_stop])

  #echo "Timestamp? ", el.hasTimestamp
  #  if el.hasTimestamp:
      # NOTE: if number of timestamp containing words not mod 2 == 0, then drop last



  discard h5f.close()

dispatch main
