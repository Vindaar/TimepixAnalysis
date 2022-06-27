import std / [sequtils, strformat, endians, algorithm, tables, options, sugar, os, strutils]
import nimhdf5
import cligen

import std / intsets


from ingrid / ingrid_types import RunTypeKind

when not defined(blosc):
  import macros
  static: error("Compilation without `blosc` support not supported! Required to " &
    "read Timepix3 data. Add `-d:blosc` compilation option.")

type
  Tpx3MetaData = object
    index_start: uint64
    index_stop: uint64
    data_length: uint32
    timestamp_start: float
    timestamp_stop: float
    scan_param_id: uint32
    discard_error: uint32
    decode_error: uint32
    trigger: float

  Tpx3Data* = object
    data_header*: uint8
    header*: uint8
    hit_index*: uint64
    x*: uint8
    y*: uint8
    TOA*: uint16
    TOT*: uint16
    EventCounter*: uint16
    HitCounter*: uint8
    FTOA*: uint8
    scan_param_id*: uint16
    chunk_start_time*: cdouble
    iTOT*: uint16
    TOA_Extension*: uint64
    TOA_Combined*: uint64

  #Tpx3Config = object
  #  configuration: cstring
  #  value: uint16

type
  ## NOTE: this data structure is not efficient to store bits of course.
  ## However, we only use it at CT to generate some lookup tables. Therefore
  ## simplicity of the implementation trumps efficiency.
  BitArray[T: static[int]] = object
    data: array[T, bool]

  ReturnData*[T: static[int]] = object
    ar*: array[T, uint64]
    len*: uint8

const headers = {0b101.uint16, 0b111, 0b110, 0b100, 0b011}

const headerMap = { "Acquisition" : 0b101.uint16,
                    "StopMatrix"  : 0b111.uint16,
                    "CTPR"        : 0b110.uint16,
                    "PCR"         : 0b100.uint16,
                    "Control"     : 0b011.uint16  }.toTable()


proc len(b: BitArray): int =
  b.data.len

proc high(b: BitArray): int =
  b.data.len - 1

template `^^`(s, i: untyped): untyped =
  (when i is BackwardsIndex: s.len - int(i) else: int(i))

proc createBitarray(size: static[int]): BitArray[size] =
  for i in 0 ..< size:
    result.data[i] = false

proc `[]=`*[T, U](b: var BitArray, inds: HSlice[T, U], val: SomeInteger) =
  let iStart = b ^^ inds.a
  let iEnd   = b ^^ inds.b
  let nInds = abs(iEnd - iStart) + 1

  if nInds > b.len:
    raise newException(IndexError, &"Slice of {inds} is out of range for BitArray of size {b.len}")
  if val.uint64 > (2 ^ nInds).uint64:
    raise newException(ValueError, &"Value of {val} is too large for {nInds} bits slice! " &
                                   &"Max size is {2 ^ nInds}")

  var m = 0
  var mval = val.uint
  var i = 0
  while mval > 0:
      b.data[i] = (mval and 1).bool
      mval = mval shr 1
      inc i
  return
  if iEnd > iStart:
    for x in iStart .. iEnd:
      let isBitOne = (mval and 1.uint).bool
      b.data[x] = if isBitOne: true else: false
      mval = mval shr 1
  else:
    for x in countdown(iStart, iEnd):
      let isBitOne = (mval and 1.uint).bool
      b.data[x] = if isBitOne: true else: false
      mval = val shr x

proc `[]=`[T: not HSlice, U: SomeInteger | bool](b: var BitArray, ind: T, val: U) =
  when val is SomeInteger:
    let boolVal = if val == 1: true else: false
  elif val is bool:
    let boolVal = val
  let i = b ^^ ind
  b.data[i] = boolVal

proc `[]`[T: BackwardsIndex | SomeInteger](b: BitArray, ind: T): uint =
  let i = b ^^ ind
  if i >= b.len:
    raise newException(IndexError , &"Index of value {i} is out of range for BitArray of size {b.len}")
  result = if b.data[i] == false: 0.uint else: 1.uint

proc `[]`[T, U](b: BitArray, inds: HSlice[T, U]): uint =
  if inds.len > b.len:
    raise newException(IndexError, &"Slice of {inds} is out of range for BitArray of size {b.len}")
  let iStart = b ^^ inds.a
  let iEnd   = b ^^ inds.b
  var m = 0
  for x in iStart .. iEnd:
    if b[x] == 1.uint:
      result += (2 ^ m).uint
    inc m

proc toValue(b: BitArray): uint = b[0 .. b.high]

proc printBytes*(ba: BitArray, asBytes = false): string =
  ## prints the BitArray as a list of individual bytes
  result = newStringOfCap(8 * ba.len + 50)
  result = "["
  let nbytes = ba.len div 8
  if asBytes == false:
    for i in 0 ..< nbytes:
      for j in 0 ..< 8:
        let ind = (i * 8 + j).int
        result.add($(ba[ind]))
      if i != nbytes - 1:
        result.add ", "
    result.add "]"
  else:
    result = $(ba.toByteList)

proc `$`(b: BitArray): string =
  b.printBytes

proc toByteList(b: BitArray): seq[uint] =
  ## returns a seq of bytes based for the given `BitArray`
  ## Notes: needs to be a BitArray of size N bytes. If not, last
  ## elements are dropped
  let nbytes = b.len div 8
  for i in 0 ..< nbytes:
    let ind = i * 8
    result.add b[ind .. ind + 7]

proc bitwordToByteSeq(val: SomeInteger, size: static[int]): seq[uint] =
  var b = createBitarray(size)
  b[0 ..< size] = val
  result = b.toByteList

proc initLfsr14Lut(): seq[uint16] =
  ## Generates a 14bit LFSR according to Manual v1.9 page 19
  result = newSeq[uint16](2^14)
  var lfsr = createBitarray(14)
  lfsr[0 .. 7] = 0xFF'u16
  lfsr[8 .. 13] = 63'u16
  var dummy = 0'u16
  for i in 0 ..< 2^14:
    result[lfsr[0 .. 13].int] = i.uint16
    dummy = lfsr[13].uint16
    lfsr[13] = lfsr[12].uint16
    lfsr[12] = lfsr[11].uint16
    lfsr[11] = lfsr[10].uint16
    lfsr[10] = lfsr[9].uint16
    lfsr[9] = lfsr[8].uint16
    lfsr[8] = lfsr[7].uint16
    lfsr[7] = lfsr[6].uint16
    lfsr[6] = lfsr[5].uint16
    lfsr[5] = lfsr[4].uint16
    lfsr[4] = lfsr[3].uint16
    lfsr[3] = lfsr[2].uint16
    lfsr[2] = lfsr[1].uint16
    lfsr[1] = lfsr[0].uint16
    lfsr[0] = (lfsr[2] xor dummy xor lfsr[12] xor lfsr[13]).uint16
  result[2 ^ 14 - 1] = 0'u16

proc initLfsr10Lut(): seq[uint16] =
  ## Generates a 10bit LFSR according to Manual v1.9 page 19
  result = newSeq[uint16](2 ^ 10)
  var lfsr = createBitarray(10)
  lfsr[0 .. 7] = 0xff'u16
  lfsr[8 .. 9] = 0b11'u16
  var dummy = 0'u16
  for i in 0 ..< 2^10:
    result[lfsr[0 .. 9].int] = i.uint16
    dummy = lfsr[9].uint16
    lfsr[9] = lfsr[8].uint16
    lfsr[8] = lfsr[7].uint16
    lfsr[7] = lfsr[6].uint16
    lfsr[6] = lfsr[5].uint16
    lfsr[5] = lfsr[4].uint16
    lfsr[4] = lfsr[3].uint16
    lfsr[3] = lfsr[2].uint16
    lfsr[2] = lfsr[1].uint16
    lfsr[1] = lfsr[0].uint16
    lfsr[0] = (lfsr[7] xor dummy).uint16
  result[2 ^ 10 - 1] = 0'u16

proc initLfsr4Lut(): seq[uint16] =
  ## Generates a 4bit LFSR according to Manual v1.9 page 19
  result = newSeq[uint16](2 ^ 4)
  var lfsr = createBitarray(4)
  lfsr[0 .. 3] = 0xF'u16
  var dummy = 0'u16
  for i in 0 ..< 2^4:
    result[lfsr[0 .. 3].int] = i.uint16
    dummy = lfsr[3].uint16
    lfsr[3] = lfsr[2].uint16
    lfsr[2] = lfsr[1].uint16
    lfsr[1] = lfsr[0].uint16
    lfsr[0] = (lfsr[3] xor dummy).uint16
  result[2 ^ 4 - 1] = 0'u16

proc initGray14Lut(): seq[uint16] =
  ## Generates a 14bit gray according to Manual v1.9 page 19
  result = newSeq[uint16](2 ^ 14)
  var i = 0
  for j in 0 ..< 2^14:
    var encodedValue = createBitArray(14) #48
    encodedValue[0 .. 13] = j.uint16 #47
    var grayDecryptV = createBitArray(14) #48
    grayDecryptV[13] = encodedValue[13] #47
    for i in countdown(12, 0): #46
      grayDecryptV[i] = (grayDecryptV[i+1] xor encodedValue[i]).uint16
    result[j] = grayDecryptV.toValue.uint16

# lookup tables for LFSR values
const Lfsr14Lut = initLfsr14Lut()
const Lfsr10Lut = initLfsr10Lut()
const Lfsr4Lut = initLfsr4Lut()
const Gray14Lut = initGray14Lut()

## XXX: replace `ToAExtension` `Option` by a static boolean. That way don't have
## to do runtime check on whether variable isSome
proc toData(x: uint64, opMode: uint8, vco = false, ToAExtension = none(uint64)): Tpx3Data =
  #echo x
  let pixel = (x shr 28) and 0b111'u64
  let super_pixel = (x shr (28 + 3)) and 0x3f
  let right_col = pixel > 3
  let eoc = (x shr (28 + 9)) and (0x7f)

  result.data_header = (x shr 47).uint8
  result.header = (x shr 44).uint8
  result.y = ((super_pixel * 4).int + (pixel.int - (if right_col: 1 else: 0) * 4)).uint8
  result.x = ((eoc * 2).int + (if right_col: 1 else: 0)).uint8
  if not vco:
    result.HitCounter = Lfsr4Lut[x and 0xF].uint8
    result.FTOA = 0'u8
  else:
    result.HitCounter = 0'u8
    result.FTOA = (x and 0xF).uint8

  proc assignToAExtension(res: var Tpx3Data, ToAExtension: Option[uint64]) {.inline.} =
    if ToA_Extension.isSome:
      res.TOA_Extension = ToA_Extension.unsafeGet and 0xFFFFFFFFFFFF'u64 # remove header marking it as timestamp
      res.TOA_Combined = (ToA_Extension.unsafeGet and 0xFFFFFFFFC000'u64) + res.TOA
    else:
      res.TOA_Extension = 0'u64
      res.TOA_Combined = 0'u64

  case opMode
  of 0b00: # ToT and ToA
    result.TOT = Lfsr10Lut[(x shr 4'u64) and 0x3ff]
    result.TOA = Gray14Lut[(x shr 14'u64) and 0x3fff]
    result.EventCounter = 0
    result.assignToAExtension(ToAExtension)
  of 0b01: # ToA
    result.TOA = Gray14Lut[(x shr 14'u64) and 0x3fff]
    result.EventCounter = 0
    result.TOT = 0
    result.assignToAExtension(ToAExtension)
  else: # Event and iToT
    result.iTOT = Lfsr14Lut[(x shr 14'u64) and 0x3fff]
    result.EventCounter = Lfsr10Lut[(x shr 4'u64) and 0x3ff]
    result.TOT = 0
    result.TOA = 0
    result.TOA_Extension = 0
    result.TOA_Combined = 0

template hasTimestamp(x: uint32 | uint64): untyped =
  when typeof(x) is uint32:
    (x and 0xF0000000'u32) shr 28 == 0b0101
  else:
    (x and 0xF000000000000'u64) shr 48 == 0b0101

proc rawDataToDut(data: openArray[uint32], chunkNr: int,
                  tstamp: var seq[uint64], indices: var seq[int], res: var seq[uint64],
                  idxToKeep: var IntSet, resultSeq: var seq[Tpx3Data]) =
  if data.len < 10: return # TODO: arbitrary
  var idx = 0
  ## it has capacity, so setlen is "free"
  tstamp.setLen(0)
  indices.setLen(0)
  res.setLen(data.len)

  var linkIdxs = newSeq[int](8)
  var linksEven = newSeqWith[bool](8, true)
  # need notion of even, odd to know if started an element or already has entry
  # essentially like `j` for timestamps, but general for all links.
  ## TODO: check for 0b0101 TOA extension packets.
  ## when encountered:
  ## start a "new packet"
  ## this will be done offline. I.e.
  ## we store all indices of the 64 bit uint words which start
  ## with 0b0101 so that at the end we can build the packages / drop everything
  ## that does not belong to anything

  var j = 0
  var lastIndex = 0
  # keep track of all indices we have touched. Then we know which
  # indices we have to copy to final result
  idxToKeep.clear()
  for i, el in data:
    let link = (el and 0xfe000000'u32) shr 25
    if el.hasTimestamp:
      let k = el and 0xFFFFFF'u32
      if j mod 2 == 0:
        tstamp.add k.uint64
        indices.add i
        lastIndex = i
      else:
        var x = tstamp[idx] shl 24 + k.uint64
        tstamp[idx] = x or (0b0101 shl 48)
        # set tstamp to at index where this was created
        res[lastIndex] = tstamp[idx]
        idxToKeep.incl lastIndex
        inc idx
      inc j
    if link in {0 .. 7}:
      var k = el and 0xFFFFFF'u32
      if linksEven[link]:
        var x: uint32
        swapEndian32(x.addr, k.addr)
        res[i] = (x.uint64 shr 8)
        # keep track of global index so that we can insert it into result once we
        # are in odd branch below
        linkIdxs[link] = i
        idxToKeep.incl i
      else:
        var x: uint32
        swapEndian32(x.addr, k.addr)
        let lkIdx = linkIdxs[link] # this is the index last modified for link# `link`
        res[lkIdx] = (x.uint64 shl 16) + res[lkIdx]
      linksEven[link] = not linksEven[link]

  ## TODO: remove everything but idxToKeep. Need to walk data again?
  ## Why? Cannot just append (using `newSeqOfCap` and `add`?) instead of
  ## accessing `i`.
  ## Do both and compare
  resultSeq.setLen(idxToKeep.card) # `.card` is maximum needed sice, reality will be less
  if indices.len > 0:
    doAssert res.len >= indices[^1]
    # TOA extension
    let indIdx = indices.toIntSet
    var tstamp: uint64
    var i = 0
    for idx in idxToKeep.toSeq.sorted:
      if idx in indIdx:
        # timestamp
        tstamp = res[idx]
      else:
        # data
        resultSeq[i] = toData(res[idx], 0b00, ToA_Extension = some(tstamp))
        inc i
    resultSeq.setLen(i) # clip length to length that is actually used
  else:
    # no TOA extension
    for i, val in idxToKeep.toSeq.sorted:
      resultSeq[i] = toData(res[val], 0b00)

  resultSeq.sort((x, y: Tpx3Data) => system.cmp(x.TOA, y.TOA))

proc rawDataToDut(data: openArray[uint32], chunkNr: int): seq[Tpx3Data] =
  var tstamp = newSeqOfCap[uint64](data.len div 2 + 1)
  # indices stores the indices at which the timestamps start in the raw data
  var indices = newSeqOfCap[int](data.len div 2 + 1)
  var res = newSeq[uint64](data.len) # div 2 + 1)
  var idxToKeep = initIntSet()
  result = newSeqOfCap[Tpx3Data](data.len)
  rawDataToDut(data, chunkNr, tstamp, indices, res, idxToKeep, result)

proc toRunPath(run: int): string = "/interpreted/run_" & $run & "/"

proc writeAttributes(h5f: H5File, runType: RunTypeKind, run, badSliceCount, badBatchCount, batchSize: int) =
  ## Writes all (at this point knowable) attributes to the `/interpreted` group and run
  ## group (if run number given).
  ## These attributes are the main ones that describe the `FileInfo` object defined in `ingrid_types`.
  var grp = h5f["/interpreted".grp_str]
  grp.attrs["TimepixVersion"] = "Timepix3"
  grp.attrs["runType"] = $runType
  grp.attrs["runFolderKind"] = "rfUnknown"
  grp.attrs["centerChip"] = 0 ## TODO: once multi links are in place, this needs to be adjusted
  #grp.attrs["centerChipName"] = chipName
  if run >= 0:
    var runGrp = h5f[toRunPath(run).grp_str]
    runGrp.attrs["numChips"] = 1 ## TODO: once multi links are in place, this needs to be adjusted
    runGrp.attrs["batchSize"] = batchSize
    runGrp.attrs["badSliceCount"] = badSliceCount # number of dropped slices due to errors
    runGrp.attrs["badBatchCount"] = badBatchCount # number of dropped batches (of `batchSize` words)
                                                  # due to decompression failure
proc processSlice(data: seq[uint32], slice: Tpx3MetaData, oldIdx, loopIdx: int): seq[Tpx3Data] =
  ## Performs processing of a single slice of data from `sliceStart` to `sliceStop`, while taking into
  ## account sanity checks and conversion of global slicing indices to `data` 'local' ones (`data` is
  ## a single `batch` instead of all data in the file).
  let
    sliceStart = slice.index_start
    sliceStop = slice.index_stop
  doAssert oldIdx.uint64 <= sliceStart and oldIdx.uint64 <= sliceStop, "uint64 underflow detected! oldIdx = " & $oldIdx &
    " sliceStart = " & $sliceStart & " sliceStop = " & $sliceStop
  let startIdx = sliceStart - oldIdx.uint64
  let stopIdx = sliceStop - oldIdx.uint64
  doAssert startIdx < (int64.high).uint64 and stopIdx < (int64.high).uint64, "int64 overflow detected"
  # can't really have `stopIdx > int64.high`, so `int` conversion is fine
  doAssert stopIdx.int - 1 <= data.len, " Stop idx: " & $stopIdx & " data.len " & $data.len
  result = rawDataToDut(toOpenArray(data, startIdx.int, stopIdx.int - 1), chunkNr = loopIdx)

proc findSliceIdx(num: int, cumNum: seq[int]): int =
  ## finds the correct index up to where to read given a desired
  ## number (upper bound) of `num` elements.
  ## Returns the *index* of the slice and not the data word index!
  ##
  ## Note: we use `upperBound` to handle the case `num == cumNum[result]`
  ## such that we get the *next* index.
  result = cumNum.lowerBound(num)
  if result == cumNum.len: # `lowerBound` returns index *after* seq, if bigger than every element
    dec result

proc readNextChunk(data: var seq[uint32], h5f: var H5File, inputDset: var H5Dataset, cumNum: seq[int],
                   dataFromIdx, dataToIdx, batchIdx, processedIdx, sliceIdx, badBatchCount: var int,
                   batchSize, numSlices: int) =
  ## Reads the next (readable) chunk of size `batchSize` from the input dataset `inputDset`.
  ## Takes into account to adjust the indices correctly.
  ##
  ## `h5f` and `inputDset` are given to reopen file in case of decompression error.
  ##
  ## If the decompression of the data fails, will jump to the next batch. This means we
  ## might drop data of the size of `batchSize` if "bad sectors" are found in the file.
  ## In the future we might want to bisect the data to find the good data inside the
  ## `batchSize` blocks.
  ##
  ## `dataFromIdx`   : index of first slice part of current `data`
  ## `dataToIdx`     : index of last slice still part of current `data`
  ## `batchIdx`      : index of `batchSize` words processed so far
  ## `processedIdx`  : number of data word (`uint32`) elements processed in total
  ## `sliceIdx`      : current data slice we are processing
  ## `batchSize`     : size of batch chunk (O(50 Mio.) data words)
  ## `badBatchCount` : number of batches (of size `batchSize`) dropped due to bad data.
  # we loop here to handle possible `blosc` decompression errors
  var success = false
  while not success:
    dataFromIdx = dataToIdx # new start is previous end
    let fromIdx = cumNum[dataFromIdx] # (we read up to `toIdx - 1` in last read)
    dataToIdx = findSliceIdx(fromIdx + batchSize, cumNum) # find index up to where to read
    let toIdx = cumNum[dataToIdx] # and the corresponding data index
    try:
      if toIdx > fromIdx: # if the last slice is empty, `toIdx == fromIdx` will hold
        let perc = (sliceIdx.float / numSlices.float) * 100.0
        echo &"[INFO]: Reading from word {fromIdx} to {toIdx} of {numSlices} slices, {perc:.2f} % processed"
        data = inputDset.read_hyperslab(uint32, @[fromIdx],
                                        count = @[toIdx - fromIdx])
        echo "[INFO]: ...reading done"
        inc batchIdx
      # regardless, we succeeded (even if we didn't have anything to do)
      success = true
    except HDF5BloscDecompressionError:
      # faled to decompress blosc data! update indices, continue and try again
      echo "[WARN]: Decompression of batch after ", batchIdx, " failed!"
      inc processedIdx, (toIdx - fromIdx) # skip `batchSize` words (mark them "processed")
      inc badBatchCount
      inc batchIdx
      # after skip, continue with next slice after one up to which we tried to read
      sliceIdx = dataToIdx + 1
      ## NOTE: in order to continue after a decompression error, we have to fully close the file
      ## and reopen. Not sure why, but if we don't any further read will also trigger a decompression
      ## error (or even a completely different "duplicate link" error).
      let fname = h5f.name
      discard h5f.close()
      h5f = H5open(fname, "r")
      inputDset = h5f[inputDset.name.dset_str]
      continue

proc parseInputFile(h5fout: H5File, # file we write to
                    fname: string, # file we will parse
                    runType: RunTypeKind,
                    run: int,
                    allowErrors: bool,
                    verbose: bool,
                    batchSize: int
                   ) =
  var h5f = H5file(fname, "r")
  let path = "raw_data"
  let meta = h5f["meta_data", Tpx3MetaData]
  # TODO: reading config still broken
  #let config = h5f["configuration/generalConfig", Tpx3Config]
  #for el in config:
  #  echo el.configuration
  #  echo el.value

  # check if run exists in output file, if not copy over `configuration` group
  let exists = toRunPath(run) in h5fout
  var dset: H5Dataset
  if not exists:
    let cfg = h5f["/configuration".grp_str]
    let status = h5f.copy(cfg, some(toRunPath(run) / "configuration"), some(h5fout))
    if not status:
      raise newException(IOError, "Could not copy over `/configuration` from " & $fname & " to " & $h5fout.name)
    let filter = H5Filter(kind: fkZlib, zlibLevel: 2)
    dset = h5fout.create_dataset(toRunPath(run) / "hit_data", 0, dtype = Tpx3Data,
                                 chunksize = @[50_000],
                                 maxshape = @[int.high], filter = filter)
  else:
    echo "[INFO]: Appending input data from ", fname, " to ", h5fout.name
    dset = h5fout[(toRunPath(run) / "hit_data").dset_str]
  var all: seq[Tpx3Data] = newSeqOfCap[Tpx3Data](batchSize)
  var inputDset = h5f[path.dset_str]
  let cumNum = meta.mapIt(it.data_length.int).cumsum
  var
    batchIdx = 0
    processedIdx = 0
    badSliceCount = 0 ## number of bad slices (size of the chunk) discarded with errors in the data
    badBatchCount = 0 ## number of bad batches (size `batchSize`) discarded as it could not be decompressed
    dataFromIdx = 0
    dataToIdx = 0
    i = 0

  # read first chunk of data
  var data: seq[uint32]
  data.readNextChunk(h5f, inputDset, cumNum, dataFromIdx, dataToIdx, batchIdx, processedIdx, i, badBatchCount, batchSize, meta.len)
  while i < meta.len:
    let slice = meta[i]
    let num = (slice.index_stop - slice.index_start).int
    if not allowErrors and
      (slice.discard_error > 0'u32 or slice.decode_error > 0'u32):
      # drop chunks with discard or decoding errors
      if verbose:
        echo "[INFO]: Skipping chunk with discard_error=", slice.discard_error, " decode_error= ", slice.decode_error, " at index ", i
      inc processedIdx, num
      inc badSliceCount
      inc i
      continue
    if num > 0: # only process if slice length > 0
      all.add processSlice(data, slice, cumNum[dataFromIdx], i)
      inc processedIdx, num.int
    if i >= dataToIdx: # if this slice is the last slice we read, write & read more data
      echo "[INFO]: Writing dataset batch: ", batchIdx, " from slice = ", dataFromIdx, " to ", dataToIdx, " processedIdx = ", processedIdx
      dset.add all # write to file
      all.setLen(0) # reset resulting seq
      data.readNextChunk(h5f, inputDset, cumNum, # and read new data
                         dataFromIdx, dataToIdx, batchIdx,
                         processedIdx, i, badBatchCount,
                         batchSize, meta.len)
    # increase slice index
    inc i
  dset.add all

  ## Write some attributes
  h5fout.writeAttributes(runType, run, badSliceCount, badBatchCount, batchSize)
  echo &"[INFO]: === Summary ==="
  echo &"\tProcessed file: {fname}, run number: {run} written to {h5fout.name}"
  echo &"\tNumber of slices processed: {i}"
  echo &"\tNumber of dropped slices: {badSliceCount}"
  echo &"\tNumber of dropped batches: {badBatchCount} of batch size: {batchSize}"
  echo "[INFO]: Closing input file ", h5f.name
  discard h5f.close()

proc checkValidTpx3H5(s: string) =
  let fname = s.extractFilename
  if not fname.startsWith("DataTake") or not fname.endsWith(".h5"):
    raise newException(IOError, "The given input file does not match the pattern `DataTake*.h5`. Please hand " &
      "a valid Tpx3 H5 file.")

proc main(path: string,
          outf: string = "/tmp/testtpx3.h5",
          runType: RunTypeKind = rtNone,
          run: int = 0,
          allowErrors = false,
          verbose = false,
          batchSize = 100_000_000
         ) =
  ## This tool is used to parse the raw Tpx3 data as it comes from the Tpx3 DAQ, sorts it from
  ## the 48-bit words into single 64-bit data chunks and writes it to an output H5 file.
  ##
  ## If `path` is the path to a single H5 data file, we perform the data parsing and output that
  ## single run into the given `outf`.
  ## If `path` is pointing to a directory containing multiple H5 files, we parse them according to
  ## their timestamp and place all of them into `outf`. They will be split into separate groups
  ## in the output of the type `/run_i` where `i` is an increasing integer, starting from `run`.
  ##
  ## For a single input file `run` is simply the run number associated with the given file. It can
  ## also be used to manually place multiple files into the same output file by calling this tool
  ## multiple times with different run numbers.
  ##
  ## If `allowErrors` is set to `true`, data chunks with errors (as noted in the `discard_error` and
  ## `decode_error` column of the `meta_data`) will be *kept*. By default those chunks are dropped
  ## completely.

  # 0. open output h5 file
  var h5fout = H5open(outf, "rw")
  # 1. determine if input is an H5 file or a directory
  let fInfo = getFileInfo(path)
  case fInfo.kind
  of pcDir:
    # iterate all files in directory that match `DataTake_*.h5`
    var runNum = run
    for file in walkFiles(path / "DataTake*.h5"):
      echo file
      h5fout.parseInputFile(file, runType, runNum, allowErrors, verbose, batchSize)
      inc runNum
  of pcFile:
    checkValidTpx3H5(path) # raises if not a valid file
    h5fout.parseInputFile(path, runType, run, allowErrors, verbose, batchSize)
  else:
    raise newException(IOError, "The given input is neither a file containing Tpx3 H5 files, nor a H5 file.")

  echo "[INFO]: Closing output file ", h5fout.name
  discard h5fout.close()


when isMainModule:
  dispatch main
