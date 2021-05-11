import std / [sequtils, strformat, endians, sets, algorithm, tables, options]
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

  proc assignToAExtension(res: var Tpx3Data, ToAExtension: Option[uint64]) =
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

proc rawDataToDut(data: seq[uint32], chunkNr: int): seq[Tpx3Data] =
  if data.len < 10: return # TODO: arbitrary
  var idx = 0
  var tstamp = newSeqOfCap[uint64](data.len div 2 + 1)
  # indices stores the indices at which the timestamps start in the raw data
  var indices = newSeqOfCap[int](data.len div 2 + 1)
  var res = newSeq[uint64](data.len) # div 2 + 1)

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
  var idxToKeep = initHashSet[int]()
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
        let lkIdx = linkIdxs[link] # this is the index last modifieyd for link# `link`
        res[lkIdx] = (x.uint64 shl 16) + res[lkIdx]
      linksEven[link] = not linksEven[link]

  ## TODO: remove everything but idxToKeep. Need to walk data again?
  ## Why? Cannot just append (using `newSeqOfCap` and `add`?) instead of
  ## accessing `i`.
  ## Do both and compare
  result = newSeq[Tpx3Data](idxToKeep.card)
  if indices.len > 0:
    doAssert res.len >= indices[^1]
    # TOA extension
    let indIdx = indices.toSet
    var tstamp: uint64
    var i = 0
    for idx in idxToKeep.toSeq.sorted:
      if idx in indIdx:
        # timestamp
        tstamp = res[idx]
      else:
        # data
        result[i] = toData(res[idx], 0b00, ToA_Extension = some(tstamp))
        inc i
    result.setLen(i)
  else:
    # no TOA extension
    for i, val in idxToKeep.toSeq.sorted:
      result[i] = toData(res[val], 0b00)
  result = result.sortedByIt(it.TOA)

proc main(fname: string) =
  let h5f = H5file(fname, "r")
  let path = "raw_data"
  let data = h5f[path, uint32]
  let meta = h5f["meta_data", Tpx3MetaData]

  #let config = h5f["configuration/generalConfig", Tpx3Config]
  #for el in config:
  #  echo el.configuration
  #  echo el.value

  var h5fout = H5File("/tmp/testtpx3.h5", "rw")

  var all: seq[Tpx3Data]
  for i, slice in meta:
    if slice.index_stop - slice.index_start > 10:
      all.add rawDataToDut(data[slice.index_start ..< slice.index_stop], chunkNr = i)

  all = all.filterIt(it.x != 0 and it.y != 0)
  let dset = h5fout.create_dataset("interpreted/hit_data_0", all.len, dtype = Tpx3Data)
  dset[dset.all] = all
  discard h5fout.close()


  discard h5f.close()

dispatch main
