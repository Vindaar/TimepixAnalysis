import strutils, strformat
import sequtils
import macros
import tables
import algorithm
import future
import os
import osproc
import re
import times
import seqmath
import memfiles

template withDocopt*(doc: string): untyped =
  ## a helper template to build a commit hash and build time
  ## aware docopt string called `doc`

  when defined(linux):
    const commitHash = staticExec("git rev-parse --short HEAD")
  else:
    const commitHash = ""

  # get date using `CompileDate` magic
  const currentDate = CompileDate & " at " & CompileTime

  const docTmplPrefix = """
  Version: $# built on: $#
  """
  # concat the user given doc to it
  const doc = docTmplPrefix % [commitHash, currentDate] & doc
  doc

macro `+`*[N, M: int](a: array[N, string], b: array[M, string]): untyped =
  ## macro to concat two const arrays `a`, `b` at compile time to return a new
  ## array
  let aImpl = a.symbol.getImpl
  let bImpl = b.symbol.getImpl
  doAssert aImpl.kind == nnkConstDef
  doAssert bImpl.kind == nnkConstDef
  var tree = nnkBracket.newTree()
  for x in aImpl[2]:
    tree.add x
  for x in bImpl[2]:
    tree.add x
  result = nnkStmtList.newTree(
    tree
  )

macro getField*(tup: typed, field: static string): untyped =
  ## Mini helper to access the field `field` of a tuple from a
  ## static string, as returned by `fieldPairs`.
  let id = ident(field)
  result = nnkDotExpr.newTree(tup, id)

proc removePref*(s: string, pref: string): string =
  result = s
  result.removePrefix(pref)

template asType*[T](s: seq[T], dtype: typedesc): untyped =
  ## convenience template to convert type of a `seq` to a different type,
  ## if possible
  mapIt(s, dtype(it))

proc traverseTree(input: NimNode): NimNode =
  # iterate children
  for i in 0 ..< input.len:
    case input[i].kind
    of nnkSym:
      # if we found a symbol, take it
      result = input[i]
    of nnkBracketExpr:
      # has more children, traverse
      result = traverseTree(input[i])
    else:
      error("Unsupported type: " & $input.kind)

proc findReplace(cImpl: NimNode, assgn: NimNode): NimNode =
  ## iterates through `cImpl`, looks for the assignment in `assgn`
  ## and if found replaces the value within. If not found appends
  ## the given statement to the NimNode
  result = copy(cImpl)
  for i in 1 ..< cImpl.len:
    let curNode = cImpl[i]
    case curNode.kind
    of nnkExprColonExpr:
      # check if field is `assgn` field
      if eqIdent(curNode[0].strVal, assgn[0].strVal):
        # replace this nodes value
        result[i][1] = assgn[1]
        # return now
        return result
    else:
      let msg = "No, shouldn't be here: " & $curNode.kind & " with repr " & curNode.treeRepr
      error(msg)
  # if we end up here we didn't find the identifier
  let newExpr = nnkExprColonExpr.newTree(
    assgn[0],
    assgn[1]
  )
  result.add newExpr

macro replace*(c: typed, x: untyped): untyped =
  ## mini dsl to replace specific fields of a given object
  ## with values given in body
  ##
  ## Example:
  ## .. code-block:
  ##   type
  ##     MyObj = object
  ##       a: int
  ##       b: seq[int]
  ##       c: bool
  ##       d: float
  ##   let obj = MyObj(a: 5,
  ##                   b: @[1, 2, 3]
  ##                   c: false)
  ##   let nObj = replace(obj):
  ##     c = true
  ##     d = 5.5
  ## # will be rewritten to
  ##   let nObj = MyObj(a: 5,
  ##                   b: @[1, 2, 3]
  ##                   c: true,
  ##                   d: 5.5)
  # in new AST `getImpl` returns symbol for which it was called
  # too
  var cImpl = c.getImpl[2]
  # now just replace the correct cImpl fields
  # by the given fields of the user
  for ch in x:
    cImpl = findReplace(cImpl, ch)
  result = cImpl
  #echo result.repr

proc readNumLinesMemFile*(ff: var MemFile, buf: var seq[string], stop: int) {.inline.} =
  ## reads memory mapped slices from file `ff`, adds to buffer and
  ## stops at line `stop`
  ## NOTE: for performance make sure thatt `buf` is a sequence of cap `stop`
  buf.setLen(stop)
  var count = 0
  var lineBuf = newStringOfCap(80)
  for slice in memSlices(ff):
    lineBuf.setLen(slice.size)
    copyMem(addr lineBuf[0], slice.data, slice.size)
    buf[count] = lineBuf
    inc count
    if count == stop:
      break

iterator memLines*(ff: var MemFile, buf: var string, start = 0, stop = -1): string {.inline.} =
  var count = 0
  for slice in memSlices(ff):
    inc count
    if count < start:
      continue
    elif count == stop:
      break
    buf.setLen(slice.size)
    copyMem(addr buf[0], slice.data, slice.size)
    yield buf

iterator linesIter*(ff: string, buf: var string, start = 0, stop = -1): string {.inline.} =
  var count = 0
  var idx = 0
  while true:
    case ff[idx]
    of '\n':
      inc count
      if count == stop:
        break
      elif count >= start:
        yield buf
      buf.setLen(0)
    else:
      buf.add ff[idx]
    inc idx

proc parseUint16*(x: string): uint16 =
  var i = 0
  var digits = 0
  while i < x.len:
    case x[i]
    of ' ': discard
    of '0' .. '9':
      let val = uint16(ord(x[i]) - ord('0'))
      result *= 10
      result += val
      inc digits
    else:
      raise newException(ValueError, "Unexpected character in uint16: " & $x[i] & ", full number: " & x)
    inc i

proc getNewBound*(ind, width, size: int, up_flag: bool = true): int {.inline.} =
  # procedure to select a new bound, either upper or lower for
  # a given index, window width and size of tensor
  # inputs:
  #    ind: the index which will be at the center of the new window
  #    width: the width of the window
  #    size: the size of the tensor
  #    up_flag: a flag to select whether we choose upper or lower bounds
  let width_half = int(width / 2)
  if up_flag == false:
    result = ind - width_half
    result = if result > 0: result else: 0
  else:
    result = ind + width_half
    result = if result < (size - 1): result else: size - 1

proc map*[S, T, U](t: OrderedTable[S, T], op: proc(k: S, v: T): U {.closure.}):
                                                                 seq[U] {.inline.} =
  ## Returns a new sequence with the results of `op` applied to every item in
  ## the OrderedTable t
  newSeq(result, t.len)
  var i = 0
  for k, v in pairs(t):
    result[i] = op(k, v)
    inc i


template delByElement*[T](a: var seq[T], p: T) =
  # template to delete the given element p from a in place
  let ind = find(a, p)
  if ind != -1:
    del(a, ind)

proc deleteIntersection*[T](a: var seq[T], b: seq[T]) =
  # procedure to delete intersection of a and b in a
  # finds all elements of b in a and deletes them, effectively
  # reducing a to the difference: a - b
  for p in b:
    delByElement(a, p)

# proc arange*(start, stop, step: int): seq[int] =
#   result = @[]
#   for i in start..<stop:
#     if (i - start) mod step == 0:
#       result.add(i)

# proc linspace*(start, stop: float, num: int): seq[float] =
#   # linspace similar to numpy's linspace
#   result = @[]
#   var step = start
#   let diff = (stop - start) / float(num)
#   if diff < 0:
#     # in case start is bigger than stop, return an empty sequence
#     return @[]
#   else:
#     for i in 0..<num:
#       result.add(step)
#       # for every element calculate new value for next iteration
#       step += diff

proc boolFromArrayOfIndices*(array: seq[int]): seq[bool] =
  result = @[]
  for i in 0..<array[array.high]:
    if i in array:
      result.add(true)
    else:
      result.add(false)

proc getSubarrayByIndices*[T](array: seq[T], inds: seq[int]): seq[T] =
  # this is the verbose function for the `[]` below. Below does not call
  # this function to avoid function calling overhead
  result = map(inds, proc(ind: int): T = array[ind])

proc `[]`*[T](a: seq[T], inds: openArray[int]): seq[T] {.inline.} =
  ## given two openArrays, return a sequence of all elements whose indices
  ## are given in 'inds'
  ## inputs:
  ##    a: seq[T] = the sequence from which we take values
  ##    inds: openArray[int] = the array which contains the indices for the
  ##         arrays, which we take from 'array'
  ## outputs:
  ##    seq[T] = a sequence of all elements s.t. array[ind] in numpy indexing
  result = map(inds, (ind: int) -> T => a[ind])

# NOTE: not as easy to implement as I first thought. Leave it for now.
# proc `[]`*[T](a: var seq[T], inds: openArray[int]) {.inline.} =
#   ## same as `[]` above, but changes a in place
#   ## given two openArrays, return a sequence of all elements whose indices
#   ## are given in 'inds'
#   ## inputs:
#   ##    a: seq[T] = the sequence from which we take values
#   ##    inds: openArray[int] = the array which contains the indices for the
#   ##         arrays, which we take from 'array'
#   ## outputs:
#   ##    seq[T] = a sequence of all elements s.t. array[ind] in numpy indexing
#   if isSorted(inds) == false:
#     let i_sorted = sort(inds)
#   result = keepIf(a, (ind: int) -> T => a[ind])


proc mean*[T](array: seq[T]): float =
  # returns the mean of the array
  let n_elements: float = float(array.len)
  result = foldl(array, a + b)
  result = result / n_elements


proc getListOfFiles*(folder: string, regex = ""): seq[string] =
  # returns a list of files from folder
  # NOTE: see a unit test for getListOfFiles for a suitable regex to
  # get a list of data*.txt files in a run folder
  result = @[]
  if existsDir(folder) == false:
    return result
  var count = 0
  let reg = re(regex)
  for file in walkDirRec(folder):
    count = count + 1
    if match(file, reg):
      result.add(file)

proc sortInodeTable*(inode_table: var OrderedTable[int, string]) =
  # this procedure sorts the given inode table by inode, to provide faster read
  # speeds. It uses in place sorting!
  # inputs:
  #   inode_table: OrderedTable[int, string] which contains a filename together with the
  #     corresponding inode

  # this one liner uses the sorted() to create a sorted copy of inode_table. It uses
  # system.cmp[int] (the int comparison procedure) to compare the inodes
  sort(inode_table, proc(x, y: (int, string)): int = system.cmp[int](x[0], y[0]))


proc createInodeTable*(list_of_files: seq[string]): OrderedTable[int, string] =
  # returns a table containing filenames and Inode number
  result = initOrderedTable[int, string]()
  # functional way to create table of inodes
  result = toOrderedTable(
    map(
      list_of_files,
      # the following line is the normal way to write the line after it
      # proc(n: string): (string, int) = (n, int(getFileInfo(n).id.file))
      # read line as: take string 'n' and make a tuple (string, int) from it
      # by applying getFileInfo to n and taking the file id
      (n: string) -> (int, string) => (int(getFileInfo(n).id.file), n)
  ))
  # and normal procedural way to create list of inodes
  # for file in list_of_files:
  #   let ino = int(getFileInfo(file).id.file)
  #   result[file] = ino

proc createSortedInodeTable*(list_of_files: seq[string]): OrderedTable[int, string] =
  # convenience wrapper around creation and sorting of the Inode table from a list of
  # files
  result = createInodeTable(list_of_files)
  sortInodeTable(result)


proc untarFile*(filepath: string, outdir = ""): string =
  # this procedure extracts the given *.tar.gz file of a run to the folder, in
  # which it is located, by making a system call to tar -xzf
  # inputs:
  #   filepath: string = the absolute path to the run file to be extracted
  # outputs:
  #   on success: the path to the extracted run folder
  #   on failure: ""

  # if successful, it will return the path to the extracted run folder
  let (dir, name_tar, ext) = splitFile(filepath)
  # the name given by splitFile only removes the leading dot. Leaves us with
  # .tar, which we need to remove
  let name = split(name_tar, ".tar")[0]
  # # given the directory, make system call to tar and extract the folder
  let outdir = if outdir.len > 0: outdir else: dir

  let cmd_tar = "tar -xzf " & filepath & " --directory " & outdir

  echo "System call to tar:\n\t", cmd_tar
  var (x, y) = execCmdEx(cmd_tar)
  if y != 0:
    echo "Warning: the extraction failed with exit status: x = ", x, " y = ", y
    result = ""
  else:
    # in this case tar returned 0 (== success)
    # now that we have extracted the folder, get list of files in run folder
    result = joinPath(outdir, name)

proc removeFolder*(folderpath: string): bool =
  # this procedure removes the folder with the given path, by making a system call to
  # rm -rf
  # WARNING: IF YOU HAND A FOLDER TO THIS FUNCTION, IT WILL BE DELETED
  # inputs:
  #   folderpath: string = full path to the folder to be deleted
  # outputs:
  #   bool = true: successfully deleted folder
  #   bool = false: either problem during deletion or action stopped by user
  let cmd_rm = "rm -rf " & folderpath
  echo "System call to rm:\n\t", cmd_rm
  echo "Are you sure you want to perform this system call? The folder WILL be deleted! (y/N)"
  let cont = stdin.readLine()
  if cont in @["y", "Y"]:
    let (x, y) = execCmdEx(cmd_rm)
    if y == 0:
      echo "... removed folder"
      result = true
    else:
      echo "Warning: something went wrong during deletion of folder. x = ", x, " y =", y
      result = false
  else:
    result = false

proc echoBenchCounted*(count: var int,
                  t: var float,
                  modby = 500, msg = " files read.") =
  inc count
  if count mod modby == 0:
    echo $count & msg
    let t1 = epochTime()
    echo "Took ", t1 - t, " s"
    t = t1

proc echoCount*(count: var int, modby = 500, msg = " files read.") =
  inc count
  if count mod modby == 0:
    echo $count & msg

proc echoCounted*(count: int, modby = 500, msg = " files read.") =
  if count mod modby == 0:
    echo $count & msg

proc getDaysHoursMinutes*(dur: Duration): string =
  ## returns a string of
  ## n days HH:MM
  ## based on a `Duration`
  let
    days = dur.inDays
    hours = dur.inHours - convert(Days, Hours, dur.inDays)
    minutes = dur.inMinutes - convert(Days, Minutes, dur.inDays) - convert(Hours, Minutes, hours)
  result = &"{days} days {hours:02}:{minutes:02}"

template getDateSyntax*(): string =
  ## returns the default syntax used when echoing a `DateTime` or `Time` object
  ## to parse a thusly created string
  # "yyyy-MM-dd'T'HH-mm-sszzz"
  "yyyy-MM-dd\'T\'HH:mm:sszzz"

# import the arraymancer related procs only if the `-d:pure` flag is not
# set. This allows us to keep this dependency out, if desired
when not defined(pure):
  import arraymancer
  proc argmin*[T](a: AnyTensor[T], axis: int): int =
    let `min` = min(a)
    for i, x in a:
      if x == `min`:
        return i[axis]

  proc argmin*[T](a: AnyTensor[T]): int =
    # argmin for 1D tensors
    let `min` = min(a)
    for i, x in a:
      if x == `min`:
        return i[0]

  proc findArgOfLocalMin*[T](view: AnyTensor[T], current_ind: int): int {.inline.} =
    # procedure to find the argument in a given view of a tensor
    result = argmin(view) + current_ind

  proc findArgOfLocalMinAlt*[T](view: AnyTensor[T], current_ind: int): int {.inline.} =
    # old very slow implementation. Why so complicated?
    # procedure to find the argument in a given view of a tensor
    var r_min = 1000.0
    for i, x in view:
      if r_min > x:
        r_min = x
        result = i[0] + current_ind

  proc findPeaks*(t: Tensor[float], width: int, skip_non_outliers: bool = true): seq[int] =
    # NOTE: trying to use a generic for this function (for any datatype in the tensor)
    # fails due to a bug in arraymancer, issue #62:
    # https://github.com/mratsim/Arraymancer/issues/62
    # this procedure searches for peaks (or dips) in the given array, by scanning
    # the array for elements, which are larger than any of the surrounding
    # 'width' elements. Thus, width defines the 'resolution' and locality of
    # the peaks for which to scan
    # TODO: rework this so that it's faster

    result = @[]
    #let t = x.toTensor()
    let size = size(t)
    let steps = int(size / width)
    # the value of the minimum in a given interval
    var int_min = 0.0

    # we define the cut value (which determines skip_non_outliers) as values, which are
    # smaller than one sigma smaller than the mean value
    var cut_value: float
    if skip_non_outliers:
      cut_value = mean(t) - std(t)

    var
      ind, u, min_ind, min_range, max_range, min_from_min: int = 0

    for i in 0 ..< steps:
      ind = i * width
      u   = ind + width - 1
      let view = t[ind .. u]
      min_ind = findArgOfLocalMin(view, ind)
      # given the current minimum, check whether we skip this element
      # since it is outside our bounds of interest anyway
      if skip_non_outliers:
        int_min = t[min_ind]
        if int_min > cut_value:
          continue

      min_range = getNewBound(min_ind, width, size, false)
      max_range = getNewBound(min_ind, width, size, true)

      min_from_min = findArgOfLocalMin(t[min_range..max_range], min_range)
      if min_ind == min_from_min:
        result.add(min_ind)

proc splitSeq*[T, U](s: seq[seq[T]], dtype: typedesc[U]): (seq[U], seq[U]) =
  ## splits a (N, 2) nested seq into two seqs
  result[0] = newSeq[dtype](s.len)
  result[1] = newSeq[dtype](s.len)
  for i in 0..s.high:
    result[0][i] = s[i][0].U
    result[1][i] = s[i][1].U

when isMainModule:
  # unit test for a regex to check for
  import re
  let test_names = @["/home/schmidt/CastData/data/2017/DataRuns/Run_84_171108-17-49/data013062.txt",
                     "/home/schmidt/CastData/data/2017/DataRuns/Run_84_171108-17-49/data015665.txt-fadc"]

  let regex = r"^/([\w-_]+/)*data\d{6}\.txt$"
  assert match(test_names[0], re(regex)) == true
  assert match(test_names[1], re(regex)) == false
  echo "All unit tests for regex passed for data*.txt files passed."
