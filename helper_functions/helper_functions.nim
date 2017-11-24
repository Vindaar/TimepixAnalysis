import strutils
import sequtils
import tables
import algorithm
import future
import os
import osproc
import re

# a simple collection of useful function for nim, mostly regarding arrays
# and sequences

# proc sum*[T](s: seq[T]): T {.inline.} =
#   # this procedure sums the given array along the given axis
#   # if T is itself e.g. a tuple, we will return a tuple, one
#   # element for each field in the tuple
#   assert s.len > 0, "Can't sum empty sequences"
#   let n_fields = len(s[0])
#   var sum_t = T
#   for p in s:
#     for i in 0..<n_fields:
#       sum_t[i] += p[i]


# template orderedTableMap[S, T](result, t: typed) =
  
#   newSeq(result, t.len)
#   var i = 0
#   for k, v in pairs(t):
#     result[i] = op(k, v)
#     inc i
  

proc map*[S, T](t: OrderedTable[S, T], op: proc(k: S, v: T): S {.closure.}):
                                                                 seq[S] {.inline.} =
  ## Returns a new sequence with the results of `op` applied to every item in
  ## the OrderedTable t
  ## This proc expects a return value of the anonymous proc of the `key` type of
  ## the OrderedTable
  newSeq(result, t.len)
  var i = 0
  for k, v in pairs(t):
    result[i] = op(k, v)
    inc i

proc map*[S, T](t: OrderedTable[S, T], op: proc(k: S, v: T): T {.closure.}):
                                                                 seq[T] {.inline.} =
  ## Returns a new sequence with the results of `op` applied to every item in
  ## the OrderedTable t
  ## This proc expects a return value of the anonymous proc of the `value` type of
  ## the OrderedTable  
  newSeq(result, t.len)
  var i = 0
  for k, v in pairs(t):
    result[i] = op(k, v)
    inc i    

template delByElement[T](a: var seq[T], p: T) =
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

proc arange*(start, stop, step: int): seq[int] = 
  result = @[]
  for i in start..<stop:
    if (i - start) mod step == 0:
      result.add(i)

proc linspace*(start, stop: float, num: int): seq[float] =
  # linspace similar to numpy's linspace
  result = @[]
  var step = start
  let diff = (stop - start) / float(num)
  if diff < 0:
    # in case start is bigger than stop, return an empty sequence
    return @[]
  else:
    for i in 0..<num:
      result.add(step)
      # for every element calculate new value for next iteration
      step += diff

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
  for file in walkDirRec(folder):
    #if file.match re(regex):
    count = count + 1
    if match(file, re(regex)):
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


proc untarFile*(filepath: string): string =
  # this procedure extracts the given *.tar.gz file of a run to the folder, in
  # which it is located, by making a system call to tar -xzf
  # inputs:
  #   filepath: string = the absolute path to the run file to be extracted
  # outputs:
  #   on success: the path to the extracted run folder
  #   on failure: nil

  # if successful, it will return the path to the extracted run folder
  let (dir, name_tar, ext) = splitFile(filepath)
  # the name given by splitFile only removes the leading dot. Leaves us with
  # .tar, which we need to remove
  let name = split(name_tar, ".tar")[0]
  # # given the directory, make system call to tar and extract the folder
  let cmd_tar = "tar -xzf " & filepath & " --directory " & dir

  echo "System call to tar:\n\t", cmd_tar
  var (x, y) = execCmdEx(cmd_tar)
  if y != 0:
    echo "Warning: the extraction failed with exit status: x = ", x, " y = ", y
    return nil
  else:
    # in this case tar returned 0 (== success)
    # now that we have extracted the folder, get list of files in run folder
    result = joinPath(dir, name)


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
  
template echoFilesCounted*(count: int) =
  inc count
  if count mod 500 == 0:
    echo count, " files read."




when isMainModule:
  # unit test for a regex to check for 
  import re
  let test_names = @["/home/schmidt/CastData/data/2017/DataRuns/Run_84_171108-17-49/data013062.txt",
                     "/home/schmidt/CastData/data/2017/DataRuns/Run_84_171108-17-49/data015665.txt-fadc"]

  let regex = r"^/([\w-_]+/)*data\d{6}\.txt$"
  assert match(test_names[0], re(regex)) == true
  assert match(test_names[1], re(regex)) == false
  echo "All unit tests for regex passed for data*.txt files passed."
