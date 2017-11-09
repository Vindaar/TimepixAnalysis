import strutils
import sequtils

# a simple collection of useful function for nim, mostly regarding arrays
# and sequences

proc arange*(start, stop, step: int): seq[int] = 
  result = @[]
  for i in start..<stop:
    if (i - start) mod step == 0:
      result.add(i)

proc boolFromArrayOfIndices*(array: seq[int]): seq[bool] = 
  result = @[]
  for i in 0..<array[array.high]:
    if i in array:
      result.add(true)
    else:
      result.add(false)

proc getSubarrayByIndices*[T](array: seq[T], inds: seq[int]): seq[T] = 
  result = map(inds, proc(ind: int): T = array[ind])

proc mean*[T](array: seq[T]): float =
  # returns the mean of the array
  let n_elements: float = float(array.len)
  result = foldl(array, a + b)
  result = result / n_elements


proc sortInodeTable*(inode_table: OrderedTable[int, string]): OrderedTable[int, string] =
  # this procedure sorts the given inode table by inode, to provide faster read
  # speeds
  # inputs:
  #   inode_table: OrderedTable[int, string] which contains a filename together with the
  #     corresponding inode

  # this one liner uses the sorted() to create a sorted copy of inode_table. It uses 
  # system.cmp[int] (the int comparison procedure) to compare the inodes
  result = sorted(inode_table, proc(x, y: (int, string)): int = system.cmp[int](x[0], y[0]))


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

