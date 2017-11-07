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

