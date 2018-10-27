import arraymancer
import ingrid/ingrid_types

proc dumpFrameToFile*[T](filepath: string, ar: Tensor[T]) =
  ## this procedure dumps a given frame (tensor ar needs to be of shape (256, 256)
  ## to a file 'filepath'
  ## inputs:
  ##   filepath: string = the file to write to
  ##   ar: Tensor[int] = a tensor of shape (256, 256) containing the data to be written
  doAssert(ar.shape == [256, 256])
  var f: File
  if open(f, filepath, fmWrite):
    for x in 0..<256:
      for y in 0..<256:
        f.write($ar[x, y] & "\t")
      f.write("\n")
    f.close()
  else:
    echo "Warning: File to dump frame data could not be opened! Does the output folder exist? Path was ", filepath

template addPixelsToOccupancy*[T](ar: Tensor[T], pixels: Pixels) =
  ## template to add pixels to occupancy by using map
  for p in pixels:
    ar[p.x, p.y] += 1#p.ch

template addPixelsToOccupancySeptem*[T](ar: var Tensor[T], pixels: Pixels, ch_num: int) =
  ## template to add pixels to occupancy by using map
  for p in pixels:
    ar[ch_num, p.x.int, p.y.int] += 1

proc createTensorFromZeroSuppressed*[T](pixels: Pixels): Tensor[T] =
  ## procedure to create a (256, 256) int array from a Pixels (seq[tuple[x, y, ch]])
  ## object
  result = zeros[T](256, 256)
  for p in pixels:
    result[p.x, p.y] = p.ch

