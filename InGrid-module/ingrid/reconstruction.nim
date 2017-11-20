# this file contains the necessary code for reconstruction of X-rays in the raw data
# consists of
#   - finding clusters in data
#   - calculating properties of Events




proc findSimpleCluster() =
  # this procedure searches for clusters based on a fixed search radius, whether
  # a pixel is found within that boundary, e.g. searchRadius = 50:
  # if within 50 pixels another pixel is found, add pixel to cluster, continue
  # inputs:
  #   -
  # ouputs:
  #   - 

  # - iterate over all hit pixels in event
  # - check next pixel, is it within search bound? yes,
  #   add pixel to hit
  var cluster: seq[tuple[x, y, ch: int]] = @[]
  var
    i = 0
    x_old = 0
    y_old = 0    
  for x, y, ch in pixels:
    if i > 0:
      x_old = pixels[i-1]
      y_old = pixels[i-1]
    else:
      cluster.add(pixels[i])
      continue
    if x_old + search_r > x and
       y_old - search_r < y:
      cluster.add(pixels[i])
    inc i
  
  

proc main() =
  discard
  
  

when isMainModule:
  main()
