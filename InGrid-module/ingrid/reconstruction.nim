# this file contains the necessary code for reconstruction of X-rays in the raw data
# consists of
#   - finding clusters in data
#   - calculating properties of Events

import tos_helper_functions
import helper_functions
import sequtils, future

type
  # Pix type which 
  Pix* = tuple[x, y: int]

  Cluster = seq[tuple[x, y, ch: int]]
  
  

proc isPixInSearchRadius(p1, p2: Pix, search_r: int): bool =
  # given two pixels, p1 and p2, we check whether p2 is within one square search
  # of p1
  # inputs:
  #   p1: Pix = pixel from which to start search
  #   p2: Pix = pixel for which to check
  #   search_r: int = search radius (square) in which to check for p2 in (p1 V search_r)
  # outpuits:
  #   bool = true if within search_r
  #          false if not
  let
    # determine boundary of search space
    right = p1.x + search_r
    left  = p1.x - search_r
    up    = p1.y + search_r
    down  = p1.y - search_r
  # NOTE: for performance we may want to use the fact that we may know that
  # p1 is either to the left (if in the same row) or below (if not in same row)
  var
    in_x: bool = false
    in_y: bool = false

  if p2.x < right and p2.x > left:
    in_x = true
  if p2.y < up and p2.y > down:
    in_y = true
  result = if in_x == true and in_y == true: true else: false

proc findSimpleCluster*(pixels: Pixels): seq[Cluster] =
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
  var
    # sequence in which to store the hits, which are part of a cluster
    c: Cluster = @[]
    # create copy of pixels so that we can remove elements from it
    raw_event = pixels
    # counter
    i = 0 
  result = @[]
  
  let
    search_r = 50
    cutoff_size = 3

  # add the first pixel of the given sequence to have a starting pixel, from which we
  # look for other pixels in the cluster
  c.add(pixels[0])
  
  while raw_event.len > 0 and i < c.len:
    let p1: Pix = (x: c[i].x, y: c[i].y)
    let t = filter(raw_event, (p: tuple[x, y, ch: int]) -> bool => isPixInSearchRadius(p1, (p.x, p.y), search_r))

    # add all found pixels to current cluster
    c = concat(c, t)
    # remove elements from t in raw_event
    deleteIntersection(raw_event, t)
    
    if i == c.len - 1 and raw_event.len > 0: 
      # if we are at the last hit pixel in the event, but raw_events is not empty,
      # this means there is more than 1 cluster in the frame. Thus, add the current
      # cluster 'c' to the seq of clusters and call this function recursively with
      # raw_events as the starting parameter
      if len(c) > cutoff_size:
        result.add(c)
      if len(raw_event) > cutoff_size:
        result.add(findSimpleCluster(raw_event))
    elif raw_event.len == 0 and len(c) > cutoff_size:
      result.add(c)
    inc i

proc main() =

  
  discard
  
  

when isMainModule:
  main()
