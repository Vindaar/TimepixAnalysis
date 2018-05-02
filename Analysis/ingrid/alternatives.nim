proc findSimpleClusterSet*(pixels: HashSet[Pix]): seq[HashSet[Pix]] = #seq[Cluster] =
  # a working implementation of cluster finder using hashsets instead of sequences
  # should be faster, due to union, exclusion etc being faster than for a sequence?
  # it's not though. So forget about it for now. Implementation is probably shitty for sets?
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
    c = initSet[Pix]()
    # create copy of pixels so that we can remove elements from it
    raw_event = pixels#.toSet
    # counter
    i = 0 
  result = @[] #new seq[Cluster]
  
  let
    search_r = 50
    cutoff_size = 5

  # add the first pixel of the given sequence to have a starting pixel, from which we
  # look for other pixels in the cluster
  var pix0: Pix
  for p in items(pixels):
    pix0 = p
    break
  c.incl(pix0)
  #raw_event.deleteIntersection(@[pixels[0]])
  raw_event.excl(pix0)
  
  while raw_event.card > 0 and i < c.card:
    var p1: Coord
    for p in c:
      p1.x = p.x
      p1.y = p.y
      break
    # alternatively:
    var t = initSet[Pix]()
    for p in items(raw_event):
      if isPixInSearchRadius(p1, (p.x, p.y), search_r) == true:
        t.incl(p)
    # add all found pixels to current cluster
    c = union(c, t)
    # remove elements from t in raw_event
    raw_event.excl(t)
    
    if i == c.card - 1 and raw_event.card > 0: 
      # if we are at the last hit pixel in the event, but raw_events is not empty,
      # this means there is more than 1 cluster in the frame. Thus, add the current
      # cluster 'c' to the seq of clusters and call this function recursively with
      # raw_events as the starting parameter
      if c.card > cutoff_size:
        result.add c
      if raw_event.card > cutoff_size:
        result.add(findSimpleClusterSet(raw_event))
    elif raw_event.card == 0 and c.card > cutoff_size:
      result.add c
    inc i
    
