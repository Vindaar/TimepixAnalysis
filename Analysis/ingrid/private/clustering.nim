#[
This file contains the code used to perform the clustering of the
`findSimpleCluster` procedure, i.e. our default clustering algorithm
(`caDefault`)
]#

import std / sets
import ingrid / [ingrid_types]

type
  PixelSearch*[T: SomePix] = object
    data*: ptr UncheckedArray[T]
    dataLen*: int
    when false: #T is Pix or T is PixTpx3:
      mask*: set[uint16]  ## This should maybe be a better structure!
    else: #elif T is PixInt:
      mask*: HashSet[int]
    #else: discard

  ClusterCandidate*[T] = object
    cid*: uint16 # id of this cluster
    data*: ptr UncheckedArray[T]
    dataLen*: int
    idxs*: seq[int] # indices in the raw event that are in this cluster
    merged*: bool # indicates cluster is merged into another cluster candidate
    mergeable*: set[uint16] # indicates the cluster IDs this one can be merged with


proc initPixelSearch[T](s: seq[T]): PixelSearch[T] =
  result = PixelSearch[T](data: cast[ptr UncheckedArray[T]](s[0].unsafeAddr), ## XXX: avoid copy
                          dataLen: s.len)
  when true:# T is PixInt:
    ## All indices as mask
    result.mask = toHashSet( toSeq( 0 ..< s.len ) )#initHashSet[int]() #
  ## no need to init `set[uint16]`

proc `[]`[T](p: PixelSearch[T], idx: int): T = p.data[idx]

func len[T](p: PixelSearch[T]): int = p.mask.card

iterator items[T](p: PixelSearch[T]): T =
  when false:
    for i in 0 ..< p.dataLen:
      when false: #T is Pix or T is PixTpx3:
        let idx = i.uint16
      else:
        let idx = i
      if idx notin p.mask:
        yield p.data[i]
  else:
    for idx in p.mask:
      yield p.data[idx]

iterator enumerate[T](p: PixelSearch[T]): (int, T) =
  when false:
    for i in 0 ..< p.dataLen:
      when false: #T is Pix or T is PixTpx3:
        let idx = i.uint16
      else:
        let idx = i
      if idx notin p.mask:
        yield (i, p.data[i])
  else:
    for idx in p.mask:
      yield (idx, p.data[idx])

proc next[T](p: var PixelSearch[T]): int =
  ## Returns the *index* of the `next` element in `p`. This essentially means the "first" element
  ## that is not yet masked.
  for i, x in enumerate(p): # items takes care of only looking at unmasked element
    when false:
      p.mask.incl i
      return i
    else:
      p.mask.excl i
      return i

template filter*[T](p: var PixelSearch[T], cond: untyped): untyped =
  ## return the indices of all those indices in `p` that pass the condition
  var res = newSeqOfCap[int](p.len div 2)
  for i, it {.inject.} in enumerate(p):
    if cond:
      res.add i
  res

proc maskIntersection[T](p: var PixelSearch[T], toMask: seq[int]) =
  ## Masks all elements in `toMask` in `p`
  when false:
    for m in toMask:
      when false: # T is Pix or T is PixTpx3:
        p.mask.incl m.uint16
      else: #elif T is PixInt:
        p.mask.incl m
  else:
    for m in toMask:
      p.mask.excl m

proc toCluster[T](pixels: seq[T], idxs: seq[int]): Cluster[T] =
  ## Converts the given pixels to a `Cluster`. For Timepix3 it also takes care
  ## of removing any pixels that are activated multiple times within a short
  ## time window (i.e. within this same ToA based cluster). We remove all
  ## those pixels that have the largest ToA value. Only the first is kept
  ## (according to the given order, which *should* {but check} be ToA ordered).
  when T is PixTpx3:
    result = newSeqOfCap[T](idxs.len)
    var pixMap = initHashSet[(uint8, uint8)]()
  else:
    result = newSeq[T](idxs.len)
  for i, idx in idxs:
    when T is PixTpx3: # in Tpx3 case filter out all pixels that appear multiple times in a
                       # short time, (multi ToA pixel filter)
      let p = pixels[idx]
      let pi = (x: p.x, y: p.y)
      if pi notin pixMap:
        pixMap.incl pi
        result.add p
    else:
      result[i] = pixels[idx]

proc isPixInSearchRadius[T: SomeInteger](p1, p2: Coord[T], search_r: int): bool =
  ## given two pixels, p1 and p2, we check whether p2 is within one square search
  ## of p1
  ## inputs:
  ##   p1: Pix = pixel from which to start search
  ##   p2: Pix = pixel for which to check
  ##   search_r: int = search radius (square) in which to check for p2 in (p1 V search_r)
  ## outpuits:
  ##   bool = true if within search_r
  ##          false if not

  # XXX: THIS searches in a ``*square*``. Add option to search in a ``*circle*``
  let
    # determine boundary of search space
    right = p1.x.int + search_r
    left  = p1.x.int - search_r
    up    = p1.y.int + search_r
    down  = p1.y.int - search_r
  # NOTE: for performance we may want to use the fact that we may know that
  # p1 is either to the left (if in the same row) or below (if not in same row)
  var
    in_x: bool = false
    in_y: bool = false

  if p2.x.int < right and p2.x.int > left:
    in_x = true
  if p2.y.int < up and p2.y.int > down:
    in_y = true
  result = if in_x == true and in_y == true: true else: false

proc initClusterCandidate[T](s: seq[T]): ClusterCandidate[T] =
  result = ClusterCandidate[T](data: cast[ptr UncheckedArray[T]](s[0].unsafeAddr),
                               dataLen: s.len,
                               idxs: newSeqOfCap[int](400))

proc add[T](c: var ClusterCandidate[T], idx: int) {.inline.} =
  c.idxs.add idx

proc add[T](c: var ClusterCandidate[T], idxs: seq[int]) {.inline.} =
  c.idxs.add idxs

proc setZero[T](c: var ClusterCandidate[T]) {.inline.} =
  c.idxs.setLen(0)

proc intersect[T](c1, c2: ClusterCandidate[T], searchRadius: int): bool =
  if c1.merged or c2.merged: return false ## shortcut if one is already merged, i.e. 'dead'
  for i in c1.idxs: # for each index in c1, check if any in reach of c2s
    for j in c2.idxs:
      let c1i = c1.data[i]
      let c2j = c2.data[j]
      if isPixInSearchRadius((x: c1i.x, y: c1i.y), (x: c2j.x, y: c2j.y), searchRadius):
        ## if *any* pixel is in radius, stop immediately
        return true

proc merge[T](c1, c2: var ClusterCandidate[T]) =
  ## Merges c2 into c1
  #doAssert not c1.merged and not c2.merged
  if c1.merged or c2.merged: return # nothing to do
  c1.add c2.idxs
  c2.merged = true

proc setMergeable[T](c1, c2: var ClusterCandidate[T]) =
  ## Assigns `c1` and `c2` to be mergeable
  c1.mergeable.incl c2.cid
  c2.mergeable.incl c1.cid

proc mergeChain[T](c: var ClusterCandidate[T], clusters: var seq[ClusterCandidate[T]],
                   ignoreIdx: set[uint16] = {}) =
  ## Merges all mergeable things into `c`
  var lastCard = 0#c.mergeable.card
  # 1. iteratively find all other cluster candidates that can be merged with this one
  #    by including all "mergeable" clusters into this clusters IDs until we the set
  #    of mergeable clusters doesn't change anymore
  while c.mergeable.card != lastCard:
    for id in c.mergeable:
      c.mergeable.incl clusters[id].mergeable
    lastCard = c.mergeable.card
  # 2. perform the merging of all clusters that are determined mergeable
  for id in c.mergeable:
    if id == c.cid: continue # don't merge itself
    var cMerge = clusters[id]
    doAssert cMerge.cid == id, "Index in clusters should be cluster id?"
    # merge this cluster
    c.merge(cMerge)
    # and set the modified back into the sequence
    ## XXX: this shouldn't be done via copying the data! It should use reference semantics
    ## probably!
    clusters[id] = cMerge


proc findSimpleCluster*[T: SomePix](pixels: seq[T], searchRadius: int): seq[Cluster[T]] =
  ## this procedure searches for clusters based on a fixed search radius, whether
  ## a pixel is found within that boundary, e.g. searchRadius = 50:
  ## if within 50 pixels another pixel is found, add pixel to cluster, continue
  ## inputs:
  ##   -
  ## ouputs:
  ##   -
  var
    # sequence in which to store the hits, which are part of a cluster
    c = initClusterCandidate[T](pixels)
    # create copy of pixels so that we can remove elements from it
    raw_event = initPixelSearch[T](pixels) #pixels
  when defined(onlySingleCluster):
    let
      cutoff_size = 1
  else:
    let
      cutoff_size = 2
  var clusterCands = newSeq[ClusterCandidate[T]]()
  var cid = 0'u16
  # 1. create proto clusters. Each cluster removes elements from `PixelSearch`. Each
  #    time we start from the next "available" pixel in `PixelSeach` and just search
  #    in the search radius.
  while raw_event.len > 0:
    let n = raw_event.next()       # get next pixel index
    let cAtIdx = raw_event[n]      # and related cluster
    c.add(n)                       # add to cluster candidate
    let t = raw_event.filter(      # get other pixel in search radius
      isPixInSearchRadius((x: cAtIdx.x, y: cAtIdx.y), (it.x, it.y), searchRadius)
    )
    maskIntersection(raw_event, t) # mask added pixels
    c.add t                        # add pixel in search radius
    c.cid = cid
    clusterCands.add c
    c.setZero()                    # reset candidate
    inc cid # cluster id
  # 2. with proto clusters in hand, check which clusters have an intersection. If so
  #    mark them as mergeable. This results in pairs of mergeable clusters
  var pairs = initHashSet[(int, int)]()
  for i, cOuter in mpairs(clusterCands):
    for j, cInner in mpairs(clusterCands):
      if i == j: continue
      if (min(i, j), max(i, j)) in pairs:
        continue # skip pairs already checked
      else:
        pairs.incl (min(i, j), max(i, j))
      if intersect(cOuter, cInner, searchRadius):
        cOuter.setMergeable cInner
  # 3. Merge all possible "chains" of clusters (A ⇔ B ⇔ C, but not A ⇔ C etc.)
  for c in mitems(clusterCands):
    c.mergeChain(clusterCands)
  # 4. finally convert the result back to real clusters as we're working with indices
  #    only here.
  ## XXX: performance improvement: *only* read and use the x/y information. Everything
  ## in addition is ballast. However, we barely touch the rest of the data, so probably
  ## irrelevant.
  for c in clusterCands:
    if not c.merged:
      if c.idxs.len > cutoffSize:
        result.add toCluster(pixels, c.idxs)

  when false:
    while raw_event.len > 0:
      let cAtIdx = raw_event[c[i]]
      let p1: Coord[type(cAtIdx.x)] = (x: cAtIdx.x, y: cAtIdx.y)
      # alternatively:
      #let t = raw_event.filterIt(isPixInSearchRadius(p1, (it.x, it.y), searchRadius))
      let t = raw_event.filter(isPixInSearchRadius(p1, (it.x, it.y), searchRadius))
      if isit:
        echo "index ", i, " rawlen : ", raw_event.len, " mask len ", raw_event.mask.len, " t len ", t.len, " at ", cAtIdx
      if t.len > 0:
        c.add t
        # remove elements from t in raw_event
        maskIntersection(raw_event, t)

      inc i

      if t.len > 0:
        c.add t
        # remove elements from t in raw_event
        #deleteIntersection(raw_event, t)
        maskIntersection(raw_event, t)
      if i == c.high or raw_event.len == 0:
        # at last index of this cluster (after possibly adding!) or
        # at last index in the event. In any case, add cluster to result
        if c.len > cutoff_size:
          result.add toCluster(pixels, c)
        # if shorter than dropoff, just ignored
        c.setLen(0)

        if raw_event.len > 0:
          ## Still indices left. Start a new cluster
          i = 0
          # add next element
          #let nextIdx = raw_event.next()
          c.add(raw_event.next()) #nextIdx
          #raw_event.maskIntersection(@[nextIdx])
          continue
      else: # for any other i we just continue
        discard
      inc i
