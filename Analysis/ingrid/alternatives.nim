proc fitPolya*(charges, counts: seq[float], chipNumber: int): FitResult =
  ## proc to fit a polya distribution to the charge values of the
  ## reconstructed run. Called if `reconstruction` ran with --only_charge.
  ## After charge calc from TOT calib, this proc calculates the gas gain
  ## for this run.
  # determine start parameters
  # estimate of 3000 or gas gain
  # TODO: test again with mpfit. Possible to get it working?
  # start parameters
  let
    # typical gain
    gain = 3500.0 
    # start parameter for p[0] (scaling) is gas gain * max count value / 2.0
    # as good guess
    scaling = max(counts) * gain / 2.0
    # factor 23.0 found by trial and error! Works 
    rms = standardDeviation(counts) / 23.0
    # combine paramters, 3rd arg from `polya.C` ROOT script by Lucian
    p = @[scaling, gain, gain * gain / (rms * rms) - 1.0]

  # start parameters used to calc counts. Useful for debugging 
  let toFit = charges.mapIt(polyaImpl(p, it.float))
  let toPlot = toFit.mapIt(if classify(it) != fcInf: it else: 0)
  let trToFit = getTrace(charges, toPlot, "polya to fit")
    
  # data trace
  let trData = getTrace(charges,
                        counts.asType(float64),
                        &"polya data {chipNumber}",
                        PlotType.Bar)

  let chErrs = mapIt(toSeq(0 .. counts.high), 0.1) #counts.mapIt(sqrt(it))
  var parConf = newSeq[mp_par](3)
  for i in 0 .. 2:
    parConf[i].deriv_debug = 1
  # parConf[0].step = 10.0
  # parConf[1].step = 50.0
  # parConf[2].step = 1.0
  parConf[0].relstep = 1e-5
  parConf[1].relstep = 1e-5
  parConf[2].relstep = 1e-5
  
  echo "p is ", p
  echo "charges is ", charges.len
  echo "counts ", counts.len
  echo "chErr ", chErrs
  plotHist(@[tr, tr2])
  let (pRes, res) = fit(polya, p, charges, counts, chErrs, bounds = parConf)
  
  # echo parameters
  echoResult(pRes, res = res)
  
  let fitted = charges.mapIt(polya(pRes, it.float))
  let tr3 = getTrace(charges, fitted, "fit result")
  
  plotHist(@[tr, tr2, tr3])

  result.x = linspace(charges[0], charges[^1], 100)
  result.y = result.x.mapIt(polyaImpl(params, it))

  let tr3 = getTrace(result.x, result.y, "fit gas gain")
  
  plotGasGaing(@[tr, tr3], chipNumber)
  
  result.pRes = pRes
  result.pErr = res.error
  result.redChiSq = res.reducedChiSq
  

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
    
