# this file contains the necessary code for reconstruction of X-rays in the raw data
# consists of
#   - finding clusters in data
#   - calculating properties of Events

# nim stdlib
import os
import sequtils, future
import threadpool
import nimnlopt
import math
import tables

# custom modules
import tos_helper_functions
import helper_functions


type
  # Coord type which contains (x, y) coordinates of a pixel
  Coord* = tuple[x, y: int]

  Cluster = seq[Pix]

const NPIX* = 256
const PITCH* = 0.055

type
  FitObject = object
    cluster: Cluster
    xy: tuple[x, y: float64]

proc excentricity(n: cuint, p: array[1, cdouble], grad: var array[1, cdouble], func_data: var pointer): cdouble {.cdecl.} =
  # this function calculates the excentricity of a found pixel cluster using nimnlopt.
  # Since no proper high level library is yet available, we need to pass a var pointer
  # of func_data, which contains the x and y arrays in which the data is stored, in
  # order to calculate the RMS variables

  # first recover the data from the pointer to func_data, by casting the
  # raw pointer to a Cluster object
  let fit = cast[FitObject](func_data)
  let c = fit.cluster
  let (x, y) = fit.xy

  var
    sum_x: cdouble = 0
    sum_y: cdouble = 0
    sum_x2: cdouble = 0
    sum_y2: cdouble = 0

  # echo "Starting new calc with param ", p[0]
  for i in 0..<len(c):
    let
      new_x = cos(p[0]) * (cdouble(c[i].x) - cdouble(x)) * PITCH - sin(p[0]) * (cdouble(c[i].y) - cdouble(y)) * PITCH
      new_y = sin(p[0]) * (cdouble(c[i].x) - cdouble(x)) * PITCH + cos(p[0]) * (cdouble(c[i].y) - cdouble(y)) * PITCH
    # echo "new_x is ", new_x, " and new_y ", new_y
    sum_x += new_x
    sum_y += new_y
    sum_x2 += (new_x * new_x)
    sum_y2 += (new_y * new_y)
  
  let
    n_elements: cdouble = cdouble(len(c))
    rms_x: cdouble = sqrt( (sum_x2 / n_elements) - (sum_x * sum_x / n_elements / n_elements))
    rms_y: cdouble = sqrt( (sum_y2 / n_elements) - (sum_y * sum_y / n_elements / n_elements))    

  #echo "sum_x2 / elements ", sum_x2 / n_elements, "sum_x * sum_x / elements / elements ", sum_x * sum_x / n_elements / n_elements
  let exc = rms_x / rms_y

  # need to check whether grad is nil. Only used for some algorithms, otherwise a
  # NULL pointer is handed in C
  if addr(grad) != nil:
    # normally we'd calculate the gradient for the current parameters, but we're
    # not going to use it. Can also remove this whole if statement
    discard

  #echo "parameters ", sum_x, " ", sum_y, " ", sum_x2, " ", sum_y2, " ", rms_x, " ", rms_y, " ", n_elements, " exc ", exc
  result = -exc
  

proc isPixInSearchRadius(p1, p2: Coord, search_r: int): bool =
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
  result = @[] #new seq[Cluster]
  
  let
    search_r = 50
    cutoff_size = 5

  # add the first pixel of the given sequence to have a starting pixel, from which we
  # look for other pixels in the cluster
  c.add(pixels[0])
  
  while raw_event.len > 0 and i < c.len:
    let p1: Coord = (x: c[i].x, y: c[i].y)
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

proc recoEvent(c: Cluster): float64 =
  let
    clustersize: int = len(c)
    (sum_x, sum_y, sumTotInCluster) = sum(c)
    # using map and sum[Pix] we can calculate sum of x^2, y^2 and x*y in one line
    (sum_x2, sum_y2, sum_xy) = sum(map(c, (p: Pix) -> Pix => (p.x * p.x,
                                                              p.y * p.y,
                                                              p.x * p.y)))
    pos_x: float64 = float64(sum_x) / float64(clustersize)
    pos_y: float64 = float64(sum_y) / float64(clustersize)
  var
    rms_x = sqrt(float64(sum_x2) / float64(clustersize) - pos_x * pos_x)
    rms_y = sqrt(float64(sum_y2) / float64(clustersize) - pos_y * pos_y)
    rotAngleEstimate = arctan( (float64(sum_xy) / float64(clustersize)) -
                               pos_x * pos_y / (rms_x * rms_x))

  if rotAngleEstimate < 0:
    echo "correcting 1"
    rotAngleEstimate += 8 * arctan(1.0)
  if rotAngleEstimate > 4 * arctan(1.0):
    echo "correcting 2"
    rotAngleEstimate -= 4 * arctan(1.0)
  elif classify(rotAngleEstimate) != fcNormal:
    echo "Rot angle estimate is NaN, vals are ", rms_x, " ", rms_y
    return NaN
  
  rms_x *= PITCH
  rms_y *= PITCH

  #echo "rot angle starting point is ", rotAngleEstimate

  var
    # define the optimizer
    opt: nlopt_opt
    # set the boundary values
    lb: array[2, cdouble] = [cdouble(-4 * arctan(1.0)), cdouble(4 * arctan(1.0))]
    # set the fit object with which we hand the necessary data to the
    # excentricity function
    fit_object: FitObject
    # the resulting fit parameter
    minf: cdouble
    p: array[1, cdouble] = [rotAngleEstimate]

  # set  the values of the fit objectn
  fit_object.cluster = c
  fit_object.xy = (x: pos_x, y: pos_y)
    
  #opt = nlopt_create(NLOPT_LN_COBYLA, 1)
  #opt = nlopt_create(NLOPT_LN_BOBYQA, 1)
  #opt = nlopt_create(NLOPT_LN_NELDERMEAD, 1)
  #opt = nlopt_create(NLOPT_LN_SBPLX, 1)
  # opt = nlopt_create(NLOPT_GN_DIRECT_L, 1)
  #opt = nlopt_create(NLOPT_GN_CRS2_LM, 1)
  #opt = nlopt_create(NLOPT_GN_ISRES, 1)
  opt = nlopt_create(NLOPT_GN_ESCH, 1)  

  # next one is useless, as dimensions needs n >= 2
  #opt = nlopt_create(NLOPT_LN_NEWUOA_BOUND, 1)
  #opt = nlopt_create(NLOPT_LN_SBPLX, 1)
  var status: nlopt_result
  status = nlopt_set_lower_bounds(opt, addr(lb[0]))
  status = nlopt_set_upper_bounds(opt, addr(lb[1]))
  if status != NLOPT_SUCCESS:
    echo "WARNING: could not set bounds!"
  # set the function to be minimized
  status = nlopt_set_min_objective(opt, cast[nlopt_func](excentricity), cast[pointer](addr fit_object))
  # set the minimization tolerance
  #status = nlopt_set_stopval(opt, cdouble(-Inf))
  status = nlopt_set_xtol_rel(opt, 1e-8)
  status = nlopt_set_ftol_rel(opt, 1e-8)
  status = nlopt_set_maxtime(opt, cdouble(1))

  var dx: cdouble = 2.0
  status = nlopt_set_initial_step(opt, addr(dx))
  
  # start minimization
  let t = nlopt_optimize(opt, addr(p[0]), addr(minf))
  if cast[int](t) < 0:
    echo t
    echo "nlopt failed!\n"
  else:
    echo t, " ;found minimum at f(", p[0], "), = ", minf #,"\n"

  nlopt_destroy(opt)

  result = float64(p[0])
    
proc reconstructSingleRun(folder: string) =
  # procedure which receives path to a run folder and reconstructs the objects
  # in that folder
  # inputs:
  #    folder: string = the run folder from which to reconstruct events

  echo "Starting to read list of files"
  let    
    files = getSortedListOfFiles(folder, EventSortType.inode)
    regex_tup = getRegexForEvents()
    f_to_read = if files.high < 30000: files.high else: 15000
    data = readListOfFiles(files[0..f_to_read], regex_tup)
  var
    min_val = 10.0
    min_seq = newSeq[seq[float64]](7)
  apply(min_seq, (x: seq[float64]) -> seq[float64] => @[])
  for e in data:
    let a: Event = (^e)[]
    let chips = a.chips
    for c in chips:
      let num: int = c.chip.number
      let cluster = findSimpleCluster(c.pixels)
      for cl in cluster:
        echo "Starting reco of ", a.evHeader["eventNumber"]
        let ob = recoEvent(cl)
        if ob < min_val:
          min_val = ob
        min_seq[num].add(ob)
  
        #sleep(100)
        #quit()
        #echo "Stuff... ", ob
  dumpRotAngle(min_seq)
  
  
proc main() =
  let args_count = paramCount()
  var folder: string
  if args_count < 1:
    echo "Please either hand a single run folder or a folder containing run folder, which to process."
    quit()
  else:
    folder = paramStr(1)
    
  # first check whether given folder is valid run folder
  let (is_run_folder, contains_run_folder) = isTosRunFolder(folder)
  echo "Is run folder       : ", is_run_folder
  echo "Contains run folder : ", contains_run_folder
  reconstructSingleRun(folder)

when isMainModule:
  main()
