# this file contains the necessary code for reconstruction of X-rays in the raw data
# consists of
#   - finding clusters in data
#   - calculating properties of Events

# nim stdlib
import os
import sequtils, future
import threadpool
import math
import tables
import docopt
import strutils
import typetraits
import re

# external modules
import nlopt
import nlopt_wrapper
import nimhdf5

# custom modules
import tos_helper_functions
import helper_functions
import ingrid_types



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


let doc = """
InGrid reconstruction and energy calibration.

Usage:
  reconstruction <HDF5file> [options]
  reconstruction -h | --help
  reconstruction --version

Options:
  --run_number <number>  Only work on this run
  --out <name>           Filename and path of output file
  -h --help              Show this help
  --version              Show version.
"""
    

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
    pos_x = float64(sum_x) / float64(clustersize)
    pos_y = float64(sum_y) / float64(clustersize)
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

  var
    # set the boundary values corresponding to range of 360 deg
    lb = (-4.0 * arctan(1.0), 4.0 * arctan(1.0))
    # set the fit object with which we hand the necessary data to the
    # excentricity function
    fit_object: FitObject
    # the resulting fit parameter
    p = @[rotAngleEstimate]

  # set  the values of the fit objectn
  fit_object.cluster = c
  fit_object.xy = (x: pos_x, y: pos_y)
    
  
  var opt = newNloptOpt("LN_COBYLA", lb)
  opt.setFunction(excentricity, fit_object)

  # set relative precisions of x and y, as well as limit max time the algorithm
  # should take to 1 second
  opt.xtol_rel = 1e-8
  opt.ftol_rel = 1e-8
  opt.maxtime  = 1.0
  opt.initial_step = 0.02

  # start minimization
  let (params, min_val) = opt.optimize(p)
  if opt.status < NLOPT_SUCCESS:
    echo opt.status
    echo "nlopt failed!\n"
  else:
    echo opt.status, " ;found minimum at f(", params[0], "), = ", min_val

  nlopt_destroy(opt.optimizer)

  result = float64(params[0])

proc readDataFromH5(h5f: var H5FileObj, group: string, run_number: int) =
  # proc to read data from the HDF5 file from `group`
  var chip_base = rawDataChipBase(run_number)
  let raw_data_basename = rawDataBase()
  for grp in keys(h5f.groups):
    if chip_base in grp:
      # now can start reading
      var group = h5f[grp.grp_str]
      # get the chip number from the dataset
      let chip_number = group.attrs["chipNumber", int]
      # given group and chip number, we can read vlen data
      let vlen = special_type(int)
      var raw_x_dset  = h5f[(grp / "raw_x").dset_str]
      var raw_y_dset  = h5f[(grp / "raw_x").dset_str]
      var raw_ch_dset = h5f[(grp / "raw_x").dset_str]      
      let raw_x  = raw_x_dset[vlen, int]
      let raw_y  = raw_y_dset[vlen, int]
      let raw_ch = raw_ch_dset[vlen, int]      

      echo raw_x.len, " ", raw_y.len, " ", raw_ch.len
      
  
  
proc reconstructSingleRun(folder: string) =
  # procedure which receives path to a run folder and reconstructs the objects
  # in that folder
  # inputs:
  #    folder: string = the run folder from which to reconstruct events

  echo "Starting to read list of files"
  let    
    files = getSortedListOfFiles(folder, EventSortType.inode, EventType.InGridType)
    regex_tup = getRegexForEvents()
    f_to_read = if files.high < 30000: files.high else: 30000
    data = readListOfInGridFiles(files[0..f_to_read], regex_tup)
  var
    min_val = 10.0
    min_seq = newSeq[seq[float64]](7)
  apply(min_seq, (x: seq[float64]) -> seq[float64] => newSeqOfCap[float64](files.high))
  for e in data:
    let
      a: Event = (^e)[]
      chips = a.chips
    for c in chips:
      let
        num: int = c.chip.number
        cluster = findSimpleCluster(c.pixels)
      for cl in cluster:
        echo "Starting reco of ", a.evHeader["eventNumber"]
        let ob = recoEvent(cl)
        if ob < min_val:
          min_val = ob
        min_seq[num].add ob
  
  dumpRotAngle(min_seq)
  
  
proc main() =

  # use the usage docstring to generate an CL argument table
  let args = docopt(doc)
  echo args
  
  let h5f_name = $args["<HDF5file>"]
  var run_number = $args["--run_number"]
  var outfile = $args["--out"]
  if run_number == "nil":
    run_number = ""
  if outfile == "nil":
    echo "Currently not implemented. What even for?"
    outfile = "None"
    

  var h5f = H5file(h5f_name, "r")
  # visit the whole file to read which groups exist
  h5f.visitFile
  let raw_data_basename = rawDataBase()
  let run_regex = re(raw_data_basename & r"(\d+)$")
  if run_number == "":
    var run: array[1, string]
    for grp in keys(h5f.groups):
      if grp.match(run_regex, run) == true:
        # now read some data. Return value will be added later
        h5f.readDataFromH5(grp, parseInt(run[0]))
  else:
    h5f.readDataFromH5(raw_data_basename / run_number, parseInt(run_number))
  
  # NOTE: there's no point for this anymore, at least not at the moment
  # # first check whether given folder is valid run folder
  # let (is_run_folder, contains_run_folder) = isTosRunFolder(folder)
  # echo "Is run folder       : ", is_run_folder
  # echo "Contains run folder : ", contains_run_folder
  # reconstructSingleRun(folder)

when isMainModule:
  main()
