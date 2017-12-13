## this script performs the raw data manipulation of a given run or list of runs
## and outputs the resulting data to a HDF5 file
## steps which are performed
## - reading all data*.txt and data*.txt-fadc files and writing them to a
##   HDF5 file, one group for each run
## - calculating the occupancy of each run
## - calculating the num_pix / event histogram
## - caluclating the FADC signal depth / event histogram
## - calculate the ToT per pixel histogram

## standard lib
import os
import re
import sequtils, future
import algorithm
import tables
import times
import threadpool
import memfiles
import strutils
import docopt

# InGrid-module
import helper_functions
import tos_helper_functions
import reconstruction

# other modules  
import nimhdf5/H5nimtypes
import arraymancer

# global experimental pragma to use parallel: statement in readRawEventData()
{.experimental.}

const FILE_BUFSIZE = 30000


let doc = """
InGrid raw data manipulation.

Usage:
  raw_data_manipulation <folder> [options]
  raw_data_manipulation -h | --help
  raw_data_manipulation --version

Options:
  -h --help           Show this help
  --version           Show version.
  --run_type <type>   Select run type (Calib | Data)
"""


proc readRawEventData(run_folder: string): seq[FlowVar[ref Event]] =
  # given a run_folder it reads all event files (data<number>.txt) and returns
  # a sequence of FlowVars of references to Events, which store the raw
  # event data
  # NOTE: this procedure does the reading of the data in parallel, thanks to
  # using spawn
  let
    regex_tup = getRegexForEvents()
    t0 = cpuTime()
  var
    count = 0
    # get a sorted list of files, sorted by inode
    inode_files: seq[string] = getSortedListOfFiles(run_folder, EventSortType.inode)
  result = @[]     
  let n_files = len(inode_files)
  
  while len(inode_files) > 0:
    # variable to set last index to read to
    var ind_high = FILE_BUFSIZE
    if len(inode_files) < FILE_BUFSIZE:
      ind_high = len(inode_files) - 1
    # read files into buffer sequence
    let buf_seq = readListOfFiles(inode_files[0..ind_high], regex_tup)
      
    echo "... removing read elements from list"
    # sequtils.delete removes the element with ind_high as well!
    inode_files.delete(0, ind_high)
    echo "... and concating buffered sequence to result"
    result = concat(result, buf_seq)
    count += FILE_BUFSIZE
  echo "All files read. Number = " & $len(result)
  echo "Compared with starting files " & $n_files

proc processRawEventData(ch: seq[FlowVar[ref Event]]): ProcessedRun =
  # procedure to process the raw data read from the event files by readRawEventData
  # inputs:
  #    ch: seq[FlowVar[ref Event]] = seq of FlowVars of references to Event objects, which
  #        each store raw data of a single event.
  #        We read data of FlowVars and store it in normal events, perform calculations
  #        to obtain ToT per pixel, number of hits and occupancies of that data
  # outputs:
  #   ProcessedRun containing:
  #    events:    seq[Event] = the raw data from the seq of FlowVars saved in seq of Events
  #    tuple of:
  #      tot:  seq[seq[int]] = ToT values for each chip of Septemboard for run
  #      hits: seq[seq[int]] = number of hits for each chip of Septemboard fo run
  #      occ:    Tensor[int] = (7, 256, 256) tensor containing occupancies of all chips for
  #        this data.

  # variable to count number of processed files
  var
    count = 0
    # store ToT data of all events for each chip
    # Note: still contains a single seq for each event! Need to concat
    # these seqs at the end
    tot_run: seq[seq[seq[int]]] = newSeq[seq[seq[int]]](7)
    # store occupancy frames for each chip
    occ = zeros[int](7, 256, 256)
    # store number of hits for each chip
    hits = newSeq[seq[int]](7)
    # initialize the events sequence of result, since we add to this sequence
    # instead of copying ch to it!
    events = newSeq[Event](len(ch))

  # initialize empty sequences. Input to anonymous function is var
  # as we change each inner sequence in place with newSeq
  apply(tot_run, (x: var seq[seq[int]]) => newSeq[seq[int]](x, 0))
  apply(hits, (x: var seq[int]) => newSeq[int](x, 0))

  echo "starting to process events..."
  count = 0
  for i in 0..<ch.high:
  #for i in 0..<10000:
    let a: Event = (^ch[i])[]
    events.add(a)
    let chips = a.chips
    for c in chips:
      let
        num = c.chip.number
        pixels = c.pixels
      addPixelsToOccupancySeptem(occ, pixels, num)
      let tot_event = pixelsToTOT(pixels)
      tot_run[num].add(tot_event)
      let n_pix = len(tot_event)
      if n_pix > 0:
        # if the compiler flag (-d:CUT_ON_CENTER) is set, we cut all events, which are
        # in the center 4.5mm^2 square of the chip
        when defined(CUT_ON_CENTER):
          if isNearCenterOfChip(pixels) == true:
            hits[num].add(n_pix)
        else:
          hits[num].add(n_pix)            
    echoFilesCounted(count)

  # using concatenation, flatten the seq[seq[int]] into a seq[int] for each chip
  # in the run (currently tot_run is a seq[seq[seq[int]]]. Convert to seq[seq[int]]
  let tot = map(tot_run, (t: seq[seq[int]]) -> seq[int] => concat(t))

  result.events = events
  result.tots = tot
  result.hits = hits
  result.occupancies = occ

proc processSingleRun(run_folder: string): ProcessedRun =
  # this procedure performs the necessary manipulations of a single
  # run. This is the main part of the raw data manipulation
  # inputs:
  #     run_folder: string = the run folder (has to be one!, check with isTosRunFolder())
  #         to be processed

  # need to:
  # - create list of all data<number>.txt files in the folder
  #   - and corresponding -fadc files
  # - read event header for each file
  # -
  # read the raw event data into a seq of FlowVars
  let ch = readRawEventData(run_folder)

  # process the data read into seq of FlowVars, save as result
  result = processRawEventData(ch)

proc writeProcessedRunToH5(h5file_id: hid_t, run: ProcessedRun) =
  # this procedure writes the data from the processed run to a HDF5
  # (opened already) given by h5file_id
  # inputs:
  #   h5file_id: hid_t = the file id (integer) pointing to the writable HDF5 file
  #   run: ProcessedRun = a tuple of the processed run data
  discard
      
proc main() =

  # use the usage docstring to generate an CL argument table
  let args = docopt(doc)
  echo args
  
  let folder = $args["<folder>"]
  var run_type = $args["--run_type"]
  if isNil(run_type) == true:
    run_type = ""
    
  # first check whether given folder is valid run folder
  let (is_run_folder, contains_run_folder) = isTosRunFolder(folder)
  echo "Is run folder       : ", is_run_folder
  echo "Contains run folder : ", contains_run_folder
  
  # if we're dealing with a run folder, go straight to processSingleRun()
  if is_run_folder == true and contains_run_folder == false:
    let r = processSingleRun(folder)

    let a = squeeze(r.occupancies[2,_,_])
    dumpFrameToFile("tmp/frame.txt", a)

    echo "occupied memory so far $# \n\n" % [$getOccupiedMem()]
    GC_fullCollect()
    echo "occupied memory after gc $#" % [$getOccupiedMem()]        

    #dumpNumberOfInstances()
    # dump sequences to file



    # var h5file_id: hid_t = 0
    # h5file_id = H5Fcreate("run_file.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)
    # writeProcessedRunToH5(h5file_id, r)
    # H5Fclose( h5file_id )
    # let events = r.events
    # var count = 0
    # for ev in events:
    #   for c in ev.chips:
    #     let t = spawn findSimpleCluster(c.pixels)
    #     # if len(t) > 1:
    #     #   let tmp = createTensorFromZeroSuppressed(c.pixels)
    #     #   echo "chip ", c.chip, " found ", len(t), " clusters"          
    #     #   dumpFrameToFile("tmp/frame.txt", tmp)
    #     #   sleep(100)
    #   echoFilesCounted(count)
    # sync()
    echo "Size of total ProcessedRun object = ", sizeof(r)
    echo "Length of tots and hits for each chip"
    # dump sequences to file
    dumpToTandHits(folder, run_type, r.tots, r.hits)
        
  elif is_run_folder == false and contains_run_folder == true:
    # in this case loop over all folder again and call processSingleRun() for each
    # run folder

    var
      # create variables to store a combined sequence for all runs
      hits_all: seq[seq[int]] = newSeq[seq[int]](7)
      tots_all: seq[seq[int]] = newSeq[seq[int]](7)
    for kind, path in walkDir(folder):
      if kind == pcDir:
        # only support run folders, not nested run folders
        echo "occupied memory before run $# \n\n" % [$getOccupiedMem()]
        let
          is_rf = if isTosRunFolder(path) == (true, false): true else: false
          r = processSingleRun(path)
        for i in 0..<len(r.hits):
          hits_all[i] = concat(hits_all[i], r.hits[i])
          tots_all[i] = concat(tots_all[i], r.tots[i])
        #echo "Size of sequences : $# / $#" % [$(sizeof(hits_all) * len(hits_all[3]) * 7), $(sizeof(tots_all))]
        echo "occupied memory so far $# \n\n" % [$getOccupiedMem()]
        GC_fullCollect()
        echo "occupied memory after gc $#" % [$getOccupiedMem()]        

        #dumpNumberOfInstances()
    # dump sequences to file
    dumpToTandHits(folder, run_type, tots_all, hits_all)
        
  elif is_run_folder == true and contains_run_folder == true:
    echo "Currently not implemented to run over nested run folders."
    quit()
  else:
    echo "No run folder found in given path."
    quit()

when isMainModule:
  main()



# proc processRawEventDataParallel(ch: seq[FlowVar[ref Event]]): ProcessedRun =
#   # procedure to process the raw data read from the event files by readRawEventData
#   # NOTE: this is the parallel version of the function with the same name
#   # inputs:
#   #    ch: seq[FlowVar[ref Event]] = seq of FlowVars of references to Event objects, which
#   #        each store raw data of a single event.
#   #        We read data of FlowVars and store it in normal events, perform calculations
#   #        to obtain ToT per pixel, number of hits and occupancies of that data
#   # outputs:
#   #   ProcessedRun containing:
#   #    events: seq[Event] = the raw data from the seq of FlowVars saved in seq of Events
#   #    tuple of:
#   #      tot: seq[seq[int]] = ToT values for each chip of Septemboard for run
#   #      hits: seq[seq[int]] = number of hits for each chip of Septemboard fo run
#   #      occ: Tensor[int] = (7, 256, 256) tensor containing occupancies of all chips for
#   #        this data.

#   # variable to count number of processed files
#   var
#     count = 0
#     # store ToT data of all events for each chip
#     # Note: still contains a single seq for each event! Need to concat
#     # these seqs at the end
#     tot_run: seq[seq[seq[int]]] = newSeq[seq[seq[int]]](7)
#     # store occupancy frames for each chip
#     occ = zeros[int](7, 256, 256)
#     #occ: seq[Tensor[int]] = newSeq[Tensor[int]](7)
#     # store number of hits for each chip
#     hits = newSeq[seq[int]](7)
#     # initialize the events sequence of result, since we add to this sequence
#     # instead of copying ch to it!
#     events = newSeq[Event](len(ch))

#   # initialize empty sequences. Input to anonymous function is var
#   # as we change each inner sequence in place with newSeq
#   apply(tot_run, (x: var seq[seq[int]]) => newSeq[seq[int]](x, 0))
#   apply(hits, (x: var seq[int]) => newSeq[int](x, 0))
#   #map(occ, (x: var Tensor[int]) => zeros[int](256, 256))
#   # for i in 0..<occ.high:
#   #   occ[i] = zeros[int](256, 256)

#   # proc processChip(c: ChipEvent, o: var Tensor[int], tot: var seq[seq[int]], hit: var seq[int]) =
#   #   echo "1"
#   #   let pixels = c.pixels
#   #   echo "2"
#   #   addPixelsToOccupancy(o, pixels)
#   #   echo "3"
#   #   let tot_event = pixelsToTOT(pixels)
#   #   echo "4"
#   #   tot.add(tot_event)
#   #   echo "5"
#   #   let n_pix = len(tot_event)
#   #   echo "6"
#   #   if n_pix > 0:
#   #     hit.add(n_pix)
#   #   echo "7"
    
#   count = 0
#   for i in 0..<10000:
#     let a: Event = (^ch[i])[]
#     events.add(a)
#     let chips = a.chips
#     #let chip3 = filter(chips) do (ch: ChipEvent) -> bool : ch.chip.number == 3
#     for c in chips:
#       #parallel:
#       let num = c.chip.number
#       #processChip(c, occ[num], tot_run[num], hits[num])
#       let pixels = c.pixels
#       #spawn addPixelsToOccupancySeptem(occ, pixels, num)
#       addPixelsToOccupancySeptem(occ, pixels, num)
#       let tot_event = pixelsToTOT(pixels)
#       tot_run[num].add(tot_event)
#       let n_pix = len(tot_event)
#       if n_pix > 0:
#         hits[num].add(n_pix)
#       #sync()
#     #sync()
#     echoFilesCounted(count)
#     #echo n_pix

#   # using concatenation, flatten the seq[seq[int]] into a seq[int] for each chip
#   # in the run (currently tot_run is a seq[seq[seq[int]]]. Convert to seq[seq[int]]
#   let tot = map(tot_run, (t: seq[seq[int]]) -> seq[int] => concat(t))

#   # the data we have now represents our needed processed data for each run
#   # (excluding FADC)
#   # events: the raw event data (header plus hits)
#   # tot_run: for each chip seq of ToT values (excl 11810) to create
#   #   ToT per pixel histogram for each chip
#   # occupancies: for each chip tensor of ints containing number of hits of
#   #   each pixel for whole run
#   # hits: for each chip seq of number of hits (excl 11810 hits) to create
#   #   Hits per event for whole run (e.g. energy calibration)
#   result.events = events
#   result.tots = tot
#   result.hits = hits
#   result.occupancies = occ
  

  
