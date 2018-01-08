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
import typetraits

# InGrid-module
import helper_functions
import tos_helper_functions
import reconstruction

# other modules  
import nimhdf5/H5nimtypes
import nimhdf5
import nimhdf5/hdf5_wrapper
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
    inode_files: seq[string] = getSortedListOfFiles(run_folder, EventSortType.fname)
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
  apply(tot_run, (x: var seq[seq[int]]) => newSeq[seq[int]](x, len(ch)))
  apply(hits, (x: var seq[int]) => newSeq[int](x, len(ch)))

  echo "starting to process events..."
  count = 0
  for i in 0..ch.high:
  #for i in 0..<10000:
    let a: Event = (^ch[i])[]
    events[i] = a
    let chips = a.chips
    for c in chips:
      let
        num = c.chip.number
        pixels = c.pixels
      addPixelsToOccupancySeptem(occ, pixels, num)
      let tot_event = pixelsToTOT(pixels)
      tot_run[num][i] = tot_event
      let n_pix = len(tot_event)
      if n_pix > 0:
        # if the compiler flag (-d:CUT_ON_CENTER) is set, we cut all events, which are
        # in the center 4.5mm^2 square of the chip
        when defined(CUT_ON_CENTER):
          if isNearCenterOfChip(pixels) == true:
            hits[num][i] = n_pix
        else:
          hits[num][i] = n_pix
        if n_pix > 4000:
          echo a.evHeader
          echo c.chip, " hits ", n_pix, " ", len(pixels)
    echoFilesCounted(count)


  # using concatenation, flatten the seq[seq[int]] into a seq[int] for each chip
  # in the run (currently tot_run is a seq[seq[seq[int]]]. Convert to seq[seq[int]]
  let tot = map(tot_run, (t: seq[seq[int]]) -> seq[int] => concat(t))

  # now fill the table in ProcessedRun tuple
  result.run_number = parseInt(events[0].evHeader["runNumber"])
  # use first event of run to fill event header. Fine, because event
  # header is contained in every file
  result.runHeader = fillRunHeader(^ch[0])
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

proc getGroupNameForRun(run_number: int): string =
  # generates the group name for a given run number
  result = "/runs/run_$#" % $run_number

proc writeProcessedRunToH5(h5f: var H5FileObj, run: ProcessedRun) =
  # this procedure writes the data from the processed run to a HDF5
  # (opened already) given by h5file_id
  # inputs:
  #   h5f: H5file = the H5 file object of the writeable HDF5 file
  #   run: ProcessedRun = a tuple of the processed run data
  const nchips = 7
  let nevents = run.events.len
  # easy to write:
  # occupancies
  # tots
  # hits
  # hard to write:
  # events
  #   for events we need to be able to write attributes (write data first, attributes later, not needed right now)
  #   but even so thanks to being zero suppressed, need vlen data
  
  # start by creating groups
  let group_name = getGroupNameForRun(run.run_number)
  var run_group = h5f.create_group(group_name)

  # TODO: write the run information into the meta data of the group
  echo "Create data to write to HDF5 file"
  let t0 = epochTime()
  # first write the raw data
  let
    ev_type = special_type(int)
    chip_group_name = group_name & "/chip_$#" #% $chp
    event_header_keys = ["eventNumber", "useHvFadc", "fadcReadout", "szint1ClockInt", "szint2ClockInt", "fadcTriggerClock"]
  var
    # create group for each chip
    chip_groups = mapIt(toSeq(0..<nchips), h5f.create_group(chip_group_name % $it))
    x_dsets  = mapIt(chip_groups, h5f.create_dataset(it.name & "/raw_x", nevents, dtype = ev_type))
    y_dsets  = mapIt(chip_groups, h5f.create_dataset(it.name & "/raw_y", nevents, dtype = ev_type))
    ch_dsets = mapIt(chip_groups, h5f.create_dataset(it.name & "/raw_ch", nevents, dtype = ev_type))

    # datasets to store the header information for each event
    evHeadersDsetTab = mapIt(event_header_keys, (it, h5f.create_dataset(group_name & "/" & it, nevents, dtype = int))).toTable
    
    evHeaders = initTable[string, seq[int]]()
    x  = newSeq[seq[seq[int]]](7)
    y  = newSeq[seq[seq[int]]](7)
    ch = newSeq[seq[seq[int]]](7)

  # prepare event header keys and value (init seqs)
  for key in event_header_keys:
    evHeaders[key] = newSeq[int](nevents)

  for i in 0..6:
    x[i]  = newSeq[seq[int]](nevents)
    y[i]  = newSeq[seq[int]](nevents)
    ch[i] = newSeq[seq[int]](nevents)
  for event in run.events:
    # given chip number, we can write the data of this event to the correct dataset
    let evNumber = parseInt(event.evHeader["eventNumber"])
    # add event header information
    for key in keys(evHeaders):
      evHeaders[key][evNumber] = parseInt(event.evHeader[key])
    # add raw chip pixel information
    for chp in event.chips:
      let
        num = chp.chip.number
      let hits = chp.pixels.len
      var
        xl: seq[int] = newSeq[int](hits)
        yl: seq[int] = newSeq[int](hits)
        chl: seq[int] = newSeq[int](hits)
      for i, p in chp.pixels:
        xl[i] = p[0]
        yl[i] = p[1]
        chl[i] = p[2]
      x[num][evNumber] = xl
      y[num][evNumber] = yl
      ch[num][evNumber] = chl

  echo "Writing all dset x data "
  let all = x_dsets[0].all
  for i in 0..x_dsets.high:
    echo "Writing dsets"
    x_dsets[i][all]  = x[i]
    y_dsets[i][all]  = y[i]
    ch_dsets[i][all] = ch[i]
  for key, dset in mpairs(evHeadersDsetTab):
    echo "Writing $# in $#" % [$key, $dset]
    dset[dset.all] = evHeaders[key]
    
  # alternative approach writing each event one after another.
  # however, this is about 9 times slower...
  # let
  #   ev_type = special_type(int)
  #   chip_group_name = group_name & "/chip_$#" #% $chp
  #   # create group for each chip 
  #   chip_groups = mapIt(toSeq(0..<nchips), h5f.create_group(chip_group_name % $it))
  # var
  #   x_dsets  = mapIt(chip_groups, h5f.create_dataset(it.name & "/raw_x", nevents, dtype = ev_type)) 
  #   y_dsets  = mapIt(chip_groups, h5f.create_dataset(it.name & "/raw_y", nevents, dtype = ev_type)) 
  #   ch_dsets = mapIt(chip_groups, h5f.create_dataset(it.name & "/raw_ch", nevents, dtype = ev_type))    
  # for event in run.events:
  #   for chp in event.chips:
  #     let
  #       num = chp.chip.number
  #     # given chip number, we can write the data of this event to the correct dataset
  #     let evNumber = parseInt(event.evHeader["eventNumber"])
  #     let hits = chp.pixels.len
  #     var
  #       xl: seq[int] = newSeq[int](hits)
  #       yl: seq[int] = newSeq[int](hits)
  #       chl: seq[int] = newSeq[int](hits)
  #     for i, p in chp.pixels:
  #       xl[i] = p[0]
  #       yl[i] = p[1]
  #       chl[i] = p[2]
  #     x_dsets[num].write(@[@[evNumber, 0]], xl)
  #     y_dsets[num].write(@[@[evNumber, 0]], yl)
  #     ch_dsets[num].write(@[@[evNumber, 0]], chl)
  echo "took a total of $# seconds" % $(epochTime() - t0)


  # create the datasets needed
  let reco_group_name = "/reconstruction/run_$#" % $run.run_number
  var reco_group = h5f.create_group(reco_group_name)

  # create hard links of header data to reco group
  for key, dset in mpairs(evHeadersDsetTab):
    echo "Writing $# in $#" % [$key, $dset]
    h5f.create_hardlink(dset.name, joinPath(reco_group_name, key))

  # now write attribute data (containing the event run header, for a start
  # NOTE: unfortunately we cannot write all of it simply using applyIt,
  # because we need to parse some numbers as ints, leave some as strings
  let asInt = ["runNumber", "runTime", "runTimeFrames", "numChips", "shutterTime",
               "runMode", "fastClock", "externalTrigger"]
  let asString = ["pathName", "dateTime", "shutterMode"]
  # write run header
  for it in asInt:
    run_group.attrs[it]  = parseInt(run.runHeader[it])
    reco_group.attrs[it] = parseInt(run.runHeader[it])
  for it in asString:
    run_group.attrs[it] = run.runHeader[it]
    reco_group.attrs[it] = run.runHeader[it]

  # write attributes for each chip
  var i = 0
  for grp in mitems(chip_groups):
    grp.attrs["chipNumber"] = run.events[0].chips[i].chip.number
    grp.attrs["chipName"]   = run.events[0].chips[i].chip.name
    inc i
  
  echo "ToTs shape is ", run.tots.shape
  echo "hits shape is ", run.hits.shape
  # into the reco group name we now write the ToT and Hits information
  # var tot_dset = h5f.create_dataset(reco_group & "/ToT")
  for chip in 0..run.tots.high:
    # since not every chip has hits on each event, we need to create one group
    # for each chip and store each chip's data in these
    let chip_group_name = reco_group_name & "/chip_$#" % $chip
    var chip_group = h5f.create_group(chip_group_name)
    # write chip name and number, taken from first event
    chip_group.attrs["chipNumber"] = run.events[0].chips[chip].chip.number
    chip_group.attrs["chipName"]   = run.events[0].chips[chip].chip.name    
    # in this chip group store:
    # - occupancy
    # - ToTs
    # - Hits
    # ... more later
    let
      tot = run.tots[chip]
      hit = run.hits[chip]
      occ = run.occupancies[chip, _, _].clone
    var
      tot_dset = h5f.create_dataset(chip_group_name & "/ToT", tot.len, int)
      hit_dset = h5f.create_dataset(chip_group_name & "/Hits", hit.len, int)
      occ_dset = h5f.create_dataset(chip_group_name & "/Occupancy", (256, 256), int)
    tot_dset[tot_dset.all] = tot
    hit_dset[hit_dset.all] = hit
    occ_dset[occ_dset.all] = occ
  
      
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
    # GC_fullCollect()
    # echo "occupied memory after gc $#" % [$getOccupiedMem()]        

    #dumpNumberOfInstances()
    # dump sequences to file

    # in order to write the processed run to file, open the HDF5 file
    var h5f = H5file("run_file.h5", "rw")

    # now write run to this file
    writeProcessedRunToH5(h5f, r)
    echo "Closing h5file with code ", h5f.close()
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
