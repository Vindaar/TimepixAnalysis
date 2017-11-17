## this script performs the raw data manipulation of a given run or list of runs
## and outputs the resulting data to a HDF5 file
## steps which are performed
## - reading all data*.txt and data*.txt-fadc files and writing them to a
##   HDF5 file, one group for each run
## - calculating the occupancy of each run
## - calculating the num_pix / event histogram
## - caluclating the FADC signal depth / event histogram
## - calculate the ToT per pixel histogram

import os
import re
import sequtils, future
import H5nimtypes

proc processSingleRun(run_folder: string, h5file_id: hid_t) =
  # this procedure performs the necessary manipulations of a single
  # run. This is the main part of the raw data manipulation
  # inputs:
  #     run_folder: string = the run folder (has to be one!, check with isTosRunFolder())
  #         to be processed
  #     h5file_id: hid_t   = the file id of the HDF5 file, to which we write the data

  # need to:
  # - create list of all data<number>.txt files in the folder
  #   - and corresponding -fadc files
  # - read event header for each file
  # - 
  

proc isTosRunFolder(folder: string): tuple[is_rf: bool, contains_rf: bool] =
  # this procedure checks whether the given folder is a valid run folder of
  # TOS
  # done by
  # - checking whether the name of the folder is a valid name for a
  #   run folder (contains Run_<number>) in the name and 
  # - checking whether folder contains data<number>.txt files
  # inputs:
  #    folder: string = the given name of the folder to check
  # outputs:
  # returns a tuple which not only says whether it is a run folder, but also
  # whether the folder itself contains a run folder
  #    tuple[bool, bool]:
  #        is_rf:       is a run folder
  #        contains_rf: contains run folders
  let run_regex = r".*Run_(\d+)_.*"
  let event_regex = r".*data\d{4,6}\.txt$"
  var matches_rf_name: bool = false
  if match(folder, re(run_regex)) == true:
    # set matches run folder flag to true, is checked when we find
    # a data<number>.txt file in the folder, so that we do not think a
    # folder with a single data<number>.txt file is a run folder
    matches_rf_name = true
    
  for kind, path in walkDir(folder):
    if kind == pcFile:
      if match(path, re(event_regex)) == true and matches_rf_name == true:
        result.is_rf = true
        # in case we found an event in the folder, we might want to stop the
        # search, in order not to waste time. Nested run folders are
        # undesireable anyways
        # for now we leave this comment here, since it may come in handy
        # break
    else:
      # else we deal with a folder. call this function recuresively
      let (is_rf, contains_rf) = isTosRunFolder(path)
      # if the underlying folder contains an event file, this folder thus
      # contains a run folder
      if is_rf == true:
        result.contains_rf = true
  
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

  # if we're dealing with a run folder, go straight to processSingleRun()
  if is_run_folder == true and contains_run_folder == false:
    processSingleRun(folder)
  elif is_run_folder == false and contains_run_folder == true:
    # in this case loop over all folder again and call processSingleRun() for each
    # run folder
    for kind, path in walkDir(folder):
      if kind == pcDir:
        # only support run folders, not nested run folders
        let is_rf = if isTosRunFolder(path) == (true, false): true else: false
        processSingleRun(path)
  elif is_run_folder == true and contains_run_folder == true:
    echo "Currently not implemented to run over nested run folders."
    quit()
  else:
    echo "No run folder found in given path."
    quit()

when isMainModule:
  main()
