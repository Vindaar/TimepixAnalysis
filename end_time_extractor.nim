# this script is used to extract the end of the run time of a given 
# run. It is run from the command line as
# ./end_time_extractor <path to run .tar.gz>
# alternatively one can also simply hand a path to a run folder

import strutils
import os
import osproc
import helper_functions
import ingrid/tos_helper_functions
import times

proc main() = 
  let args_count = paramCount()
  if args_count < 1:
    echo "Please hand a .tar.gz file containing a compressed run or a run folder"
    quit()

  let input_file = paramStr(1)
  
  # first check whether the input really is a .tar.gz file
  let is_tar = ".tar.gz" in input_file

  # we need a variable for the run folder in order to deal with both cases
  # .tar.gz and run folder as input
  var run_folder = ""
  
  if is_tar:
    # in this case we need to extract the file to a temp directory
    run_folder = untarFile(input_file)
    if run_folder == nil:
      echo "Warning: Could not untar the run folder successfully. Exiting now."
      quit()
  else:
    run_folder = paramStr(1)

  # now we have a run folder, which we can work on
  let regex = r"^/([\w-_]+/)*data\d{6}\.txt$"
  let files = getListOfFiles(run_folder, regex)
  if len(files) == 0:
    echo "Either folder ", run_folder, " does not exist, or it is empty."
    quit()
  echo len(files)

  let rt_info = getRunTimeInfo(files)
  
  let parsed_first = formatAsOrgDate(rt_info.t_start)
  let parsed_last  = formatAsOrgDate(rt_info.t_end)      
  
  echo "Start of run:  <", parsed_first, ">"
  echo "End of run:    <", parsed_last, ">"
  echo "Length of run: ", rt_info.t_length.days, " days ", rt_info.t_length.hours, ":", rt_info.t_length.minutes
  
  if is_tar:
    # in this case we need to remove the temp files again
    # now that we have all information we needed from the run, we can delete the folder again
    let removed = removeFolder(run_folder)
    if removed == true:
      echo "Successfully removed all temporary files."
  

when isMainModule:
  main()


