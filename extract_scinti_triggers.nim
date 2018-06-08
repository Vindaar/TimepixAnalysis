# this script is used to extract the number of scintillator triggers in a given
# event set

import os
#import ingrid/tos_helper_functions
import tables
import sequtils
import re
import helper_functions
import strutils

proc readEventHeader*(filepath: string): Table[string, string] =
  # this procedure reads a whole event header and returns 
  # a table containing the data, where the key is the key from the data file
  result = initTable[string, string]()
  
  # we define a regex for the header of the file
  let regex = r"^\#\# (.*):\s(.*)"
  var matches: array[2, string]
  for line in lines filepath:
    if line.match(re(regex), matches):
      # get rid of whitespace and add to result
      let key = strip(matches[0])
      let val = strip(matches[1])
      result[key] = val

proc main() =
  let args_count = paramCount()
  if args_count < 1:
    echo "Please hand a run folder"
    quit()

  var scint1_hits = initTable[string, int]()
  var scint2_hits = initTable[string, int]()

  var input_folder = paramStr(1)

  # first check whether the input really is a .tar.gz file
  let is_tar = ".tar.gz" in input_folder

  if is_tar:
    # in this case we need to extract the file to a temp directory
    input_folder = untarFile(input_folder)
    if input_folder == nil:
      echo "Warning: Could not untar the run folder successfully. Exiting now."
      quit()

  # first check whether the input really is a valid folder
  if existsDir(input_folder) == true:
    # get the list of files in the folder
    #"/data/schmidt/data/2017/DataRuns/Run_84_171108-17-49/data001101.txt"
    let files = getListOfFiles(input_folder, r"^/?([\w-_]+/)*data\d{4,6}\.txt$")
    var inode_tab = createInodeTable(files)
    sortInodeTable(inode_tab)
    
    var count = 0
    for tup in pairs(inode_tab):
      let file = tup[1]

      let t = readEventHeader(file)
      let scint1 = parseInt(t["szint1ClockInt"])
      let scint2 = parseInt(t["szint2ClockInt"])
      let fadc_triggered = if parseInt(t["fadcReadout"]) == 1: true else: false
      # make sure we only read the scintillator counters, in case the fadc was 
      # actually read out. Otherwise it does not make sense (scintis not read out)
      # and the src/waitconditions bug causes overcounting
      if fadc_triggered:
        if scint1 != 0:
          scint1_hits[file] = scint1
        if scint2 != 0:
          scint2_hits[file] = scint2
      if count mod 500 == 0:
        echo count, " files read. Scint counters: 1 = ", len(scint1_hits), "; 2 = ", len(scint2_hits)
      count = count + 1
  else:
    echo "Input folder does not exist. Exiting..."
    quit()

  # all done, print some output
  echo "Reading of all files in folder ", input_folder, " finished."
  echo "\t Scint1     = ", len(scint1_hits)
  echo "\t Scint2     = ", len(scint2_hits)

  proc min_of_table(tab: Table[string, int]): tuple[min_val: int, file_min: string] =
    var min_val = 9999
    var file_min = ""
    for pair in pairs(tab):
      let val = pair[1]
      if val < min_val:
        min_val = val
        file_min = pair[0]
    result = (min_val, file_min)

  
  let unequal1 = toSeq(values(scint1_hits)).filterIt(it != 0 and it != 4095)
  let unequal2 = toSeq(values(scint2_hits)).filterIt(it != 0 and it != 4095)  

  echo "The values unequal to 0 for 1 ", unequal1
  echo "The values unequal to 0 for 2 ", unequal2

  echo "Number of unequal values for 1 ", unequal1.len
  echo "Number of unequal values for 2 ", unequal2.len
  
  let min_tup1 = min_of_table(scint1_hits)
  let min_tup2 = min_of_table(scint2_hits)

  echo "\t Scint1_min = ", min_tup1[0], " in file ", min_tup1[1]
  echo "\t Scint2_min = ", min_tup2[0], " in file ", min_tup2[1]

  # clean up after us, if desired
  if is_tar:
    # in this case we need to remove the temp files again
    # now that we have all information we needed from the run, we can delete the folder again
    let removed = removeFolder(input_folder)
    if removed == true:
      echo "Successfully removed all temporary files."
    

when isMainModule:
  main()
