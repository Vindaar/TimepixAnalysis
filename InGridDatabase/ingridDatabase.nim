import os
import docopt
import arraymancer
import re
import strutils

const doc = """
The InGrid database management tool.

Usage:
  ingridDatabase (--add=FOLDER | --modify | --delete=CHIPNAME) [options]

Options:
  --add=FOLDER        Add the chip contained in the given FOLDER 
                      to the database
  --modify            start a CLI interface to allow viewing the files
                      content and modify attributes / delete elements
  --delete=CHIPNAME   Delete contents of CHIPNAME from database
  -h, --help          Show this help
  --version           Show the version number

Documentation:
  This tool and library aims to provide two things. Running it as a main 
  module, it is a tool to add / manage the InGrid database HDF5 file, which
  stores information like FSR, thresholds, TOT calibration etc. for each chip
  with the corresponding name and location (e.g. Septem H).
  As a library it provides convenience functions to access the database to
  read and write to it, e.g.
  getTot(chipName)
  to get the ToT values of that chip.
"""

type
  ChipName = object
    col: char
    row: int
    wafer: int
    
  Chip = object
    name: ChipName
    info: Table[string, string]

  Tot = object
    pulses: seq[int]
    mean: seq[float]
    std: seq[float]

  SCurve = object
    thl: seq[int]
    hits: seq[int]
    
  SCurveSeq = object
    files: seq[string]
    curves: seq[SCurve]

  Threshold = Tensor[int]
  ThresholdMeans = Tensor[int]

  FSR = Table[string, int]

const db = "ingridDatabase.h5"
const
  ChipInfo = "chipInfo.txt"
  ChipNameReg = re(r"chipName:")
  FsrPrefix = "fsr"
  #FsrReg = re(FsrPrefix & r"([0-9])\.txt")
  TotPrefix = "TOTCalib"
  #TotReg = re(TotPrefix & r"([0-9])\.txt")
  SCurvePrefix = "voltage_"
  #SCurveReg = re(SCurvePrefix & r"([0-9]+)\.txt")

proc parseChipInfo(filename: string): Chip =
  ## parses the `chipInfo` file and returns a chip object from that
  let info = readFile(filename).splitLines
  assert match(info[0], ChipNameReg), "The chipInfo file needs to contain the " &
    "`chipName: ` as the first line"


proc addChip(folder: string) =
  ## handles adding the data for a chip from a folder
  ## it parses all available data in the folder and adds it to the correct group
  ## in the database file
  ## the given folder needs to conform to a few things to be parsed properly:
  ## Folder structure:\
  ## - chipInfo.txt
  ##   Contains:
  ##     - chipName: <name as given in TOS>
  ##     - arbitrary lines, which are individually added as string attributes for
  ##       this chips group. If lines have a colon, the string before that is the
  ##       name of the attribute
  ## Optional: if the following exists, it will be added
  ## - fsr?.txt
  ##   The FSR file of that chip
  ## - SCurve/
  ##   A folder containing the SCurves for that chip
  ## - TOTCalib?.txt
  ##   the TOT calibrration of that chip
  assert dirExists(folder), "Please hand an existing folder!"
  assert existsFile joinPath(folder, ChipInfo), "The given folder needs to contain a `chipInfo.txt`!"

  # read chip info file
  let t = parseChipInfo joinPath(folder, ChipInfo)
  

proc main() =

  let args = docopt(doc)
  echo args

  let
    addStr = $args["--add"]
    modifyStr = $args["--modify"]
    deleteStr = $args["--delete"]

  if modifyStr != "nil":
    assert false, "Modify support is not implmented yet."
  elif deleteStr != "nil":
    assert false, "Delete support is not implmented yet."
  else:
    # call proc to add data from folder to database
    addChip(addStr)  

when isMainModule:
  main()
    
