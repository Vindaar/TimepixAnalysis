import os
import docopt
import arraymancer
import re
import strutils, strformat, sequtils
import ingrid/tos_helpers
import helpers/utils
import zero_functional
import seqmath
import nimhdf5
import times

import ingridDatabase/databaseDefinitions
import ingridDatabase/databaseRead
import ingridDatabase/databaseWrite
import ingridDatabase/databaseUtils

# cannot import `Chip` due to clash with this modules `Chip`
import ingrid/ingrid_types except Chip
import ingrid/calibration

const doc = """
The InGrid database management tool.

Usage:
  ingridDatabase (--add=FOLDER | --calibrate=TYPE | --modify | --delete=CHIPNAME) [options]

Options:
  --add=FOLDER        Add the chip contained in the given FOLDER 
                      to the database
  --calibrate=TYPE    Perform given calibration (TOT / SCurve) and save
                      fit parameters as attributes. Calibrates all chips
                      with suitable data.                      
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

proc calibrateType(tcKind: TypeCalibrate) =
  ## calibrates the given type and stores the fit results as attributes
  ## of each chip with valid data for calibration
  # open H5 file with write access
  var h5f = H5file(dbPath, "rw")
  case tcKind
  of TotCalibrate:
    # perform TOT calibration for all chips, which contain a TOT
    # calibration
    # - get groups in root of file == chips in file
    # - get TOT object of each chip
    # - perform Fit as in plotCalibration
    # - take `mpfit` result and save fit parameters as attrs
    for group in h5f.items(depth = 1):
      # group name is the chip name
      let chip = group.name
      echo "Getting TOT for chip ", chip
      let tot = getTotCalib(chip)
      let totCalib = fitToTCalib(tot, 0.0)
      # given fit result write attributes:
      h5f.writeTotCalibAttrs(chip, totCalib)
    
  of SCurveCalibrate:
    discard

  let err = h5f.close()
  if err < 0:
    echo "Could not properly close database! err = ", err


proc parseCalibrateType(typeStr: string): TypeCalibrate =
  ## checks the given `typeStr` whether it's a valid type to
  ## calibrate for. If so returns that type, else raises an
  ## exception
  case typeStr.toLower
  of $TotCalibrate:
    result = TotCalibrate
  of $SCurveCalibrate:
    result = ScurveCalibrate
  else:
    raise newException(KeyError, "Given type is not a valid type to calibrate " &
      """for. Valid types are ["TOT", "SCurve"] (case insensitive)""")
    
proc parseChipInfo(filename: string): Chip =
  ## parses the `chipInfo` file and returns a chip object from that
  let info = readFile(filename).splitLines
  assert match(info[0], ChipNameLineReg), "The chipInfo file needs to contain the " &
    "`chipName: ` as the first line"
  let chipStr = info[0].split(':')[1].strip
  # fill the `name` field
  result.name = parseChipName(chipStr)
  # parse the optional content of further notes on the chip
  # e.g. board, chip number on that board etc.
  result.info = initTable[string, string]()
  for line in info[1 .. info.high]:
    if line.len > 0:
      if ':' in line:
        let
          lineSplit = line.split(':')
          n = lineSplit[0].strip
          c = lineSplit[1].strip
        result.info[n] = c
      else:
        # support for lines without a key
        result.info[line.strip] = ""
    
proc parseFsr(filename: string): FSR =
  ## parses the contents of a (potential) given FSR file
  ## uses regex `FsrContentReg` to parse the file
  echo "Filename is ", filename
  if existsFile(filename) == true:
    let fLines = readFile(filename).splitLines
    var dacVal: array[2, string]
    result = initTable[string, int]()
    for line in fLines:
      if match(line.strip, FsrContentReg, dacVal):
        result[dacVal[0]] = dacVal[1].parseInt

proc parseScurveFolder(folder: string): SCurveSeq =
  ## parses all SCurve files (matching `SCurveReg`) in the folder
  ## and returns an `SCurveSeq` from it
  var scurveMatch: array[1, string]
  echo folder
  for f in walkFiles(joinPath(folder, SCurvePattern)):
    if match(f, SCurveReg, scurveMatch):
      result.files.add f
      let voltage = scurveMatch[0].parseInt
      let scurve = readScurveVoltageFile(f)
      result.curves.add scurve

proc parseTotFile(filename: string): Tot =
  ## checks for the existence of a TOT calibration file, reads it
  ## and returns a `Tot` object
  if existsFile(filename) == true:
    let (chip, tot) = readToTFile(filename,
                                  startRead = StartTotRead,
                                  totPrefix = TotPrefix)
    result.pulses = tot.pulses
    result.mean = tot.mean
    result.std = tot.std

proc parseThresholdFile(filename: string): Threshold =
  ## parses a Threshold(Means) file and returns it
  echo filename
  let flat = readFile(filename).splitLines -->
    map(it.splitWhitespace) -->
    flatten() -->
    map(it.parseInt)
  result = flat.toTensor().reshape([256, 256])

proc addChip(folder: string) =
  ## handles adding the data for a chip from a folder
  ## it parses all available data in the folder and adds it to the correct group
  ## in the database file
  ## the given folder needs to conform to a few things to be parsed properly:
  ## Folder structure:
  ## - chipInfo.txt
  ##   Contains:
  ##     - chipName: <name as given in TOS>
  ##     - arbitrary lines, which are individually added as string attributes for
  ##       this chips group. If lines have a colon, the string before that is the
  ##       name of the attribute
  ## Optional: if the following exists, it will be added
  ## - fsr?.txt
  ##   The FSR file of that chip
  ## - threshold?.txt
  ##   The threshold equalization of that chip
  ## - thresholdMeans?.txt                         
  ##   The mean values of threshold equalization of that chip
  ## - SCurve/
  ##   A folder containing the SCurves for that chip
  ## - TOTCalib?.txt
  ##   the TOT calibrration of that chip
  assert dirExists(folder), "Please hand an existing folder!"
  assert existsFile joinPath(folder, ChipInfo), "The given folder needs to contain a `chipInfo.txt`!"

  # read chip info file
  let chip = parseChipInfo joinPath(folder, ChipInfo)
  echo "Chip is ", chip

  # after chip check for fsr
  var fsrFileName = ""
  for f in walkFiles(joinPath(folder, FsrPattern)):
    echo f
    fsrFileName = f
  let fsr = parseFsr fsrFileName
  if fsr.len > 0:
    echo fsr
  
  # check for SCurves
  var scurves: SCurveSeq
  if dirExists(joinPath(folder, SCurveFolder)) == true:
    scurves = parseScurveFolder(joinPath(folder, SCurveFolder))

  # check for TOT
  var tot: Tot
  for f in walkFiles(joinPath(folder, TotPattern)):
    tot = parseTotFile(f)

  # check for threshold / threshold means
  var
    threshold: Threshold
    thresholdMeans: ThresholdMeans
  for f in walkFiles(joinPath(folder, ThresholdPattern)):
    threshold = parseThresholdFile(f)
  #for f in walkFiles(joinPath(folder, ThresholdMeansPattern)):
    #thresholdMeans = parseThresholdFile(f)
  #  assert false, "Threshold means writing not yet implemented!"

  # given all data, add it to the H5 file
  addChipToH5(chip, fsr, scurves, tot, threshold, thresholdMeans)

proc main() =

  let args = docopt(doc)
  echo args

  let
    addStr = $args["--add"]
    calibrateStr = $args["--calibrate"]
    modifyStr = $args["--modify"]
    deleteStr = $args["--delete"]

  if modifyStr != "false":
    assert false, "Modify support is not implmented yet."
  elif deleteStr != "nil":
    assert false, "Delete support is not implmented yet."
  elif calibrateStr != "nil":
    calibrateType calibrateStr.parseCalibrateType
  else:
    # call proc to add data from folder to database
    addChip(addStr)  

when isMainModule:
  main()
    
