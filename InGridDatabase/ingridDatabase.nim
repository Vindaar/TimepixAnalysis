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


  TotType = object
    pulses: int
    mean: float
    std: float
    
  Tot = object
    pulses: seq[int]
    mean: seq[float]
    std: seq[float]

  SCurve = object
    name: string
    voltage: int
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
  FsrPrefix = "fsr"
  FsrPattern = "fsr*.txt"
  TotPrefix = "TOTCalib"
  TotPattern = TotPrefix & r"*\.txt"
  ThresholdPrefix = "threshold"
  ThresholdPattern = r"threshold[0-9]\.txt"
  ThresholdMeansPrefix = ThresholdPrefix & "Means"
  ThresholdMeansPattern = ThresholdPrefix & r"Means*\.txt"
  SCurvePrefix = "voltage_"
  SCurveRegPrefix = r".*" & SCurvePrefix
  SCurvePattern = r"voltage_*\.txt"
  SCurveFolder = "SCurve/"
  StartTotRead = 20.0
let
  ChipNameLineReg = re(r"chipName:")
  ChipNameReg = re(r"([A-Z])\s*([0-9]+)\s*W\s*([0-9]{2})")
  FsrReg = re(FsrPrefix & r"([0-9])\.txt")
  FsrContentReg = re(r"(\w+)\s([0-9]+)")
  TotReg = re(TotPrefix & r"([0-9])\.txt")
  SCurveReg = re(SCurveRegPrefix & r"([0-9]+)\.txt")

proc `$`(chip: ChipName): string =
  result = $chip.col & $chip.row & " W" & $chip.wafer

proc parseChipInfo(filename: string): Chip =
  ## parses the `chipInfo` file and returns a chip object from that
  let info = readFile(filename).splitLines
  assert match(info[0], ChipNameLineReg), "The chipInfo file needs to contain the " &
    "`chipName: ` as the first line"
  let chipStr = info[0].split(':')[1].strip
  var mChip: array[3, string]
  # parse the chip name
  if match(chipStr, ChipNameReg, mChip) == true:
    result.name.col   = mChip[0][0]
    result.name.row   = mChip[1].parseInt
    result.name.wafer = mChip[2].parseInt
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
      let (_, thlFloat, countFloat) = readScurveVoltageFile(f)
      let scurve = SCurve(name: (SCurvePrefix & $voltage),
                          voltage: voltage,
                          thl: thlFloat.asType(int),
                          hits: countFloat.asType(int))
      result.curves.add scurve

proc parseTotFile(filename: string): Tot =
  ## checks for the existence of a TOT calibration file, reads it
  ## and returns a `Tot` object
  if existsFile(filename) == true:
    let (chip, pulses, mean, std) = readToTFile(filename,
                                                startRead = StartTotRead,
                                                totPrefix = TotPrefix)
    
    result.pulses = pulses.asType(int)
    result.mean = mean
    result.std = std

proc parseThresholdFile(filename: string): Threshold =
  ## parses a Threshold(Means) file and returns it
  echo filename
  let flat = readFile(filename).splitLines -->
    map(it.splitWhitespace) -->
    flatten() -->
    map(it.parseInt)
  result = flat.toTensor().reshape([256, 256])

#template chipGroup(chipName: string): string =
iterator pairs(t: FSR): (string, int) =
  for key, value in t:
    yield (key, value)

proc writeThreshold(h5f: var H5FileObj, threshold: Threshold, chipGroupName: string) =
  var thresholdDset = h5f.create_dataset(joinPath(chipGroupName, ThresholdPrefix),
                                          (256, 256),
                                          dtype = int)
  thresholdDset[thresholdDset.all] = threshold.data.reshape([256, 256])
  

proc addChipToH5(chip: Chip,
                 fsr: FSR,
                 scurves: SCurveSeq,
                 tot: Tot,
                 threshold: Threshold,
                 thresholdMeans: ThresholdMeans) =
  ## adds the given chip to the InGrid database H5 file

  var h5f = H5File(db, "rw")
  var chipGroup = h5f.create_group($chip.name)
  # add FSR, chipInfo to chipGroup attributes
  for key, value in chip.info:
    echo &"Appending {key} with {value}"
    chipGroup.attrs[key] = if value.len > 0: value else: "nil"
  for dac, value in fsr:
    chipGroup.attrs[dac] = value

  # SCurve groups
  if scurves.files.len > 0:
    var scurveGroup = h5f.create_group(joinPath(chipGroup.name, SCurveFolder))
    for i, f in scurves.files:
      # for each file write a dataset
      let curve = scurves.curves[i]
      var scurveDset = h5f.create_dataset(joinPath(scurveGroup.name,
                                                   curve.name),
                                          (curve.thl.len, 2),
                                          dtype = int)
      # reshape the data to be two columns of [thl, hits] pairs and write
      scurveDset[scurveDset.all] = zip(curve.thl, curve.hits).mapIt(@[it[0], it[1]])

  if tot.pulses.len > 0:
    # TODO: replace the TOT write by a compound data type using TotType
    # that allows us to easily name the columns too!
    echo sizeof(TotType)
    var totDset = h5f.create_dataset(joinPath(chipGroup.name, TotPrefix),
                                     (tot.pulses.len, 3),
                                     dtype = float)
    # sort of ugly conversion to 3 columns, using double zip
    # since we don't have a zip for more than 2 seqs
    totDset[totDset.all] = zip(zip(tot.pulses.asType(float),
                               tot.mean),
                               tot.std).mapIt(@[it[0][0],
                                                it[0][1],
                                                it[1]])

  if threshold.shape == @[256, 256]:
    h5f.writeThreshold(threshold, chipGroup.name)

  #if thresholdMeans.shape == @[256, 256]:
  #  h5f.writeThreshold(thresholdMeans, chipGroup.name)    
  
  let err = h5f.close()

        
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
    modifyStr = $args["--modify"]
    deleteStr = $args["--delete"]

  if modifyStr != "false":
    assert false, "Modify support is not implmented yet."
  elif deleteStr != "nil":
    assert false, "Delete support is not implmented yet."
  else:
    # call proc to add data from folder to database
    addChip(addStr)  

when isMainModule:
  main()
    
