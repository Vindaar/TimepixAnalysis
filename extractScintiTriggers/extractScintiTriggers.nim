# this script is used to extract the number of scintillator triggers in a given
# event set

import os
import ingrid/tos_helpers
import helpers/utils
import tables
import sequtils
import re
import strutils, strformat
import seqmath
import docopt
import algorithm
import loopfusion
import plotly
import nimhdf5

const doc = """
A simple tool to extract scintillator information from a run.

Usage:
  extractScintiTriggers <runFolder> [--outfile=FILE | --outfolder=FOLDER | --infolder=FOLDER]
  extractScintiTriggers <h5file> (--runNumber=NUMBER | --allRuns) [--outfile=FILE | --outfolder=FOLDER]

Options:
  --outfolder=FOLDER        Copy files of events with scintillator triggers
                            != 0 and != 4095 to this folder. If no argument
                            is given, a folder is created of the same name as
                            the run.
  --infolder=FOLDER         If an analysis was already done with this tool using the
                            outfolder option, hand that folder using infolder,
                            to only produce the plot.
  --runNumber=NUMBER        (h5file): If the input is a h5 file, hand the run for which
                            to extract scinti information
  --allRuns                 (h5file): use this flag to concat all runs in the h5 file.
  --outfile=FILE            If set we write the filenames of all interesting events
                            (see --outfolder) to this file [default: 123]
  -h --help                 Show this help
  --version                 Show version

Documentation:
  This tool extracts the number of scintillator triggers in a given run.
"""

# TODO: plot both scintillators in one window!

proc readEventHeader*(filepath: string): Table[string, string] =
  ## this procedure reads a whole event header and returns
  ## a table containing the data, where the key is the key from the data file
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

proc readRunFolder*(runFolder: string): (Table[string, int], Table[string, int]) =

  result[0] = initTable[string, int]()
  result[1] = initTable[string, int]()

  # first check whether the input really is a valid folder
  if existsDir(runFolder) == true:
    # get the list of files in the folder
    #"/data/schmidt/data/2017/DataRuns/Run_84_171108-17-49/data001101.txt"
    let files = getListOfFiles(runFolder, r"^/?([\w-_]+/)*data\d{4,6}\.txt$")
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
          result[0][file] = scint1
        if scint2 != 0:
          result[1][file] = scint2
      if count mod 500 == 0:
        echo count, " files read. Scint counters: 1 = ", len(result[0]), "; 2 = ", len(result[1])
      count = count + 1
  else:
    echo "Input folder does not exist. Exiting..."
    quit()

proc readScintFromH5(h5file: string, runNumber = 0, allRuns = false): (seq[int64], seq[int64]) =
  var h5f = H5file(h5file, "r")

  if allRuns == true:
    for num, grp in runs(h5f, rawDataBase()):
      var
        scint1_dset = h5f[(grp / "szint1ClockInt").dset_str]
        scint2_dset = h5f[(grp / "szint2ClockInt").dset_str]
      result[0] = concat(result[0], scint1_dset[int64])
      result[1] = concat(result[1], scint2_dset[int64])
  else:
    let grp = rawDataBase()
    var
      scint1_dset = h5f[(grp & $runNumber / "szint1ClockInt").dset_str]
      scint2_dset = h5f[(grp & $runNumber / "szint2ClockInt").dset_str]
    result[0] = scint1_dset[int64]
    result[1] = scint2_dset[int64]

  result[0] = filterIt(result[0], it != 0 and it != 4095)
  result[1] = filterIt(result[1], it != 0 and it != 4095)

proc getScintTrace[T](hist: seq[T], title = ""): Trace[T] =
  ## returns histogram trace suitable to plot scintillator triggers
  result = Trace[T](`type`: PlotType.Histogram,
                    bins: (0.0, 75.0), binSize: 1.0, histNorm: Percent)
  # filter out clock cycles larger 300 and assign to `Trace`
  result.xs = filterIt(hist, it < 300)
  result.name = title

proc plotHist*[T](traces: seq[Trace[T]]) =
  ## given a seq of scintillator counts, plot them as a histogram
  let
    goldenMean = (sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
    figWidth = 1200.0                     # width in inches
    figHeight = figWidth * goldenMean     # height in inches

  let
    layout = Layout(title: "Clock cycles since last Scinti trigger",
                    width: figWidth.int, height: figHeight.int,
                    xaxis: Axis(title: "# Clock cycles / 25 ns"),
                    yaxis: Axis(title: "# events"),
                    barmode: BarMode.Overlay,
                    autosize: false)
    p = Plot[int64](layout: layout, traces: traces)
  let filename = "scinti_triggers.svg"
  p.show(filename)

proc workOnRunFolder(rf: string): (seq[int64], seq[int64]) =
  var runFolder = rf

  let args = docopt(doc)

  var
    outfolder = $args["--outfolder"]
    outfile = ""
  if outfolder == "":
    outfolder = "true"
  elif outfolder == "nil":
    outfolder = ""

  if $args["--outfile"] != "nil":
    outfile = $args["--outfile"]

  # first check whether the input really is a .tar.gz file
  let is_tar = ".tar.gz" in runFolder

  if is_tar:
    # in this case we need to extract the file to a temp directory
    runFolder = untarFile(runFolder)
    if runFolder == nil:
      echo "Warning: Could not untar the run folder successfully. Exiting now."
      quit()

  # read the actual data
  let (scint1_hits, scint2_hits) = readRunFolder(runFolder)

  # all done, print some output
  echo "Reading of all files in folder ", runFolder, " finished."
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


  let unequal1_p = toSeq(pairs(scint1_hits)).filterIt(it[1] != 0 and it[1] != 4095)
  let unequal2_p = toSeq(pairs(scint2_hits)).filterIt(it[1] != 0 and it[1] != 4095)

  let unequal1 = unequal1_p.mapIt(it[1])
  let unequal2 = unequal2_p.mapIt(it[1])

  let unequal1_f = unequal1_p.mapIt(it[0])
  let unequal2_f = unequal2_p.mapIt(it[0])

  echo "The values unequal to 0 for 1 ", unequal1
  echo "The values unequal to 0 for 2 ", unequal2

  echo &"\t with a mean value of {unequal1.mapIt(it.float32).mean} for 1"
  echo &"\t with a mean value of {unequal2.mapIt(it.float32).mean} for 2"

  echo "Number of unequal values for 1 ", unequal1.len
  echo "Number of unequal values for 2 ", unequal2.len

  let min_tup1 = min_of_table(scint1_hits)
  let min_tup2 = min_of_table(scint2_hits)

  echo "\t Scint1_min = ", min_tup1[0], " in file ", min_tup1[1]
  echo "\t Scint2_min = ", min_tup2[0], " in file ", min_tup2[1]


  let files = concat(unequal1_f, unequal2_f).sorted(system.cmp[string])
  let filesFadc = mapIt(files, it & "-fadc")

  if outfolder.len > 0:
    # get filenames and write to file
    if outfolder == "true":
      # create the correct folder name
      outfolder = extractFilename(runFolder)
      echo "Creating the folder ", outfolder
      if existsOrCreateDir(outfolder):
        echo "Folder already existed"

    forZip f in files, ffadc in filesFadc:
      let fd = joinPath(outfolder, extractFilename(f))
      let fdFadc = joinPath(outfolder, extractFilename(ffadc))
      echo &"Copying file {f} to {fd}"
      copyFile(f, fd)
      copyFile(ffadc, fdFadc)

    # now write the scintillator trigger data to the folder as well
    var outScint1 = open("scint1.txt", fmWrite)
    defer: outScint1.close()
    var outScint2 = open("scint2.txt", fmWrite)
    defer: outScint2.close()
    outScint1.write(unequal1)
    outScint2.write(unequal2)

  if outfile.len > 0:
    var outf = open(outfile, fmWrite)
    defer: outf.close()
    forZip f in files, ffadc in filesFadc:
      outf.write(f & "\n")
      outf.write(ffadc & "\n")

  # clean up after us, if desired
  if is_tar:
    # in this case we need to remove the temp files again
    # now that we have all information we needed from the run, we can delete the folder again
    let removed = removeFolder(runFolder)
    if removed == true:
      echo "Successfully removed all temporary files."


  result[0] = unequal1.mapIt(it.int64)
  result[1] = unequal2.mapIt(it.int64)


proc workOnH5File(h5file: string): (seq[int64], seq[int64]) =

  let args = docopt(doc)

  let runNumberStr = $args["--runNumber"]
  var runNumber = 0
  if runNumberStr != "nil":
    runNumber = parseInt(runNumberStr)
  let allRunsStr = $args["--allRuns"]
  var allRuns = false
  if allRunsStr != "nil":
    allRuns = true

  result = readScintFromH5(h5file, runNumber, allRuns)

proc main() =

  let args = docopt(doc)
  echo args

  let infolder = $args["--infolder"]

  var useH5file = false

  var runFolder = ""
  var h5file = ""
  if $args["<runFolder>"] != "nil":
    runFolder = $args["<runFolder>"]
    useH5file = false
  else:
    h5file = $args["<h5file>"]
    useH5file = true

  var unequal1: seq[int64]
  var unequal2: seq[int64]
  if useH5file == false:
    (unequal1, unequal2) = workOnRunFolder(runFolder)
  else:
    (unequal1, unequal2) = workOnH5File(h5file)

  # now plot the scinti events
  let traceScint1 = getScintTrace(unequal1, "SiPM")
  let traceScint2 = getScintTrace(unequal2, "Veto scinti")
  plotHist(@[traceScint2, traceScint1])

when isMainModule:
  main()
