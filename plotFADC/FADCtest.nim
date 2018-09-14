import nimhdf5, os, ingrid/tos_helpers, sequtils, strutils
import sets
import plotly
import tables
import re
import docopt


when defined(linux):
  const commitHash = staticExec("git rev-parse --short HEAD")
  const currentDate = staticExec("date")
else:
  const commitHash = ""
  const currentDate = ""

const docTmpl = """
InGrid/FADC reading and cuting.

Usage:
 FADCtest <HDF5file> [options]
 FADCtest <HDF5file> --runNumber <runnum> [options]
 FADCtest <HDF5file> --runNumber <runnum> --chipNumber <chipnum> [options]
 FADCtest <HDF5file> --chipNumber <chipnum> [options]

Options:
 -h --help show this text
 --runNumber      choose the run you would like to look at
 --chipNumber     choose the chip you would like to look at
 if not selected the program chooses the last values it gets from the file

 """

const doc = docTmpl % [commitHash, currentDate]

proc readandcut*(h5f: var H5FileObj, runNumber: int, chip: int) =

  ##get the ingrid data

  var reco_group = rawDatabase() & $runNumber
  var reco_group_chip = recoDataChipBase(runNumber) & $chip
  var group = h5f[reco_group.grp_str]
  var group_chip = h5f[reco_group_chip.grp_str]
  #echo groups

  var evNumbersDset = h5f[(group_chip.name / "eventNumber").dset_str]
  let evNumbers = evNumbersDset[int64]

  let chip_number = group_chip.attrs["chipNumber", int]
  var pos_x_dset = h5f[(group_chip.name / "centerX").dset_str]
  let pos_x = pos_x_dset[float64]
  var eccentricity_dset = h5f[(group_chip.name / "eccentricity").dset_str]
  let eccentricity = eccentricity_dset[float64]

  ##get the fadc data

  var fadc_eventNumber_name = eventNumberBasename(runNumber)
  var fadc_evNumbersDset = h5f[fadc_eventNumber_name.dset_str]
  let fadc_evNumbers = fadc_evNumbersDset[int64]
  #echo fadc_evNumbers

  var riseStart_name = riseStartBasename(runNumber)
  var fadc_riseStartDset = h5f[riseStart_name.dset_str]
  let fadc_riseStart = fadc_riseStartDset[uint16]

  var fallStop_name = fallStopBasename(runNumber)
  var fadc_fallStopDset = h5f[fallStop_name.dset_str]
  let fadc_fallStop = fadc_fallStopDset[uint16]

  var riseTime_name = riseTimeBasename(runNumber)
  var fadc_riseTimeDset = h5f[riseTime_name.dset_str]
  let fadc_riseTime = fadc_riseTimeDset[uint16]

  var fallTime_name = fallTimeBasename(runNumber)
  var fadc_fallTimeDset = h5f[fallTime_name.dset_str]
  let fadc_fallTime = fadc_fallTimeDset[uint16]
  #echo fadc_fallTime

  ##process data and make cuts

  let zipData = zip(evNumbers, eccentricity)
  #echo zipData

  let zipData_fadc = zip(fadc_evNumbers, fadc_fallTime)
  let cutData_fadc = zipData_fadc.filterIt(it[1] > 500'u16)
  let evcut_fadc = cutData_fadc.mapIt(it[0])
  let  evcut_fadc_set = toSet(evcut_fadc)

  let cutData = zipData.filterIt(it[0] in evcut_fadc_set)
  #echo cutData
  let cut = cutData.mapIt(it[1])
  #echo cut

  let
    d = Trace[float](`type`: PlotType.Histogram,
                     bins: (0.0, 15.0), binSize: 0.5)
  d.xs = cut
  let
    layout = Layout(title: "histogram in automatic range with specific bin range and size",
                    width: 1200, height: 800,
                    xaxis: Axis(title:"values"),
                    yaxis: Axis(title: "counts"),
                    autosize: false)
    p = Plot[float](layout: layout, traces: @[d])
  p.show()


proc main() =

  let args = docopt(doc)
  echo args

  let  h5file = $args["<HDF5file>"]
  var
    runNumflag = $args["--runNumber"]
    runNum = $args["<runnum>"]
    chipNumflag = $args["--chipNumber"]
    chipNum = $args["<chipnum>"]

  var h5f = H5file(h5file, "rw")
  h5f.visit_file
 # echo h5f
  var runN: int
  var runNumint: int
  var chipNumint: int
 #  const chip = 1
 # var chipNumber: int
  var dataBasename = recoBase()

  ##first idea to read name, run and chip from the shell
  ##now achieved with docopt
  #if paramCount() > 0:
   # h5file=paramStr(1)
   # runNum=paramStr(2)
   # chipNum=paramStr(3)
    #echo  h5file
   # echo  chipNum
  #else:
   # echo "missing file"
    #quit()

 # var runNumint = runNum.parseInt
 # var chipNumint = chipNum.parseInt

  ##get the runNumber and the chip number

  let groups = toSeq(keys(h5f.groups))
  let runRegex = re(data_basename & r"(\d+)$")

  var run: array[1, string]
  for grp in groups:
    if grp.match(runRegex, run) == true:
      # now read some data. Return value will be added later
      echo (run[0], grp)
  var runNumber = run[0].parseInt
  #echo groups

  var chipBasename = recoDataChipBase(runNumber)
  let chipRegex = re(chipBasename & r"(\d+)$")
  var chop: array[1, string]
  for grp in groups:
    if grp.match(chipRegex, chop) == true:
      echo(chop[0], grp)
  var chip = chop[0].parseInt
  #echo chip

  if runNum == "nil":
    runNumint = run[0].parseInt
  else:
    runNumint = runNum.parseInt
  if chipNum == "nil":
    chipNumint = chop[0].parseInt
  else:
    chipNumint = chipNum.parseInt
  readandcut(h5f, runNumint, chipNumint)

  discard h5f.close()


when isMainModule:
  main()
