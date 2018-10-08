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
InGrid/FADC reading and cutting.

Usage:
 plotFADC <HDF5file> [options]
 plotFADC <HDF5file> --runNumber <runnum> [options]
 plotFADC <HDF5file> --runNumber <runnum> --chipNumber <chipnum> [options]
 plotFADC <HDF5file> --chipNumber <chipnum> [options]
 plotFADC <HDF5file> --evParams  [options]
 plotFADC <HDF5file> --choose <chparams> [options]

Options:
 --runNumber      choose the run you would like to look at
 --chipNumber     choose the chip you would like to look at
 --evParams       shows a list of Parameters to choose from
 --chParams       choose the Parameter from the list
 -h --help        show this text

 """

const doc = docTmpl % [commitHash, currentDate]

proc readandcut*[T](h5f: var H5FileObj, runNumber: int, chip: int, chParams: string): seq[T] =

  ##get the ingrid data

  var reco_group = rawDatabase() & $runNumber
  var reco_group_chip = recoDataChipBase(runNumber) & $chip
  var group = h5f[reco_group.grp_str]
  var group_chip = h5f[reco_group_chip.grp_str]
  let chip_number = group_chip.attrs["chipNumber", int]
  let evNumbers = h5f[group_chip.name / "eventNumber", int64]
  #let eccentricity = h5f[group_chip.name / "eccentricity", float64]
  let choosenParams = h5f[group_chip.name / chParams, float64]

  ##get the fadc data

  var fadc_eventNumber_name = eventNumberBasename(runNumber)
  let fadc_evNumbers = h5f[fadc_eventNumber_name, int64]
  var fallTimename = fallTimeBasename(runNumber)
  let fadc_fallTime = h5f[fallTimename, uint16]

  ##process data and make cuts

  let zipData = zip(evNumbers, choosenParams)
  let zipData_fadc = zip(fadc_evNumbers, fadc_fallTime)
  let cutData_fadc = zipData_fadc.filterIt(it[1] > 500'u16)
  let evcut_fadc = cutData_fadc.mapIt(it[0])
  let evcut_fadc_set = toSet(evcut_fadc)

  let cutData = zipData.filterIt(it[0] in evcut_fadc_set)
  let cuta = cutData.mapIt(it[0])
  let cutb = cutData.mapIt(it[1])
  result = cutb

proc plotcuts*[T](cuts: seq[T], chParams: string) =
  ##plot the results
  let
    d = Trace[float](`type`: PlotType.Histogram,
                    bins: (0.0, 60.0), binSize: 0.01 ) ##somehow change bin and binsize in respect to the data
  d.xs = cuts

  let
    layout = Layout(title: "distribution histogram of the choosen parameter",
                    width: 1200, height: 800,
                    xaxis: Axis(title: chParams),
                    yaxis: Axis(title:"counts"),
                    autosize: false)
    p = Plot[float](layout: layout, traces: @[d])
  p.show()
#  p.saveImage("test.pdf")


proc main() =

  let args = docopt(doc)
  let  h5file = $args["<HDF5file>"]
  var
    runNumflag = $args["--runNumber"]
    runNum = $args["<runnum>"]
    chipNumflag = $args["--chipNumber"]
    chipNum = $args["<chipnum>"]
    evParamsflag = $args["--evParams"]
    chParamsflag = $args["--choose"]
    chParams = $args["<chparams>"]

  var h5f = H5file(h5file, "rw")
  h5f.visit_file
  var runNumint: int
  var chipNumint: int
  var cuts: seq[float]
  var chParamstring: string
  var dataBasename = recoBase()

  if $args["--choose"] == "true":
    chParamstring = chParams
  else:
    chParamstring = "eccentricity"
    echo "When no Ingrid Parameter was choosen, the eccentricity will be used"

  ##get the runNumber and the center chip number

  let groups = toSeq(keys(h5f.groups))
  var rawG = h5f["runs".grp_str]
  let center_chip = rawG.attrs["centerChip", int]
  if chipNum == "nil":
    chipNumint = center_chip
  else:
    chipNumint = chipNum.parseInt

  let runRegex = re(data_basename & r"(\d+)$")
  var run: array[1, string]

  if runNum == "nil":
    for grp in groups:
      if grp.match(runRegex, run) == true:
        #echo (run[0], grp)
        var runNumber = run[0].parseInt
        runNumint = run[0].parseInt
        cuts.add readandcut[float](h5f, runNumint, chipNumint, chParamstring)
  else:
    runNumint = runNum.parseInt
    cuts = readandcut[float](h5f, runNumint, chipNumint, chParamstring)

  if $args["--evParams"] ==  "true":
    echo getFloatGeometryNames()
  else:
    plotcuts(cuts, chParamstring)

  discard h5f.close()


when isMainModule:
  main()
