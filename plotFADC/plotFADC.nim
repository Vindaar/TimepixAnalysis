import nimhdf5, os, ingrid/tos_helpers, sequtils, strutils, math
import sets
import plotly
import tables
import re
import docopt

let doc = """
InGrid/FADC reading and cutting.


Usage:
  plotFADC <HDF5file> [options]
  plotFADC <HDF5file> [--runNumber <runnum>] [--chipNumber <chipnum>] [options]
  plotFADC <HDF5file> --evParams  [options]
  plotFADC <HDF5file> [--choose <chparams>] [--cutFADC_low <fadccutlow>] [--cutFADC_high <fadccuthigh>] [options]
  plotFADC <HDF5file>  [--choose <chparams>] [--plot_low <plotlow>] [--plot_high <plothigh>] [options]
  plotFADC -h | --help


Options:
  --runNumber <runnum>          choose the run you would like to look at
  --chipNumber <chipnum>        choose the chip you would like to look at
  --evParams                    shows a list of possible InGrid Parameters
  --choose <chparams>           choose one InGrid Parameter from the list
  --cutFADC_low <fadccutlow>    choose a low limit for cuts
  --cutFADC_high <fadccuthigh>  choose a high limit for cuts
  --plot_low <plotlow>          choose a low limit for plot range
  --plot_high <plothigh>        choose a high limit for plot range
  -h --help                     show this text
 """

proc readandcut*[T](h5f: var H5FileObj, runNumber: int, chip: int, chParams: string, cutlow: int, cuthigh: int): seq[T] =

  ##get the ingrid data

  var reco_group = rawDatabase() & $runNumber
  var reco_group_chip = recoDataChipBase(runNumber) & $chip
  var group = h5f[reco_group.grp_str]
  var group_chip = h5f[reco_group_chip.grp_str]
  let chip_number = group_chip.attrs["chipNumber", int]
  let evNumbers = h5f[group_chip.name / "eventNumber", int64]
  let choosenParams = h5f[group_chip.name / chParams, float64]
  let a = choosenParams.filterIt(it < 1000)
  let b = a.filterIt(abs(it).classify != fcInf) ## filter for fcInf, to guarantee that all Parameters can be plotted


  ##get the fadc data

  var fadc_eventNumber_name = eventNumberBasename(runNumber)
  let fadc_evNumbers = h5f[fadc_eventNumber_name, int64]
  var fallTimename = fallTimeBasename(runNumber)
  let fadc_fallTime = h5f[fallTimename, uint16]

  ##process data and make cuts

  let zipData = zip(evNumbers, b)
  let zipData_fadc = zip(fadc_evNumbers, fadc_fallTime)
  let cutData_fadc = zipData_fadc.filterIt(it[1] > cutlow.uint16 and it[1] < cuthigh.uint16)

  let evcut_fadc = cutData_fadc.mapIt(it[0])
  let evcut_fadc_set = toSet(evcut_fadc)
  let cutData = zipData.filterIt(it[0] in evcut_fadc_set)
  let cuta = cutData.mapIt(it[0])
  let cutb = cutData.mapIt(it[1])

  result = cutb

proc plotcuts*[T](cuts: seq[T], chParams: string, maxcut: float) =
  ##plot the results
  let
    d = Trace[float](`type`: PlotType.Histogram,
                     bins: (0.0, maxcut), binSize: 0.1 ) ##somehow change bin and binsize in respect to the data

  d.xs = cuts
  let
    layout = Layout(title: "distribution histogram"  ,
                    width: 1200, height: 800,
                    xaxis: Axis(title: chParams),
                    yaxis: Axis(title:"counts"),
                    autosize: false)
    p = Plot[float](layout: layout, traces: @[d])
#  p.show()
#  p.saveImage("chParams.pdf")


proc main() =

  let args = docopt(doc)
  when not defined(release):
    echo args
  let  h5file = $args["<HDF5file>"]
  var
    runNum = $args["--runNumber"]
    chipNum = $args["--chipNumber"]
    chParams = $args["--choose"]
    cutFADCnum_low = $args["--cutFADC_low"]
    cutFADCnum_high = $args["--cutFADC_high"]
    plothigh_num = $args["--plot_high"]
    plotlow_num = $args["--plot_low"]

  var h5f = H5file(h5file, "rw")
  h5f.visit_file
  var runNumint: int
  var chipNumint: int
  var cuts: seq[float]
  var chParamstring: string
  var cutlow: int
  var cuthigh: int
  var dataBasename = recoBase()
  var maxcut: int
  var maxcutfloat: float

  if chParams != "nil":
    chParamstring = chParams
  else:
    chParamstring = "eccentricity"
    echo "When no InGrid Parameter is set, the eccentricity will be used"

  if cutFADCnum_low != "nil":
    cutlow = cutFADCnum_low.parseInt
  else:
    cutlow = 1
    echo "Standard low cut limit is 1"

  if cutFADCnum_high != "nil":
    cuthigh = cutFADCnum_high.parseInt
  else:
    cuthigh = 1000
    echo "Standard high cut limit is 1000"
  #echo cuthigh
  ##get the runNumber and the center chip number

  let groups = toSeq(keys(h5f.groups))
  var rawG = h5f["runs".grp_str]
  let center_chip = rawG.attrs["centerChip", int]
  if chipNum == "nil":
    chipNumint = center_chip
    echo "When no chip is set, the center chip will be used"
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
        cuts.add readandcut[float](h5f, runNumint, chipNumint, chParamstring, cutlow, cuthigh)
  else:
    runNumint = runNum.parseInt
    cuts = readandcut[float](h5f, runNumint, chipNumint, chParamstring, cutlow, cuthigh)

  if plothigh_num != "nil":
     maxcutfloat = plothigh_num.parseFloat
  else:
     maxcutfloat = cuts.max ##find the maximum to define the upper plotting limit
 # echo plothigh_num

  if $args["--evParams"] ==  "true":
    echo getFloatDsetNames()
  else:
    plotcuts(cuts, chParamstring, maxcutfloat)

  discard h5f.close()


when isMainModule:
  main()
