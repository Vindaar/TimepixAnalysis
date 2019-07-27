## Simple script to compare the calculation of Marlin to create the
## XrayReferenceDataSet.h5 with the calculation introduced into the
## cdl_spectrum_creation.nim from the calibration-cdl.h5.

import nimhdf5
import plotly
import helpers / utils
import docopt, sequtils

const docStr = """
Usage:
  compareXrayRefMarlinTPA <RefMarlin> <RefTPA> [options]
  compareXrayRefMarlinTPA -h | --help
  compareXrayRefMarlinTPA --version

Options:
  -h, --help   Show this help
  --version    Show the version number
"""
const doc = withDocopt(docStr)


proc main =
  let args = docopt(doc)
  let marlin = $args["<RefMarlin>"]
  let tpa = $args["<RefTPA>"]

  # walk RefTPA, check if dataset found in RefMarlin, then read
  # both datasets and compare
  var h5tpa = H5file(tpa, "r")
  var h5marlin = H5file(marlin, "r")

  for grp in h5tpa:
    var mgrp = grp
    for dset in mgrp:
      if dset.name in h5marlin:
        # read both
        let
          dataTpa = h5tpa.readAs(dset.name, float).reshape2D(dset.shape).split(SplitSeq.Seq2Col)
          dataMarlin = h5marlin.readAs(dset.name, float).reshape2D(dset.shape).split(SplitSeq.Seq2Col)
        #echo dataTpa
        #echo dataMarlin
        barPlot(dataTpa[0], dataTpa[1]).show()
        barPlot(dataMarlin[0], dataMarlin[1]).show()


when isMainModule:
  main()
