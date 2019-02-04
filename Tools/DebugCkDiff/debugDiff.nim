import nimhdf5, zero_functional
import helpers / utils
import ingrid / [tos_helpers, ingrid_types]
import sequtils, strutils, plotly, strformat, os, sets, typetraits
import docopt

const docStr = """
A tool to debug the difference between Christoph's Marlin analysis and
TimepixAnalysis.

Usage:
  debugDiff <h5file> (--extract | --compare <h5file2>) [options]

Options:
  --extract              If given we will extract the raw data from the CDL file
  --compare <h5file2>    If this flag is given, we will compare the CDL file with
                         our reconstruction of it.
  -h, --help             Show this help
  -v, --version          Show the version number

The tool reads each group from the given H5 file, which has to be the
`calibration-cdl.h5`, because it contains the raw data for Christoph's
analysis as well as the calculated geometrical values. After reading
the raw data values it calculates all other properties and compares them
to those calc'd by Marlin.
"""
const doc = withDocopt(docStr)

proc writeAttrs(h5f: var H5FileObj,
                baseName, groupName: string,
                runNumber: int) =
  let chipGroup = h5f.create_group(&"{baseName}/run_{runNumber}/chip_0")
  let recoGroup = h5f.create_group(&"{baseName}")
  let runGroup = h5f.create_group(&"{baseName}/run_{runNumber}")
  runGroup.attrs["CDL_name"] = groupName
  recoGroup.attrs["centerChip"] = 0
  recoGroup.attrs["centerChipName"] = OldChipName
  runGroup.attrs["centerChip"] = 0
  runGroup.attrs["centerChipName"] = OldChipName
  runGroup.attrs["runNumber"] = runNumber
  runGroup.attrs["numChips"] = 1
  recoGroup.attrs["runFolderKind"] = $rfOldTos
  recoGroup.attrs["runRunType"] = $rtCalibration
  chipGroup.attrs["chipName"] = OldChipName
  chipGroup.attrs["chipNumber"] = 0

proc extractFromGroup(h5f: var H5FileObj, groupName: string) =
  ## read X, Y and charge from this group and write to a new H5 file matching data
  ## needed for reconstruction
  ## NOTE: runNumber is a dummy number assigned from iteration
  let tab = { "XCoordinatesVector" : "raw_x",
              "YCoordinatesVector" : "raw_y",
              "ChargeValuesVector" : "raw_ch",
              "EventNumber" : "eventNumber",
              "Timestamp" : "timestamp" }.toTable()

  var h5out = H5File("raw_from_cdl.h5", "rw")

  # determine run number from run number dataset
  let runNumberDset = h5f[(groupName / "RunNumber").dset_str]
  let runNumber = runNumberDset[0, float32].round.int

  let chipGroup = h5out.create_group(&"/runs/run_{runNumber}/chip_0")
  h5out.writeAttrs("/runs", groupName, runNumber)
  h5out.writeAttrs("/reconstruction", groupName, runNumber)

  let filter = H5Filter(kind: fkZlib, zlibLevel: 4)
  var grp = h5f[groupName.grpStr]

  let
    # create datatypes for variable length data
    ev_type_xy = special_type(uint8)
    ev_type_ch = special_type(uint16)

  template readWrite(groupRead, groupWrite: string,
                     dtypeI: typedesc,
                     dtypeF: untyped): untyped =
    when dtypeI is uint8:
      let data = h5f[groupRead, ev_type_xy, dtypeI]
    elif dtypeI is uint16:
      let data = h5f[groupRead, ev_type_xy, dtypeI]
    else:
      let data = h5f[groupRead, dtypeI]
    let dset = h5out.create_dataset(groupWrite,
                                    (data.len, 1),
                                    dtype = dtypeF,
                                    chunksize = @[data.len, 1],
                                    maxshape = @[data.len, 1],
                                    filter = filter)
    when dtypeF is int64:
      dset[dset.all] = data.mapIt(it.round.dtypeF)
    else:
      dset[dset.all] = data

  for initial, final in tab:
    if "Coordinates" in initial:
      readWrite(groupName / initial, chipGroup.name / final, uint8, ev_type_xy)
    elif "Charge" in initial:
      readWrite(groupName / initial, chipGroup.name / final, uint16, ev_type_ch)
    else:
      readWrite(groupName / initial, chipGroup.name.parentDir / final, float32, int64)
  discard h5out.close()

proc compareRun(h5M, h5Tp: var H5FileObj, grpM, grpTp: string) =
  ## creates all plots for the given run
  let tab = { # "NumberOfPixels"              : "hits",
              "Width"                       : "width",
              "SkewnessTransverse"          : "skewnessTransverse",
              "KurtosisLongitudinal"        : "kurtosisLongitudinal",
              "Length"                      : "length",
              "FractionWithinRmsTransverse" : "fractionInTransverseRms",
              "PositionX"                   : "centerX",
              "RmsLongitudinal"             : "rmsLongitudinal",
              "KurtosisTransverse"          : "kurtosisTransverse",
              #"EnergyFromCharge"            : "energyFromCharge",
              "RotationAngle"               : "rotationAngle",
              "SkewnessLongitudinal"        : "skewnessLongitudinal",
              "RmsTransverse"               : "rmsTransverse",
              "PositionY"                   : "centerY",
              "Excentricity"                : "eccentricity" }.toTable()
               # "TotalCharge"                 : "totalCharge"
  let
    mEvNum = h5M[grpM / "EventNumber", float32].mapIt(it.round.int)
    mEvSet = mEvNum.toSet
    tpEvNum = h5Tp[grpTp / "eventNumber", int]
    tpEvSet = tpEvNum.toSet
    utype = special_type(uint8)
    runPath = grpTp.replace("/reconstruction", "/runs")
  # sanity check if raw_data reconstructed such that only 1 cluster w/ any number of pixels
  doAssert mEvSet == tpEvSet, "Sanity check failed! Reconstruct with 1 cluster only with cutoff_size == 1"
  #let
    #rawEvNum = h5Tp[runPath.replace("/chip_0", "") / "eventNumber", int]
    #xdata = h5Tp[runPath / "raw_x", utype, uint8]
    #rawHits = xdata.mapIt(it.len)
    #underFive = rawHits.filterIt(it <= 5).len
  #doAssert tpEvNum.len == rawEvNum.len - underFive + 1, " no tp: " & $tpEvNum.len & " while " & $rawEvNum.len & " - " & $underFive & " = " & $(rawEvNum.len - underFive)

  #doAssert mEvNum.len == tpEvNum.len, " Number of events differs! M: " & $mEvNum.len & " Tp: " & $tpEvNum.len
  #for i in 0 .. mEvNum.high:
  #  doAssert mEvNum[i] == tpEvNum[i], " event number differs! M: " & $mEvNum[i] & " Tp: " & $tpEvNum[i]
  #doAssert mEvNum == tpEvNum, "Event numbers differ!"
  let layout = Layout(title: grpM,
                      width: 1920,
                      height: 1080)
  var grid = createGrid(numPlots = tab.len, layout = layout)
  var i = 0
  for marKey, tpKey in tab:
    let
      mData = h5M[grpM / marKey, float32].mapIt(it.float)
      mFilter = zip(mEvNum, mData) --> filter(it[0] in tpEvSet) --> map(it[1])
      tpData = h5Tp[grpTp / tpKey, float64]
      tpFilter = zip(tpEvNum, tpData) --> filter(it[0] in mEvSet) --> map(it[1])
    doAssert tpFilter.len == mFilter.len, " smaller " & $tpFilter.len & " while M " & $mFilter.len
    let diff = zip(mFilter, tpFilter) --> map(it[1] - it[0]) --> filter(abs(it) > 0.01 and abs(it) < 1e5)
    let plt = scatterPlot(x = (0 .. diff.high).toSeq.mapIt(it.float), y = diff)
      .width(1920)
      .height(1400)
      .xlabel(&"{tpKey} ratio: {diff.len.float / tpFilter.len.float:.2e}")
    echo type(plt).name
    grid[i] = plt
    inc i
  grid.show()

proc compare(h5Marlin, h5TpAnalysis: string) =
  ## compares the Marlin analysis with the TimepixAnalysis by creating
  ## plots of the difference between the Marlin properties and the TP
  ## properties
  echo "h5f ", h5TpAnalysis
  var
    h5M = H5file(h5Marlin, "r")
    h5Tp = H5file(h5TpAnalysis, "r")

  for run, grp in runs(h5Tp):
    echo "!"
    let tpGroupRun = h5Tp[grp.grp_str]
    let mGroupName = tpGroupRun.attrs["CDL_name", string]
    echo "Cm ", mGroupName, " and ", grp
    compareRun(h5M, h5Tp, mGroupName, grp / "chip_0")
    #compareRun((file: h5M, grp: mGroupName),
    #           (file: h5M, grp: mGroupName))
    echo "?"


proc main =
  let args = docopt(doc)
  echo args
  let h5file = $args["<h5file>"]
  let extract = args["--extract"].toBool
  let compareTo = $args["--compare"]

  if extract:
    var h5f = H5file(h5file, "r")
    for grp in h5f:
      echo "Reading group ", grp.name
      h5f.extractFromGroup(grp.name)
    discard h5f.close()
  else:
    doAssert compareTo != "nil"
    compare(h5file, compareTo)


when isMainModule:
  main()
