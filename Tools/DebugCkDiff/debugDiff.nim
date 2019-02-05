import nimhdf5, zero_functional
import helpers / utils
import ingrid / [tos_helpers, ingrid_types]
import sequtils, strutils, plotly, strformat, os, sets, typetraits, stats
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
            #"TotalCharge"                 : "totalCharge",
            "Excentricity"                : "eccentricity" }.toTable()

proc compareRun(h5M, h5Tp: var H5FileObj, grpM, grpTp: string) =
  ## creates all plots for the given run
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
    #if tpKey == "rotationAngle":
    #  let diff = zip(mEvNum, mFilter, tpFilter) -->
    #    map((it[2] - it[1], it[1], it[2])) -->
    #    filter(abs(it[0]) > 3.0) -->
    #    map((it[1], it[2]))
    #  echo "Diffs are: "
    #  echo diff
    #  let x = diff --> map(it[0])
    #  let y = diff --> map(it[1])
    #  scatterPlot(x, y).show()

    let cutVal = percentile(tpData, 95) * 0.01
    var passesCut: proc(x: float): bool
    if tpKey == "totalCharge":
      echo (zip(mFilter, tpFilter) --> map(it[1] - it[0]))[0 .. 200]
      echo percentile(tpData, 95)
      passesCut = (
        proc(x: float): bool =
          abs(x) > cutVal and abs(x) < 1e4
      )
    elif tpKey == "eccentricity":
      passesCut = (
        proc(x: float): bool =
          abs(x) > cutVal and abs(x) < 1e4
      )
    else:
      passesCut = (
        proc(x: float): bool =
          abs(x) > cutVal and abs(x) < 1e9
      )
    let diff = zip(mFilter, tpFilter) -->
      map(it[1] - it[0]) -->
      filter(it.passesCut)

    let plt = scatterPlot(x = (0 .. diff.high).toSeq.mapIt(it.float), y = diff)
      .width(1920)
      .height(1400)
      .xlabel(&"{tpKey} ratio: {diff.len.float / tpFilter.len.float:.2e}")
    echo type(plt).name
    grid[i] = plt
    inc i
  grid.show(grpM.strip(chars = {'/'}) & ".svg")

proc compareToFile(h5Marlin, h5TpAnalysis: string) =
  ## compares the Marlin analysis with the TimepixAnalysis by creating
  ## plots of the difference between the Marlin properties and the TP
  ## properties
  echo "h5f ", h5TpAnalysis
  var
    h5M = H5file(h5Marlin, "r")
    h5Tp = H5file(h5TpAnalysis, "r")

  for run, grp in runs(h5Tp):
    let tpGroupRun = h5Tp[grp.grp_str]
    let mGroupName = tpGroupRun.attrs["CDL_name", string]
    compareRun(h5M, h5Tp, mGroupName, grp / "chip_0")

proc readAllRuns(h5f: var H5FileObj,
                 dset: string,
                 evNums: OrderedSet[int64],
                 runs: OrderedSet[int]): (seq[int64], seq[float]) =
  ## reads all runs given by the `HashSet` contained in the file for the
  ## given dataset, filtered to the event numbers given by `evNums`
  ## `evNums` should thus be from the other file (calibration-cdl.h5)
  var retData: seq[float]
  var retEvs: seq[int64]
  for r in runs:
    echo "Reading ", (recoPath(r, 0).string / dset)
    let dset = h5f[(recoPath(r, 0).string / dset).dset_str]
    let convert = dset.convertType(float64)
    let data = dset.convert
    let evNumsTarget = h5f[recoPath(r, 0).string / "eventNumber", int64]
    let filtered = zip(data, evNumsTarget) -->
      filter(it[1] in evNums)
    # add events filtered by event number
    retData.add (filtered --> map(it[0]))
    retEvs.add (filtered --> map(it[1]) --> to(seq[int64]))
  result = (retEvs, retData)

proc compareToPath(h5Marlin, tpAnaPath: string) =
  ## create comparison plots for ``calibration-cdl.h5`` vs reconstructed
  ## data from raw data
  var h5M = H5File(h5Marlin, "r")
  for grp in h5M:
    var h5Tp = H5file(tpAnaPath / grp.name & ".h5", "r")
    let
      # filter out run 8 because it doesn't exist in raw data
      runNumberAll = h5M[grp.name / "RunNumber", float32]
      evNumberAll = h5M[grp.name / "EventNumber", float32]
      filteredEvRun = zip(evNumberAll, runNumberAll) -->
        filter(it[1] != 8) -->
        map((it[0].round.int64, it[1].round.int))
      runNumber = filteredEvRun --> map(it[1])
      evNumber = filteredEvRun --> map(it[0])
      runNumSet = runNumber.toOrderedSet
    for mKey, tpKey in tab:
      let (tpEvs, tpData) = readAllRuns(h5Tp, tpKey, evNumber.toOrderedSet, runNumSet)
      doAssert tpEvs == evNumber, "Not the same event numbers! " & $tpEvs[0 .. 100] & "\n" &
        $evNumber[0 .. 100]


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
    if compareTo.endsWith(".h5"):
      compareToFile(h5file, compareTo)
    else:
      compareToPath(h5file, compareTo)


when isMainModule:
  main()
