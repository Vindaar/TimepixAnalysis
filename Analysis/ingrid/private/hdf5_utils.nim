import ../ingrid_types
import helpers/utils
import os except FileInfo
import strutils, times, strformat, sequtils, tables, re, algorithm, sets, strscans
import nimhdf5, arraymancer
import macros
import pure

#####################################################
# Procs describing the data layout in the HDF5 file #
#####################################################

const PlotDirPrefixAttr* = "plotDirPrefix"
const PlotDirRawPrefixAttr* = "plotDirRawPrefix"
const XrayReferenceDsets* = {
  igHits,
  igLengthDivRadius,
  igSkewnessLongitudinal,
  igSkewnessTransverse,
  igRmsTransverse,
  igEccentricity,
  igHits,
  igKurtosisLongitudinal,
  igKurtosisTransverse,
  igLength,
  igNumClusters,
  igFractionInHalfRadius,
  igRadiusDivRmsTrans,
  igBalance,
  igWidth,
  igRmsLongitudinal,
  igLengthDivRmsTrans,
  igRotationAngle,
  igEnergyFromCharge,
  igLikelihood,
  igRadius,
  igFractionInTransverseRms,
  igTotalCharge}

proc toDset*(igKind: InGridDsetKind, frameworkKind: FrameworkKind = fkTpa): string =
  ## Converts a InGridDsetKind enum element to the string representation given
  ## the framework the data was created with
  case igKind
  of igInvalid: doAssert false, "Invalid dataset kind!"
  of igCenterX:
    case frameworkKind
    of fkTpa: result = "centerX"
    of fkMarlin: result = "PositionX"
  of igCenterY:
    case frameworkKind
    of fkTpa: result = "centerY"
    of fkMarlin: result = "PositionY"
  of igEventNumber:
    case frameworkKind
    of fkTpa: result = "eventNumber"
    of fkMarlin: result = "EventNumber"
  of igHits:
    case frameworkKind
    of fkTpa: result = "hits"
    of fkMarlin: result = "NumberOfPixels"
  of igEccentricity:
    case frameworkKind
    of fkTpa: result = "eccentricity"
    of fkMarlin: result = "Excentricity"
  of igSkewnessLongitudinal:
    case frameworkKind
    of fkTpa: result = "skewnessLongitudinal"
    of fkMarlin: result = "SkewnessLongitudinal"
  of igSkewnessTransverse:
    case frameworkKind
    of fkTpa: result = "skewnessTransverse"
    of fkMarlin: result = "SkewnessTransverse"
  of igRmsTransverse:
    case frameworkKind
    of fkTpa: result = "rmsTransverse"
    of fkMarlin: result = "RmsTransverse"
  of igRmsLongitudinal:
    case frameworkKind
    of fkTpa: result = "rmsLongitudinal"
    of fkMarlin: result = "RmsLongitudinal"
  of igKurtosisLongitudinal:
    case frameworkKind
    of fkTpa: result = "kurtosisLongitudinal"
    of fkMarlin: result = "KurtosisLongitudinal"
  of igKurtosisTransverse:
    case frameworkKind
    of fkTpa: result = "kurtosisTransverse"
    of fkMarlin: result = "KurtosisLongitudinal"
  of igLength:
    case frameworkKind
    of fkTpa: result = "length"
    of fkMarlin: result = "Length"
  of igWidth:
    case frameworkKind
    of fkTpa: result = "width"
    of fkMarlin: result = "Width"
  of igLengthDivRmsTrans:
    case frameworkKind
    of fkTpa: result = "lengthDivRmsTrans"
    of fkMarlin: doAssert false, "`igLengthDivRmsTrans` does not exist as a dataset in Marlin. Is calculated!"
  of igRotationAngle:
    case frameworkKind
    of fkTpa: result = "rotationAngle"
    of fkMarlin: result = "RotationAngle"
  of igEnergyFromCharge:
    case frameworkKind
    of fkTpa: result = "energyFromCharge"
    of fkMarlin: result = "EnergyFromCharge"
  of igEnergyFromPixel:
    case frameworkKind
    of fkTpa: result = "energyFromPixel"
    of fkMarlin: result = "EnergyFromPixel"
  of igLikelihood:
    case frameworkKind
    of fkTpa: result = "likelihood"
    of fkMarlin: result = "LikelihoodMarlin"
  of igFractionInTransverseRms:
    case frameworkKind
    of fkTpa: result = "fractionInTransverseRms"
    of fkMarlin: result = "FractionWithinRmsTransverse"
  of igTotalCharge:
    case frameworkKind
    of fkTpa: result = "totalCharge"
    of fkMarlin: result = "TotalCharge"
  of igNumClusters, igFractionInHalfRadius, igRadiusDivRmsTrans,
     igRadius, igBalance, igLengthDivRadius:
    doAssert false, "Only exists in 2014 XrayReferenceFile.h5: " & $igKind

proc toIngridDset*(dset: string): InGridDsetKind =
  ## Converts a InGridDsetKind enum element to the string representation given
  ## the framework the data was created with
  if dset in ["hits", "NumberOfPixels"]: result = igHits
  elif dset == "centerX": result = igCenterX
  elif dset == "centerY": result = igCenterY
  elif dset in ["eventNumber", "EventNumber"]:  result = igEventNumber
  elif dset in ["eccentricity", "Excentricity"] : result = igEccentricity
  elif dset in ["skewnessLongitudinal", "SkewnessLongitudinal"]: result = igSkewnessLongitudinal
  elif dset in ["skewnessTransverse", "SkewnessTransverse"]: result = igSkewnessTransverse
  elif dset in ["rmsLongitudinal", "RmsLongitudinal"]: result = igRmsLongitudinal
  elif dset in ["rmsTransverse", "RmsTransverse"]: result = igRmsTransverse
  elif dset in ["kurtosisLongitudinal", "KurtosisLongitudinal"]: result = igKurtosisLongitudinal
  elif dset in ["kurtosisTransverse", "KurtosisTransverse"]: result = igKurtosisTransverse
  elif dset in ["length", "Length"]: result = igLength
  elif dset in ["width", "Width"]: result = igWidth
  elif dset == "lengthDivRmsTrans": result = igLengthDivRmsTrans
  elif dset in ["rotationAngle", "RotationAngle"]: result = igRotationAngle
  elif dset in ["energyFromCharge", "EnergyFromCharge"]: result = igEnergyFromCharge
  elif dset in ["energyFromPixel", "EnergyFromPixel"]: result = igEnergyFromPixel
  elif dset in ["likelihood", "LikelihoodMarlin"]: result = igLikelihood
  elif dset in ["fractionInTransverseRms", "FractionWithinRmsTransverse"]: result = igFractionInTransverseRms
  elif dset in ["totalCharge", "TotalCharge"]: result = igTotalCharge
  elif dset == "xrayperevent": result = igNumClusters
  elif dset == "fractionwithin0.5radius": result = igFractionInHalfRadius
  elif dset == "radiusdivbyrmsy": result = igRadiusDivRmsTrans
  elif dset == "radius": result = igRadius
  elif dset == "balance": result = igBalance
  elif dset == "lengthdivbyradius": result = igLengthDivRadius
  else: result = igInvalid

func cdlPath*(tfKindStr: string, year = "2019"): string =
  var myr = ""
  case year
  of "2014", "2015":
    myr = "apr2014"
  of "2018", "2019":
    myr = "feb2019"
  result = &"calibration-cdl-{myr}-{tfKindStr}"

func cdlGroupName*(tfKindStr, year, dset: string): string =
  let dsetName = dset.extractFilename
  result = &"{cdlPath(tfKindStr, year)}/{dsetName}"

func cdlToXrayBinning2014Map(): Table[InGridDsetKind, tuple[bins: int, min, max: float]] =
  ## Maps the names of the `XrayReferenceDataSet.h5` (2014) to the
  ## number of bins and min, max values that must be given to the histogram function
  ## to arrive at the result from the `calibration-cdl.h5` (2014) file.
  result = { igLengthDivRadius : (bins: 100, min: 0.9950000047683716, max: 1.985000014305115),
             igSkewnessLongitudinal : (bins: 100, min: -5.050000190734863, max: 4.849999904632568),
             igSkewnessTransverse : (bins: 100, min: -5.050000190734863, max: 4.849999904632568),
             igRmsTransverse : (bins: 150, min: -0.01666666753590107, max: 4.949999809265137),
             igEccentricity : (bins: 150, min: 0.9700000286102295, max: 9.909999847412109),
             igHits : (bins: 250, min: -0.5, max: 497.5),
             igKurtosisLongitudinal : (bins: 100, min: -5.050000190734863, max: 4.849999904632568),
             igKurtosisTransverse : (bins: 100, min: -5.050000190734863, max: 4.849999904632568),
             igLength : (bins: 200, min: -0.05000000074505806, max: 19.85000038146973),
             igNumClusters : (bins: 6, min: -1.0, max: 4.0),
             igFractionInHalfRadius : (bins: 100, min: -0.004999999888241291, max: 0.9850000143051147),
             igRadiusDivRmsTrans : (bins: 100, min: -0.05000000074505806, max: 9.850000381469727),
             igBalance : (bins: 500, min: -0.0003999999898951501, max: 0.3987999856472015),
             igWidth : (bins: 100, min: -0.05000000074505806, max: 9.850000381469727),
             igRmsLongitudinal : (bins: 150, min: -0.01666666753590107, max: 4.949999809265137),
             igLengthDivRmsTrans : (bins: 150, min: -0.1000000014901161, max: 29.70000076293945),
             igRotationAngle : (bins: 100, min: -0.0157079640775919, max: 3.094468832015991),
             igEnergyFromCharge : (bins: 100, min: -0.05000000074505806, max: 9.850000381469727),
             igLikelihood : (bins: 200, min: -40.125, max: 9.625),
             igRadius : (bins: 100, min: -0.02500000037252903, max: 4.925000190734863),
             igFractionInTransverseRms : (bins: 100, min: -0.004999999888241291, max: 0.9850000143051147),
             igTotalCharge : (bins: 200, min: -6250.0, max: 2481250.0) }.toTable

func cdlToXrayBinning2014*(name: string): tuple[bins: int, min, max: float] =
  const map = cdlToXrayBinning2014Map()
  let nameIgKind = name.toIngridDset
  if nameIgKind in map:
    result = map[nameIgKind]

func cdlToXrayBinning2018Map(): Table[InGridDsetKind, tuple[bins: int, min, max: float]] =
  ## Maps the names of the `XrayReferenceDataSet.h5` (2019) to the
  ## number of bins and min, max values that must be given to the histogram function
  ## to arrive at the result equivalent to the `calibration-cdl.h5` (2019) file.
  result = { igSkewnessLongitudinal : (bins: 100, min: -5.050000190734863, max: 4.849999904632568),
             igSkewnessTransverse : (bins: 100, min: -5.050000190734863, max: 4.849999904632568),
             igRmsTransverse : (bins: 150, min: -0.01666666753590107, max: 4.949999809265137),
             igEccentricity : (bins: 150, min: 0.9700000286102295, max: 9.909999847412109),
             igHits : (bins: 250, min: -0.5, max: 497.5),
             igKurtosisLongitudinal : (bins: 100, min: -5.050000190734863, max: 4.849999904632568),
             igKurtosisTransverse : (bins: 100, min: -5.050000190734863, max: 4.849999904632568),
             igLength : (bins: 200, min: -0.05000000074505806, max: 19.85000038146973),
             igWidth : (bins: 100, min: -0.05000000074505806, max: 9.850000381469727),
             igRmsLongitudinal : (bins: 150, min: -0.01666666753590107, max: 4.949999809265137),
             igLengthDivRmsTrans : (bins: 150, min: -0.1000000014901161, max: 29.70000076293945),
             igRotationAngle : (bins: 100, min: -0.0157079640775919, max: 3.094468832015991),
             igEnergyFromCharge : (bins: 100, min: -0.05000000074505806, max: 9.850000381469727),
             igLikelihood : (bins: 200, min: -40.125, max: 9.625),
             igFractionInTransverseRms : (bins: 100, min: -0.004999999888241291, max: 0.9850000143051147),
             igTotalCharge : (bins: 200, min: -6250.0, max: 2481250.0) }.toTable

func cdlToXrayBinning2018*(name: string): tuple[bins: int, min, max: float] =
  const map = cdlToXrayBinning2018Map()
  let nameIgKind = name.toIngridDset
  if nameIgKind in map:
    result = map[nameIgKind]

func inCdl2018*(name: string): bool =
  const dsetSet = [ "skewnessLongitudinal",
                    "skewnessTransverse",
                    "rmsTransverse",
                    "eccentricity",
                    "hits",
                    "kurtosisLongitudinal",
                    "kurtosisTransverse",
                    "length",
                    "width",
                    "rmsLongitudinal",
                    "lengthDivRmsTrans",
                    "rotationAngle",
                    "energyFromCharge",
                    "likelihood",
                    "fractionInTransverseRms",
                    "totalCharge" ].toHashSet
  result = name in dsetSet

func inCdl2014*(name: string): bool =
  const dsetSet = [ "EnergyFromCharge",
                    "SkewnessLongitudinal",
                    "RotationAngle",
                    "Radius",
                    "Width",
                    "Length",
                    "KurtosisLongitudinal",
                    "Excentricity",
                    "KurtosisTransverse",
                    "RmsLongitudinal",
                    "RmsTransverse",
                    "SkewnessTransverse",
                    "LikelihoodMarlin",
                    "NumberOfPixels",
                    "TotalCharge",
                    "FractionWithinRmsTransverse" ].toHashSet
  result = name in dsetSet

func cdlToXray2014Map(): Table[string, string] =
  ## Maps the datasets from the `calibration-cdl.h5` (2014) file to the
  ## `XrayReferenceDataSet.h5` (2014) file.
  ## The latter is derived from the former via the charge cuts via:
  ## `func getEnergyBinMinMaxVals*(): Table[string, Cuts]`
  ## applied to the datasets.
  ## The result is binned via the binning described in cdlToXrayBinning2014()
  ## The following datasets are unaccounted for in the map, since they must
  ## be calculated from the existing other datasets:
  ## - "xrayperevent"
  ## - "lengthdivbyradius"
  ## - "lengthdivbyrmsy" # <- exists in TPA
  ## - "fractionwithin0.5radius"
  ## - "radiusdivbyrmsy"
  ## - "balance"
  result = { "EnergyFromCharge" : "energy",
             "SkewnessLongitudinal" : "skewnessl",
             "RotationAngle" : "rotAngle",
             "Radius" : "radius",
             "PositionX" : "",
             "PositionY" : "",
             "Width" : "width",
             "Length" : "length",
             "KurtosisLongitudinal" : "kurtosisl",
             "RunType" : "",
             "EnergyFromPixels" : "",
             "Excentricity" : "excentricity",
             "KurtosisTransverse" : "kurtosist",
             "RmsLongitudinal" : "rmsx",
             "RmsTransverse" : "rmsy",
             "SkewnessTransverse" : "skewnesst",
             "EventNumber" : "",
             "LikelihoodMarlin" : "likelihood",
             "NumberOfPixels" : "pixels",
             "TotalCharge" : "charge",
             "FractionWithinRmsTransverse" : "fractionwithinrmsy",
             "Timestamp" : "",
             "RunNumber" : ""}.toTable

# mapping of 2014 calibration-cdl data to XrayReference data
proc cdlToXray2014*(name: string): string =
  const map = cdlToXray2014Map()
  result = map[name]

proc getFloatGeometryNames*(): array[12, string] =
  ## returns all dataset names in the H5 output file, which are members
  ## of a `ClusterGeometry` object
  result = ["rmsLongitudinal", "rmsTransverse", "skewnessLongitudinal", "skewnessTransverse",
            "kurtosisLongitudinal", "kurtosisTransverse", "eccentricity", "rotationAngle",
            "length", "width", "fractionInTransverseRms", "lengthDivRmsTrans"]

proc getFloatToANames*(): array[5, string] =
  ## returns all dataset names in the H5 output file, which are members
  ## of a `ToAGeometry` object
  result = ["toaLength", "toaMean", "toaRms", "toaSkewness", "toaKurtosis"]

proc getUint16ToANames*(): array[1, string] =
  ## returns all dataset names in the H5 output file, which are members
  ## of a `ToAGeometry` object
  result = ["toaMin"]

proc getFloatClusterNames*(): array[2, string] =
  ## returns all dataset names in the H5 output file, which are members of
  ## a `ClusterObject` object
  result = ["centerX", "centerY"]

proc getFloatDsetNames*(): array[14, string] =
  ## returns the names of all datasets in the H5 output file, which appear as
  ## (N, 1) data columns. Combination of two above procs
  # need to define consts of arrays to use `+` macro
  const
    float_geo = getFloatGeometryNames()
    float_obj = getFloatClusterNames()
  result = float_geo + float_obj

proc getIntClusterNames*(): array[2, string] =
  ## returns names of datasets in H5 output file, which are integer datasets and
  ## members of a `ClusterObject`
  result = ["hits", "sumTot"]

proc getIntDsetNames*(): array[3, string] =
  ## returns all names of integer dataset for the H5 output file, which appear
  ## as (N, 1) data columns
  # need to define consts of arrays to use `+` macro
  const
    int_cluster = getIntClusterNames()
    int_dset = ["eventNumber"]
  result = int_cluster + int_dset

template rawGroupGrpStr*(): grp_str =
  "/runs/".grp_str

template recoGroupGrpStr*(): grp_str =
  "/reconstruction/".grp_str

template likelihoodGroupGrpStr*(): grp_str =
  "/likelihood/".grp_str

template tpx3InterpGroupGrpStr*(): grp_str =
  "/interpreted/".grp_str

template rawDataBase*(): string =
  "/runs/run_"

template recoBase*(): string =
  "/reconstruction/run_"

template likelihoodBase*(): string =
  "/likelihood/run_"

template tpx3Base*(): string =
  "/interpreted/run_"

template rawDataChipBase*(runNumber: int): string =
  "/runs/run_$#/chip_" % $runNumber # & "$#"

template recoDataChipBase*(runNumber: int): string =
  "/reconstruction/run_$#/chip_" % $runNumber # & "$#"

template tpx3DataChipBase*(runNumber: int): string =
  "/interpreted/run_$#/chip_" % $runNumber # & "$#"

template likelihoodDataChipBase*(runNumber: int): string =
  "/likelihood/run_$#/chip_" % $runNumber # & "$#"

proc rawPath*(runNumber, chipNumber: int): grp_str {.inline.} =
  result = (rawDataChipBase(runNumber) & $chipNumber).grp_str

proc recoPath*(runNumber, chipNumber: int): grp_str {.inline.} =
  result = (recoDataChipBase(runNumber) & $chipNumber).grp_str

proc tpx3Path*(runNumber, chipNumber: int): grp_str {.inline.} =
  result = (tpx3DataChipBase(runNumber) & $chipNumber).grp_str

proc likelihoodPath*(runNumber, chipNumber: int): grp_str {.inline.} =
  result = (likelihoodDataChipBase(runNumber) & $chipNumber).grp_str

proc getGroupNameRaw*(runNumber: int): string =
  # generates the group name for a given run number
  result = rawDataBase() & $runNumber

proc recoRunGrpStr*(runNumber: int): grp_str {.inline.} =
  result = (recoBase() & $runNumber).grp_str

proc getGroupNameReco*(runNumber: int): string =
  # generates the reconstrution group name for a given run number
  result = recoBase() & $runNumber

proc fadcRawPath*(runNumber: int): string {.inline.} =
  result = getGroupNameRaw(runNumber) / "fadc"

proc fadcRecoPath*(runNumber: int): string {.inline.} =
  result = getGroupNameReco(runNumber) / "fadc"

template rawFadcBasename*(runNumber: int): string =
  fadcRawPath(runNumber) / "raw_fadc"

template trigRecBasename*(runNumber: int): string =
  fadcRawPath(runNumber) / "trigger_record"

template noiseBasename*(runNumber: int): string =
  fadcRecoPath(runNumber) / "noisy"

template minvalsBasename*(runNumber: int): string =
  fadcRecoPath(runNumber) / "minvals"

template pedestalBasename*(runNumber: int): string =
  fadcRecoPath(runNumber) / "pedestalRun"

template fadcDataBasename*(runNumber: int): string =
  fadcRecoPath(runNumber) / "fadc_data"

template fadcBaselineBasename*(runNumber: int): string =
  fadcRecoPath(runNumber) / "baseline"

template argMinvalBasename*(runNumber: int): string =
  fadcRecoPath(runNumber) / "argMinval"

template riseStartBasename*(runNumber: int): string =
  fadcRecoPath(runNumber) / "riseStart"

template fallStopBasename*(runNumber: int): string =
  fadcRecoPath(runNumber) / "fallStop"

template riseTimeBasename*(runNumber: int): string =
  fadcRecoPath(runNumber) / "riseTime"

template fallTimeBasename*(runNumber: int): string =
  fadcRecoPath(runNumber) / "fallTime"

template eventNumberBasenameReco*(runNumber: int): string =
  fadcRecoPath(runNumber) / "eventNumber"

template eventNumberBasenameRaw*(runNumber: int): string =
  fadcRawPath(runNumber) / "eventNumber"

proc hasDset*(h5f: H5File, runNumber, chipNumber: int, dset: string): bool =
  ## returns `true` if the given run and chip has the given `dset`
  let path = recoDataChipBase(runNumber) & $chipNumber / dset
  result = if path in h5f: true else: false

proc hasTotalChargeDset*(h5f: H5File, runNumber, chipNumber: int):
                           bool {.inline.} =
  ## returns `true` if the given run and chip has the
  ## `totalCharge`
  ## dataset
  result = h5f.hasDset(runNumber, chipNumber, "totalCharge")

proc hasFadc*(h5f: H5File, runNumber: int): bool =
  ## Returns `true` if the given H5 file has the `fadc` group in the given run number
  result = fadcRecoPath(runNumber) in h5f

proc hasRawRun*(h5f: H5File, runNumber: int): bool =
  ## checks for the existence of the given run in the file
  let path = rawDataBase & $runNumber
  # this will be the implementation once we reran the whole analysis...
  if path in h5f:
    # check if attribute `done` on run
    var grp = h5f[path.grp_str]
    if "rawDataFinished" in grp.attrs:
      let isDone = grp.attrs["rawDataFinished", string]
      if isDone == "true":
        result = true
      else:
        result = false
    else:
      result = false
  else:
    result = false
  when false:
    # old implementation
    result = if path in h5f: true else: false

proc runFinished*(h5f: H5File, runNumber: int) =
  ## writes the `rawDataFinished` attribute to the run with
  ## `runNumber`
  let path = rawDataBase() & $runNumber
  var grp = h5f[path.grp_str]
  grp.attrs["rawDataFinished"] = "true"

proc getCenterChip*(h5f: H5File, runNumber: int): int =
  ## reads the `centerChip` attribute from the run group corresponding to
  ## `runNumber`
  result = h5f[recoRunGrpStr(runNumber)].attrs["centerChip", int]


################################################################################
##################### HDF5 related helper functions ############################
################################################################################

proc getTrackingEvents*(h5f: H5File, group: H5Group, num_tracking: int = -1, tracking = true): seq[int] =
  ## given a `group` in a `h5f`, filter out all indices, which are part of
  ## a tracking (`tracking == true`) or not part of a tracking (`tracking == false`)
  ## NOTE: the indices of the whole timestamp array correspond to the event
  ## numbers of a run, since this is a 1:1 mapping.
  ## NOTE2: This only works as long as the raw input does ``actually`` contain
  ## all files! While this was true for Virtex 6 TOS files in 2017/18, it does
  ## not at all hold for V6 2014/15 and SRS TOS files!. Thus get the event numbers
  ## and map to them.
  result = @[]
  # attributes of this group
  var attrs = group.attrs
  try:
    # try except for check of num_trackings
    let ntrackings = attrs["num_trackings", int]
    const date_syntax = getDateSyntax()

    var
      tr_starts = newSeq[DateTime](ntrackings)
      tr_stops  = newSeq[DateTime](ntrackings)
    for i in 0 ..< ntrackings:
      # get start and stop time of each tracking
      tr_starts[i] = attrs[&"tracking_start_{i}", string].parse(date_syntax)
      tr_stops[i]  = attrs[&"tracking_stop_{i}", string].parse(date_syntax)
    # get the timestamp of all events
    let
      evNums = h5f[group.name / "eventNumber", int64].asType(int)
      tstamp = h5f[group.name / "timestamp", int64]
      # start and stop in seconds
      tr_starts_s = mapIt(tr_starts, it.toTime.toUnix)
      tr_stops_s  = mapIt(tr_stops,  it.toTime.toUnix)
    # filter out all indices of timestamps, which lie inside the tracking

    # first get all indices of all trackings in a seq[seq[int]]
    var allTrackingInds: seq[seq[int]] = @[]
    for k in 0 ..< ntrackings:
      allTrackingInds.add filterIt(toSeq(0 ..< tstamp.len)) do:
        tstamp[it] > tr_starts_s[k] and tstamp[it] < tr_stops_s[k]
    if tracking == true:
      if num_tracking >= 0:
        # simply get the correct indices from allTrackingInds
        result = allTrackingInds[num_tracking].mapIt(evNums[it])
      else:
        # flatten the allTrackingInds nested seq and return
        result = flatten(allTrackingInds).mapIt(evNums[it])
    else:
      # all outside trackings are simply the indices, which are not part of a flattened
      # allTrackingInds
      let allTrackingsFlat = flatten(allTrackingInds)
      # and now filter all indices not part of flattened index
      result = toSeq(0 ..< tstamp.len)
        .filterIt(it notin allTrackingsFlat)
        .mapIt(evNums[it])
  except KeyError:
    # in this case there is no tracking information. Keep all indices
    echo &"No tracking information in {group.name} found, use all clusters"
    result = @[]


proc filterTrackingEvents*[T: SomeInteger](cluster_events: seq[T], eventsInTracking: seq[int]): seq[int] =
  ## filters out all event numbers of a reconstructed run for one chip
  ## Need to remove all indices, which are within the tracking indices, but for which
  ## no cluster is found in the datasets, so that we can only read the clusters, which
  ## happened during (or outside) of a tracking
  ## inputs:
  ##   `cluster_events`: all events for one reconstructed chip
  ##   `tracking_inds`: the event numbers which are part (or not) of a tracking
  # set result to the indices of tracking (non tracking), i.e.
  # all allowed events
  result = toSeq(0 .. cluster_events.high)
  if eventsInTracking.len == 0:
    # in case we are handed an empty seq, we simply use all cluster events
    discard
  else:
    # now given all indices describing the cluster event numbers, filter out
    # those, which are ``not`` part of `eventsInTracking``
    result = toSeq(0 ..< cluster_events.len)
      .filterIt(cluster_events[it].int in eventsInTracking)

proc filterTrackingEvents*(h5f: H5File, group: H5Group, tracking_inds: seq[int]): seq[int] =
  ## wrapper around the above proc, which reads the data about which events are allowed
  ## by itself
  ## inputs:
  ##   `h5f`: H5file from which to read the data
  ##   `group`: H5Group object of the specific chip, which contains the clustes
  ##   `tracking_inds`: the indices which are part (or not) of a tracking
  let
    # event numbers of clusters of this chip
    evNumbers = h5f[(group.name / "eventNumber").dset_str][int64]
  result = filterTrackingEvents(evNumbers, tracking_inds)

proc removePrefix(s: string, prefix: string): string =
  result = s
  result.removePrefix(prefix)

iterator runs*(h5f: H5File, data_basename = recoBase()): (int, string) =
  ## simple iterator, which yields the run number and group name of runs in the file.
  ## If reco is true (default) we yield reconstruction groups, else raw groups
  ## Iterator saves the state of `h5f` during the first call to this iterator! If
  ## additional groups are added while iterating, they will be ignored.
  if h5f.visited == false:
    h5f.visit_file

  # get a seq of all groups currently in the H5 file
  # only want to iterate over groups existing at the time, when
  # this proc is being called.
  # If we just use the `keys` iterator for `h5f.groups` we end up
  # skipping runs randomly, since we insert new groups, changing the
  # iterator while iterating. Bad! Solves issue #8
  let groups = toSeq(keys(h5f.groups))
  var runNumber: int
  for grp in groups:
    if grp.startsWith(data_basename) and
       grp.removePrefix(data_basename).scanf("$i$.", runNumber):
      # now read some data. Return value will be added later
      yield (runNumber, grp)

iterator chipGroups*(h5f: H5File, data_basename = recoBase()): (int, int, string) =
  ## simple iterator, which yields the run number and chip group name of runs in the file.
  ## If reco is true (default) we yield reconstruction groups, else raw groups
  ## Iterator saves the state of `h5f` during the first call to this iterator! If
  ## additional groups are added while iterating, they will be ignored.
  if h5f.visited == false:
    h5f.visit_file

  let groups = toSeq(keys(h5f.groups))
  var
    runNumber = -1 # if run number is already part of `data_basename` we return -1
    chipNumber: int
  # determine whether user gave us a `base` path (i.e. ending on `run_`)
  let isBasePath = data_basename.endsWith("run_")
  if isBasePath:
    for grp in groups:
      if grp.startsWith(data_basename) and
         grp.removePrefix(data_basename).scanf("$i/chip_$i$.", runNumber, chipNumber):
        yield (runNumber, chipNumber, grp)
  else:
    # this must already contain a run number!
    let nameFormat = data_basename
    doAssert nameFormat[^1] in Digits, "Given path must be a run group, i.e. with a digit for a run number."
    for grp in groups:
      if grp.startsWith(nameFormat) and
         grp.removePrefix(nameFormat).scanf("/chip_$i$.", chipNumber):
        yield (runNumber, chipNumber, grp)


iterator dsets*(h5f: H5File,
                dsetName: string,
                dtype: typedesc,
                chipNumber: int,
                dataBasename = recoBase()): (int, seq[dtype]) =
  ## given a dataset name, its corresponding datatype (needed to define the return type)
  ## and a chip number, read all datasets of all runs stored in the given H5 file.
  ## Choose a base location, by default reconstruction group
  ## NOTE: this cannot yield any datatypes with variable length data!

  if h5f.visited == false:
    h5f.visit_file

  let runChipName = joinPath(dataBasename, r"run_(\d+)/chip_" & $chipNumber)
  let dsetPath = joinPath(runChipName, dsetName)
  let dsetLocationReg = re(dsetPath)
  var runNumber = newSeq[string](1)
  let dsets = toSeq(keys(h5f.datasets))
  for dset in dsets:
    if match(dset, dsetLocationReg, runNumber):
      # found a matching dataset, yield the run number as well as the actual
      # data
      var mdset = h5f[dsetPath.dset_str]
      yield (runNumber[0].parseInt, mdset[dtype])

proc getRunInfo*(path: string): RunInfo =
  ## wrapper around the above proc if only the path to the run is known
  let regex = r"^/([\w-_]+/)*data\d{6,9}\.txt$"
  let fadcRegex = r"^/([\w-_]+/)*data\d{6,9}\.txt-fadc$"
  let (_, runNumber, rfKind, _) = isTosRunFolder(path)
  let files = getListOfFiles(path, regex)
  let fadcFiles = getListOfFiles(path, fadcRegex)
  if files.len > 0:
    result.timeInfo = getRunTimeInfo(files)
  result.runNumber = runNumber
  result.rfKind = rfKind
  result.runType = rtNone
  result.path = path
  result.nEvents = files.len
  result.nFadcEvents = fadcFiles.len

proc timepixVersion*(h5f: H5File): TimepixVersion =
  result = Timepix1
  template tryGroup(grp: untyped): untyped =
    if grp.string in h5f:
      let group = h5f[grp]
      if "TimepixVersion" in group.attrs:
        result = parseEnum[TimepixVersion](group.attrs["TimepixVersion", string])
  tryGroup(rawGroupGrpStr())
  tryGroup(recoGroupGrpStr())
  tryGroup(likelihoodGroupGrpStr())
  tryGroup(tpx3InterpGroupGrpStr())

proc dataPath*(fileInfo: FileInfo, run, chip: int): grp_str =
  ## Returns the path to the data based on the `FileInfo` for a given run
  ## and chip number.
  case fileInfo.tpaFileKind
  of tpkRawData:    result = rawPath(run, chip)
  of tpkReco:       result = recoPath(run, chip)
  of tpkLogL:       result = likelihoodPath(run, chip)
  of tpkTpx3Interp: result = tpx3Path(run, chip)
  of tpkCDL:        result = fileInfo.cdlGroup.grp_str
  of tpkTpx3Raw:
    doAssert false, "Invalid file kind to plot data from! " & $tpkTpx3Raw

proc dataBase*(fileInfo: FileInfo): string =
  ## Returns the base path to the data based on the `FileInfo`
  case fileInfo.tpaFileKind
  of tpkRawData:    result = rawDataBase()
  of tpkReco:       result = recoBase()
  of tpkLogL:       result = likelihoodBase()
  of tpkTpx3Interp: result = tpx3Base()
  of tpkCDL:        result = fileInfo.cdlGroup
  of tpkTpx3Raw:
    doAssert false, "Invalid file kind to plot data from! " & $tpkTpx3Raw

proc getFileInfo*(h5f: H5File): FileInfo =
  ## returns a set of all run numbers in the given file
  # 1. determine the `TpaFileKind` (determines `baseGroup`)
  result.name = h5f.name
  result.tpaFileKind =
    if "/runs" in h5f: tpkRawData
    elif "/reconstruction" in h5f: tpkReco
    elif "/likelihood" in h5f: tpkLogL
    elif "meta_data" in h5f and "raw_data" in h5f: tpkTpx3Raw
    elif "/interpreted" in h5f: tpkTpx3Interp
    elif cdlPath("C-EPIC-0.6kV") in h5f: tpkCDL # attempt to read a group from a CDL file
    else:
      raise newException(IOError, "The input file " & $h5f.name & " is not a valid TPA H5 file.")
  let baseGroup = ($result.tpaFileKind).grp_str
  # 2. read information based on type of file
  case result.tpaFileKind
  of tpkTpx3Raw:
    result.timepix = Timepix3
  of tpkCDL:
    result.timepix = Timepix1
    result.runType = rtCalibration # of a special kind...
    result.centerChip = 3
    result.runs = @[0] # dummy run
    result.chips = @[3] # dummy chip
  else:
    var readAux = false
    # get reconstruction group
    let group = h5f[baseGroup]
    if "runType" in group.attrs:
      result.runType = parseEnum[RunTypeKind](group.attrs["runType", string], rtNone)
    if "runFolderKind" in group.attrs:
      result.rfKind = parseEnum[RunFolderKind](group.attrs["runFolderKind", string],
                                               rfUnknown)
    if "centerChip" in group.attrs:
      result.centerChip = group.attrs["centerChip", int]
    if "centerChipName" in group.attrs:
      result.centerChipName = group.attrs["centerChipName", string]

    for runNumber, group in runs(h5f, data_basename = baseGroup.string / "run_"):
      result.runs.add runNumber
      if not readAux:
        let grp = h5f[group.grp_str]
        let nChips = grp.attrs["numChips", int]
        result.chips = toSeq(0 ..< nChips)
        readAux = true
    # set timepix version
    result.timepix = h5f.timepixVersion()
    # sort the run numbers
    result.runs.sort

proc getFileInfo*(fname: string): FileInfo =
  ## Returns the `FileInfo` for an input file if it is a valid TPA H5 file.
  ## Raises an `IOError` if the input file is not a H5 file.
  if not fname.endsWith(".h5"):
    raise newException(IOError, "Input file `" & $fname & "` is not an HDF5 file. Cannot " &
      "return the `FileInfo`.")
  withH5(fname, "r"):
    result = getFileInfo(h5f)

proc parseTracking(grp: H5Group, idx: int): RunTimeInfo =
  let
    start = grp.attrs[&"tracking_start_{idx}", string]
    stop = grp.attrs[&"tracking_stop_{idx}", string]
  result = RunTimeInfo(t_start: start.parseTime("yyyy-MM-dd\'T\'HH:mm:sszzz", utc()),
                       t_end: stop.parseTime("yyyy-MM-dd\'T\'HH:mm:sszzz", utc()))
  # t_end for a tracking does not take into account length of an event, as it is decoupled
  # from actual events
  result.t_length = result.t_end - result.t_start

proc getExtendedRunInfo*(h5f: H5File, runNumber: int,
                         runType: RunTypeKind,
                         rfKind: RunFolderKind = rfNewTos): ExtendedRunInfo =
  ## reads the extended run info from a H5 file for `runNumber`
  result.runNumber = runNumber
  let
    grp = h5f[(recoBase() & $runNumber).grp_str]
    tstamp = h5f[grp.name / "timestamp", int64]
    evDuration = h5f[grp.name / "eventDuration", float64]
    nEvents = h5f[(grp.name / "eventNumber").dset_str].shape[0]
  var nFadcEvents = 0
  if minvalsBasename(runNumber) in h5f:
    nFadcEvents = h5f[minvalsBasename(runNumber).dset_str].shape[0]

  result.activeTime = initDuration(seconds = evDuration.foldl(a + b).round.int)
  result.nEvents = nEvents
  result.nFadcEvents = nFadcEvents
  var tInfo = RunTimeInfo(t_start: tstamp[0].fromUnix,
                          t_end: (tstamp[^1].float + evDuration[^1]).fromUnixFloat) # end is last tstamp + duration of that event!
  tInfo.t_length = tInfo.t_end - tInfo.t_start
  result.timeInfo = tInfo

  # calculate background time from trackings
  var trackingDuration = initDuration()

  var numTracking = 0
  if "num_trackings" in grp.attrs:
    numTracking = grp.attrs["num_trackings", int]
  for i in 0 ..< numTracking:
    let tracking = parseTracking(grp, i)
    result.trackings.add tracking
    trackingDuration = trackingDuration + tracking.t_length

  result.trackingDuration = trackingDuration
  result.nonTrackingDuration = result.timeInfo.t_length - trackingDuration

  result.rfKind = rfKind
  result.runType = runType

proc readFadcFromH5*(h5f: H5File, runNumber: int): ProcessedFadcRun =
  ## proc to read data from the HDF5 file from `group`
  ## returns the chip number and a sequence containing the pixel data for this
  ## event and its event number
  let fadcGroup = fadcRawPath(runNumber)
  doAssert fadcGroup in h5f
  let group = h5f[fadcGroup.grp_str]
  let evNumbers = h5f[group.name / "eventNumber", int]
  let trigRec = h5f[group.name / "trigger_record", int]
  let dset = h5f[(group.name / "raw_fadc").dset_str]
  var fadcData = newTensorUninit[uint16](dset.shape)
  dset.read(fadcData.toUnsafeView())
  let settings = FadcSettings(isValid: true,
                              channelMask: group.attrs["channel_mask", int],
                              frequency: group.attrs["frequency", int],
                              postTrig: group.attrs["posttrig", int],
                              preTrig: group.attrs["pretrig", int],
                              nChannels: group.attrs["n_channels", int],
                              samplingMode: group.attrs["sampling_mode", int],
                              pedestalRun: group.attrs["pedestal_run", int] == 1)
  result = ProcessedFadcRun(
    settings: settings,
    rawFadcData: fadcData,
    trigRecs: trigRec,
    eventNumber: evNumbers
  )

proc read[T](h5f: H5File, path: string, _: typedesc[T]): Tensor[T] =
  let dset = h5f[path.dset_str]
  result = newTensorUninit[T](dset.shape)
  dset.read(result.toUnsafeView())

proc readRecoFadc*(h5f: H5File, runNumber: int): RecoFadc =
  ## Returns the `RecoFadc` data, i.e. all 1D data computed in `calcRiseAndFallTime`
  ## for a given run.
  let fadcGroup = fadcRecoPath(runNumber)
  doAssert fadcGroup in h5f
  let group = h5f[fadcGroup.grp_str]
  result = RecoFadc(
    baseline: h5f.read(fadcBaselineBasename(runNumber), float),
    xmin: h5f.read(argMinvalBasename(runNumber), uint16),
    riseStart: h5f.read(riseStartBasename(runNumber), uint16),
    fallStop: h5f.read(fallStopBasename(runNumber), uint16),
    riseTime: h5f.read(riseTimeBasename(runNumber), uint16),
    fallTime: h5f.read(fallTimeBasename(runNumber), uint16),
  )

proc readRecoFadcRun*(h5f: H5File, runNumber: int): ReconstructedFadcRun =
  ## Returns the `RecoFadc` data, i.e. all 1D data computed in `calcRiseAndFallTime`
  ## for a given run.
  let fadcGroup = fadcRecoPath(runNumber)
  doAssert fadcGroup in h5f
  let group = h5f[fadcGroup.grp_str]
  result = ReconstructedFadcRun(
    eventNumber: h5f[group.name / "eventNumber", int],
    fadcData: h5f.read(fadcDataBasename(runNumber), float),
    noisy: h5f[noiseBasename(runNumber), int],
    minVals: h5f[minvalsBasename(runNumber), float]
  )

proc genPlotDirname*(h5f: H5File, outpath: string, attrName: string): string =
  ## generates a unique name for the directory in which all plots for this H5
  ## file will be created.
  # first check whether this file already has such a name stored in its attributes
  if attrName in h5f.attrs:
    # nothing to do
    result = h5f.attrs[attrName, string]
  else:
    # generate a new one. base filename w/o file extension and current date
    let (_, name, _) = splitFile(h5f.name)
    let timeStr = format(now(), "yyyy-MM-dd'_'HH-mm-ss")
    result = outpath / name & "_" & timeStr
    h5f.attrs[attrName] = result

proc toString[N](a: array[N, char]): string =
  ## Returns all non `\0` characters as a string
  for c in a:
    if c != '\0':
      result.add c
    else:
      break

proc readTpx3CompoundKeyValue*[T; U](h5f: H5File, path: string): U =
  ## Reads the `run_config` dataset in the `/configuration` group and returns
  ## it turned into a `Table[string, string]`
  let data = h5f[path, T]
  for el in data:
    let attr = el.attribute.toString
    when typeof(el.value) is array[128, char]:
      let val = el.value.toString
    elif typeof(el.value) is uint16:
      let val = el.value
    for field, fval in fieldPairs(result):
      if field.nimIdentNormalize == attr.nimIdentNormalize:
        when typeof(fval) is int:
          fval = parseInt(val)
        elif typeof(fval) is char:
          fval = val[0]
        else:
          fval = val

proc readTpx3RunConfig*(h5f: H5File, run: int): Tpx3RunConfig =
  ## Reads the `run_config` dataset in the `/configuration` group and returns
  ## it turned into a `Table[string, string]`
  let runPath = tpx3Base() & $run
  let cfgPath = runPath / "configuration/run_config"
  result = readTpx3CompoundKeyValue[Tpx3RunConfigRaw, Tpx3RunConfig](h5f, cfgPath)

proc readTpx3Dacs*(h5f: H5File, path: string): Tpx3Dacs =
  ## Reads the `run_config` dataset in the `/configuration` group and returns
  ## it turned into a `Table[string, string]`
  result = readTpx3CompoundKeyValue[Tpx3DacsRaw, Tpx3Dacs](h5f, path)

proc timeFromTpx3RunConfig*(tpx3: Tpx3RunConfig): DateTime =
  ## Parses the date / time from the `runName` and returns it as a `DateTime`
  # either two options is possible. lower case is older file version
  let tstr = tpx3.runName.removePrefix("data_take_").removePrefix("DataTake_")
  result = parse(tstr, "YYYY-MM-dd'_'HH-mm-ss")

proc chipNameFromTpx3RunConfig*(tpx3: Tpx3RunConfig): string =
  ## Parses the date / time from the `runName` and returns it as a `DateTime`
  result = &"{tpx3.chipX}{tpx3.chipY:<2}W{tpx3.chipWafer}"


type
  Tpx3RunConfigToT = object
    scan_id: string
    run_name: string
    software_version: string
    board_name: string
    firmware_version: int
    chip_wafer: int
    chip_x: char
    chip_y: int
    VTP_fine_start: int
    VTP_fine_stop: int
    mask_step: int
    tp_period: int
    thrfile: string
    maskfile: string

type
  Tpx3TotMean = object
    tot: float32
    totError: float32

import unchained
proc readToTFileTpx3*(filename: string,
                      startRead = 0.0,
                      totPrefix = "TOTCalib"): (int, Tot) =
  # 1. read chip number ⇐ does not exist in H5 file at the moment for Tpx3
  # 2. read VTP values and use them to compute voltages

  var h5f = H5open(filename, "r")
  let tpx3Dacs = h5f.readTpx3Dacs("/configuration/dacs")

  let vtpCoarse = tpx3Dacs.VTP_coarse

  # 3. read start and stop VTP_fine values via runConfig
  let runConfig = readTpx3CompoundKeyValue[Tpx3RunConfigRaw, Tpx3RunConfigToT](h5f, "/configuration/run_config")
  let vtpStart = runConfig.VTP_fine_start
  let vtpStop = runConfig.VTP_fine_stop
  # 3. read ToT & error values

  proc toMilliVolt(vtpCoarse, vtpFine: int): MilliVolt =
    result = (-5.0.mV * vtpCoarse + 2.5.mV * vtpFine).to(mV)

  let data = h5f["interpreted/mean_curve", Tpx3TotMean]
  let vtps = arange(vtpStart, vtpStop)
  doAssert vtps.len == data.len
  var tot: ToT
  for i in 0 ..< data.len:
    if data[i].totError > 0.0:
      tot.mean.add data[i].tot
      tot.std.add data[i].totError
      tot.pulses.add toMilliVolt(vtpCoarse.int, vtps[i]).float
  result = (0, tot) # chip assume 0 for now

when isMainModule:
  assert combineRawBasenameToT(0, 1) == "/runs/combined/ToT_0_1"
  assert combineRecoBasenameToT(0, 1) == "/reconstruction/combined/ToT_0_1"
