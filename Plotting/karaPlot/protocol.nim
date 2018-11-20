import karax / kbase
export kbase
import sets, json, strutils

import ingrid / ingrid_types

type
  Messages* = enum
    Connected = "Connected"
    DataStart = "DataStart"
    DataStop = "DataStop"
    Request = "Request"

  DataPacket* = object
    recvData*: bool
    done*: bool
    data*: kstring



  DataPacket* = object
    header*: PacketHeader
    payload*: kstring

# payload size of the data packets. A header is added to this, which
# is why it's not close to 32768
const FakeFrameSize* = 32000

const InGridFnameTemplate* = "$1_run$2_chip$3_$4"
const InGridTitleTemplate* = "Dataset: $1 for run $2, chip $3 in range: $4"
const FadcFnameTemplate* = "fadc_$1_run$2_$3"
const FadcTitleTemplate* = "Dataset: $1 for run $2, fadc in range: $3"
const PolyaFnameTemplate* = "polya_run$1_chip$2"
const PolyaTitleTemplate* = "Polya for run $1 of chip $2"
const CombPolyaFnameTemplate* = "combined_polya_run$1"
const CombPolyaTitleTemplate* = "Polyas for all chips of run $1"
const OccupancyFnameTemplate* = "occupancy_run$1_chip$2_$3"
const OccupancyTitleTemplate* = "Occupancy of chip $1 for run $2, $3"
const OccClusterFnameTemplate* = "occupancy_clusters_run$1_chip$2_$3"
const OccClusterTitleTemplate* = "Occupancy of cluster centers for run $1, chip $2, $3"
const FeSpecFnameTemplate* = "fe_pixel_spectrum_run$1"
const FeSpecChargeFnameTemplate* = "fe_charge_spectrum_run$1"
const FeSpecTitleTemplate* = "Fe pixel spectrum of run $1 for chip $2"
const FeSpecChargeTitleTemplate* = "Fe charge spectrum of run $1 for chip $2"
const EnergyCalibFnameTemplate* = "fe_energy_calib_run$1"
const PhotoVsTimeFnameTemplate* = "photopeak_vs_time_runs$1"
const PhotoVsTimeTitleTemplate* = "Photopeak pixel peak position of Fe spectra (runs $1) vs time"
const PhotoPixDivChVsTimeFnameTemplate* = "photopeak_pix_div_charge_pos_vs_time_runs$1"
const PhotoPixDivChVsTimeTitleTemplate* = "Photopeak pix / charge of Fe spectra (runs $1) vs time"
const InGridEventTitleTemplate* = "InGrid event for run $1, chip $2, event index $3"
const InGridEventFnameTemplate* = "ingrid_event_run$1_chip$2_event$3"

proc initDataPacket*(): DataPacket =
  result.recvData = false
  result.done = false
  result.data = kstring""

# PlotDescriptor object
# storing all needed information to create a specific plot
func `%`*(pd: PlotDescriptor): JsonNode =
  ## serialize `PlotDescriptor` to Json
  result = newJObject()
  result["runType"] = % pd.runType
  result["name"] = % pd.name
  result["runs"] = % pd.runs
  result["chip"] = % pd.chip
  result["xlabel"] = % pd.xlabel
  result["ylabel"] = % pd.ylabel
  result["title"] = % pd.title
  result["plotKind"] = % pd.plotKind
  case pd.plotKind
  of pkInGridDset, pkFadcDset:
    result["range"] = %* { "low": % pd.range[0],
                           "high": % pd.range[1],
                           "name": % pd.range[2] }
  of pkOccupancy, pkOccCluster:
    result["clampKind"] = % pd.clampKind
    case pd.clampKind
    of ckAbsolute:
      result["clampA"] = % pd.clampA
    of ckQuantile:
      result["clampQ"] = % pd.clampQ
    else: discard
  of pkCombPolya:
    result["chipsCP"] = % pd.chipsCP
  else: discard

func `$`*(pd: PlotDescriptor): string =
  ## use JSON representation as string
  result = (% pd).pretty

func parsePd*(pd: JsonNode): PlotDescriptor =
  ## parse a PlotDescriptor stored as JsonNode back to an object
  result.runType = parseEnum[RunTypeKind](pd["runType"].getStr, rtNone)
  result.name = pd["name"].getStr
  result.runs = to(pd["runs"], seq[int])
  result.chip = pd["chip"].getInt
  result.xlabel = pd["xlabel"].getStr
  result.ylabel = pd["ylabel"].getStr
  result.title = pd["title"].getStr
  result.plotKind = parseEnum[PlotKind](pd["plotKind"].getStr)
  case result.plotKind
  of pkInGridDset, pkFadcDset:
    result.range = (low: pd["range"]["low"].getFloat,
                    high: pd["range"]["high"].getFloat,
                    name: pd["range"]["name"].getStr)
  of pkOccupancy, pkOccCluster:
    result.clampKind = parseEnum[ClampKind](pd["clampKind"].getStr)
    case result.clampKind
    of ckAbsolute:
      result.clampA = pd["clampA"].getFloat
    of ckQuantile:
      result.clampQ = pd["clampQ"].getFloat
    else: discard
  of pkCombPolya:
    result.chipsCP = to(pd["chipsCP"], seq[int])
  else: discard
func cKindStr(pd: PlotDescriptor, sep: string): string =
  case pd.clampKind
  of ckAbsolute:
    result = &"{pd.clampKind}{sep}{pd.clampA}"
  of ckQuantile:
    result = &"{pd.clampKind}{sep}{pd.clampQ}"
  of ckFullRange:
    result = &"{pd.clampKind}"

func `%%`[T](formatstr: string, a: openArray[T]): string =
  var b: seq[string]
  for x in a:
    b.add $x
  result = formatstr % b

proc buildOutfile*(pd: PlotDescriptor): kstring =
  var name = ""
  var runsStr = ""
  if pd.runs.len > 3:
    runsStr = &"{pd.runs.min}_{pd.runs.max}"
  else:
    runsStr = pd.runs.foldl($a & " " & $b, "").strip(chars = {' '})
  case pd.plotKind
  of pkInGridDset:
    name = InGridFnameTemplate %% [pd.name,
                                   runsStr,
                                   $pd.chip,
                                   $pd.range[2]]
  of pkFadcDset:
    name = FadcFnameTemplate %% [pd.name,
                                 runsStr,
                                 $pd.range[2]]
  of pkOccupancy:
    let clampStr = cKindStr(pd, "_")
    name = OccupancyFnameTemplate %% [runsStr,
                                      $pd.chip,
                                      clampStr]
  of pkOccCluster:
    let clampStr = cKindStr(pd, "_")
    name = OccClusterFnameTemplate %% [runsStr,
                                      $pd.chip,
                                      clampStr]
  of pkPolya:
    name = PolyaFnameTemplate %% [runsStr,
                                 $pd.chip]
  of pkCombPolya:
    name = CombPolyaFnameTemplate %% [runsStr]
  of pkFeSpec:
    name = FeSpecFnameTemplate %% [runsStr]
  of pkFeSpecCharge:
    name = FeSpecChargeFnameTemplate %% [runsStr]
  of pkFeVsTime:
    name = PhotoVsTimeFnameTemplate %% [runsStr]
  of pkFePixDivChVsTime:
    name = PhotoPixDivChVsTimeFnameTemplate %% [runsStr]
  of pkInGridEvent:
    name = InGridEventFnameTemplate %% [runsStr,
                                       $pd.chip,
                                       $pd.event]
  else:
    discard
  echo "Result ", $name
  result = &"figs/{name}.svg"

proc buildTitle*(pd: PlotDescriptor): kstring =
  var runsStr = ""
  if pd.runs.len > 3:
    runsStr = &"{pd.runs.min}..{pd.runs.max}"
  else:
    runsStr = pd.runs.foldl($a & " " & $b, "").strip(chars = {' '})
  case pd.plotKind
  of pkInGridDset:
    result = InGridTitleTemplate %% [pd.name,
                                    runsStr,
                                    $pd.chip,
                                    $pd.range[2]]
  of pkFadcDset:
    result = FadcTitleTemplate %% [pd.name,
                                  runsStr,
                                  $pd.range[2]]
  of pkOccupancy:
    let clampStr = cKindStr(pd, "@")
    result = OccupancyTitleTemplate %% [runsStr,
                                       $pd.chip,
                                       clampStr]
  of pkOccCluster:
    let clampStr = cKindStr(pd, "@")
    result = OccClusterTitleTemplate %% [runsStr,
                                       $pd.chip,
                                       clampStr]
  of pkPolya:
    result = PolyaTitleTemplate %% [runsStr,
                                   $pd.chip]
  of pkCombPolya:
    result = CombPolyaTitleTemplate %% [runsStr]
  of pkFeSpec:
    result = FeSpecTitleTemplate %% [runsStr,
                                    $pd.chip]
  of pkFeSpecCharge:
    result = FeSpecChargeTitleTemplate %% [runsStr,
                                          $pd.chip]
  of pkFeVsTime:
    result = PhotoVsTimeTitleTemplate %% [runsStr]
  of pkFePixDivChVsTime:
    result = PhotoPixDivChVsTimeTitleTemplate %% [runsStr]
  of pkInGridEvent:
    result = InGridEventTitleTemplate %% [runsStr,
                                         $pd.chip,
                                         $pd.event]
  else:
    discard
