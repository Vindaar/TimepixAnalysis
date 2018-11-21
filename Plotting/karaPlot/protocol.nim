import karax / kbase
export kbase
import sets, strutils, strformat, sequtils, sugar

when defined(js):
  import karax / jjson
  import components / utils
else:
  import json

import common_types
import ingrid / ingrid_types

type
  Messages* = enum
    MessageUnknown = "Error"
    Connected = "Connected"
    DataStart = "DataStart"
    DataStop = "DataStop"
    DataSingle = "DataSingle"
    Data = "Data"
    Request = "Request"

  PacketKind* = enum
    PacketKindUnknown, Descriptors, Plots

  PacketHeader* = object
    msg*: Messages
    kind*: PacketKind
    recvData*: bool
    done*: bool

  DataPacket* = object
    header*: PacketHeader
    payload*: kstring

  # a simple object storing the runs, chips etc. from a given
  # H5 file
  FileInfo* = object
    runs*: seq[int]
    chips*: seq[int]
    runType*: RunTypeKind
    rfKind*: RunFolderKind
    centerChip*: int
    centerChipName*: string
    hasFadc*: bool # reads if FADC group available
    # TODO: move the following to a CONFIG object
    plotlySaveSvg*: bool
    # NOTE: add other flags for other optional plots?
    # if e.g. FeSpec not available yet, we can just call the
    # procedure to create it for us



# payload size of the data packets. A header is added to this, which
# is why it's not close to 32768
const FakeFrameSize* = 32000

const InGridFnameTemplate* = "$1_run$2_chip$3_$4_binSize$5_binRange$6_$7"
const InGridTitleTemplate* = "Dataset: $1 for run $2, chip $3 in range: $4"
const FadcFnameTemplate* = "fadc_$1_run$2_$3_binSize$4_binRange$5_$6"
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

#proc initDataPacket*(): DataPacket =
#  result.recvData = false
#  result.done = false
#  result.data = kstring""

proc initDataPacket*(): DataPacket =
  result.header = PacketHeader(msg: MessageUnknown,
                               kind: PacketKindUnknown,
                               recvData: false,
                               done: false)
  result.payload = ""

proc initDataPacket*(kind: PacketKind,
                     header: Messages,
                     payload: kstring = "",
                     recvData: bool = false,
                     done: bool = false): DataPacket =
  result.header = PacketHeader(kind: kind,
                               msg: header,
                               recvData: recvData,
                               done: done)
  result.payload = payload

#proc getStartPacket*(kind: PacketKind): JsonNode =
#  ## returns the data start packet
#  result = initDataPacket(kind,
#                          header = Messages.DataStart)
#
#proc getStopPacket*(kind: PacketKind): JsonNode =
#  ## returns the data start packet
#  result = initDataPacket(kind,
#                          header = Messages.DataStop)

func `%`*(header: PacketHeader): JsonNode =
  ## serializes a `PacketHeader`
  result = newJObject()
  result[kstring"kind"] = % header.kind
  result[kstring"msg"] = % header.msg
  # bools need to be stored as strings for JS interop
  result[kstring"recvData"] = % $(header.recvData)
  result[kstring"done"] = % $(header.done)

func `%`*(packet: DataPacket): JsonNode =
  ## serialize `packet` to send via websocket
  result = newJObject()
  result[kstring"header"] = % packet.header
  result[kstring"payload"] = % packet.payload

proc `$`*(packet: DataPacket): kstring =
  result = (% packet).pretty

when not defined(js):
  proc asData*(packet: DataPacket): kstring =
    let node = % packet
    toUgly(result, node)

proc parsePacketHeader(header: JsonNode): PacketHeader =
  ## unmarshals a header as a JsonNode to a `PacketHeader`
  result.msg = parseEnum[Messages](header["msg"].getStr,
                                   Messages.MessageUnknown)
  result.kind = parseEnum[PacketKind](header["kind"].getStr,
                                      PacketKind.PacketKindUnknown)
  result.recvData = parseBool($(header["recvData"].getStr))
  result.done = parseBool($(header["done"].getStr))

proc parseDataPacket*(packet: kstring): DataPacket =
  ## unmarshals the data packet string
  # first create a JsonNode from the packet
  let pJson = parseJson(packet)
  result.header = parsePacketHeader(pJson[kstring"header"])
  result.payload = pJson[kstring"payload"].getStr

#proc parseHeaderPacket*(packet: kstring): (Messages, PacketKind) =
#  ## parses the header packet and returns the tuple of `Messages`
#  ## and `PacketKind` to tell whether data follows or is finished
#  ## and how to interpret it
#  echo "Parsing packet"
#  let pJson = parseJson(packet)
#  echo "Done ", pJson.pretty
#  echo pJson[kstring"MessageKind"].getStr
#  result[0] = parseEnum[Messages](pJson[kstring"MessageKind"].getStr,
#                                  Messages.MessageUnknown)
#  result[1] = parseEnum[PacketKind](pJson[kstring"PacketKind"].getStr,
#                                    PacketKind.PacketKindUnknown)

# PlotDescriptor object
# storing all needed information to create a specific plot

func `%`*(pd: PlotDescriptor): JsonNode =
  ## serialize `PlotDescriptor` to Json
  result = newJObject()
  result[kstring"runType"] = % $pd.runType
  result[kstring"name"] = % pd.name
  result[kstring"runs"] = % pd.runs
  result[kstring"chip"] = % pd.chip
  result[kstring"xlabel"] = % pd.xlabel
  result[kstring"ylabel"] = % pd.ylabel
  result[kstring"title"] = % pd.title
  result[kstring"plotKind"] = % $pd.plotKind
  case pd.plotKind
  of pkInGridDset, pkFadcDset:
    result[kstring"range"] = %* { kstring"low": % ($pd.range[0]),
                                  kstring"high": % ($pd.range[1]),
                                  kstring"name": % pd.range[2] }
    result[kstring"binSize"] = % ($pd.binSize)
    # binRange is a float, but we should never encounter `Inf`, thus keep as float
    result[kstring"binRange"] = %* { kstring"low": % ($pd.binRange[0]),
                                     kstring"high": % ($pd.binRange[1]) }
  of pkOccupancy, pkOccCluster:
    result[kstring"clampKind"] = % $pd.clampKind
    case pd.clampKind
    of ckAbsolute:
      result[kstring"clampA"] = % ($pd.clampA)
    of ckQuantile:
      result[kstring"clampQ"] = % ($pd.clampQ)
    else: discard
  of pkCombPolya:
    result[kstring"chipsCP"] = % pd.chipsCP
  else: discard

func `$`*(pd: PlotDescriptor): kstring =
  ## use JSON representation as string
  result = (% pd).pretty

func parsePd*(pd: JsonNode): PlotDescriptor =
  ## parse a PlotDescriptor stored as JsonNode back to an object
  result.runType = parseEnum[RunTypeKind](pd["runType"].getStr, rtNone)
  result.name = pd[kstring"name"].getStr
  result.runs = to(pd[kstring"runs"], seq[int])
  result.chip = pd[kstring"chip"].getInt
  result.xlabel = pd[kstring"xlabel"].getStr
  result.ylabel = pd[kstring"ylabel"].getStr
  result.title = pd[kstring"title"].getStr
  result.plotKind = parseEnum[PlotKind](pd[kstring"plotKind"].getStr,
                                        pkInGridDset)
  case result.plotKind
  of pkInGridDset, pkFadcDset:
    result.range = (low: pd[kstring"range"][kstring"low"].getStr.parseFloat,
                    high: pd[kstring"range"][kstring"high"].getStr.parseFloat,
                    name: pd[kstring"range"][kstring"name"].getStr)
    result.binSize = pd[kstring"binSize"].getStr.parseFloat
    result.binRange = (low: pd[kstring"binRange"]["low"].getStr.parseFloat,
                       high: pd[kstring"binRange"]["high"].getStr.parseFloat)
  of pkOccupancy, pkOccCluster:
    result.clampKind = parseEnum[ClampKind](pd[kstring"clampKind"].getStr,
                                            ckFullRange)
    case result.clampKind
    of ckAbsolute:
      result.clampA = pd[kstring"clampA"].getStr.parseFloat
    of ckQuantile:
      result.clampQ = pd[kstring"clampQ"].getStr.parseFloat
    else: discard
  of pkCombPolya:
    result.chipsCP = to(pd[kstring"chipsCP"], seq[int])
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
                                   $pd.range[2],
                                   $pd.binSize,
                                   $pd.binRange[0],
                                   $pd.binRange[1]]
  of pkFadcDset:
    name = FadcFnameTemplate %% [pd.name,
                                 runsStr,
                                 $pd.range[2],
                                 $pd.binSize,
                                 $pd.binRange[0],
                                 $pd.binRange[1]]
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
