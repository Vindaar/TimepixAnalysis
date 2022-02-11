import karax / kbase
export kbase
import sets, strutils, strformat, sequtils, sugar, times
from os import splitFile

when defined(js):
  import karax / jjson
  import components / utils
else:
  import json

import common_types
import ingrid / ingrid_types

type
  Messages* {.pure.} = enum
    MessageUnknown = "Error"
    Connected = "Connected"
    DataStart = "DataStart"
    DataStop = "DataStop"
    DataSingle = "DataSingle"
    Data = "Data"
    Request = "Request"
    # TODO: add message for "close nominally?"

  RequestKind* = enum
    rqPing     # just a ping, server/client will return pong
    rqFileInfo # request FileInfo of current file
    rqPlotDescriptors # request all PlotDescriptors the server believes
                      # available. Note: this might not be every possible
                      # plot, because the config.toml might limit the
                      # created PDs
    rqPlot # request a specific plot. Payload needs to be a
           # serialized PlotDescriptor

  PacketKind* {.pure.} = enum
    # answer packet server -> after a request
    PacketKindUnknown, Descriptors, FileInfo, Plots, Request

  PacketHeader* = object
    # A real Data packet typically server -> client
    msg*: Messages
    recvData*: bool
    done*: bool

  DataPacket* = object
    payload*: kstring
    case kind*: PacketKind
    of PacketKind.Request:
      # a Request packet typically client -> server
      reqKind*: RequestKind
    of PacketKind.Descriptors, PacketKind.Plots, PacketKind.FileInfo:
      # A real Data packet typically server -> client
      header*: PacketHeader
    of PacketKind.PacketKindUnknown: discard

# payload size of the data packets. A header is added to this, which
# is why it's not close to 32768
const FakeFrameSize* = 32000

const InGridFnameTemplate* = "$1_run$2_chip$3_$4_binSize$5_binRange$6_$7_chipRegion_$8"
const InGridTitleTemplate* = "Dataset: $1 for run $2, chip $3 in range: $4, chipRegion: $5"
const AnyScatterFnameTemplate* = "$1_run$2_chip$3_x_$4_y_$5_z_$6"
const AnyScatterTitleTemplate* = "Dataset: $4 / $5 by $6, for run $2, chip $3"
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
const PhotoVsTimeFnameTemplate* = "photopeak_vs_time_runs$1_$2"
const PhotoChargeVsTimeFnameTemplate* = "photopeak_charge_vs_time_runs$1_$2"
const PhotoVsTimeTitleTemplate* = "Photopeak pixel peak position of Fe spectra (runs $1) vs time"
const PhotoPixDivChVsTimeFnameTemplate* = "photopeak_pix_div_charge_pos_vs_time_runs$1"
const PhotoPixDivChVsTimeTitleTemplate* = "Photopeak pix / charge of Fe spectra (runs $1) vs time"
const PhotoDivEscapeFnameTemplate* = "photo_div_escape_pos_vs_time_runs$1"
const PhotoDivEscapeTitleTemplate* = "Ratio of photo and escape peak in charge (runs $1) vs time"
const InGridEventTitleTemplate* = "InGrid event for run $1, chip $2, event index $3"
const InGridEventFnameTemplate* = "ingrid_event_run$1_chip$2_event$3"
const FadcEventTitleTemplate* = "FADC event for run $1, event index $2"
const FadcEventFnameTemplate* = "fadc_event_run$1_event$2"
const OuterChipFnameTemplate* = "outer_chip_$1"
const ToTPerPixelFnameTemplate* = "tot_per_pixel_run$1_chip$2"
const ToTPerPixelTitleTemplate* = "ToT histogram for run $1 of chip $2"

#proc initDataPacket*(): DataPacket =
#  result.recvData = false
#  result.done = false
#  result.data = kstring""

proc initDataPacket*(kind: PacketKind = PacketKindUnknown,
                     rqKind: RequestKind = rqPing): DataPacket =
  result.kind = kind
  case result.kind
  of PacketKind.Descriptors, PacketKind.Plots, PacketKind.FileInfo:
    result.header = PacketHeader(msg: MessageUnknown,
                                 recvData: false,
                                 done: false)
  of PacketKind.Request:
    result.reqKind = rqKind
  of PacketKind.PacketKindUnknown: discard
  result.payload = ""

proc initDataPacket*(kind: PacketKind,
                     header: Messages,
                     reqKind: RequestKind = rqPing,
                     payload: kstring = "",
                     recvData: bool = false,
                     done: bool = false): DataPacket =
  result.kind = kind
  case result.kind
  of PacketKind.Descriptors, PacketKind.Plots, PacketKind.FileInfo:
    result.header = PacketHeader(msg: header,
                                 recvData: recvData,
                                 done: done)
  of PacketKind.Request:
    result.reqKind = reqKind
  of PacketKind.PacketKindUnknown: discard
  result.payload = payload

type
  DataPacketStorage* = object
    # object to store one DataPacket of each kind
    pltPacket*: DataPacket
    descPacket*: DataPacket
    reqPacket*: DataPacket
    fiPacket*: DataPacket

proc initDataPacketStorage*(): DataPacketStorage =
  result.pltPacket = initDataPacket(kind = PacketKind.Plots)
  result.descPacket = initDataPacket(kind = PacketKind.Descriptors)
  result.reqPacket = initDataPacket(kind = PacketKind.Request)
  result.fiPacket = initDataPacket(kind = PacketKind.FileInfo)

func `%`*(header: PacketHeader): JsonNode =
  ## serializes a `PacketHeader`
  result = newJObject()
  result[kstring"msg"] = % $(header.msg)
  # bools need to be stored as strings for JS interop
  result[kstring"recvData"] = % $(header.recvData)
  result[kstring"done"] = % $(header.done)

func `%`*(packet: DataPacket): JsonNode =
  ## serialize `packet` to send via websocket
  result = newJObject()
  result[kstring"kind"] = % $packet.kind
  case packet.kind
  of PacketKind.Descriptors, PacketKind.Plots, PacketKind.FileInfo:
    result[kstring"header"] = % packet.header
  of PacketKind.Request:
    result[kstring"reqKind"] = % $packet.reqKind
  of PacketKind.PacketKindUnknown: discard
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
  result.recvData = parseBool($(header["recvData"].getStr))
  result.done = parseBool($(header["done"].getStr))

proc parseDataPacket*(packet: kstring): DataPacket =
  ## unmarshals the data packet string
  # first create a JsonNode from the packet
  echo "Parsing packet ", packet
  let pJson = parseJson(packet)
  result.kind = parseEnum[PacketKind](pJson["kind"].getStr,
                                      PacketKind.PacketKindUnknown)
  case result.kind
  of PacketKind.Descriptors, PacketKind.Plots, PacketKind.FileInfo:
    result.header = parsePacketHeader(pJson[kstring"header"])
  of PacketKind.Request:
    result.reqKind = parseEnum[RequestKind](pJson["reqKind"].getStr,
                                            rqPing)
  of PacketKind.PacketKindUnknown: discard
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

func `%`*(c: GenericCut): JsonNode =
  result = %* { kstring"dset": % ($c.dset),
                kstring"lower": % ($c.lower),
                kstring"upper": % ($c.upper) }

func `%`*(s: DataSelector): JsonNode =
  result = newJObject()
  result[kstring"region"] = % $(s.region)
  result[kstring"cuts"] = % $(s.cuts)
  result[kstring"idxs"] = % $(s.idxs)
  result[kstring"applyAll"] = % $(s.applyAll)

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
    result[kstring"selector"] = % pd.selector
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

func parseGenericCut*(j: JsonNode): GenericCut =
  result = (dset: j[kstring"dset"].getStr,
            lower: j[kstring"lower"].getStr.parseFloat,
            upper: j[kstring"upper"].getStr.parseFloat)

func parseSelector*(j: JsonNode): DataSelector =
  result = DataSelector(region: parseEnum[ChipRegion](j[kstring"region"].getStr, crAll),
                        applyAll: parseBool(j[kstring"applyAll"].getStr))
  for jch in j[kstring"cuts"]:
    result.cuts.add parseGenericCut(jch)
  for idx in j[kstring"idxs"]:
    result.idxs.add parseInt(idx.getStr)

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
  result.selector = parseSelector(pd[kstring"selector"])
  case result.plotKind
  of pkInGridDset, pkFadcDset:
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

func `%`*(fileInfo: FileInfo): JsonNode =
  ## serializes a `FileInfo` object to a `JsonNode`
  result = newJObject()
  result["runs"] = % fileInfo.runs
  result["chips"] = % fileInfo.chips
  result["runType"] = % $fileInfo.runType
  result["rfKind"] = % $fileInfo.rfKind
  result["centerChip"] = % fileInfo.centerChip
  result["centerChipName"] = % fileInfo.centerChipName
  result["hasFadc"] = % $fileInfo.hasFadc

func `$`*(fileInfo: FileInfo): kstring =
  result = (% fileInfo).pretty

func parseFileInfo*(fi: JsonNode): FileInfo =
  ## parse a FileInfo stored as JsonNode back to an object
  result.runs = to(fi[kstring"runs"], seq[int])
  result.chips = to(fi[kstring"chips"], seq[int])
  result.runType = parseEnum[RunTypeKind](fi["runType"].getStr, rtNone)
  result.rfKind = parseEnum[RunFolderKind](fi["rfKind"].getStr, rfUnknown)
  result.centerChip = fi["centerChip"].getInt
  result.centerChipName = fi["centerChipName"].getStr
  result.hasFadc = parseBool($(fi["hasFadc"].getStr))

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

proc getXlabel*(pKind: PlotKind): string =
  ## returns the correct xlabel for a given plot kind
  discard

proc getYlabel*(pKind: PlotKind): string =
  ## returns the correct ylabel for a given plot kind
  discard

proc getPlotName*(pKind: PlotKind): string =
  ## returns the correct name for a given plot kind, if possible
  ## For pkInGrid or pkFadc this is not a unique name, because in those
  ## plots it's the name of the dataset we display.
  discard

proc getRunsStr*(runs: seq[int]): kstring =
  ## Returns a usable string of the given runs in the seq to use for
  ## a file name. If less than 4 runs, then simply all runs separated
  ## by underscores.
  ## If more than 3 runs, simply show a range of runs.
  if runs.len > 3:
    result = &"{runs.min}_{runs.max}"
  else:
    result = runs.foldl($a & " " & $b, "").strip(chars = {' '})

proc toFilename(s: DataSelector): string =
  let sCuts = s.cuts.mapIt(&"{it.dset}_{it.lower}_{it.upper}").join("_")
  let numIdxs = if s.idxs.len > 0: &"numIdxs_{s.idxs.len}" else: ""
  let sApplyAll = if s.cuts.len > 0: &"applyAll_{s.applyAll}" else: ""
  let sRegion = &"region_{s.region}"
  result = @[sRegion, sCuts, sApplyAll, numIdxs].join("_")

proc buildOutfile*(pd: PlotDescriptor, prefix, filetype: string): kstring =
  var name = ""
  let runsStr = getRunsStr(pd.runs)
  case pd.plotKind
  of pkInGridDset:
    name = InGridFnameTemplate %% [pd.name,
                                   runsStr,
                                   $pd.chip,
                                   $pd.binSize,
                                   $pd.binRange[0],
                                   $pd.binRange[1]]
  of pkAnyScatter:
    name = AnyScatterFnameTemplate %% [pd.name,
                                       runsStr,
                                       $pd.chip,
                                       $pd.x,
                                       $pd.y,
                                       $pd.color]
  of pkFadcDset:
    name = FadcFnameTemplate %% [pd.name,
                                 runsStr,
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
    name = PhotoVsTimeFnameTemplate %% [runsStr,
                                        $pd.splitBySec]
  of pkFeChVsTime:
    name = PhotoChargeVsTimeFnameTemplate %% [runsStr,
                                              $pd.splitBySec]
  of pkFePixDivChVsTime:
    name = PhotoPixDivChVsTimeFnameTemplate %% [runsStr]
  of pkFePhotoDivEscape:
    name = PhotoDivEscapeFnameTemplate %% [runsStr]
  of pkInGridEvent:
    name = InGridEventFnameTemplate %% [runsStr,
                                        $pd.chip,
                                        $pd.event]
  of pkFadcEvent:
    name = FadcEventFnameTemplate %% [runsStr,
                                      $pd.event]

  of pkSubPlots:
    for p in pd.plots:
      name &= pd.buildOutfile("", "")
  of pkOuterChips:
    name &= OuterChipFnameTemplate %% [$pd.runType]
  of pkToTPerPixel:
    name = ToTPerPixelFnameTemplate %% [runsStr,
                                        $pd.chip]
  else:
    discard
  # now add generic selector
  name.add "_" & toFilename(pd.selector)
  result = &"{prefix}/{name}.{filetype}"

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
                                     $pd.chip]
  of pkAnyScatter:
    result = AnyScatterTitleTemplate %% [pd.name,
                                       runsStr,
                                       $pd.chip,
                                       $pd.x,
                                       $pd.y,
                                       $pd.color]
  of pkFadcDset:
    result = FadcTitleTemplate %% [pd.name,
                                  runsStr]

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
  of pkFePhotoDivEscape:
    result = PhotoDivEscapeTitleTemplate %% [runsStr]
  of pkInGridEvent:
    result = InGridEventTitleTemplate %% [runsStr,
                                          $pd.chip,
                                          $pd.event]
  of pkFadcEvent:
    result = FadcEventTitleTemplate %% [runsStr,
                                      $pd.event]
  of pkToTPerPixel:
    result = ToTPerPixelTitleTemplate %% [runsStr,
                                          $pd.chip]
  of pkSubPlots:
    result = ""
    for p in pd.plots:
      result &= pd.buildTitle()
  else:
    discard
