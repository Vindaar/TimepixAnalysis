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

  # enum listing all available `plot types` we can produce
  PlotKind* = enum
    pkInGridDset           # histogram InGrid property
    pkFadcDset             # histogram FADC property
    pkPolya                # InGrid polya distribution
    pkCombPolya            # combined polya of all chips
    pkOccupancy            # Occupancy of InGrid chip
    pkOccCluster           # Occupancy of clusters of InGrid chip
    pkFeSpec               # Fe pixel (or different) spectrum
    pkEnergyCalib          # Energy calibration from Fe pixel spectrum
    pkFeSpecCharge         # Fe charge (or different) spectrum
    pkEnergyCalibCharge    # Energy calibration from Fe charge spectrum
    pkFeVsTime             # Evolution of Fe pix peak location vs tim
    pkFePixDivChVsTime     # Evolution of Fe (pix peak / charge peak) location vs tim
    pkInGridEvent          # Individual InGrid event
    pkFadcEvent            # Individual FADC event
    pkCalibRandom          # ? to be filled for different calibration plots
    pkAnyScatter           # Scatter plot of some x vs. some y
    pkMultiDset            # Plot of multiple histograms. Will be removed and replaced
                           # by redesign of `createPlot`
    pkInGridCluster        # superseeded by pkInGridEvent?

  ClampKind* = enum
    ckFullRange, ckAbsolute, ckQuantile

  DataKind* = enum
    dkInGrid, dkFadc

  CutRange* = tuple[low, high: float, name: string]

  PlotDescriptor* = object
    runType*: RunTypeKind
    name*: string
    runs*: seq[int]
    chip*: int
    xlabel*: string
    ylabel*: string
    title*: string
    # bKind: BackendKind <- to know which backend to use for interactive plot creation
    case plotKind*: PlotKind
    of pkInGridDset, pkFadcDset:
      range*: CutRange
    of pkAnyScatter:
      # read any dataset as X and plot it against Y
      x*: string
      y*: string
    of pkMultiDset:
      # histogram of all these datasets in one
      names*: seq[string]
    of pkInGridCluster:
      eventNum*: int
    of pkOccupancy, pkOccCluster:
      case clampKind*: ClampKind
      of ckAbsolute:
        # absolute clamp tp `clampA`
        clampA*: float
      of ckQuantile:
        # clamp to `clampQ` quantile
        clampQ*: float
      of ckFullRange:
        # no field for ckFullRange
        discard
    of pkCombPolya:
      chipsCP*: seq[int]
    of pkInGridEvent:
      events*: OrderedSet[int] # events to plot (indices at the moment, not event numbers)
      event*: int # the current event being plotted
    else:
      discard



const FakeFrameSize* = 32760

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
