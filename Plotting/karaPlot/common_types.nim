import karax / kbase
import sets, tables

import ingrid / ingrid_types

type
  ## enum listing all available `plot types` we can produce
  PlotKind* = enum
    pkInGridDset           ## histogram InGrid property
    pkFadcDset             ## histogram FADC property
    pkPolya                ## InGrid polya distribution
    pkCombPolya            ## combined polya of all chips
    pkOccupancy            ## Occupancy of InGrid chip
    pkOccCluster           ## Occupancy of clusters of InGrid chip
    pkFeSpec               ## Fe pixel (or different) spectrum
    pkEnergyCalib          ## Energy calibration from Fe pixel spectrum
    pkFeSpecCharge         ## Fe charge (or different) spectrum
    pkEnergyCalibCharge    ## Energy calibration from Fe charge spectrum
    pkFeVsTime             ## Evolution of Fe pix peak location vs time
    pkFeChVsTime           ## Evolution of Fe charge peak location vs time
    pkFePixDivChVsTime     ## Evolution of Fe (pix peak / charge peak) location vs time"
    pkFePhotoDivEscape     ## Evolution of Fe photo peak / escape peak location vs time
    pkInGridEvent          ## Individual InGrid event
    pkFadcEvent            ## Individual FADC event
    pkCalibRandom          ## ? to be filled for different calibration plots
    pkAnyScatter           ## Scatter plot of some x vs. some y
    pkMultiDset            ## Plot of multiple histograms. Will be removed and replaced
                           ## by redesign of `createPlot`
    pkSubPlots             ## several subplots in one plot
    pkInGridCluster        ## superseeded by pkInGridEvent?
    pkOuterChips           ## histogram of # hits of outer chips
    pkToTPerPixel          ## histogram of ToT values. Essentially Polya but with raw ToT
                           ## Jochen always calls it ToT/Pixel...

  ClampKind* = enum
    ckFullRange, ckAbsolute, ckQuantile

  DataKind* = enum
    dkInGrid, dkFadc

  CutRange* = tuple[low, high: float, name: kstring]

  Domain* = tuple
    left, bottom, width, height: float

  PlotDescriptor* = object
    runType*: RunTypeKind
    name*: kstring
    runs*: seq[int]
    chip*: int
    isCenterChip*: bool
    xlabel*: kstring
    ylabel*: kstring
    title*: kstring
    # bKind: BackendKind <- to know which backend to use for interactive plot creation
    case plotKind*: PlotKind
    of pkInGridDset, pkFadcDset, pkToTPerPixel:
      range*: CutRange
      cutRegion*: ChipRegion
      # optional fields for bin size and range
      binSize*: float
      binRange*: tuple[low, high: float]
    of pkAnyScatter:
      # read any dataset as X and plot it against Y
      x*: kstring
      y*: kstring
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
    of pkInGridEvent, pkFadcEvent:
      # events*: OrderedSet[int] # events to plot (indices at the moment, not event numbers)
      event*: int # the current event being plotted
    of pkFeVsTime, pkFeChVsTime, pkFePixDivChVsTime:
      # If unequal to 0 will create the plot not just split by runs, but rather split the
      # calib data for each run in pieces of `splitBySec` seconds of time slices.
      splitBySec*: int
      # allowed divergence of last slice's length in percent
      lastSliceError*: float
      # if splitBySec doesn't fit into splitBySec within `lastSliceError` decide if to drop
      # that slice or keep it
      dropLastSlice*: bool
    of pkSubPlots:
      # a way to combine several plots into a single plot of subplots
      plots*: seq[PlotDescriptor]
      domain*: seq[Domain] # relative location within [0, 1] of the
                           # plot canvas for each subplot
    of pkOuterChips:
      outerChips*: seq[int] # seq of all chips considered "outer"
      rangeCenter*: CutRange
    else:
      discard
