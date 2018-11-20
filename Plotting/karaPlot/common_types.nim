import karax / kbase
import sets, tables

import ingrid / ingrid_types

type
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

  CutRange* = tuple[low, high: float, name: kstring]

  PlotDescriptor* = object
    runType*: RunTypeKind
    name*: kstring
    runs*: seq[int]
    chip*: int
    xlabel*: kstring
    ylabel*: kstring
    title*: kstring
    # bKind: BackendKind <- to know which backend to use for interactive plot creation
    case plotKind*: PlotKind
    of pkInGridDset, pkFadcDset:
      range*: CutRange
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
    of pkInGridEvent:
      events*: OrderedSet[int] # events to plot (indices at the moment, not event numbers)
      event*: int # the current event being plotted
    else:
      discard
