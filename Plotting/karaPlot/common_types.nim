import karax / kbase
import sets, tables

import ingrid / ingrid_types

when defined(js):
  # just some opaque thing
  type
    H5File* = object
    Config* = object
    PlotV* = object
else:
  from nimhdf5 import H5File
  from nimpy import PyObject
  import plotly
  import ggplotnim
  import ingrid / tos_helpers
  type
    CustomPlotKind* = enum
      cpScatter, cpHistogram #, ...
    ## For now: simple object storing x/y datasets to plot against in scatter
    ## From it a correct `PlotDescriptor` will be generated once the file is known
    CustomPlot* = object
      kind*: CustomPlotKind
      x*: kstring
      y*: kstring     # not used for histogram
      color*: kstring

    OutputFiletypeKind* = enum
      ofUnknown,
      ofJson = "json"
      ofOrg = "org"

    PlotFiletypeKind* = enum
      ftSvg = "svg"
      ftPdf = "pdf"
      ftPng = "png"

    ConfigFlagKind* = enum
      cfNone, cfFadc, cfInGrid, cfOccupancy, cfPolya, cfFeSpectrum, cfTotPerPixel, cfProvideServer, cfShow,
      cfApplyAllCuts,  # if true, will apply all cuts to all datasets. If not given (default) a cut on
                       # a single dataset will only apply on that dataset
      cfCutFePeak,     # if true and input is calibration, will create plots cut to Photopeak & Escape peak
      cfCompiledCustom # if given will create the compiled custom plots (i.e. in file `moreCustomPlots`).
                       # A change requires a recompilation!

    ## A generic cut on input data using dset & low / high values
    GenericCut* = object
      applyFile*: seq[string] ## apply this cut to all files in this seq (all if empty)
      applyDset*: seq[string] ## apply this cut when reading all datasets in this seq (all if empty)
      dset*: string
      min*: float
      max*: float
    ChipCoord* = range[0 .. 255]
    MaskRegion* = object
      applyFile*: seq[string] ## apply this cut to all files in this seq (all if empty)
      applyDset*: seq[string] ## apply this cut when reading all datasets in this seq (all if empty)
      x*: tuple[min, max: ChipCoord]
      y*: tuple[min, max: ChipCoord]

    Config* = object
      flags*: set[ConfigFlagKind]
      chips*: set[uint16]
      runs*: set[uint16]
      outputType*: OutputFiletypeKind # the file format to use for the file containing all plots
      fileType*: PlotFiletypeKind
      ingridDsets*: set[IngridDsetKind]
      fadcDsets*: seq[string] # currently don't have an enum for them
      cuts*: seq[GenericCut] ## Used to fill the `DataSelector`
      maskRegion*: seq[MaskRegion] ## Used to fill the `DataSelector`
      region*: ChipRegion    ## From input to preselect a region
      idxs*: seq[int]        ## Indices to read. Negative indices are interpreted as seen from the end of dset
      plotlySaveSvg*: bool
      customPlots*: seq[CustomPlot]
      cdlGroup*: string ## The (optional) CDL dataset (group, target/filter kind) from which
                        ## to read data if an input file is a CDL file (`calibration-cdl*.h5`)

    PlottingBackendKind* = enum
      bNone, bMpl, bPlotly, bGgPlot

    # variant object for the layout combining both
    # TODO: make generic or always use float?
    PlotV* = object
      annotations*: seq[string]
      invalid*: bool
      case kind*: PlottingBackendKind
      of bMpl:
        # what needs to go here?
        plt*: PyObject
        fig*: PyObject
        ax*: PyObject
      of bPlotly:
        plLayout*: Layout
        plPlot*: Plot[float]
        plPlotJson*: PlotJson
      of bGgPlot:
        pltGg*: GgPlot
        width*: float
        height*: float
        theme*: Theme
      else: discard


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
    pkCustomPlot           ## Custom plot of: scatter: some x vs. some y, histogram, ...
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

  ## Helper that is used to store combined cuts as well as a region
  DataSelector* = object
    region*: ChipRegion ## Region defaults to full chip
    cuts*: seq[GenericCut]
    maskRegion*: seq[MaskRegion]
    idxs*:  seq[int] ## Optional indices. If given we extract *these* indices from the
                    ## resulting data after cuts are applied. Used to implement `head`, `tail`
    applyAll*: bool ## set to indicate to apply all cuts instead of matching dataset names

  Domain* = tuple
    left, bottom, width, height: float

  ## The type corresponding to the procedures that
  PlotHandlerProc = proc(h5f: H5File, fileInfo: FileInfo, pd: PlotDescriptor, config: Config): (string, PlotV)

  PlotDescriptor* = object
    runType*: RunTypeKind
    name*: kstring
    runs*: seq[int]
    chip*: int
    isCenterChip*: bool
    xlabel*: kstring
    ylabel*: kstring
    title*: kstring
    selector*: DataSelector
    # bKind: BackendKind <- to know which backend to use for interactive plot creation
    case plotKind*: PlotKind
    of pkInGridDset, pkFadcDset, pkToTPerPixel:
      # optional fields for bin size and range
      binSize*: float
      binRange*: tuple[low, high: float]
    of pkCustomPlot:
      # read any dataset as X and plot it against Y
      customPlot*: CustomPlot
      processData*: PlotHandlerProc
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
    else:
      discard

proc `==`*(s1, s2: DataSelector): bool =
  result = s1.region == s2.region and
     s1.cuts == s2.cuts and
     s1.maskRegion == s2.maskRegion and
     s1.idxs == s2.idxs and
     s1.applyAll == s2.applyAll

proc `==`*(p1, p2: PlotDescriptor): bool =
  ## `runType` is ignored to allow to compare rtBackground with rtCalibration files!
  #if p1.runType != p2.runType: result = false
  if p1.plotKind != p2.plotKind: result = false
  else:
    result = true # start with true
    template cmpField(f: untyped): untyped =
      result = result and (p1.f == p2.f)
    cmpField(name)
    #cmpField(runs) # `runs` is *not* compared as they don't describe the *kind* of plot
    cmpField(isCenterChip)
    if not p1.isCenterChip and not p2.isCenterChip: # only compare the chip numbers if they are not
      cmpField(chip)                                # center chips. Else
    #cmpField(xlabel) # labels irrelevant
    #cmpField(ylabel)
    #cmpField(title) # title irrelevant

    cmpField(selector)
    case p1.plotKind
    of pkInGridDset, pkFadcDset, pkToTPerPixel:
      cmpField(binSize)
      cmpField(binRange)
    of pkCustomPlot:
      cmpField(customPlot)
    of pkMultiDset:
      cmpField(names)
    of pkInGridCluster:
      cmpField(eventNum)
    of pkOccupancy, pkOccCluster:
      if p1.clampKind != p2.clampKind: result = false
      else:
        case p1.clampKind
        of ckAbsolute:
          cmpField(clampA)
        of ckQuantile:
          cmpField(clampQ)
        else: discard
    of pkCombPolya:
      cmpField(chipsCP)
    of pkInGridEvent, pkFadcEvent:
      cmpField(event)
    of pkFeVsTime, pkFeChVsTime, pkFePixDivChVsTime:
      cmpField(splitBySec)
      cmpField(lastSliceError)
      cmpField(dropLastSlice)
    of pkSubPlots:
      cmpField(plots)
      cmpField(domain)
    of pkOuterChips:
      cmpField(outerChips)
    else:
      discard
