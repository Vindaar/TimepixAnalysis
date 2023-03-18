# module which contains the used type definitions in the InGrid module
when not defined(js):
  import times
import tables, strformat, memfiles, options
import karax / kbase

type
  EventHeader* = Table[string, string]
  ChipHeader*  = Table[string, string]
  Pix*         = tuple[x, y: uint8, ch: uint16]
  PixTpx3*     = tuple[x, y: uint8, ch, toa: uint16, toaCombined: uint64]
  # Integer based pixels are used for full Septemboard frames, due to 3x256 pixels per direction
  PixInt*      = tuple[x, y: int, ch: int]
  PixIntTpx3*  = tuple[x, y, ch: int, toa: uint16, toaCombined: uint64]
  SomePix*     = Pix | PixInt | PixTpx3 | PixIntTpx3
  Pixels*      = seq[Pix]
  PixelsInt*   = seq[PixInt]

  # Coord type which contains (x, y) coordinates of a pixel
  Coord*[T] = tuple[x, y: T]
  # cluster object
  Cluster*[T: SomePix] = seq[T]

  Chip* = tuple[name: string, number: int]

  TimepixVersion* = enum
    Timepix1, Timepix3

  ChipEvent* = object
    chip*: Chip
    pixels*: Pixels
    case version*: TimepixVersion
    of Timepix1: discard
    of Timepix3:
      toa*: seq[uint16]         # Time of Arrival (ToA) in "local" time
      toaCombined*: seq[uint64] # ToA extended by "run wide" counter

  ProtoFile* = object
    name*: string
    fileData*: string

  Event* = object
    isValid*: bool
    evHeader*: Table[string, string]
    chips*: seq[ChipEvent]
    nChips*: int
    # time the shutter was open in seconds
    length*: float

  Run* = object
    # a run stores raw run data of a sequence of events
    # and information about the chips
    events*: seq[Event]
    chips*: seq[Chip]
    runNumber*: int
    runHeader*: Table[string, string]

  # define a distinct `OldEvent` to differentiate the (in principle) not
  # different old TOS storage format
  OldEvent* = Event
  # alias for SRS events for the Event. Only difference is meta information
  SrsEvent* = Event

  HistTuple* = tuple[bins: seq[float64], hist: seq[float64]]

  #############################
  # Calibration related types #
  #############################

  Tot* = object
    pulses*: seq[float]
    mean*: seq[float]
    std*: seq[float]

  SCurve* = object
    name*: string
    voltage*: int
    thl*: seq[int]
    hits*: seq[float]

  SCurveSeq* = object
    files*: seq[string]
    curves*: seq[SCurve]

  FSR* = Table[string, int]

  ######################
  # File related enums #
  ######################

  EventSortType* = enum
    fname, inode

  EventType* = enum
    FadcType, InGridType

  RunTypeKind* = enum
    rtNone, rtCalibration, rtBackground, rtXrayFinger

  RunFolderKind* = enum
    ## TODO: add `rfTpx3Daq` and replace current usage of `rfUnknown`
    rfNewTos, rfOldTos, rfSrsTos, rfUnknown

  # The enum which determines how the gas gain vs energy calbration is done
  # and which values to use for the energy calibration
  GasGainVsChargeCalibKind* = enum
    gcNone = ""
    gcMean = "Mean"
    gcIndividualFits = "Individual"

  ## selection of the cluster finding algorithm
  ClusteringAlgorithm* = enum
    caDefault = "default" ## our default dumb r pix radius algorithm
    caDBSCAN = "dbscan"   ## DBSCAN as the clustering algo

  ## Type that stores the data files (the file names!) of a single run
  DataFiles* = object
    files*: seq[string] ## path to the file including filename
    eventNumbers*: seq[int] ## evnet numbers of all the files
    kind*: EventType ## Type of event (might be a temperature
    rfKind*: RunFolderKind
    temperatureLog*: string ## Optional temperature log file, if available.

  ## Entry of a temperature log. It's a simple object like this to write it to H5
  ## as a compound dataset
  TemperatureLogEntry* = object
    imb*: float
    septem*: float
    timestamp*: int64

  TemperatureLog* = seq[TemperatureLogEntry]

when not defined(js):
  type
    # an object, which stores information about a run's start, end and length
    RunTimeInfo* = object
      t_start*: Time
      t_end*: Time
      t_length*: Duration # total duration of the run / tracking, from first to last event!

    # an object which stores general information about a run
    RunInfo* = object
      timeInfo*: RunTimeInfo
      runNumber*: int
      rfKind*: RunFolderKind
      runType*: RunTypeKind
      path*: string
      nEvents*: int
      nFadcEvents*: int

    # extension of the above read from H5 file, including effective
    # run times, possible trackings etc
    ExtendedRunInfo* = object
      timeInfo*: RunTimeInfo
      runNumber*: int
      rfKind*: RunFolderKind
      runType*: RunTypeKind
      nEvents*: int
      nFadcEvents*: int
      totalTime*: Duration # total time from beginning of each run to endy
      activeTime*: Duration # total time shutter was open
      activeRatio*: float # ratio of total time the shutter was open
      # reuse `RunTimeInfo` to store possible tracking starts / ends
      trackings*: seq[RunTimeInfo]
      nonTrackingDuration*: Duration
      activeNonTrackingTime*: Duration # total time shutter was open during no tracking
      trackingDuration*: Duration
      activeTrackingTime*: Duration # total time shutter was open during tracking

  ################################
  # Reconstruction related types #
  ################################
type

  ## The base element for a single event to be reconstructed
  RecoInputEvent*[T: SomePix] = tuple
    pixels: seq[T] # SomePix
    eventNumber: int
    toa: seq[uint16]
    toaCombined: seq[uint64]
  RecoInputData*[T: SomePix] = seq[RecoInputEvent[T]]

  # object which stores the geometry information of a single
  # `ClusterObject`
  ClusterGeometry* = object
    rmsLongitudinal*: float
    rmsTransverse*: float
    eccentricity*: float
    rotationAngle*: float
    skewnessLongitudinal*: float
    skewnessTransverse*: float
    kurtosisLongitudinal*: float
    kurtosisTransverse*: float
    length*: float
    width*: float
    fractionInTransverseRms*: float
    lengthDivRmsTrans*: float

  ## An object that stores information that describes the "geometry" of the
  ## ToA data in an object
  ToAGeometry* = object
    toaLength*: float   # length in ToA ## XXX: ideally this would be a small integer type, i.e. int32
    toaMin* : uint16    # the minimal ToA value found before subtraction
    toaMean*: float     # mean ToA value of the ToA distribution of the cluster
    toaRms* : float     # RMS of the ToA distribution of the cluster
    toaSkewness*: float # skewness of the ToA distribution of the cluster
    toaKurtosis*: float # kurtosis of the ToA distribution of the cluster

  # object which stores a single `Cluster` in combination with information
  # about itself, e.g. energy, geometry etc.
  ClusterObject*[T: SomePix] = object
    data*: Cluster[T]
    hits*: int
    centerX*: float
    centerY*: float
    # total tot in the whole cluster
    sumTot*: int
    energy*: float
    geometry*: ClusterGeometry
    case version*: TimepixVersion
    of Timepix1: discard # nothing extra
    of Timepix3:
      toa*: seq[uint16]
      toaCombined*: seq[uint64]
      toaGeometry*: ToAGeometry

  # object which stores information about a reconstructed event, i.e.
  # split into different clusters and information about it, chip and
  # event number (run number is left out, because it will be stored in
  # the group of a run anyways)
  RecoEvent*[T: SomePix] = object
    cluster*: seq[ClusterObject[T]]
    event_number*: int
    chip_number*: int

  ################################
  #### Analysis related types ####
  ################################

  ChipRegion* = enum
    crGold, crSilver, crBronze, crAll

  CutsKind* = enum
    ckReference, ckXray

  # variant object to store cut values to either get
  # the reference spectra from the Xray spectra or
  # to build the Xray spectra from the raw calibration-cdl.h5
  # data
  # NOTE: in principle one would combine the variant
  # object into one! This is only not done, to reflect the
  # 2 stage process described in Christoph's PhD thesis
  Cuts* = object
    minRms*: float
    maxRms*: float
    maxLength*: float
    minPix*: float
    case kind*: CutsKind
    of ckReference:
      minCharge*: float
      maxCharge*: float
    of ckXray:
      maxEccentricity*: float
      # we also cut to the silver region
      cutTo*: ChipRegion

  # object to store region based cuts (gold, silver, bronze region)
  CutsRegion* = object
    xMin*: float
    xMax*: float
    yMin*: float
    yMax*: float
    radius*: float

  # type to store results of fitting with mpfit / NLopt / Python
  FitResult* = object
    xMin*: float # the actual fit range minimum
    xMax*: float # the actual fit range maximum
    x*: seq[float]
    y*: seq[float]
    pRes*: seq[float]
    pErr*: seq[float]
    redChiSq*: float
    resText*: string

  TpaFileKind* = enum
    tpkRawData = "/runs"           # output of `raw_data_manipulation`
    tpkReco = "/reconstruction"    # output of `reconstruction`
    tpkLogL = "/likelihood"        # output of `likelihood`
    tpkTpx3Raw = ""                # raw Tpx3 data straight from the Tpx3 DAQ. Has no base group
    tpkTpx3Interp = "/interpreted" # interpreted raw data, output of `parse_raw_tpx3`
    tpkCDL = "CDL"                 # CDL data files (reconstructed by TPA)

  # a simple object storing the runs, chips etc. from a given
  # H5 file
  FileInfo* = object
    name*: string # the name of the file
    runs*: seq[int]
    chips*: seq[int]
    runType*: RunTypeKind
    rfKind*: RunFolderKind
    centerChip*: int
    centerChipName*: kstring
    hasFadc*: bool # reads if FADC group available
    timepix*: TimepixVersion
    case tpaFileKind*: TpaFileKind
    of tpkCDL:
      cdlGroup*: string ## An (optional) selector of a specific CDL group (target / filter type)
                        ## if the input file corresponds to a CDL file. In that case the data
                        ## path returned for this file will be from that dataset.
    else: discard


  ## Tpx3Data is the compound data type stored in the "interpreted" branch of the
  ## H5 files after running "analyse_data.py" from `tpx3-daq`.
  ## IMPORTANT: At the moment the order ``matters``! It must be the same order
  ## as in the file, because the data is parsed into this structure from binary
  ## using the data sizes of the types as offsets.
  Tpx3Data* = object
    data_header*: uint8
    header*: uint8
    hit_index*: uint64
    x*: uint8
    y*: uint8
    TOA*: uint16
    TOT*: uint16
    EventCounter*: uint16
    HitCounter*: uint8
    FTOA*: uint8
    scan_param_id*: uint16
    chunk_start_time*: cdouble
    iTOT*: uint16
    TOA_Extension*: uint64
    TOA_Combined*: uint64

  ## The `run_config` dataset in the Tpx3 H5 files is stored as a compound
  ## type of this
  Tpx3RunConfigRaw* = object
    ## IMPORTANT: I don't know if these sizes may ever change!
    attribute*: array[64, char]
    value*: array[128, char]


  ## The above merged into a sane typed object
  Tpx3RunConfig* = object
    scanId*: string ## why is this a string?
    runName*: string
    softwareVersion*: string
    boardName*: string
    firmwareVersion*: int
    chipWafer*: int   ## The number of the wafer
    chipX*: char      ## The 'X' coordinate on the wafer, a character
    chipY*: int       ## The 'Y' coordinate on the wafer, a number
    thrFile*: string  ## Path to the threshold file
    maskFile*: string ## Path to the mask file

  Tpx3DacsRaw* = object
    ## IMPORTANT: I don't know if these sizes may ever change!
    attribute*: array[64, char]
    value*: uint16

  Tpx3Dacs* = object
    Ibias_Preamp_ON*:   uint16
    Ibias_Preamp_OFF*:  uint16
    VPreamp_NCAS*:      uint16
    Ibias_Ikrum*:       uint16
    Vfbk*:              uint16
    Vthreshold_fine*:   uint16
    Vthreshold_coarse*: uint16
    Ibias_DiscS1_ON*:   uint16
    Ibias_DiscS1_OFF*:  uint16
    Ibias_DiscS2_ON*:   uint16
    Ibias_DiscS2_OFF*:  uint16
    Ibias_PixelDAC*:    uint16
    Ibias_TPbufferIn*:  uint16
    Ibias_TPbufferOut*: uint16
    VTP_coarse*:        uint16
    VTP_fine*:          uint16
    Ibias_CP_PLL*:      uint16
    PLL_Vcntrl*:        uint16
    Sense_DAC*:         uint16

  ## GasGainIntervalData stores the information about binning the data for the
  ## calculation of the gas gain in a run by a fixed time interval, e.g.
  ## 100 minutes
  GasGainIntervalData* = object
    idx*: int # index
    interval*: float # interval length in minutes
    minInterval*: float # minimum length of an interval in min (end of run)
    tStart*: int # timestamp of interval start
    tStop*: int # timestamp of interval stop

  # simple object duplicating fields of `GasGainIntervalData`, because this type will
  # be stored as a composite type in the H5 file
  GasGainIntervalResult* = object
    idx*: int
    interval*: float
    minInterval*: float # minimum length of an interval in min (end of run)
    tStart*: int
    tStop*: int
    tLength*: float # actual length of the slice in min: `(tStop - tStart) / 60`
    sliceStart*: int # start and stop of the indices which belong to this slice
    sliceStop*: int
    sliceStartEvNum*: int # start and stop of the indices which belong to this slice
    sliceStopEvNum*: int
    numClusters*: int # the actual number of clusters in this slice
    N*: float # amplitude (?) of the polya fit, parameter 0
    G_fit*: float # gas gain as determined from fit, parameter 1
    G_fitmean*: float # gas gain as determined from mean of fit "distribution"
    G*: float # fit parameter as determined from mean of data histogram
    theta*: float # theta parameter of polya fit, parameter 2
    redChiSq*: float # reduced χ² of the polya fit of this slice

  FeSpecFitData* = object
    hist*: seq[float]
    binning*: seq[float]
    idx_kalpha*: int
    idx_sigma*: int
    kalpha*: float
    sigma_kalpha*: float
    pRes*: seq[float]
    pErr*: seq[float]
    xFit*: seq[float]
    yFit*: seq[float]
    chiSq*: float
    nDof*: int

  EnergyCalibFitData* = object
    energies*: seq[float]
    peaks*: seq[float]
    peaksErr*: seq[float]
    pRes*: seq[float]
    pErr*: seq[float]
    xFit*: seq[float]
    yFit*: seq[float]
    aInv*: float
    aInvErr*: float
    chiSq*: float
    nDof*: int

  ## Object which stores the technique used to morph between different CDL reference spectra.
  ## Set in `config.toml`.
  MorphingKind* = enum
    mkNone = "None"
    mkLinear = "Linear"

  YearKind* = enum
    yr2014 = "2014"
    yr2018 = "2018"

  ## Types used in the `cdl_spectrum_creation`

  TargetKind* = enum
    tEmpty = ""
    tCu = "Cu"
    tMn = "Mn"
    tTi = "Ti"
    tAg = "Ag"
    tAl = "Al"
    tC = "C"

  FilterKind* = enum
    fEmpty = ""
    fEpic = "EPIC"
    fCr = "Cr"
    fNi = "Ni"
    fAg = "Ag"
    fAl = "Al"
    fTi = "Ti"

  TargetFilterKind* = enum
    tfCuNi15 = "Cu-Ni-15kV"
    tfMnCr12 = "Mn-Cr-12kV"
    tfTiTi9 = "Ti-Ti-9kV"
    tfAgAg6 = "Ag-Ag-6kV"
    tfAlAl4 = "Al-Al-4kV"
    tfCuEpic2 = "Cu-EPIC-2kV"
    tfCuEpic0_9 = "Cu-EPIC-0.9kV"
    tfCEpic0_6 =  "C-EPIC-0.6kV"

  CdlRun* = object
    number*: int
    runType*: RunTypeKind
    hasFadc*: bool
    target*: TargetKind
    filter*: FilterKind
    hv*: float

  ## This stupidly named object stores the parameters for the stretching of the
  ## CDL data to suit the detector's data based on a set of 55Fe data.
  LineParams* = tuple[m, b: float]
  CdlStretch* = object
    fe55*: string ## file name used to fill this object
    eccEsc*: float ## eccentricity value determined best fit for the escape peak
    eccPho*: float ## eccentricity value determined best fit for the photo peak
    ## maps each property to its slope / intercept pair for the minima / maxima of
    ## that property so that we can compute the correct needed parameters for the
    ## stretching on the fly for any energy
    lines*: Table[string, tuple[mins, maxs: LineParams]]

    #fns: seq[FormulaNode]
    # something something min max of 55fe data?

  InGridDsetKind* = enum
    igInvalid, # invalid dataset
    igCenterX,
    igCenterY,
    igHits, # hits, NumberOfPixels
    igEventNumber,
    igEccentricity, # eccentricity, excentricity
    igSkewnessLongitudinal, # skewnessLongitudinal, SkewnessLongitudinal
    igSkewnessTransverse, # skewnesssTransverse, SkewnessTransverse
    igKurtosisLongitudinal, # kurtosisLongitudinal, KurtosisLongitudinal
    igKurtosisTransverse, # kurtosisTransverse, KurtosisTransverse
    igLength, # length, Length
    igWidth, # width, Width
    igRmsLongitudinal, # rmsLongitudinal, RmsLongitudinal
    igRmsTransverse, # rmsTransverse, RmsTransverse
    igLengthDivRmsTrans, # lengthDivRmsTrans,
    igRotationAngle, # rotationAngle, RotationAngle
    igEnergyFromCharge, # energyFromCharge, Energy
    igEnergyFromPixel, # energyFromPixel, Energy
    igLikelihood, # likelihood, LikelihoodMarlin
    igFractionInTransverseRms, # fractionInTransverseRms, FractionWithinRmsTransverse
    igTotalCharge, # totalCharge, TotalCharge
    igNumClusters, # number of clusters in a single event; Marlin: xrayperevent
    igFractionInHalfRadius, # fraction of pixels in half the radius
    igRadiusDivRmsTrans, # "radiusdivbyrmsy"
    igRadius, # radius of cluster
    igBalance, # balance ?
    igLengthDivRadius

  FrameworkKind* = enum
    fkTpa, fkMarlin

  ## A generic cut on input data using dset & low / high values
  GenericCut* = object
    applyFile*: seq[string] ## apply this cut to all files in this seq (all if empty)
    applyDset*: seq[string] ## apply this cut when reading all datasets in this seq (all if empty)
    dset*: string
    isFadc*: bool # true if dataset is an FADC dataset. ``Must`` contain `"fadc/*"` prefix!
    min*: float
    max*: float
    inverted*: bool ## If true the cut is inverted, i.e. we remove everything *in* the cut

  ## A helper distinct version of an unchecked array. Mainly used in the septem veto related
  ## handling of all chip data
  DataView*[T] = distinct ptr UncheckedArray[T]

## Accessor and converter of a `ptr` to a `DataView` (all we need)
func `[]`*[T](dv: DataView[T], idx: int): T = cast[ptr UncheckedArray[T]](dv)[idx]
#proc `=copy`[T](dv1, dv2: DataView[T]) {.error: "Copying a data view is forbidden!".}
func toDataView*[T](p: ptr seq[T]): DataView[T] = DataView[T](cast[ptr UncheckedArray[T]](p[0].addr))


const TosDateString* = "yyyy-MM-dd'.'hh:mm:ss"

# and some general InGrid related constants
const NPIX* = 256
const PITCH* = 0.055
const TimepixSize* = NPIX * PITCH

const SrsRunIncomplete* = "incomplete"
const SrsRunIncompleteMsg* = "This run does not contain a run.txt and so " &
  "is incomplete!"
const SrsNoChipId* = "ChipIDMissing"
const SrsNoChipIdMsg* = "The chip IDs are missing from the run.txt. Old format!"
const SrsDefaultChipName* = "SRS Chip"

const TpaIngridDsetKinds* = block:
  var dsets = newSeq[IngridDsetKind]()
  for dset in IngridDsetKind:
    if dset notin { igInvalid, igNumClusters, igLikelihood,
                    igFractionInHalfRadius, igRadiusDivRmsTrans, igRadius, igBalance, igLengthDivRadius }:
      dsets.add dset
  dsets


# the following will not be available, if the `-d:pure` flag is set,
# to allow importing the rest of the types, without a `arraymancer`
# dependency
when not defined(pure) and not defined(js):
  import arraymancer, datamancer
  type
    Threshold* = Tensor[int]
    ThresholdMeans* = Tensor[int]

    ToTCut* = object
      low*: int ## *exclusive* low cut `(x < low) -> remove`
      high*: int ## *exclusive* high cut `(x > high) -> remove`
      rmLow*: int ## number of elements removed due to being too low
      rmHigh*: int ## number of elements removed due to being too high

    # process events stores all data for septemboard
    # of a given run
    ProcessedRun* = object
      # indicates which Timepix was/were used to take this run
      timepix*: TimepixVersion
      # just the number of chips in the run
      nChips*: int
      # the chips as (name, number) tuples
      chips*: seq[Chip]
      # run number
      runNumber*: int
      # table containing run header ([General] in data file)
      runHeader*: Table[string, string]
      # event which stores raw data
      events*: seq[Event]
      # time the shutter was open in seconds, one value for each
      # event
      length*: seq[float]
      # tots = ToT per pixel of whole run
      tots*: seq[seq[uint16]]
      # hits = num hits per event of whole run
      hits*: seq[seq[uint16]]
      # occupancies = occupancies of each chip for run
      occupancies*: Tensor[int64]
      # number of pixels removed due to ToT cut
      totCut*: ToTCut

    ##############
    # FADC types #
    ##############
    FadcSettings* = object of RootObj
      isValid*: bool
      postTrig*: int
      preTrig*: int
      bitMode14*: bool
      nChannels*: int
      channelMask*: int
      frequency*: int
      samplingMode*: int
      pedestalRun*: bool

    # object to save FADC data from file into
    # inherits from FadcObject, only adds a sequence
    # to store the data
    FadcFile* = object of FadcSettings
      data*: seq[uint16]
      eventNumber*: int
      trigRec*: int

    # object to store actual FADC data, which is
    # used (ch0 already extracted)
    # instead of a sequence for the data, we store the
    # converted data in an arraymancer tensor
    FadcData* = object of FadcSettings
      # will be a 2560 element tensor
      data*: Tensor[float]
      trigRec*: int

    ## Stores all the 1D information of the FADC data in a single
    ## run, i.e. everything that is computed related to the fall
    ## and rise time.
    ##
    ## NOTE: First of all the order of the fields is important. And it is
    ## important that it _only_ contains 1D data. We use `fieldPairs` to
    ## simplify dealing with the object (serializing to H5).
    RecoFadc* = object
      baseline*: Tensor[float] # all tensors 1D, 1 element per FADC event
      argMinval*: Tensor[uint16]
      riseStart*: Tensor[uint16]
      fallStop*: Tensor[uint16]
      riseTime*: Tensor[uint16]
      fallTime*: Tensor[uint16]
      skewness*: Tensor[float] # the skewness of the full FADC event, useful to detect noise events
      noisy*: Tensor[int]      # flag which says whether event was noisy
      minVal*: Tensor[float]  # minimum values of events (voltage of dips)

    ## This object stores the FADC data of a (possibly partial) run. That means
    ## each field contains N (= number of FADC events in a run) or a subset of
    ## that elements. The raw FADC data tensor is a `[N, 2560 * 4]` sized tensor,
    ## where each "row" is one FADC event of raw data.
    ProcessedFadcRun* = object
      settings*: FadcSettings
      # raw fadc data (optional)
      rawFadcData*: Tensor[uint16]
      # trigger record times, stored
      trigRecs*: seq[int]
      # register of minimum value
      minRegs*: seq[int]
      # eventNumber for FADC
      eventNumber*: seq[int]

    ReconstructedFadcRun* = object
      eventNumber*: seq[int]
      # processed and converted FADC data
      fadcData*: Tensor[float]


    CutMethodKind* = enum
      cmLnLCut, cmNnCut

    NeuralNetCutKind* = enum
      nkNone = "" # none selected
      nkGlobal = "global"
      nkLocal = "local"
      nkInterpolated = "interpolated"

    CutValueInterpolator* = object
      case cutMethod*: CutMethodKind
      of cmLnLCut:
        case morphKind*: MorphingKind
        of mkNone:
          ## just the regular table containing a cut value for each target/filter combination
          lnLCutTab*: OrderedTable[string, float]
        of mkLinear:
          ## two tensors: `cutEnergies` stores the energy values we used to
          ## compute different morphed distributions. `cutValues` stores the cut
          ## value of each distribution. So for index `idx` the cut value of energy
          ## `cutEnergies[idx]` is given as `cutValue[idx]`.
          lnLCutEnergies*: seq[float] # is a `seq` to use `lowerBound`
          lnLCutValues*: Tensor[float]
      of cmNnCut: # neural network cut method (only supported on cpp backend!)
        case nnCutKind*: NeuralNetCutKind
        of nkGlobal: cut*: float # global fixed efficiency
        of nkLocal: nnCutTab*: OrderedTable[string, float] # local (per target!) fixed efficiency
        of nkInterpolated: discard # liner interpolation between `nkLocal` targets
        else: discard

    ## Explanation of the 3 different line veto kinds. Terminology used here:
    ## 'OC': original cluster. This is the cluster that was identified to pass the
    ##    lnL cut on the center chip.
    ## 'HLC': Hypothetical Larger Cluster, a cluster found using the full Septemboard
    ##    event, which ``includes`` the OC in it.
    ## 1. 'regular' line veto. *Every* cluster checks the line to the center
    ##    cluster. Without septem veto this includes HLC checking OC.
    ## 2. 'regular without HLC' line veto: Lines check the OC, but the HLC is
    ##    explicitly *not* considered.
    ## 3. 'checking the HLC' line veto: In this case *all* clusters check the
    ##    center of the HLC.
    ## (From experience I would argue 3 in particular is bad and 1 is essentially a
    ## bad version of the septem veto + the line veto. Best to just use 2)
    LineVetoKind* = enum
      lvNone, # the default, if none given (only for initialization)
      lvRegular,
      lvRegularNoHLC,
      lvCheckHLC

    ## These veto related types are here as one of them uses a `Tensor` and they are
    ## anyhow only intended for usage in `likelihood.nim`.

    ################################################
    ############## FADC veto related ###############
    ################################################

    ## Helper object to store the cuts used in the FADC veto
    FadcCuts* = tuple
      active: bool # only a valid cut if `active`, used to indicate we had data for this FADC setting
      riseLow: float
      riseHigh: float
      fallLow: float
      fallHigh: float
      skewness: float

    ## Helper enum to map different FADC settings to a common string.
    ## The run numbers matching each setting are defined by the `toFadcSettings`
    ## helper in `fadc_utils.nim`.
    FadcSetting* = enum
      fs1 = "Setting 1"
      fs2 = "Setting 2"
      fs3 = "Setting 3"
      fs4 = "Setting 4"

    ################################################
    ############# Septem veto related ##############
    ################################################

    ## Stores all the data of all chips needed to reconstruct the septemboard events in the context
    ## of the septem veto
    AllChipData* = object
      x*: seq[seq[seq[uint8]]]
      y*: seq[seq[seq[uint8]]]
      ToT*: seq[seq[seq[uint16]]]
      charge*: seq[seq[seq[float]]]

    ## Stores information about a reconstructed septem event
    SeptemFrame* = object
      pixels*: PixelsInt        ## pure pixel data of all pixels in the septem frame (zero suppressed)
      centerCluster*: PixelsInt ## pure pixel data of the original center cluster that passes lnL
      charge*: Tensor[float]    ## charge values of all pixels in a [768,768] tensor
      centerEvIdx*: int         ## index of the center event
      numRecoPixels*: int       ## number of pixels ``after`` reconstruction. To truncate `pixels` after recoEvent

    ## Stores information about the clusters that pass the lnL cut on the center chip
    ## for which the septem veto is applied
    CenterChipData* = object
      lhoodCenter*: seq[float]
      energies*: seq[float]
      energyCenter*: seq[float]
      cXCenter*: seq[float]
      cYCenter*: seq[float]
      hitsCenter*: seq[int]
      rmsTCenter*: seq[float]
      rmsLCenter*: seq[float]

    ## Equivalent of the above, but for a single event.
    CenterClusterData* = object
      lhood*: float
      energy*: float
      energyCenter*: float
      cX*: float
      cY*: float
      hits*: int
      rmsT*: float
      rmsL*: float

    ## The name might not be the most descriptive: This object stores information about different
    ## geometric properties of the separate clusters found in the septem event. This includes
    ## information about the original "center" event (the one that passed the logL cut on the center
    ## chip) (`center` in the field name) and sequences for the centers of each other cluster
    ## as well as lines that go through the centers of each cluster along the long axis.
    SeptemEventGeometry* = object
      lines*: seq[tuple[m, b: float]]   ## lines along the long axis through the center of each cluster
      centers*: seq[tuple[x, y: float]] ## centers of each found cluster in the septem event
      xCenter*: int ## x center pixel of the original "center event" (that passed logL)
      yCenter*: int ## y center pixel of the original "center event" (that passed logL)
      centerRadius*: float ## "radius" of the original "center event" (based on 3 * (RMS_T + RMS_L)/2 in pixel)

    LogLFlagKind* = enum
      # vetoes
      fkTracking, fkLogL, fkMLP, fkConvNet, fkFadc, fkScinti, fkSeptem, fkLineVeto, fkAggressive,
      # other options
      fkRocCurve, fkComputeLogL, fkPlotLogL, fkPlotSeptem,
      fkEstRandomCoinc, # used to estimate the random coincidence of the septem & line veto
      fkEstRandomFixedEvent, # use a fixed center cluster and only vary around the outer ring
      fkReadOnly # makes the input file read only

    ## All the veto related setting parameters for easy serialization
    VetoSettings* = object
      # NN cut method
      useNeuralNetworkCut*: bool
      nnSignalEff*: float = 0.95
      nnModelPath*: string
      nnCutKind*: NeuralNetCutKind = nkLocal
      # lnL settings
      useLnLCut*: bool
      signalEfficiency*: float = 0.8 # the signal efficiency that defines the cuts
      # Septem & line veto related
      clusterAlgo*: ClusteringAlgorithm = caDBSCAN
      searchRadius*: int = 50 # for caDefault the search radius in septem events
      dbscanEpsilon*: float = 65.0 # for caDBSCAN the epsilon
      centerChip*: int = 3 # center chip on the detector
      numChips*: int = 7 # number of chips on the detector
      useRealLayout*: bool # whether to use Septemboard layout with gaps between chips or not
      # Line veto related
      lineVetoKind*: LineVetoKind  # which kind of line veto to use (see enum def for explanation)
      eccLineVetoCut*: float = 1.0 # eccentricity value that needs to be exceeded by a cluster to participate
                                   # in line veto
      # FADC veto related
      calibFile*: string # path to the calibration file from which to read FADC
                         # data to deduce rise/fall time cuts
      vetoPercentile*: float = 0.99 # percentile to use for rise/fall time distribution
      fadcScaleCutoff*: float = 1.45 # scale factor applied to the `maximum` of the rise/fall time distribution
                                     # to define the ``upper hard cutoff` on the data (due to double X-ray events
                                     # defining a _very long_ tail.
      fadcVetoes*: array[FadcSetting, FadcCuts] # computed `FadcCuts` based on `calibFile`. Only those entries
                                                # contained in `calibFile` will be filled!

    ## Helper object to store configuration parameters and relevant data pieces
    ## that are used recurringly in `likelihood.nim`.
    LikelihoodContext* = object
      cdlFile*: string ## the CDL file used
      year*: YearKind
      region*: ChipRegion
      morph*: MorphingKind
      energyDset*: InGridDsetKind
      timepix*: TimepixVersion
      stretch*: Option[CdlStretch]
      refSetTuple*: tuple[ecc, ldivRms, fracRms: Table[string, HistTuple]]
      numMorphedEnergies*: int = 1000
      refDf*: DataFrame
      refDfEnergy*: seq[float]
      # general
      useTeX*: bool # whether to generate TikZ plots or cairo
      #when defined(cpp):
      flags*: set[LogLFlagKind]
      vetoCfg*: VetoSettings

proc initFeSpecData*(hist: seq[float],
                     binning: seq[float],
                     idx_kalpha: int,
                     idx_sigma: int,
                     pRes: seq[float],
                     pErr: seq[float],
                     xFit: seq[float],
                     yFit: seq[float],
                     chiSq: float,
                     nDof: int): FeSpecFitData =
  result = FeSpecFitData(hist: hist,
                         binning: binning,
                         idx_kalpha: idx_kalpha,
                         kalpha: pRes[idx_kalpha],
                         idx_sigma: idx_sigma,
                         sigma_kalpha: pRes[idx_sigma],
                         pRes: pRes,
                         pErr: pErr,
                         xFit: xFit,
                         yFit: yFit,
                         chiSq: chiSq,
                         nDof: nDof)

proc initEnergyCalibData*(energies: seq[float],
                          peaks: seq[float],
                          peaksErr: seq[float],
                          pRes: seq[float],
                          pErr: seq[float],
                          xFit: seq[float],
                          yFit: seq[float],
                          chiSq: float,
                          nDof: int): EnergyCalibFitData =
  let aInv = 1.0 / pRes[0] * 1000
  result = EnergyCalibFitData(energies: energies,
                              peaks: peaks,
                              peaksErr: peaksErr,
                              pRes: pRes,
                              pErr: pErr,
                              xFit: xFit,
                              yFit: yFit,
                              aInv: aInv,
                              aInvErr: aInv * pErr[0] / pRes[0],
                              chiSq: chiSq,
                              nDof: nDof)

proc initInterval*(idx: int, interval, minInterval: float,
                   tStart = 0.0, tStop = 0.0): GasGainIntervalData =
  result = GasGainIntervalData(idx: idx,
                               interval: interval,
                               minInterval: minInterval,
                               tStart: tStart.int,
                               tStop: tStop.int)

proc `$`*(g: GasGainIntervalData): string =
  let tStamp = $(g.tStart.fromUnix.utc)
  result = &"idx = {g.idx}, interval = {g.interval} min, start = {tStamp}"

proc toAttrPrefix*(g: GasGainIntervalData): string =
  result = "interval_" & $g.idx & "_"

proc toGainAttr*(g: GasGainIntervalData): string =
  result = g.toAttrPrefix & "G"

proc toSliceStartAttr*(g: GasGainIntervalData): string =
  result = g.toAttrPrefix & "tstart"

proc toSliceStopAttr*(g: GasGainIntervalData): string =
  result = g.toAttrPrefix & "tstop"

proc toPathSuffix*(g: GasGainIntervalData): string =
  result = &"_{g.idx}_{g.interval.int}_min_{g.tstart}"

proc toDsetSuffix*(g: GasGainIntervalData): string =
  result = &"_{g.idx}_{g.interval.int}"

proc initCutValueInterpolator*(kind: MorphingKind): CutValueInterpolator =
  case kind
  of mkNone:
    result = CutValueInterpolator(cutMethod: cmLnLCut, morphKind: mkNone)
    result.lnLCutTab = initOrderedTable[string, float]()
  of mkLinear:
    result = CutValueInterpolator(cutMethod: cmLnLCut, morphKind: mkLinear)

proc initCutValueInterpolator*(kind: NeuralNetCutKind): CutValueInterpolator =
  case kind
  of nkNone: doAssert false, "Should not happen"
  of nkGlobal:
    result = CutValueInterpolator(cutMethod: cmNnCut, nnCutKind: nkGlobal)
  of nkLocal:
    result = CutValueInterpolator(cutMethod: cmNnCut, nnCutKind: nkLocal)
    result.nnCutTab = initOrderedTable[string, float]()
  of nkInterpolated:
    result = CutValueInterpolator(cutMethod: cmNnCut, nnCutKind: nkInterpolated)

proc initCutValueInterpolator*(tab: OrderedTable[string, float],
                               cutMethod = cmLnLCut): CutValueInterpolator =
  case cutMethod
  of cmLnLCut:
    result = initCutValueInterpolator(mkNone)
    result.lnLCutTab = tab
  of cmNnCut:
    result = initCutValueInterpolator(nkLocal)
    result.nnCutTab = tab

proc initCutValueInterpolator*(energies: seq[float]): CutValueInterpolator =
  result = initCutValueInterpolator(mkLinear)
  result.lnLCutEnergies = energies
  result.lnLCutValues = zeros[float](energies.len.int)

proc `[]`*(cv: CutValueInterpolator, s: string): float =
  case cv.cutMethod
  of cmLnLCut: result = cv.lnLCutTab[s]
  of cmNnCut: result = cv.nnCutTab[s]

proc `[]=`*(cv: var CutValueInterpolator, s: string, val: float) =
  doAssert (cv.cutMethod == cmLnLCut and cv.morphKind == mkNone) or
    (cv.cutMethod == cmNnCut and cv.nnCutKind == nkLocal)
  case cv.cutMethod
  of cmLnLCut: cv.lnLCutTab[s] = val
  of cmNnCut: cv.nnCutTab[s] = val

proc `==`*[T: SomePix](c1, c2: ClusterObject[T]): bool =
  result = true
  result = result and c1.data     == c2.data
  result = result and c1.hits     == c2.hits
  result = result and c1.centerX  == c2.centerX
  result = result and c1.centerY  == c2.centerY
  result = result and c1.sumTot   == c2.sumTot
  result = result and c1.energy   == c2.energy
  result = result and c1.energy   == c2.energy
  result = result and c1.geometry == c2.geometry
  result = result and c1.version  == c2.version
  if result and c1.version == Timepix3:
    result = result and c1.toa    == c2.toa
    result = result and c1.toaCombined == c2.toaCombined

proc `==`*[T: SomePix](r1, r2: RecoEvent[T]): bool =
  result = true
  for name, f1, f2 in fieldPairs(r1, r2):
    result = result and f1 == f2

when false:
  ## XXX:  write some macro code to clean up the whole int / float dataset debacle
  macro objFields*(t: typed): untyped =
    ## Given the type `t`, generate code returning names of each field
    ## that is of a primitive type (int + float)
    let typ = t.getType[1].getImpl
    doAssert typ.kind == nnkTypeDef
    let recLs = typ[2][2] # objTy, recList

    echo typ.treerepr
