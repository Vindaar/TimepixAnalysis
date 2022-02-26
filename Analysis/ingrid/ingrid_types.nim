# module which contains the used type definitions in the InGrid module
when not defined(js):
  import times
import tables, strformat, memfiles
import karax / kbase

type
  EventHeader* = Table[string, string]
  ChipHeader*  = Table[string, string]
  Pix*         = tuple[x, y: uint8, ch: uint16]
  PixTpx3*     = tuple[x, y: uint8, ch, toa: uint16, toaCombined: uint64]
  # Integer based pixels are used for full Septemboard frames, due to 3x256 pixels per direction
  PixInt*      = tuple[x, y: int, ch: int]
  SomePix*     = Pix | PixInt | PixTpx3
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

  histTuple* = tuple[bins: seq[float64], hist: seq[float64]]

  #############################
  # Calibration related types #
  #############################

  Tot* = object
    pulses*: seq[int]
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

when not defined(js):
  type
    # an object, which stores information about a run's start, end and length
    RunTimeInfo* = object
      t_start*: Time
      t_end*: Time
      t_length*: Duration

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
      activeTime*: Duration # total time shutter was open
      # reuse `RunTimeInfo` to store possible tracking starts / ends
      trackings*: seq[RunTimeInfo]
      nonTrackingDuration*: Duration
      trackingDuration*: Duration

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
    toaLength*: float # length in ToA ## XXX: ideally this would be a small integer type, i.e. int32
    toaMean*: float   # mean ToA value
    toaRms* : float   # RMS of ToA values
    toaMin* : uint16  # the minimal ToA value found before subtraction

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

  ##############
  # FADC types #
  ##############

  FadcObject* = object of RootObj
    isValid*: bool
    postTrig*: int
    preTrig*: int
    trigRec*: int
    bitMode14*: bool
    nChannels*: int
    channelMask*: int
    frequency*: int
    samplingMode*: int
    pedestalRun*: bool


  # object to save FADC data from file into
  # inherits from FadcObject, only adds a sequence
  # to store the data
  FadcFile* = object of FadcObject
    data*: seq[uint16]
    eventNumber*: int

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

  # a simple object storing the runs, chips etc. from a given
  # H5 file
  FileInfo* = object
    runs*: seq[int]
    chips*: seq[int]
    runType*: RunTypeKind
    rfKind*: RunFolderKind
    centerChip*: int
    centerChipName*: kstring
    hasFadc*: bool # reads if FADC group available
    timepix*: TimepixVersion

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


const TosDateString* = "yyyy-MM-dd'.'hh:mm:ss"

# and some general InGrid related constants
const NPIX* = 256
const PITCH* = 0.055

const SrsRunIncomplete* = "incomplete"
const SrsRunIncompleteMsg* = "This run does not contain a run.txt and so " &
  "is incomplete!"
const SrsNoChipId* = "ChipIDMissing"
const SrsNoChipIdMsg* = "The chip IDs are missing from the run.txt. Old format!"
const SrsDefaultChipName* = "SRS Chip"

# the following will not be available, if the `-d:pure` flag is set,
# to allow importing the rest of the types, without a `arraymancer`
# dependency
when not defined(pure) and not defined(js):
  import arraymancer
  type
    Threshold* = Tensor[int]
    ThresholdMeans* = Tensor[int]

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

    # object to store actual FADC data, which is
    # used (ch0 already extracted)
    # instead of a sequence for the data, we store the
    # converted data in an arraymancer tensor
    FadcData* = object of FadcObject
      # will be a 2560 element tensor
      data*: Tensor[float]

    ProcessedFadcData* = tuple[
      # raw fadc data
      rawFadcData: seq[seq[uint16]],
      # processed and converted FADC data
      fadcData: Tensor[float],
      # trigger record times, stored
      trigRecs: seq[int],
      # flag which says whether event was noisy
      noisy: seq[int],
      # minimum values of events (voltage of dips)
      minVals: seq[float],
      # register of minimum value
      minRegs: seq[int],
      #eventNumber for FADC
      eventNumber: seq[int],
    ]

    CutValueInterpolator* = object
      case kind*: MorphingKind
      of mkNone:
        ## just the regular table containing a cut value for each target/filter combination
        cutTab*: Table[string, float]
      of mkLinear:
        ## two tensors: `cutEnergies` stores the energy values we used to
        ## compute different morphed distributions. `cutValues` stores the cut
        ## value of each distribution. So for index `idx` the cut value of energy
        ## `cutEnergies[idx]` is given as `cutValue[idx]`.
        cutEnergies*: seq[float] # is a `seq` to use `lowerBound`
        cutValues*: Tensor[float]

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
    result = CutValueInterpolator(kind: mkNone)
    result.cutTab = initTable[string, float]()
  of mkLinear:
    result = CutValueInterpolator(kind: mkLinear)

proc initCutValueInterpolator*(tab: Table[string, float]): CutValueInterpolator =
  result = initCutValueInterpolator(mkNone)
  result.cutTab = tab

proc initCutValueInterpolator*(energies: seq[float]): CutValueInterpolator =
  result = initCutValueInterpolator(mkLinear)
  result.cutEnergies = energies
  result.cutValues = zeros[float](energies.len.int)

proc `[]`*(cv: CutValueInterpolator, s: string): float =
  result = cv.cutTab[s]

proc `[]=`*(cv: var CutValueInterpolator, s: string, val: float) =
  doAssert cv.kind == mkNone
  cv.cutTab[s] = val

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
