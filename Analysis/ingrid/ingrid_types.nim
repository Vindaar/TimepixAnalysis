# module which contains the used type definitions in the InGrid module
when not defined(js):
  import times
import tables, strformat
import karax / kbase

type
  EventHeader* = Table[string, string]
  ChipHeader*  = Table[string, string]
  Pix*         = tuple[x, y: uint8, ch: uint16]
  # Integer based pixels are used for full Septemboard frames, due to 3x256 pixels per direction
  PixInt*      = tuple[x, y: int, ch: int]
  SomePix*     = Pix | PixInt
  Pixels*      = seq[Pix]
  PixelsInt*   = seq[PixInt]

  # Coord type which contains (x, y) coordinates of a pixel
  Coord*[T] = tuple[x, y: T]
  # cluster object
  Cluster*[T: SomePix] = seq[T]

  Chip* = tuple[name: string, number: int]

  ChipEvent* = object
    chip*: Chip
    pixels*: Pixels

  Event* = object
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
    # TODO: move the following to a CONFIG object
    plotlySaveSvg*: bool
    # NOTE: add other flags for other optional plots?
    # if e.g. FeSpec not available yet, we can just call the
    # procedure to create it for us

type
  ## GasGainIntervalData stores the information about binning the data for the
  ## calculation of the gas gain in a run by a fixed time interval, e.g.
  ## 100 minutes
  GasGainIntervalData* = object
    idx*: int # index
    interval*: float # interval length in minutes
    tStart*: int # timestamp of interval start
    tStop*: int # timestamp of interval stop

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
    ProcessedRun* = tuple[
      # just the number of chips in the run
      nChips: int,
      # the chips as (name, number) tuples
      chips: seq[Chip],
      # run number
      runNumber: int,
      # table containing run header ([General] in data file)
      runHeader: Table[string, string],
      # event which stores raw data
      events: seq[Event],
      # time the shutter was open in seconds, one value for each
      # event
      length: seq[float],
      # tots = ToT per pixel of whole run
      tots: seq[seq[uint16]],
      # hits = num hits per event of whole run
      hits: seq[seq[uint16]],
      # occupancies = occupancies of each chip for run
      occupancies: Tensor[int64]
      #occupancies: seq[Tensor[int]]
    ]

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

proc initInterval*(idx: int, interval: float,
                   tStart = 0.0, tStop = 0.0): GasGainIntervalData =
  result = GasGainIntervalData(idx: idx,
                               interval: interval,
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
