# module which contains the used type definitions in the InGrid module
import times
import tables
import arraymancer

type
  # an object, which stores information about a run's start, end and length
  RunTimeInfo* = object
    t_start*: Time
    t_end*: Time
    t_length*: Duration

  EventHeader* = Table[string, string]
  ChipHeader*  = Table[string, string]
  Pix*         = tuple[x, y: uint8, ch: uint16]
  Pixels*      = seq[Pix]

  # Coord type which contains (x, y) coordinates of a pixel
  Coord* = tuple[x, y: uint8]
  # cluster object
  Cluster* = seq[Pix]

  Pixels_prot = object#Table[string, seq[int]]
    # x:  seq[int]
    # y:  seq[int]
    # ch: seq[int]
    x:  seq[uint8]
    y:  seq[uint8]
    ch: seq[uint16]

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

  # define a distinct `OldEvent` to differentiate the (in principle) not
  # different old TOS storage format
  OldEvent* = Event

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

  Threshold* = Tensor[int]
  ThresholdMeans* = Tensor[int]

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
    rfNewTos, rfOldTos

  ################################
  # Reconstruction related types #
  ################################

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
  ClusterObject* = object
    data*: Cluster
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
  RecoEvent* = object
    cluster*: seq[ClusterObject]
    event_number*: int
    chip_number*: int

  # process events stores all data for septemboard
  # of a given run
  ProcessedRun* = tuple[
    # just the number of chips in the run
    nChips: int,
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

  ################################
  #### Analysis related types ####
  ################################

  ChipRegion* = enum
    crGold, crSilver, crBronze, crAll

  Cuts* = object
    minCharge*: float
    maxCharge*: float
    minRms*: float
    maxRms*: float
    maxLength*: float
    minPix*: float
  # object to store region based cuts (gold, silver, bronze region)
  CutsRegion* = object
    xMin*: float
    xMax*: float
    yMin*: float
    yMax*: float
    radius*: float

# and some general InGrid related constants
const NPIX* = 256
const PITCH* = 0.055
