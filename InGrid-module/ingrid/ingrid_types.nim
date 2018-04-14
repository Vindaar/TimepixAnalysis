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
  Pix*         = tuple[x, y, ch: int]
  Pixels*      = seq[Pix]

  # Coord type which contains (x, y) coordinates of a pixel
  Coord* = tuple[x, y: int]
  # cluster object
  Cluster* = seq[Pix]  

  Pixels_prot = object#Table[string, seq[int]]
    x:  seq[int]
    y:  seq[int]
    ch: seq[int]

  Chip* = tuple[name: string, number: int]

  ChipEvent* = object
    chip*: Chip
    pixels*: Pixels
    
  Event* = object
    evHeader*: Table[string, string]
    chips*: seq[ChipEvent]
    # time the shutter was open in seconds
    length*: float

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
    # run number
    run_number: int,
    # table containing run header ([General] in data file)
    runHeader: Table[string, string],
    # event which stores raw data    
    events: seq[Event],
    # time the shutter was open in seconds, one value for each
    # event
    length: seq[float],
    # tots = ToT per pixel of whole run
    tots: seq[seq[int]],
    # hits = num hits per event of whole run
    hits: seq[seq[int]],
    # occupancies = occupancies of each chip for run
    occupancies: Tensor[int]
    #occupancies: seq[Tensor[int]]
  ]

  EventSortType* = enum
    fname, inode

  EventType* = enum
    FadcType, InGridType

  ##############
  # FADC types #
  ##############

  FadcObject* = object of RootObj
    posttrig*: int
    pretrig*: int
    trigrec*: int
    bit_mode14*: bool
    n_channels*: int
    channel_mask*: int
    frequency*: int
    sampling_mode*: int
    pedestal_run*: bool
    
    
  # object to save FADC data from file into
  # inherits from FadcObject, only adds a sequence
  # to store the data
  FadcFile* = object of FadcObject
    data*: seq[uint16]

  # object to store actual FADC data, which is
  # used (ch0 already extracted)
  # instead of a sequence for the data, we store the
  # converted data in an arraymancer tensor
  FadcData* = object of FadcObject
    # will be a 2560 element tensor
    data*: Tensor[float]
    
  ProcessedFadcData* = tuple[
    # raw fadc data
    raw_fadc_data: seq[seq[uint16]],
    # processed and converted FADC data
    fadc_data: Tensor[float],
    # trigger record times, stored 
    trigrecs: seq[int],
    # flag which says whether event was noisy
    noisy: seq[int],
    # minimum values of events (voltage of dips)
    minvals: seq[float]
    # more?
  ]
    

# and some general InGrid related constants
const NPIX* = 256
const PITCH* = 0.055

