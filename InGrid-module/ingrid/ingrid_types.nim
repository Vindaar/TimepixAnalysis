# module which contains the used type definitions in the InGrid module
import times
import tables
import arraymancer

type 
  # an object, which stores information about a run's start, end and length
  RunTimeInfo* = object
    t_start*: Time
    t_end*: Time
    t_length*: TimeInterval

  EventHeader* = Table[string, string]
  ChipHeader*  = Table[string, string]
  Pix*         = tuple[x, y, ch: int]
  Pixels*      = seq[Pix]

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

  # process events stores all data for septemboard
  # of a given run
  ProcessedRun* = tuple[
    # run number
    run_number: int,
    # table containing run header ([General] in data file)
    runHeader: Table[string, string],
    # event which stores raw data    
    events: seq[Event],
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
    
  # object to save FADC data from file into
  FadcFile* = object
    vals*: seq[float]
    posttrig*: int 
    trigrec*: int
    bit_mode14*: bool

  # object to store actual FADC data, which is
  # used (ch0 already extracted)
  FadcData* = object
    data*: seq[float]
    posttrig*: int
    trigrec*: int
