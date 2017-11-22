import os, ospaths
import strutils
import times
import algorithm
import re
import tables
import memfiles
import sequtils, future

# other modules
import arraymancer

type 
  # an object, which stores information about a run's start, end and length
  RunTimeInfo* = object
    t_start*: Time
    t_end*: Time
    t_length*: TimeInterval

  EventHeader* = Table[string, string]
  ChipHeader*  = Table[string, string]
  Pixels*      = seq[tuple[x, y, ch: int]]

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
    
    

proc parseTOSDateString*(date_str: string): Time = 
  # function receives a string from a date time from TOS and creates
  # a Time object from it
  let date = toTime(parse(date_str, "yyyy-MM-dd.hh:mm:ss"))
  return date

proc readDateFromEvent*(filepath: string): string = 
  # procedure opens the given file and returns the dateTime value
  # (TOS syntaxed date string)
  for line in lines filepath:
    if "dateTime" in line:
      # if dateTime is found, split, assign to string and break from while
      let line_seq = split(line, " ")
      result = line_seq[high(line_seq)]
      break
    else:
      continue

proc formatAsOrgDate*(t: Time, org_format = "yyyy-MM-dd ddd H:mm"): string =
  # this procedure formats the given Time object as an org-mode date
  # by first converting it to a TimeInfo object
  let ti = getLocalTime(t)
  result = format(ti, org_format)

proc getTimeFromEvent*(file: string): Time = 
  # this procedure returns the time info from an event,
  # given the filename by parsing the TOS date syntax
  let date_str = readDateFromEvent(file)
  result = parseTOSDateString(date_str)

proc readListOfFilesAndGetTimes*(path: string, list_of_files: seq[string]): seq[Time] = 
  # function analogues to walkRunFolderAndGetTimes except we read all files
  # given by list_of_files
  var i: int = 0
  result = @[]

  for file in list_of_files:
    let filename = joinPath(path, file)
    let date_str = readDateFromEvent(filename)
    if date_str isnot "":
      result.add(parseTOSDateString(date_str))
    if i mod 500 == 0:
      echo i, " files read."
    i += 1    
  return result

proc walkRunFolderAndGetTimes*(folder: string): seq[Time] = 
  var i: int = 0
  result = @[]

  # for file in walkDir(folder):
  #   let path = file.path
  for file in walkFiles(joinPath(folder, "data*.txt-fadc")):
    let path = file
    if "data" in path and "fadc" in path:
      let filename = strip(path, false, true, {'-', 'f', 'a', 'd', 'c'})
      let date_str = readDateFromEvent(filename)
      if date_str isnot "":
        result.add(parseTOSDateString(date_str))
      if i mod 500 == 0:
        echo i, " files read."
      i += 1
  return result

proc writeDateSeqToFile*(date_seq: seq[Time]): void = 
  let outfile = "times.txt"
  var f: File
  if open(f, outfile, fmWrite):
    f.write("# date\t unix time\n")
    for time in date_seq:
      let sec: string = formatFloat(toSeconds(time))
      #let t: string = format(getLocalTime(time), "yyyy-MM-dd hh:mm:ss")
      #let str = join(@[t, sec, "\n"], sep = "\t")
      let str = join(@[$(time), sec, "\n"], sep = "\t")
      f.write(str)
  f.close()


proc getRunTimeInfo*(run_files: seq[string]): RunTimeInfo =
  # this procdure creates a RunTimeInfo object from a given list of files
  # in a run folder (not sorted). The list is sorted and then the first and
  # last event are read, converted to Time objects and the length of the run 
  # is calculated from the difference between both events
  result = RunTimeInfo()

  # sort list of files
  let sorted_files = sorted(run_files, system.cmp[string])
  let first_file = sorted_files[0]
  let last_file  = sorted_files[^1]

  let time_first = getTimeFromEvent(first_file)
  let time_last  = getTimeFromEvent(last_file)

  let run_length = initInterval(seconds=int(time_last - time_first))
  
  result.t_start = time_first
  result.t_end = time_last 
  result.t_length = run_length


proc getRegexForEvents*(): tuple[header, chips, pixels: string] = 
  # this procedure returns a tuple of the regexes used to
  # read different parts of the events
  # outputs:
  #     tuple[header, chips, pixels: string] =
  #       header: the regex to read the event header
  #       chips: the regex to read the chip headers
  #       pixels: the regex to read the pixel data
  result.header = r"^\#{2}\s(\w+):\s+\b(\S*)\b"
  result.chips  = r"^\#\s(\w+):\s+\b(\w\ \d{1,2} W\d{2}|\S*)\b"
  result.pixels = r"^(\d+)\s+(\d+)\s+(\d+)$"

proc readEventHeader*(filepath: string): Table[string, string] =
  # this procedure reads a whole event header and returns 
  # a table containing the data, where the key is the key from the data file
  result = initTable[string, string]()
  
  # we define a regex for the header of the file
  let regex = r"^\#{2}\s(\w+):\s+\b(\S*)\b"
  var matches: array[2, string]
  for line in lines filepath:
    if line.match(re(regex), matches):
      # get rid of whitespace and add to result
      let key = matches[0]
      let val = matches[1]
      result[key] = val

# template createTuple(ln: int) =
  
##   for i in 0..<ln:


proc readMemFilesIntoBuffer*(list_of_files: seq[string]): seq[seq[string]] =
  # procedure which reads a list of files via memory mapping and returns
  # the thus read data as a sequence of MemFiles
  # inputs:
  #    list_of_files: seq[string] = seq of strings containing the filenames to be read
  # outputs:
  #    seq[var MemFile] = seq of memmapped files with data to work on
  result = @[]
  
  for f in list_of_files:
    var ff: MemFile = memfiles.open(f, mode = fmRead, mappedSize = -1)
    var dat: seq[string] = @[]
    for l in lines(f):
      dat.add(l)
    ff.close()
    result.add(dat)

proc processEventWithRegex*(data: seq[string], regex: tuple[header, chips, pixels: string]): ref Event =
  # this template is used to create the needed functions to
  # - read the event header
  # - read the chip information for an event
  # - read the pixel data in an event
  # either as single functions or as a combination of
  # all
  # when type(regex) == string:
  #   # in this case we deal with a single regex string
  #   result = newTable[string, string]()
  #   for line in lines filepath:
  #     if line.match(re(regex), matches):
  #       # get rid of whitespace and add to result
  #       let key = matches[0]
  #       let val = matches[1]
  #       result[key] = val
  # regex[0] == header
  # regex[1] == chips
  # regex[2] == pixels
  # in this case we have a sequence of strings
  var
    # variable to count already read pixels
    pix_counter = 0
    # variables for regex
    h_matches: array[2, string]
    p_matches: array[3, string]
    # variables to read data into
    e_header = initTable[string, string]()
    c_header = initTable[string, string]()
    pixels: Pixels = @[]
    # variable to store resulting chip events
    chips: seq[ChipEvent] = @[]
  result = new Event

  for line in data:
    if match(line, re(regex[0]), h_matches) == true:
      e_header[h_matches[0]] = h_matches[1]
    elif match(line, re(regex[1]), h_matches) == true:
      # in case we match a chip header, pix_counter to 0,
      # need to make sure it is 0, once we reach the first
      # hit pixel
      pix_counter = 0      
      c_header[h_matches[0]] = h_matches[1]
    elif match(line, re(regex[2]), p_matches) == true:
      # in this case we have matched a pixel hit line
      # get number of hits to process
      let
        x: int  = parseInt(p_matches[0])
        y: int  = parseInt(p_matches[1])
        ch: int = parseInt(p_matches[2])
      pixels.add((x, y, ch))
      # after adding pixel, increase pixel counter so that we know when we're
      # reading the last hit of a given chip
      inc pix_counter
      if pix_counter == parseInt(c_header["numHits"]):
        # now we are reading the last hit, process chip header and pixels
        var ch_event = ChipEvent()
        ch_event.chip = (c_header["chipName"], parseInt(c_header["chipNumber"]))
        ch_event.pixels = pixels
        # add  the chip event object to the sequence
        chips.add(ch_event)
        # reset the pixels sequence
        pixels = @[]

  result.evHeader = e_header
  result.chips = chips

  

#template readEventWithRegex(filepath, regex: string): typed =
proc readEventWithRegex*(filepath: string, regex: tuple[header, chips, pixels: string]): ref Event =
  # this procedure reads the lines from a given event and then calls processEventWithRegex
  # to process the file. Returns a ref to an Event object
  # inputs:
  #    filepath: string = the path to the file to be read
  #    regex: tuple[...] = tuple of 3 regex's, one for each part of an event file
  # outputs:
  #    ref Event: reference to an event object
  var f = memfiles.open(filepath, mode = fmRead, mappedSize = -1)
  var data: seq[string] = @[]
  for line in lines(f):
    data.add(line)
  result = processEventWithRegex(data, regex)
  


proc pixelsToTOT*(pixels: Pixels): seq[int] {.inline.} =
  # extracts all charge values (ToT values) of a given pixels object (all hits
  # of a single chip in a given event, filters the full values of 11810
  # inputs:
  #    pixels: Pixels: input pixel sequence of tuples containing hits of one event
  # output:
  #    seq[int]: all ToT values smaller ToT
  result = filter(map(pixels, (p: tuple[x, y, ch: int]) -> int => p.ch),
                  (ch: int) -> bool => ch < 11810)
  
template addPixelsToOccupancy*(ar: Tensor[int], pixels: Pixels) =
  # template to add pixels to occupancy by using map
  #map(pixels, (p: tuple[x, y, c: int]) => ar[p.x, p.y] = p.c)
  for p in pixels:
    ar[p.x, p.y] += 1#p.ch

template addPixelsToOccupancySeptem*(ar: var Tensor[int], pixels: Pixels, ch_num: int) =
  # template to add pixels to occupancy by using map
  #map(pixels, (p: tuple[x, y, c: int]) => ar[p.x, p.y] = p.c)
  #proc(ar: var Tensor[int], pixels: Pixels, ch_num: int) =
  for p in pixels:
    ar[ch_num, p.x, p.y] += 1#p.ch

proc createTensorFromZeroSuppressed*(pixels: Pixels): Tensor[int] =
  # procedure to create a (256, 256) int array from a Pixels (seq[tuple[x, y, ch]])
  # object
  result = zeros[int](256, 256)
  for p in pixels:
    result[p.x, p.y] = p.ch
  

proc rawEventManipulation(filepath: string, regex: tuple[header, chips, pixels: string]): ref Event =
  # this procedure performs all processing of a single event, from a raw data event to
  # a processedEvent
  discard
#template 


proc dumpFrameToFile*(filepath: string, ar: Tensor[int]) =
  # this procedure dumps a given frame (tensor ar needs to be of shape (256, 256)
  # to a file 'filepath'
  # inputs:
  #   filepath: string = the file to write to
  #   ar: Tensor[int] = a tensor of shape (256, 256) containing the data to be written
  doAssert(ar.shape == [256, 256])
  var f: File
  if open(f, filepath, fmWrite):
    for x in 0..<256:
      for y in 0..<256:
        f.write($ar[x, y] & "\t")
      f.write("\n")
  f.close()
  
  
