import os, ospaths
import strutils
import times
import algorithm
import re
import tables
import memfiles
import sequtils, future
import threadpool
import math

# custom modules
import helper_functions
import ingrid_types
import nimhdf5

# other modules
import arraymancer

import macros

{.deadCodeElim: on.}

proc sum*(c: seq[Pix]): Pix {.inline.} =
  # this procedure sums the sequence of pixels such that it returns
  # a tuple of (sum_x, sum_y, sum_charge)
  # if T is itself e.g. a tuple, we will return a tuple, one
  # element for each field in the tuple
  assert c.len > 0, "Can't sum empty sequences"
  for p in c:
    result.x  += p.x
    result.y  += p.y
    result.ch += p.ch

proc sum2*(c: seq[Pix]): Pix {.inline.} =
  # this procedure sums the squares of the pixels in the sequence
  assert c.len > 0, "Can't sum empty sequences"
  for p in c:
    result.x += p.x * p.x
    result.y += p.y * p.y
    result.ch += p.ch * p.ch
    
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

proc initIntervalOld*(milliseconds, seconds, minutes, hours, days, months,
                   years: int = 0): TimeInterval =
  ## creates a new ``TimeInterval``.
  ##
  ## You can also use the convenience procedures called ``milliseconds``,
  ## ``seconds``, ``minutes``, ``hours``, ``days``, ``months``, and ``years``.
  ##
  ## Example:
  ##
  ## .. code-block:: nim
  ##
  ##     let day = initInterval(hours=24)
  ##     let tomorrow = getTime() + day
  ##     echo(tomorrow)
  var carryO = 0
  result.milliseconds = `mod`(milliseconds, 1000)
  carryO = `div`(milliseconds, 1000)
  result.seconds = `mod`(carryO + seconds, 60)
  carryO = `div`(carryO + seconds, 60)
  result.minutes = `mod`(carryO + minutes, 60)
  carryO = `div`(carryO + minutes, 60)
  result.hours = `mod`(carryO + hours, 24)
  carryO = `div`(carryO + hours, 24)
  result.days = carryO + days

  result.months = `mod`(months, 12)
  carryO = `div`(months, 12)
  result.years = carryO + years


proc getRunTimeInfo*(run_files: seq[string]): RunTimeInfo =
  # this procdure creates a RunTimeInfo object from a given list of files
  # in a run folder (not sorted). The list is sorted and then the first and
  # last event are read, converted to Time objects and the length of the run 
  # is calculated from the difference between both events
  # sort list of files
  result = RunTimeInfo()
  
  let
    sorted_files = sorted(run_files, system.cmp[string])
    first_file = sorted_files[0]
    last_file  = sorted_files[^1]
    # and times 
    time_first = getTimeFromEvent(first_file)
    time_last  = getTimeFromEvent(last_file)
    # calc run length
    run_length = time_last - time_first
  echo "Time first is $# and time last is $#" % [$time_first, $time_last]
  echo "Time difference in seconds $#" % $((time_last - time_first).seconds)
  echo "Time difference start end $#" % $(run_length)
  
  result.t_start = time_first
  result.t_end = time_last 
  result.t_length = run_length

proc parseShutterMode*(mode: string): int =
  # proc to parse the TOS shutter mode selection
  if mode == "verylong" or mode == "vl":
    result = 2
  elif mode == "long" or mode == "l":
    result = 1
  else:
    result = 0

proc calcLength*(event: Event): float =
  # calculates the event length in seconds, taking into account
  # whether the FADC triggered or not and the shutter opening time
  # TODO: is it really * 46 or not? Kind of important...
  # 2, 13 comes out correctly with 46!
  let
    time = parseFloat(event.evHeader["shutterTime"])
    mode = float(parseShutterMode(event.evHeader["shutterMode"]))
    fadc_triggered = parseInt(event.evHeader["fadcReadout"])
  if unlikely(fadc_triggered == 1):
    let fadc_clock_triggered = parseFloat(event.evHeader["fadcTriggerClock"])
    # in case FADC triggered it's simply the number of clockcycles the shutter was open
    # multiplied by the clock speed numbers
    result = fadc_clock_triggered * 46'f / 40'f / 1_000_000'f
  else:
    result = pow(256'f, mode) * 46'f * time / 40'f / 1_000_000'f

proc calcLength*(event: Event, time, mode: float): float =
  # overload of above, with time and mode already read
  # calculates the event length in seconds, taking into account
  # whether the FADC triggered or not and the shutter opening time
  let fadc_triggered = parseInt(event.evHeader["fadcReadout"])
  if unlikely(fadc_triggered == 1):
    let fadc_clock_triggered = parseFloat(event.evHeader["fadcTriggerClock"])
    # in case FADC triggered it's simply the number of clockcycles the shutter was open
    # divided by the clock speed
    result = fadc_clock_triggered / 40'f / 1_000_000'f
  else:
    result = pow(256'f, mode) * 46'f * time / 40'f / 1_000_000'f

proc getRegexForEvents*(): tuple[header, chips, pixels: string] = 
  # this procedure returns a tuple of the regexes used to
  # read different parts of the events
  # outputs:
  #     tuple[header, chips, pixels: string] =
  #       header: the regex to read the event header
  #       chips: the regex to read the chip headers
  #       pixels: the regex to read the pixel data
  result.header = r"^\#{2}\s(\w+):\s+(-?\b\S*)\b"
  result.chips  = r"^\#\s(\w+):\s+\b(\w\ \d{1,2} W\d{2}|\S*)\b"
  result.pixels = r"^(\d+)\s+(\d+)\s+(\d+)$"

proc readEventHeader*(filepath: string): Table[string, string] =
  # this procedure reads a whole event header and returns 
  # a table containing the data, where the key is the key from the data file
  result = initTable[string, string]()
  
  # we define a regex for the header of the file
  let regex = re(r"^\#{2}\s(\w+):\s+\b(\S*)\b")
  var matches: array[2, string]
  for line in lines filepath:
    if line.match(regex, matches):
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
  result = newSeqOfCap[seq[string]](len(list_of_files))

  var
    ff: MemFile
    # reserver enough space for most events, only for large events do we have to reserve
    # more space
    dat: seq[string] = @[] #newSeqOfCap[string](100)

  echo "free memory ", getFreeMem()
  echo "occ memory ", getOccupiedMem()    
  for f in list_of_files:
    ff = memfiles.open(f, mode = fmRead, mappedSize = -1)
    for l in lines(ff):
      dat.add(l)
    ff.close()
    result.add(dat)
    dat.setLen(0)
  echo "free memory ", getFreeMem()
  echo "occ memory ", getOccupiedMem()
  
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
    # create a sequence large enough to hold most events, so that only for very large
    # events we need to resize the sequence
    pixels: Pixels = newSeqOfCap[Pix](400)
    # variable to store resulting chip events
    chips: seq[ChipEvent] = newSeqOfCap[ChipEvent](7)
    # variable to determine, when we are reading the last pixel of one chip
    # important, because we cannot trust numHits in chip header, due to zero
    # suppression!
    pix_to_read: int = 0
  result = new Event
  let
    regex_h = re(regex[0])
    regex_c = re(regex[1])
    regex_p = re(regex[2])

  for line in data:
    if match(line, regex_h, h_matches) == true:
      e_header[h_matches[0]] = h_matches[1]
    elif match(line, regex_c, h_matches) == true:
      # in case we match a chip header, pix_counter to 0,
      # need to make sure it is 0, once we reach the first
      # hit pixel
      pix_counter = 0
      c_header[h_matches[0]] = h_matches[1]
      if h_matches[0] == "numHits":
        let nhits = parseInt(h_matches[1])
        pix_to_read = if nhits < 4096: nhits else: 4095
        if pix_to_read == 0:
          var ch_event = ChipEvent()
          ch_event.chip = (c_header["chipName"], parseInt(c_header["chipNumber"]))
          ch_event.pixels = pixels
          # add  the chip event object to the sequence
          chips.add(ch_event)
    elif match(line, regex_p, p_matches) == true:
      # in this case we have matched a pixel hit line
      # get number of hits to process
      let
        x: int  = parseInt(p_matches[0])
        y: int  = parseInt(p_matches[1])
        ch: int = parseInt(p_matches[2])
      # if the compiler flag (-d:REMOVE_FULL_PIX) is set, we cut all pixels, which have
      # ToT values == 11810, the max ToT value
      when defined(REMOVE_FULL_PIX):
        if ch != 11810:
          pixels.add((x, y, ch))
      else:
        pixels.add((x, y, ch))        
      # after adding pixel, increase pixel counter so that we know when we're
      # reading the last hit of a given chip
      inc pix_counter
      if pix_counter == pix_to_read:
        # now we are reading the last hit, process chip header and pixels
        var ch_event = ChipEvent()
        ch_event.chip = (c_header["chipName"], parseInt(c_header["chipNumber"]))
        ch_event.pixels = pixels
        # add  the chip event object to the sequence
        chips.add(ch_event)
        # reset the pixels sequence
        pixels = newSeqOfCap[Pix](400)
        pix_to_read = 0

  # finally add the timestamp from the dateTime to the table as well
  e_header["timestamp"] = $(int(parseTOSDateString(e_header["dateTime"]).toSeconds))
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
  else:
    echo "Warning: File to dump frame data could not be opened! Does the output folder exist? Path was ", filepath

template writeSeqsToFile[T](filename, header, element_header: string, collection: seq[seq[T]]) =
  # a template to construct different write file procedures with different headers
  var f: File
  let n_el = len(collection)
  for index in 0..<n_el:
    # use index to determine filename
    let fname = replacef(filename, re"(.*)\.txt", by = "$1_" & $index & ".txt")
    let n_in_el = len(collection[index])
    if open(f, fname, fmWrite):
      f.write("# " & header & "for element " & element_header & $index & "\n")
      for i in 0..<n_in_el:
        f.write($collection[index][i] & "\n")
    else:
      echo "Warning: Could not open file " & filename & " to write."

proc writeHitsFile[T](filename: string, collection: seq[seq[T]]) =
  # procedure to write the number of hits in a given run to file by
  # calling writeSeqsToFile template
  let header = "Hits per Event"
  let element_header = "Chip "
  writeSeqsToFile(filename, header, element_header, collection)

proc writeToTFile[T](filename: string, collection: seq[seq[T]]) =
  # procedure to write the ToT values per pixel in a given run to file by
  # calling writeSeqsToFile template
  let header = "ToT per pixel"
  let element_header = "Chip "
  writeSeqsToFile(filename, header, element_header, collection)


proc writeRotAngleFile[T](filename: string, collection: seq[seq[T]]) =
  # procedure to write the ToT values per pixel in a given run to file by
  # calling writeSeqsToFile template
  let header = "Rotation angle"
  let element_header = "Chip "
  writeSeqsToFile(filename, header, element_header, collection)
  

proc dumpToTandHits*(name, run_type: string, tots, hits: seq[seq[int]]) =
  # this procedure dumps the ToT and Hits sequences to .txt files
  # in the out/ folder
  for i in 0..<7:
    echo "Chip " & $i & ":"
    echo "\tToT : " & $len(tots[i]) & "\t" & $len(tots)
    echo "\tHits: " & $len(hits[i])
  let
    outfile = "out/$#_$#_$#.txt"
    hitfile = outfile % ["hits", run_type, extractFilename(name)]
    totfile = outfile % ["tots", run_type, extractFilename(name)]
  writeHitsFile(hitfile, hits)
  writeToTFile(totfile, tots)

proc dumpRotAngle*(angles: seq[seq[float64]]) =
  # this procedure dumps the ToT and Hits sequences to .txt files
  # in the out/ folder
  for i in 0..<7:
    echo "Chip " & $i & ":"
    echo "\tRotAngle : " & $len(angles[i]) & "\t" & $len(angles)
  writeRotAngleFile("out/rotAngle.txt", angles)
  

proc getSortedListOfFiles*(run_folder: string, sort_type: EventSortType, event_type: EventType): seq[string] =
  # this procedure returns a sorted list of event files from a
  # run folder. The returned files are sorted by the event number
  # inputs:
  #    run_folder: string = folder from which to read filenames
  #    sort_type: EventSortType = enum which decides whether we sort by
  #               inode or filename
  #    event_type: EventType = enum which decides which files we read
  # outputs:
  #    seq[string] = sequence containing event filnames, index 0 corresponds
  #                  to event data000000.txt
  var event_regex = ""
  case event_type:
  of EventType.InGridType:
    event_regex = r".*data\d{4,6}\.txt$"
  of EventType.FadcType:
    event_regex = r".*data\d{4,6}\.txt-fadc$"    
  # get the list of files from this run folder and sort it
  case sort_type
  of fname:
    result = sorted(getListOfFiles(run_folder, event_regex),
                    (x: string, y: string) -> int => system.cmp[string](x, y))
  of inode:
    result = map(createSortedInodeTable(getListOfFiles(run_folder, event_regex)),
                 (i: int, name: string) -> string => name)

proc readListOfFiles*[T](list_of_files: seq[string],
                         regex: tuple[header, chips, pixels: string] = ("", "", "")): seq[FlowVar[ref T]] = #{.inline.} =
  ## As this is to be called from a function specifying the datatype, see the calling functions
  ## for descriptions of input and output
  let nfiles = len(list_of_files)
  echo "Reading files into buffer from " & $0 & " to " & $(nfiles - 1)
  # seq of lines from memmapped files
  let mmfiles = readMemFilesIntoBuffer(list_of_files)
  echo "...done reading"

  # create a buffer sequence, into which we store the results processed
  # in parallel (cannot add to the result seq with arbitrary indexes)
  # need ind_high + 1, since newSeq creates a seq with as many elements, while
  # the slicing syntax a[0..10] includes (!) the last element, thus this slice
  # has 11 elements
  result = newSeq[FlowVar[ref T]](nfiles)

  parallel:
    var f_count = 0
    for i, s in mmfiles:
      # loop over each file and call work on data function
      if i < len(result):
        when T is Event:
          result[i] = spawn processEventWithRegex(s, regex)
        elif T is FadcFile:
          result[i] = spawn readFadcFile(s)
      echoFilesCounted(f_count)
  sync()
  

# set experimental pragma to enable parallel: block
{.experimental.}
proc readListOfInGridFiles*(list_of_files: seq[string],
                           regex_tup: tuple[header, chips, pixels: string]): seq[FlowVar[ref Event]] =
  # this procedure receives a list of files, reads them into memory (as a buffer)
  # and processes the content into a seq of ref Events
  # inputs:
  #    list_of_files: seq[string] = a seq of filenames, which are to be read in one go
  #    regex_tup: tuple[...] = a tuple of the different regexes needed to read the different
  #                            parts of a file
  # outputs:
  #    seq[FlowVar[ref Event]] = a seq of flow vars pointing to events, since we read
  #                              in parallel
  result = readListOfFiles[Event](list_of_files, regex_tup)
  

proc isTosRunFolder*(folder: string): tuple[is_rf: bool, contains_rf: bool] =
  # this procedure checks whether the given folder is a valid run folder of
  # TOS
  # done by
  # - checking whether the name of the folder is a valid name for a
  #   run folder (contains Run_<number>) in the name and 
  # - checking whether folder contains data<number>.txt files
  # inputs:
  #    folder: string = the given name of the folder to check
  # outputs:
  # returns a tuple which not only says whether it is a run folder, but also
  # whether the folder itself contains a run folder
  #    tuple[bool, bool]:
  #        is_rf:       is a run folder
  #        contains_rf: contains run folders
  let run_regex = re(r".*Run_(\d+)_.*")
  let event_regex = re(r".*data\d{4,6}\.txt$")
  var matches_rf_name: bool = false
  if match(folder, run_regex) == true:
    # set matches run folder flag to true, is checked when we find
    # a data<number>.txt file in the folder, so that we do not think a
    # folder with a single data<number>.txt file is a run folder
    matches_rf_name = true
    
  for kind, path in walkDir(folder):
    if kind == pcFile:
      if match(path, event_regex) == true and matches_rf_name == true:
        result.is_rf = true
        # in case we found an event in the folder, we might want to stop the
        # search, in order not to waste time. Nested run folders are
        # undesireable anyways
        # for now we leave this comment here, since it may come in handy
        # break
    else:
      # else we deal with a folder. call this function recuresively
      let (is_rf, contains_rf) = isTosRunFolder(path)
      # if the underlying folder contains an event file, this folder thus
      # contains a run folder
      if is_rf == true:
        result.contains_rf = true

# proc sum*[T: tuple](s: seq[T]): T {.inline.} =
#   # this procedure sums the given array along the given axis
#   # if T is itself e.g. a tuple, we will return a tuple, one
#   # element for each field in the tuple
#   assert s.len > 0, "Can't sum empty sequences"
#   var sum_t: T
#   for p in s:
#     for n, f in fieldPairs(p):
#       sum_t[f] += p[n]
  
proc calcCentroidOfEvent*(pix: Pixels): tuple[x, y: float] =
  # proc to calc centroid of the given pixels
  # inputs:
  #    pixels object (seq[tuple[x, y, ch: int]]) containing raw event
  # outputs:
  #    tuple[x, y: int]: tuple containing centroid x and y position
  # let x = map(pix, (p: tuple[x, y, ch: int]) -> int => p.x)
  # let y = map(pix, (p: tuple[x, y, ch: int]) -> int => p.y)
  # let sum_x = foldl(x, a + b)
  # let sum_y = foldl(y, a + b)
  var
    sum_x: int = 0
    sum_y: int = 0
  for p in pix:
    sum_x += p.x
    sum_y += p.y
  #let (sum_x, sum_y, sum_ch) = sum(pix)
  result.x = float(sum_x) / float(len(pix))
  result.y = float(sum_y) / float(len(pix))
  
  
proc isNearCenterOfChip*(pix: Pixels): bool =
  # proc to check whether event is located around center of chip
  # inputs:
  #    pixels object (seq[tuple[x, y, ch: int]]) containing raw event
  # outputs:
  #    true if within 4.5mm center square, false otherwise
  let (center_x, center_y) = calcCentroidOfEvent(pix)
  # pitch in um 
  let pitch = 0.05
  let n_pix_to_bound = 2.25 / pitch
  # center pixel is (127, 127)
  let center_pix = 127'f
  var
    in_x = false
    in_y = false
  if center_x > (center_pix - n_pix_to_bound) and center_x < (center_pix + n_pix_to_bound):
    in_x = true
  if center_y > (center_pix - n_pix_to_bound) and center_y < (center_pix + n_pix_to_bound):  
    in_y = true
  if in_x == true and in_y == true:
    result = true
  else:
    result = false

# template which calculates euclidean distance between 2 points    
template distance*(x, y): float = sqrt(x * x + y * y)

# template which returns pitch converted positions on chip pixel values
# to mm from center of chip
# constants are:
# const NPIX = 256
# const PITCH = 0.0055 (see ingrid_types)
template applyPitchConversion*[T: (float | int)](x, y: T): (float, float) =
  # template which returns the converted positions on a Timepix
  # pixel position --> position from center in mm
  ((float(NPIX) - float(x) - 0.5) * PITCH, (float(y) + 0.5) * PITCH)


proc fillRunHeader*(event: ref Event): Table[string, string] =
  result = initTable[string, string]()
  # run number
  result["runNumber"] = event.evHeader["runNumber"]
  # run time (TOS command, n_shutters or time)
  result["runTimeFrames"] = event.evHeader["runTimeFrames"]
  # run time, setting in `run_time_frames`
  result["runTime"] = event.evHeader["runTime"]
  # path of run folder
  result["pathName"] = event.evHeader["pathName"]
  # date of run start
  result["dateTime"] = event.evHeader["dateTime"]
  # number of chips in run
  result["numChips"] = event.evHeader["numChips"]
  # shutter mode (standard, long, very long)
  result["shutterMode"] = event.evHeader["shutterMode"]
  # time in shutter mode (0-255)
  result["shutterTime"] = event.evHeader["shutterTime"]
  # run mode (ToT or ToA)
  result["runMode"] = event.evHeader["runMode"]
  # clock rate, true: 80MHz, false: 40MHz
  result["fastClock"] = event.evHeader["fastClock"]
  # trigger type (internal = 0, external = 1)
  result["externalTrigger"] = event.evHeader["externalTrigger"]


#####################################################
# Procs describing the data layout in the HDF5 file #
#####################################################

template rawDataBase*(): string =
  "/runs/run_"

template rawDataChipBase*(run_number: int): string =
  "/runs/run_$#/chip_" % $run_number # & "$#"

template recoDataChipBase*(run_number: int): string =
  "/reconstruction/run_$#/chip_" % $run_number # & "$#"  

proc getGroupNameForRun*(run_number: int): string =
  # generates the group name for a given run number
  result = rawDataBase() & "$#" % $run_number

template recoBase*(): string =
  "/reconstruction/run_"

proc getRecoNameForRun*(run_number: int): string =
  # generates the reconstrution group name for a given run number
  result = recoBase() & "$#" % $run_number

proc getRawCombineName*(): string =
  # generates the base path for the combine folder
  result = "/runs/combined/"

proc getRecoCombineName*(): string =
  # generates the base path for the combine folder
  result = "/reconstruction/combined/"


macro createCombineTemplates(name, datatype: string): typed =
  # creates a template, which returns a basename of the type
  # combineBasename`name`(chip_number, run_number): string =
  #   "/runs/combined/`name`_$#_$#" % [$chip_number, $run_number]
  var source = ""
  case datatype.strVal
  of "runs":
    source &= "template combineRawBasename" & name.strVal & "*(chip_number, run_number: int): string =\n"
  of "reconstruction":
    source &= "template combineRecoBasename" & name.strVal & "*(chip_number, run_number: int): string =\n"
  else:
    discard
  case name.strVal
  of "Hits", "ToT":
    source &= "  \"/" & datatype.strVal & "/combined/" & name.strVal
  else:
    source &= "  \"/" & datatype.strVal.toLowerAscii & "/combined/" & name.strVal.toLowerAscii
  source &= "_$#_$#\" % [$chip_number, $run_number]"
  result = parseStmt(source)
  echo toStrLit(result)

createCombineTemplates("ToT", "runs")
createCombineTemplates("Hits", "runs")
createCombineTemplates("ToT", "reconstruction")
createCombineTemplates("Hits", "reconstruction")
# these don't work, since we use the name once in upper and once in lower case
#createCombineTemplates("Fadc", "reconstruction")
#createCombineTemplates("Noisy", "reconstruction")
#createCombineTemplates("Minvals", "reconstruction")
#createCombineTemplates("Noisy", "reconstruction")


# template combineRawBasenameToT*(chip_number, run_number: int): string =
#   "/runs/combined/ToT_$#_$#" % [$chip_number, $run_number]

# template combineRawBasenameHits*(chip_number, run_number: int): string =
#   "/runs/combined/Hits_$#_$#" % [$chip_number, $run_number]

# template combineRecoBasenameToT*(chip_number, run_number: int): string =
#   "/reconstruction/combined/ToT_$#_$#" % [$chip_number, $run_number]

# template combineRecoBasenameHits*(chip_number, run_number: int): string =
#   "/reconstruction/combined/Hits_$#_$#" % [$chip_number, $run_number]

template combineRecoBasenameFadc*(): string =
  "/reconstruction/combined/fadc/"
  
template combineRecoBasenameNoisy*(run_number: int): string =
  "/reconstruction/combined/fadc/noisy_$#" % [$run_number]

template combineRecoBasenameMinvals*(run_number: int): string =
  "/reconstruction/combined/fadc/minvals_$#" % [$run_number]

template noiseBasename*(run_number: int): string =
  getRecoNameForRun(run_number) / "fadc/noisy"

template minvalsBasename*(run_number: int): string =
  getRecoNameForRun(run_number) / "fadc/minvals"

template rawFadcBasename*(run_number: int): string =
  getGroupNameForRun(run_number) / "fadc/raw_fadc"

template trigrecBasename*(run_number: int): string =
  getGroupNameForRun(run_number) / "fadc/trigger_record"
  
template fadcDataBasename*(run_number: int): string =
  getRecoNameForRun(run_number) / "fadc/fadc_data"

template fadcBaselineBasename*(run_number: int): string =
  getRecoNameForRun(run_number) / "fadc/baseline"

template argMinvalBasename*(run_number: int): string =
  getRecoNameForRun(run_number) / "fadc/argMinval"  

template riseStartBasename*(run_number: int): string =
  getRecoNameForRun(run_number) / "fadc/riseStart"

template fallStopBasename*(run_number: int): string =
  getRecoNameForRun(run_number) / "fadc/fallStop"

template riseTimeBasename*(run_number: int): string =
  getRecoNameForRun(run_number) / "fadc/riseTime"

template fallTimeBasename*(run_number: int): string =
  getRecoNameForRun(run_number) / "fadc/fallTime"




################################################################################
##################### HDF5 related helper functions ############################
################################################################################  

iterator runs*(h5f: var H5FileObj, reco = true): (string, string) =
  # simple iterator, which yields the group name of runs
  # in the file. If reco is true (default) we yield
  # reconstruction groups, else raw grous
  var data_basename: string = ""
  if reco == true:
    data_basename = recoBase()
  else:
    data_basename = rawDataBase()

  if h5f.visited == false:
    h5f.visit_file
    
  let run_regex = re(data_basename & r"(\d+)$")
  var run: array[1, string]
  var reco_run: seq[FlowVar[ref RecoEvent]] = @[]
  for grp in keys(h5f.groups):
    if grp.match(run_regex, run) == true:
      # now read some data. Return value will be added later
      #let run_number = parseInt(run[0])
      yield (run[0], grp)

  
  


when isMainModule:

  assert combineBasenameToT(0, 1) == "/runs/combined/ToT_0_1"
  assert combineRecoBasenameToT(0, 1) == "/reconstruction/combined/ToT_0_1"

  echo "All tests passed!"




  
