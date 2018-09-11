import os, ospaths
import strutils, strformat, strscans
import times
import algorithm
import re
import tables
import memfiles
import sequtils, future
import threadpool
import math
import streams, parsecsv
import sets

import typetraits

# cus modules
import helpers/utils
import ingrid_types
import nimhdf5
import seqmath

# other modules
import arraymancer
import loopfusion
import zero_functional

import macros

{.deadCodeElim: on.}

const
  # some helper constants
  StartTot* = 20.0
  # constant regex for InGrid type events
  eventRegexInGrid = r".*data\d{4,6}(_1_[0-9]+)?.*\.txt$"

  # default chip names, shutter modes and shutter times
  OldShutterMode = "verylong"
  OldShutterTime = "13"
  OldChipName = "Christophs"
  OldTosRunDescriptorPrefix* = r".*\/(\d{1,3})-"

proc readToTFile*(filename: string,
                  startRead = 0.0,
                  totPrefix = "TOTCalib"): (int, Tot) =
  ## reads the given TOT file and returns a tuple of seqs containing
  ## the chip number, pulse heights, mean and std values
  let
    dataLines = readFile(filename).splitLines.filterIt(it.len > 0)
  # get the TOTCalib filename prefix and use its length as a search index
  # to find the chip number
  let jumpTo = totPrefix.len
  var chip = 0
  try:
    chip = ($filename.extractFilename[jumpTo]).parseInt
  except ValueError:
    # if we can't extract the chip number from the file, ignore it
    discard
  result[0] = chip
  # create seqs for each column
  var
    pulses: seq[float]
    mean: seq[float]
    std: seq[float]
  try:
    pulses = dataLines.mapIt(it.splitWhitespace[1].parseFloat)
    mean   = dataLines.mapIt(it.splitWhitespace[5].parseFloat)
    # convert RMS (that's the value in the column) to one standard deviation by
    # STD = RMS / sqrt( 4 * 256 * 256 ) = RMS / 512
    std    = dataLines.mapIt(it.splitWhitespace[7].parseFloat / 512.0)
  except IndexError:
    # in this case we're *probably* reading a file, which does not contain any alphabetical
    # characters, so try 0, 1, 2 as indices
    pulses = dataLines.mapIt(it.splitWhitespace[0].parseFloat)
    mean   = dataLines.mapIt(it.splitWhitespace[1].parseFloat)
    # convert RMS (that's the value in the column) to one standard deviation by
    # STD = RMS / sqrt( 4 * 256 * 256 ) = RMS / 512
    std    = dataLines.mapIt(it.splitWhitespace[2].parseFloat / 512.0)

  # get the number of TOT calibration "starts", i.e. 20mV is the starting
  # pulse height, so search for number of these
  var startTot = 0.0
  if startRead > 0.0:
    startTot = startRead
  else:
    startTot = StartTot
  let nstarts = pulses.filterIt(it == startTot).len

  let lastInd = pulses.len - pulses.reversed.find(startToT) - 2
  if lastInd > 0:
    # if there is only a single StartToT value, lastInd will be -1
    pulses.delete(0, lastInd)
    mean.delete(0, lastInd)
    std.delete(0, lastInd)

  # filter out elements with std == 0.0
  let nonZero = zip(std, pulses, mean) --> filter(it[0] > 0.0)
  # see zips above for indices
  result[1].pulses = nonZero.mapIt(it[1].int)
  result[1].mean = nonZero.mapIt(it[2])
  result[1].std = nonZero.mapIt(it[0])

proc readScurveVoltageFile*(filename: string): SCurve =
  ## reads an SCurve file and returns an SCurve object
  let file = filename.expandTilde
  # - read file as string
  # - split all lines after header at \n
  # - filter lines with no content
  # - create tuple of (THL, Counts) for each line
  let dataTuple = readFile(file).splitLines[2..^1].filterIt(it.len > 0).mapIt(
    ((it.split('\t')[0].parseInt, it.split('\t')[1].parseFloat))
  )
  result.name = filename
  result.voltage = file.extractFilename.strip(chars = {'a'..'z', '_', '.'}).parseInt
  result.thl = dataTuple.mapIt(it[0])
  result.hits = dataTuple.mapIt(it[1])

proc sum*(c: seq[Pix]): (int, int, int) {.inline.} =
  ## this procedure sums the sequence of pixels such that it returns
  ## a tuple of (sum_x, sum_y, sum_charge)
  ## Need to return `int`, because the sum may be much larger than
  ## what fits into `uint8` or `uint16` as per `Pix`
  assert c.len > 0, "Can't sum empty sequences"
  for p in c:
    result[0] += p.x.int
    result[1] += p.y.int
    result[2] += p.ch.int

proc sum*[T: SomeInteger](c: seq[(T, T, T)]): (int, int, int) {.inline.} =
  ## same as `sum[Pix]` above, but taking in any integer
  assert c.len > 0, "Can't sum empty sequences"
  for p in c:
    result[0] += p[0].int
    result[1] += p[1].int
    result[2] += p[2].int

proc sum2*(c: seq[Pix]): Pix {.inline, deprecated.} =
  ## this procedure sums the squares of the pixels in the sequence
  ## NOTE: this proc is used nowhere!
  assert c.len > 0, "Can't sum empty sequences"
  for p in c:
    result.x += p.x * p.x
    result.y += p.y * p.y
    result.ch += p.ch * p.ch

proc parseTOSDateString*(date_str: string): Time =
  ## function receives a string from a date time from TOS and creates
  ## a Time object from it
  let date = toTime(parse(date_str, "yyyy-MM-dd'.'hh:mm:ss"))
  return date

proc parseRunType*(runType: string): RunTypeKind =
  ## given a string describing a run type, return the correct
  ## `RunTypeKind`
  if runType.normalize in ["calibration", "calib", "c"]:
    result = rtCalibration
  elif runType.normalize in ["background", "back", "b"]:
    result = rtBackground
  elif runType.normalize in ["xrayfinger", "xray", "x"]:
    result = rtXrayFinger
  else:
    result = rtNone

proc readDateFromEvent*(filepath: string): string =
  ## procedure opens the given file and returns the dateTime value
  ## (TOS syntaxed date string)
  for line in lines filepath:
    if "dateTime" in line:
      # if dateTime is found, split, assign to string and break from while
      let line_seq = split(line, " ")
      result = line_seq[high(line_seq)]
      break
    else:
      continue

proc getFilenameFromEventNumber*[T: SomeInteger](evNumber: T): string =
  ## returns the correct TOS filename for the given event number
  result = &"data{evNumber:06}.txt"

proc formatAsOrgDate*(t: Time, org_format = "yyyy-MM-dd ddd H:mm"): string =
  ## this procedure formats the given Time object as an org-mode date
  ## by first converting it to a TimeInfo object
  let ti = getLocalTime(t)
  result = format(ti, org_format)

proc getTimeFromEvent*(file: string): Time =
  ## this procedure returns the time info from an event,
  ## given the filename by parsing the TOS date syntax
  let date_str = readDateFromEvent(file)
  result = parseTOSDateString(date_str)

proc readListOfFilesAndGetTimes*(path: string, list_of_files: seq[string]): seq[Time] =
  ## function analogues to walkRunFolderAndGetTimes except we read all files
  ## given by list_of_files
  var i: int = 0
  result = @[]

  for file in list_of_files:
    let filename = joinPath(path, file)
    let date_str = readDateFromEvent(filename)
    if date_str != "":
      result.add(parseTOSDateString(date_str))
    if i mod 500 == 0:
      echoFilesCounted(i, msg = " read files and parsed times.")
    i += 1
  return result

proc walkRunFolderAndGetTimes*(folder: string): seq[Time] =
  var i: int = 0
  result = @[]

  for file in walkFiles(joinPath(folder, "data*.txt-fadc")):
    let path = file
    if "data" in path and "fadc" in path:
      let filename = strip(path, false, true, {'-', 'f', 'a', 'd', 'c'})
      let date_str = readDateFromEvent(filename)
      if date_str != "":
        result.add(parseTOSDateString(date_str))
      if i mod 500 == 0:
        echoFilesCounted(i, msg = " files walked and parsed times.")
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
  ## this procdure creates a RunTimeInfo object from a given list of files
  ## in a run folder (not sorted). The list is sorted and then the first and
  ## last event are read, converted to Time objects and the length of the run
  ## is calculated from the difference between both events
  ## sort list of files
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
  ## proc to parse the TOS shutter mode selection
  if mode == "verylong" or mode == "vl":
    result = 2
  elif mode == "long" or mode == "l":
    result = 1
  else:
    result = 0

proc calcLength*(event: Event): float =
  ## calculates the event length in seconds, taking into account
  ## whether the FADC triggered or not and the shutter opening time
  ## TODO: is it really * 46 or not? Kind of important...
  ## 2, 13 comes out correctly with 46!
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
  ## overload of above, with time and mode already read
  ## calculates the event length in seconds, taking into account
  ## whether the FADC triggered or not and the shutter opening time
  var fadc_triggered = 0
  if event.evHeader.hasKey("fadcReadout"):
    fadc_triggered = parseInt(event.evHeader["fadcReadout"])

  if unlikely(fadc_triggered == 1):
    let fadc_clock_triggered = parseFloat(event.evHeader["fadcTriggerClock"])
    # in case FADC triggered it's simply the number of clockcycles the shutter was open
    # divided by the clock speed
    result = fadc_clock_triggered / 40'f / 1_000_000'f
  else:
    result = pow(256'f, mode) * 46'f * time / 40'f / 1_000_000'f

proc getRegexForEvents*(): tuple[header, chips, pixels: string] =
  ## this procedure returns a tuple of the regexes used to
  ## read different parts of the events
  ## outputs:
  ##     tuple[header, chips, pixels: string] =
  ##       header: the regex to read the event header
  ##       chips: the regex to read the chip headers
  ##       pixels: the regex to read the pixel data
  result.header = r"^\#{2}\s(\w+):\s+(-?\b\S*)\b"
  result.chips  = r"^\#\s(\w+):\s+\b(\w\ \d{1,2} W\d{2}|\S*)\b"
  result.pixels = r"^(\d+)\s+(\d+)\s+(\d+)\s?$"

proc readEventHeader*(filepath: string): Table[string, string] =
  ## this procedure reads a whole event header and returns
  ## a table containing the data, where the key is the key from the data file
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

proc readMemFilesIntoBuffer*(list_of_files: seq[string]): seq[seq[string]] =
  ## procedure which reads a list of files via memory mapping and returns
  ## the thus read data as a sequence of MemFiles
  ## inputs:
  ##    list_of_files: seq[string] = seq of strings containing the filenames to be read
  ## outputs:
  ##    seq[seq[string]] = seq containing a seq of data for each file. Zeroth element
  ##      always contains the filename
  result = newSeqOfCap[seq[string]](len(list_of_files))

  var
    ff: MemFile
    # reserver enough space for most events, only for large events do we have to reserve
    # more space
    dat: seq[string] = @[] #newSeqOfCap[string](100)

  echo "free memory ", getFreeMem()
  echo "occ memory ", getOccupiedMem()

  # TODO: is it the smartest way to open and read the memfiles directly?
  # I guess there's a reason why we do not return a seq of memory mapped files
  # anymore
  var lineBuf = newStringOfCap(80)
  for f in list_of_files:
    ff = memfiles.open(f, mode = fmRead, mappedSize = -1)
    dat.add f
    for _ in lines(ff, lineBuf):
      dat.add lineBuf
    ff.close()
    result.add dat
    dat.setLen(0)
  echo "free memory ", getFreeMem()
  echo "occ memory ", getOccupiedMem()

proc processEventWithRegex*(data: seq[string],
                            regex: tuple[header, chips, pixels: Regex]):
                              ref Event {.deprecated.}=
  ## This proc is deprecated, use `processEventWithScanf` instead,
  ## it's a lot faster!
  ## this template is used to create the needed functions to
  ## - read the event header
  ## - read the chip information for an event
  ## - read the pixel data in an event
  ## either as single functions or as a combination of
  ## all
  ## when type(regex) == string:
  ##   # in this case we deal with a single regex string
  ##   result = newTable[string, string]()
  ##   for line in lines filepath:
  ##     if line.match(re(regex), matches):
  ##       # get rid of whitespace and add to result
  ##       let key = matches[0]
  ##       let val = matches[1]
  ##       result[key] = val
  ## regex[0] == header
  ## regex[1] == chips
  ## regex[2] == pixels
  ## in this case we have a sequence of strings
  var
    # variable to count already read pixels
    pix_counter = 0
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

  for line in data[1 .. ^1]:
    # NOTE: match using matching template
    if line =~ regex.header:
      # Event header
      e_header[matches[0]] = matches[1]
    elif line =~ regex.chips:
      # Chip header
      # set pix_counter to 0, need to make sure it is 0,
      # once we reach the first
      # hit pixel
      pix_counter = 0
      c_header[matches[0]] = matches[1]
      if matches[0] == "numHits":
        let nhits = parseInt(matches[1])
        pix_to_read = if nhits < 4096: nhits else: 4095
        if pix_to_read == 0:
          var ch_event = ChipEvent()
          ch_event.chip = (c_header["chipName"], parseInt(c_header["chipNumber"]))
          ch_event.pixels = pixels
          # add the chip event object to the sequence
          chips.add(ch_event)
    elif line =~ regex.pixels:
      # in this case we have matched a pixel hit line
      # get number of hits to process
      let
        x: uint8  = parseInt(matches[0]).uint8
        y: uint8  = parseInt(matches[1]).uint8
        ch: uint16 = parseInt(matches[2]).uint16
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
  result.nChips = chips.len

proc processEventWithScanf*(data: seq[string]): ref Event =
  ## Reads a current TOS event file using the strscans.scanf
  ## macro. It
  ## - reads the event header
  ## - reads the chip information for an event
  ## - reads the pixel data in an event
  ## able to read arbitrary many chips per file
  var
    # variable to count already read pixels
    pix_counter = 0
    # variable to count line in chip header (always 3 lines),
    # 3rd line is important for number of pixels
    cHeaderCount = 0
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
  # get filename as 0th element of data seq
  let filename = data[0]

  var
    # variables to use for `scanf` matching
    keyMatch: string
    valMatch: string
    x: int
    y: int
    ch: int

  const
    # the constants to match against using `scanf`
    headerMatch = "$*:$s${matchNonSpace}"
    pixMatch = "$i$s$i$s$i"

  ##############################
  # some helper procs for matching
  ##############################

  proc matchNonSpace(input: string, strVal: var string, start: int): int =
    ## proc for `scanf` macro to many any ascii character
    var i = 0
    while i + start < input.len:
      strVal.add input[i + start]
      inc i
    result = i

  proc matchEventHeader(line: string,
                        e_header: var Table[string, string],
                        keyMatch, valMatch: var string) {.inline.} =
    ## proc performing the match of `Event Header` part using strscans.scanf.
    ## These lines start with `## `.
    ## The given string has already been stripped of the prefix
    if likely(line[0]!= '[') and
       likely(scanf(line, headerMatch, keyMatch, valMatch)):
      e_header[keyMatch] = valMatch
      keyMatch.setLen(0)
      valMatch.setLen(0)

  proc matchChipHeader(line: string,
                       c_header: var Table[string, string],
                       pixels: var seq[Pix],
                       keyMatch, valMatch: var string,
                       pix_counter, cHeaderCount, pix_to_read: var int,
                       chips: var seq[ChipEvent],
                       filepath: string) {.inline.} =
    ## performs the parsing of the chip header using strscans.scanf.
    ## These line start with `# `. The input has already been stripped
    ## of the prefix.
    pix_counter = 0
    # start from 2 to skip ``# ``
    if likely(scanf(line, headerMatch, keyMatch, valMatch)):
      c_header[keyMatch] = valMatch
      if cHeaderCount == 2:
        # we're in the 3rd line of chip header containing numHits
        # reset cHeaderCount
        cHeaderCount = 0
        assert keyMatch == "numHits", "File seems broken: " & filepath
        let nhits = valMatch.parseInt
        pix_to_read = if nhits < 4096: nhits else: 4095
        if pix_to_read == 0:
          var ch_event = ChipEvent()
          ch_event.chip = (c_header["chipName"], parseInt(c_header["chipNumber"]))
          ch_event.pixels = pixels
          # add the chip event object to the sequence
          chips.add(ch_event)
      else:
        # only increase cHeader count if it wasn't 2
        inc cHeaderCount
      # increase chip header count
      keyMatch.setLen(0)
      valMatch.setLen(0)
    else:
      raise newException(IOError, "ChipHeader: This shouldn't happen. File is broken!" &
        " filename: " & filepath)

  proc matchPixels(line: string,
                   c_header: var Table[string, string],
                   pixels: var seq[Pix],
                   x, y, ch: var int,
                   pix_counter: var int,
                   pix_to_read: var int,
                   chips: var seq[ChipEvent],
                   filename: string) {.inline.} =
    ## matches the lines containing pixels using strscans.scanf macro
    if likely(scanf(line, pixMatch, x, y, ch)):
      # if the compiler flag (-d:REMOVE_FULL_PIX) is set, we cut all pixels, which have
      # ToT values == 11810, the max ToT value
      when defined(REMOVE_FULL_PIX):
        if ch != 11810:
          pixels.add((x.uint8, y.uint8, ch.uint16))
      else:
        pixels.add((x.uint8, y.uint8, ch.uint16))

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
    else:
      raise newException(IOError, "PixMatch: This shouldn't happen. File is broken!" &
          " filename: " & filename)

  ##############################
  # perform the matching using helpers
  ##############################

  # first parse the event header (first 19 lines)
  # we can skip line `data[0 .. 1]`, since that:
  # - data[0] == filename
  # - data[1] == `[General]` header
  for line in data[2 .. 19]:
    # event header
    # start from 3 to skip ``## ``
    line[3 .. line.high].matchEventHeader(e_header,
                                          keyMatch,
                                          valMatch)

  # parse the rest of the file
  for line in data[20 .. ^1]:
    # NOTE: match using matching template
    case line[0]
    of '#':
      # Chip header
      # set pix_counter to 0, need to make sure it is 0,
      # once we reach the first
      # hit pixel
      line[2 .. line.high].matchChipHeader(c_header,
                                           pixels,
                                           keyMatch,
                                           valMatch,
                                           pix_counter,
                                           cHeaderCount,
                                           pix_to_read,
                                           chips,
                                           filename)
    else:
      # in this case we have matched a pixel hit line
      # get number of hits to process
      # match the line with scanf
      line.matchPixels(c_header,
                       pixels,
                       x, y, ch,
                       pix_counter,
                       pix_to_read,
                       chips,
                       filename)

  # finally add the timestamp from the dateTime to the table as well
  e_header["timestamp"] = $(int(parseTOSDateString(e_header["dateTime"]).toSeconds))

  result.evHeader = e_header
  result.chips = chips
  result.nChips = chips.len

proc addOldHeaderKeys(e_header, c_header: var Table[string, string],
                      eventNumber, chipNumber, pix_counter: int,
                      filepath: string) {.inline.} =
  ## inline adds the keys to tables of event header and chip header
  e_header["eventNumber"] = $eventNumber
  # TODO: in this case simply use the "old default constants"
  e_header["shutterMode"] = OldShutterMode
  e_header["shutterTime"] = OldShutterTime
  # only support 1 chip for old storage format anyways
  e_header["numChips"] = $1
  c_header["numHits"] = $pix_counter
  c_header["chipName"] = OldChipName
  # subtract 1, because in the past TOS was 1 indexed
  c_header["chipNumber"] = $(chipNumber - 1)
  let (head, tail) = filepath.splitPath
  e_header["pathName"] = head

proc processOldEventWithRegex*(data: seq[string],
                               regPixels: Regex):
                                 ref OldEvent {.deprecated.} =
  ## This proc is deprecated, use `processEventWithScanf` instead,
  ## it's a lot faster!  ## TODO: update doc!
  # in this case we have a sequence of strings
  var
    # variable to count already read pixels
    pix_counter = 0
    # variables to read data into
    e_header = initTable[string, string]()
    c_header = initTable[string, string]()
    # create a sequence large enough to hold most events, so that only for very large
    # events we need to resize the sequence
    pixels: Pixels = newSeqOfCap[Pix](400)
    # variable to store resulting chip events
    chips: seq[ChipEvent] = newSeqOfCap[ChipEvent](1)
    # variable to determine, when we are reading the last pixel of one chip
    # important, because we cannot trust numHits in chip header, due to zero
    # suppression! We init by `-1` to deal with old TOS format, in which there
    # is no header
    pix_to_read: int = 0
  result = new OldEvent

  let filepath = data[0]
  # check for run folder kind by looking at first line of file
  assert data[1][0] != '#', "This is not a valid rfOldTos file (old storage format)!" & data[1]

  for line in data[1 .. ^1]:
    # match with template
    if line =~ regPixels:
      # in this case we have matched a pixel hit line
      # get number of hits to process
      let
        x  = parseInt(matches[0]).uint8
        y  = parseInt(matches[1]).uint8
        ch = parseInt(matches[2]).uint16
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

  # in that case we're reading an ``old event files (!)``. Get the event number
  # from the filename
  let evNumberRegex = re".*data(\d{4,6})_(\d)_.*"
  var evNumChipNumStr: array[2, string]
  if match(filepath, evNumberRegex, evNumChipNumStr) == true:
    addOldHeaderKeys(e_header,
                     c_header,
                     evNumChipNumStr[0].parseInt,
                     evNumChipNumStr[1].parseInt,
                     pix_counter,
                     filepath)

  # now we are reading the last hit, process chip header and pixels
  var ch_event = ChipEvent()
  ch_event.chip = (c_header["chipName"], c_header["chipNumber"].parseInt)
  ch_event.pixels = pixels
  # add  the chip event object to the sequence
  chips.add(ch_event)

  result.evHeader = e_header
  result.chips = chips
  result.nChips = chips.len

proc processOldEventScanf*(data: seq[string]): ref OldEvent =
# proc processOldEventScanf*(tup: tuple[fname: string, dataStream: MemMapFileStream]): ref OldEvent =
  ## Reads an Old TOS zero suppressed data file using the strscans.scanf
  ## macro
  ## - read the pixel data in an event
  ## the event and chip header are added from sane defaults for the
  ## 2014/15 data (defined as constans at top of `tos_helpers.nim`).

  # in this case we have a sequence of strings
  var
    # variable to count already read pixels
    pix_counter = 0
    # variables to read data into
    e_header = initTable[string, string]()
    c_header = initTable[string, string]()
    # create a sequence large enough to hold most events, so that only for very large
    # events we need to resize the sequence
    pixels: Pixels = newSeqOfCap[Pix](400)
    # variable to store resulting chip events
    chips: seq[ChipEvent] = newSeqOfCap[ChipEvent](1)
    # variable to determine, when we are reading the last pixel of one chip
    # important, because we cannot trust numHits in chip header, due to zero
    # suppression! We init by `-1` to deal with old TOS format, in which there
    # is no header
    pix_to_read: int = 0
  result = new OldEvent

  let filepath = data[0]
  # check for run folder kind by looking at first line of file
  assert data[1][0] != '#', "This is not a valid rfOldTos file (old storage format)!" & data[1]

  var
    x: int
    y: int
    ch: int

  for line in data[1 .. ^1]:
    # match with scanf
    if scanf(line, "$i$s$i$s$i", x, y, ch):
      # if the compiler flag (-d:REMOVE_FULL_PIX) is set, we cut all pixels, which have
      # ToT values == 11810, the max ToT value
      when defined(REMOVE_FULL_PIX):
        if ch != 11810:
          pixels.add((x.uint8, y.uint8, ch.uint16))
      else:
        pixels.add((x.uint8, y.uint8, ch.uint16))
      # after adding pixel, increase pixel counter so that we know when we're
      # reading the last hit of a given chip
      inc pix_counter

  # in that case we're reading an ``old event files (!)``. Get the event number
  # from the filename
  var
    evNumber: string
    chipNumber: string
    dummy: string
  # need to use `$*` to parse, because if we use $i for integers, we end up
  # parsing the `_` tokens as well, making scanf fail
  # TODO: implement custom parser proc
  if scanf(filepath, r"$*data$*_$*_", dummy, evNumber, chipNumber):
    addOldHeaderKeys(e_header,
                     c_header,
                     evNumber.parseInt,
                     chipNumber.parseInt,
                     pix_counter,
                     filepath)

  # now we are reading the last hit, process chip header and pixels
  var ch_event = ChipEvent()
  ch_event.chip = (c_header["chipName"], c_header["chipNumber"].parseInt)
  ch_event.pixels = pixels
  # add  the chip event object to the sequence
  chips.add(ch_event)

  result.evHeader = e_header
  result.chips = chips
  result.nChips = chips.len


proc processEventWrapper(data: seq[string],
                         regex: tuple[header, chips, pixels: Regex],
                         rfKind: RunFolderKind): ref Event {.inline.} =
  ## wrapper around both process event procs, which determines which one to call
  ## based on the run folder kind. Need a wrapper, due to usage of spawn in
  ## caling prof `readListOfFiles`.
  case rfKind
  of rfNewTos:
    #result = processEventWithRegex(data, regex)
    result = processEventWithScanf(data)
  of rfOldTos:
    # result = processOldEventWithRegex(data, regex[2])
    result = processOldEventScanf(data)

#template readEventWithRegex(filepath, regex: string): typed =
proc readEventWithRegex*(filepath: string, regex: tuple[header, chips, pixels: Regex]): ref Event =
  ## this procedure reads the lines from a given event and then calls processEventWithRegex
  ## to process the file. Returns a ref to an Event object
  ## inputs:
  ##    filepath: string = the path to the file to be read
  ##    regex: tuple[...] = tuple of 3 regex's, one for each part of an event file
  ## outputs:
  ##    ref Event: reference to an event object
  var f = memfiles.open(filepath, mode = fmRead, mappedSize = -1)
  var data: seq[string] = @[]
  for line in lines(f):
    data.add(line)
  result = processEventWithRegex(data, regex)

proc pixelsToTOT*(pixels: Pixels): seq[uint16] {.inline.} =
  ## extracts all charge values (ToT values) of a given pixels object (all hits
  ## of a single chip in a given event, filters the full values of 11810
  ## inputs:
  ##    pixels: Pixels: input pixel sequence of tuples containing hits of one event
  ## output:
  ##    seq[int]: all ToT values smaller ToT
  result = filter(map(pixels, (p: tuple[x, y: uint8, ch: uint16]) -> uint16 => p.ch),
                  (ch: uint16) -> bool => ch < 11810'u16)

template addPixelsToOccupancy*[T](ar: Tensor[T], pixels: Pixels) =
  ## template to add pixels to occupancy by using map
  for p in pixels:
    ar[p.x, p.y] += 1#p.ch

template addPixelsToOccupancySeptem*[T](ar: var Tensor[T], pixels: Pixels, ch_num: int) =
  ## template to add pixels to occupancy by using map
  for p in pixels:
    ar[ch_num, p.x.int, p.y.int] += 1#p.ch

proc createTensorFromZeroSuppressed*[T](pixels: Pixels): Tensor[T] =
  ## procedure to create a (256, 256) int array from a Pixels (seq[tuple[x, y, ch]])
  ## object
  result = zeros[T](256, 256)
  for p in pixels:
    result[p.x, p.y] = p.ch

proc dumpFrameToFile*[T](filepath: string, ar: Tensor[T]) =
  ## this procedure dumps a given frame (tensor ar needs to be of shape (256, 256)
  ## to a file 'filepath'
  ## inputs:
  ##   filepath: string = the file to write to
  ##   ar: Tensor[int] = a tensor of shape (256, 256) containing the data to be written
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
  ## a template to construct different write file procedures with different headers
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
  ## procedure to write the number of hits in a given run to file by
  ## calling writeSeqsToFile template
  let header = "Hits per Event"
  let element_header = "Chip "
  writeSeqsToFile(filename, header, element_header, collection)

proc writeToTFile[T](filename: string, collection: seq[seq[T]]) =
  ## procedure to write the ToT values per pixel in a given run to file by
  ## calling writeSeqsToFile template
  let header = "ToT per pixel"
  let element_header = "Chip "
  writeSeqsToFile(filename, header, element_header, collection)


proc writeRotAngleFile[T](filename: string, collection: seq[seq[T]]) =
  ## procedure to write the ToT values per pixel in a given run to file by
  ## calling writeSeqsToFile template
  let header = "Rotation angle"
  let element_header = "Chip "
  writeSeqsToFile(filename, header, element_header, collection)


proc dumpToTandHits*(name, run_type: string, tots, hits: seq[seq[int]]) =
  ## this procedure dumps the ToT and Hits sequences to .txt files
  ## in the out/ folder
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
  ## this procedure dumps the ToT and Hits sequences to .txt files
  ## in the out/ folder
  for i in 0..<7:
    echo "Chip " & $i & ":"
    echo "\tRotAngle : " & $len(angles[i]) & "\t" & $len(angles)
  writeRotAngleFile("out/rotAngle.txt", angles)

proc sortByInode*(listOfFiles: seq[string]): seq[string] =
  ## given a list of files, sorts them by Inode using `createSortedInodeTable`
  result = createSortedInodeTable(listOfFiles).map(
    (i: int, name: string) -> string => name
  )

proc getSortedListOfFiles*(run_folder: string, sort_type: EventSortType, event_type: EventType): seq[string] =
  ## this procedure returns a sorted list of event files from a
  ## run folder. The returned files are sorted by the event number
  ## inputs:
  ##    run_folder: string = folder from which to read filenames
  ##    sort_type: EventSortType = enum which decides whether we sort by
  ##               inode or filename
  ##    event_type: EventType = enum which decides which files we read
  ## outputs:
  ##    seq[string] = sequence containing event filnames, index 0 corresponds
  ##                  to event data000000.txt
  var eventRegex = ""
  case event_type:
  of EventType.InGridType:
    eventRegex = eventRegexInGrid#r".*data\d{4,6}\.txt$"
  of EventType.FadcType:
    eventRegex = r".*data\d{4,6}\.txt-fadc$"
  # get the list of files from this run folder and sort it
  case sort_type
  of fname:
    result = sorted(getListOfFiles(run_folder, eventRegex),
                    (x: string, y: string) -> int => system.cmp[string](x, y))
  of inode:
    result = sortByInode(getListOfFiles(run_folder, eventRegex))

proc readListOfFiles*[T](list_of_files: seq[string],
                         regex: tuple[header, chips, pixels: string] = ("", "", "")):
                           seq[FlowVar[ref T]] = #{.inline.} =
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
  var regex_tup = (re(regex[0]), re(regex[1]), re(regex[2]))

  when T is Event:
    # determine the run folder kind
    let rfKind = if mmfiles[0][1][0] == '#': rfNewTos else: rfOldTos
    echo "rf kind is ", rfKind

  parallel:
    var f_count = 0
    for i, s in mmfiles:
      # loop over each file and call work on data function
      if i < len(result):
        when T is Event:
          result[i] = spawn processEventWrapper(s, regex_tup, rfKind)
        elif T is FadcFile:
          result[i] = spawn readFadcFile(s)
      echoFilesCounted(f_count)
  sync()

# set experimental pragma to enable parallel: block
{.experimental.}
proc readListOfInGridFiles*(list_of_files: seq[string],
                            regex_tup: tuple[header, chips, pixels: string]):
                              seq[FlowVar[ref Event]] =
  ## this procedure receives a list of files, reads them into memory (as a buffer)
  ## and processes the content into a seq of ref Events
  ## inputs:
  ##    list_of_files: seq[string] = a seq of filenames, which are to be read in one go
  ##    regex_tup: tuple[...] = a tuple of the different regexes needed to read the different
  ##                            parts of a file
  ## outputs:
  ##    seq[FlowVar[ref Event]] = a seq of flow vars pointing to events, since we read
  ##                              in parallel
  result = readListOfFiles[Event](list_of_files, regex_tup)

proc readListOfOldInGridFiles*(list_of_files: seq[string],
                            regex_tup: tuple[header, chips, pixels: string]):
                              seq[FlowVar[ref Event]] =
  ## see documentation of above
  # since `OldEvent` is simply an alias for us, it can be returned as
  # an `Event` as well
  result = readListOfFiles[OldEvent](list_of_files, regex_tup)

proc isTosRunFolder*(folder: string):
  tuple[is_rf: bool, runNumber: int, rfKind: RunFolderKind, contains_rf: bool] =
  ## this procedure checks whether the given folder is a valid run folder of
  ## TOS
  ## done by
  ## - checking whether the name of the folder is a valid name for a
  ##   run folder (contains Run_<number>) in the name and
  ## - checking whether folder contains data<number>.txt files
  ## inputs:
  ##    folder: string = the given name of the folder to check
  ## outputs:
  ## returns a tuple which not only says whether it is a run folder, but also
  ## whether the folder itself contains a run folder
  ##    tuple[(bool, int), bool]:
  ##        (is_rf, runNumber): is a run folder, its Number
  ##        contains_rf:        contains run folders

  let runRegex = re(r".*Run_(\d+)_.*")
  # else check for old `Chistoph style` runs
  let oldRunRegexStr = r".*Run\d{6}_\d{2}-\d{2}-\d{2}.*"
  let oldRunRegex = re(oldRunRegexStr)
  # TODO: check whether following works
  const oldTosRunDescriptorPrefix = OldTosRunDescriptorPrefix
  var oldTosRunDescriptor: Regex

  let eventRegex = re(eventRegexInGrid) #r".*data\d{4,6}(_1_[0-9]+)?.*\.txt$")
  var matches_rf_name: bool = false
  var runNumber: array[1, string]
  if match(folder, runRegex, runNumber) == true:
    # set matches run folder flag to true, is checked when we find
    # a data<number>.txt file in the folder, so that we do not think a
    # folder with a single data<number>.txt file is a run folder
    matches_rf_name = true
    result.runNumber = runNumber[0].parseInt
    result.rfKind = rfNewTos
  elif match(folder, oldRunRegex) == true:
    matches_rf_name = true
    # in case of the old tos, extract the run number from the folder name
    result.rfKind = rfOldTos
    let (head, tail) = folder.splitPath
    oldTosRunDescriptor = re(oldTosRunDescriptorPrefix & oldRunRegexStr)
    if folder =~ oldTosRunDescriptor:
      result.runNumber = matches[0].parseInt


  for kind, path in walkDir(folder):
    if kind == pcFile:
      if match(path, eventRegex) == true and matches_rf_name == true:
        result.is_rf = true
        # found an event file -> is run folder -> break
        break
    else:
      # else we deal with a folder. call this function recursively
      let (is_rf, runNumber, rfKind, contains_rf) = isTosRunFolder(path)
      # if the underlying folder contains an event file, this folder thus
      # contains a run folder
      if is_rf == true:
        result.contains_rf = true

proc getOldRunInformation*(folder: string, runNumber: int, rfKind: RunFolderKind):
  (int, int, int, int, int) =
  ## given a TOS raw data run folder of kind `rfOldTos`, parse information based
  ## on the `*.dat` file contained in it
  case rfKind
  of rfOldTos:
    const oldTosRunDescriptorPrefix = OldTosRunDescriptorPrefix
    let
      (head, tail) = folder.splitPath
      datFile = joinPath(folder, $runNumber & "-" & tail & ".dat")
      lines = readFile(datFile).splitLines
      # now just parse the files correctly. Everything in last column, except
      # runTime
      totalEvents = lines[0].splitWhitespace[^1].parseInt
      numEvents = lines[1].splitWhitespace[^1].parseInt
      startTime = lines[2].splitWhitespace[^1].parseInt
      stopTime = lines[3].splitWhitespace[^1].parseInt
      runTime = lines[4].splitWhitespace[^2].parseInt
    result = (totalEvents, numEvents, startTime, stopTime, runTime)
  else: discard

proc parseOldTosRunlist*(path: string, rtKind: RunTypeKind): set[uint16] =
  ## parses the run list and returns a set of integers, corresponding to
  ## the valid run numbers of the `RunTypeKind` given
  var s = newFileStream(path, fmRead)
  defer: s.close()
  var typeStr: string
  case rtKind
  of rtCalibration:
    typeStr = "C"
  of rtBackground:
    typeStr = "B"
  of rtXrayFinger:
    typeStr = "X"
  else:
    return {}

  var csv: CsvParser
  open(csv, s, path)
  while readRow(csv):
    if csv.row[2] == typeStr:
      result.incl csv.row[0].strip.parseInt.uint16
  csv.close()

proc extractRunNumber*(runFolder: string): int =
  ## given a valid TOS run folder, extract its run Number
  let (is_rf, runNumber, rfKind, contains_rf) = isTosRunFolder(runFolder)
  result = runNumber

proc extractRunFolderKind*(runFolder: string): RunFolderKind =
  ## given a valid TOS run folder, extract its run Number
  let (is_rf, runNumber, rfKind, contains_rf) = isTosRunFolder(runFolder)
  result = rfKind


################################################################################
############# Geometry calculation related procs ###############################
################################################################################

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
  ## proc to calc centroid of the given pixels
  ## inputs:
  ##    pixels object (seq[tuple[x, y, ch: int]]) containing raw event
  ## outputs:
  ##    tuple[x, y: int]: tuple containing centroid x and y position
  ## let x = map(pix, (p: tuple[x, y, ch: int]) -> int => p.x)
  ## let y = map(pix, (p: tuple[x, y, ch: int]) -> int => p.y)
  ## let sum_x = foldl(x, a + b)
  ## let sum_y = foldl(y, a + b)
  var
    sum_x: int = 0
    sum_y: int = 0
  for p in pix:
    sum_x += p.x.int
    sum_y += p.y.int
  #let (sum_x, sum_y, sum_ch) = sum(pix)
  result.x = float(sum_x) / float(len(pix))
  result.y = float(sum_y) / float(len(pix))


proc isNearCenterOfChip*(pix: Pixels): bool =
  ## proc to check whether event is located around center of chip
  ## inputs:
  ##    pixels object (seq[tuple[x, y, ch: int]]) containing raw event
  ## outputs:
  ##    true if within 4.5mm center square, false otherwise
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
template applyPitchConversion*[T: (float | SomeInteger)](x, y: T): (float, float) =
  ## template which returns the converted positions on a Timepix
  ## pixel position --> position from center in mm
  ((float(NPIX) - float(x) - 0.5) * PITCH, (float(y) + 0.5) * PITCH)

proc fillRunHeader*(event: Event): Table[string, string] =
  ## TODO: what's the reason for this functions existence?
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
#### Procs specifically realted to hardware #########
#####################################################

proc getSeptemHChip*(chipNumber: int): string =
  ## returns the name of a given SeptemH chip
  const names = ["E6 W69",
                 "K6 W69",
                 "H9 W69",
                 "H10 W69",
                 "G10 W69",
                 "D9 W69",
                 "L8 W69"]
  result = names[chipNumber]

#####################################################
# Procs describing the data layout in the HDF5 file #
#####################################################

proc getFloatGeometryNames*(): array[12, string] =
  ## returns all dataset names in the H5 output file, which are members
  ## of a `ClusterGeometry` object
  result = ["rmsLongitudinal", "rmsTransverse", "skewnessLongitudinal", "skewnessTransverse",
            "kurtosisLongitudinal", "kurtosisTransverse", "eccentricity", "rotationAngle",
            "length", "width", "fractionInTransverseRms", "lengthDivRmsTrans"]

proc getFloatClusterNames*(): array[2, string] =
  ## returns all dataset names in the H5 output file, which are members of
  ## a `ClusterObject` object
  result = ["centerX", "centerY"]

proc getFloatDsetNames*(): array[14, string] =
  ## returns the names of all datasets in the H5 output file, which appear as
  ## (N, 1) data columns. Combination of two above procs
  # need to define consts of arrays to use `+` macro
  const
    float_geo = getFloatGeometryNames()
    float_obj = getFloatClusterNames()
  result = float_geo + float_obj

proc getIntClusterNames*(): array[2, string] =
  ## returns names of datasets in H5 output file, which are integer datasets and
  ## members of a `ClusterObject`
  result = ["hits", "sumTot"]

proc getIntDsetNames*(): array[3, string] =
  ## returns all names of integer dataset for the H5 output file, which appear
  ## as (N, 1) data columns
  # need to define consts of arrays to use `+` macro
  const
    int_cluster = getIntClusterNames()
    int_dset = ["eventNumber"]
  result = int_cluster + int_dset

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

template likelihoodBase*(): string =
  "/likelihood/run_"

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
  ## creates a template, which returns a basename of the type
  ## combineBasename`name`(chip_number, run_number): string =
  ##   "/runs/combined/`name`_$#_$#" % [$chip_number, $run_number]
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

template trigRecBasename*(run_number: int): string =
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

template eventNumberBasename*(run_number: int): string =
  getRecoNameForRun(run_number) / "fadc/eventNumber"


################################################################################
##################### procs related to X-ray reference datasets ################
################################################################################

template cdlPrefix*(): string =
  ## return the prefix of the group names in the `calibration-cdl.h5` file
  ## part of the names of the reference names below
  "calibration-cdl-apr2014-"

proc getXrayRefTable*(): Table[int, string] =
  ## returns a table mapping the different energy bins to the correct
  ## datasets in the X-ray reference file
  # NOTE: we could also simply store this in a seq...
  result = { 0: "C-EPIC-0.6kV",
             1: "Cu-EPIC-0.9kV",
             2: "Cu-EPIC-2kV",
             3: "Al-Al-4kV",
             4: "Ag-Ag-6kV",
             5: "Ti-Ti-9kV",
             6: "Mn-Cr-12kV",
             7: "Cu-Ni-15kV" }.toTable()

proc getEnergyBinning(): seq[float] =
  ## returns the binning of the energy (upper range for each bin)
  ## as a sequence of floats
  result = @[0.4, 0.7, 1.2, 2.1, 3.2, 4.9, 6.9, Inf]

proc toRefDset*(energy: float): string =
  ## returns the correct X-ray reference table for a given
  ## `energy`
  # define xray table as global to only initialize it once
  const
    xray_table = getXrayRefTable()
    binning = getEnergyBinning()
  let ind = binning.lowerBound(energy)
  result = xray_table[ind]

proc getEnergyBinMinMaxVals*(): Table[string, Cuts] =
  ## returns a table of Cuts objects, one for each energy bin
  let
    range0 = Cuts(minCharge: 0.0,
                  maxCharge: 5e4,
                  minRms: 0.1,
                  maxRms: 20.0,
                  maxLength: 6.0,
                  minPix: 3)
    range1 = Cuts(minCharge: 3.0e4,
                  maxCharge: 8.0e4,
                  minRms: 0.0,
                  maxRms: 1.1,
                  maxLength: 6.0,
                  minPix: 3)
    range2 = Cuts(minCharge: 7.0e4,
                  maxCharge: 1.3e5,
                  minRms: 0.0,
                  maxRms: 1.1,
                  maxLength: 7.0,
                  minPix: 3)
    range3 = Cuts(minCharge: 9.0e4,
                  maxCharge: 2.1e5,
                  minRms: 0.0,
                  maxRms: 1.1,
                  maxLength: 7.0,
                  minPix: 3)
    range4 = Cuts(minCharge: 2.0e5,
                  maxCharge: 4.0e5,
                  minRms: 0.0,
                  maxRms: 1.1,
                  maxLength: 7.0,
                  minPix: 3)
    range5 = Cuts(minCharge: 2.9e5,
                  maxCharge: 5.5e5,
                  minRms: 0.0,
                  maxRms: 1.1,
                  maxLength: 7.0,
                  minPix: 3)
    range6 = Cuts(minCharge: 3.5e5,
                  maxCharge: 6.0e5,
                  minRms: 0.0,
                  maxRms: 1.1,
                  maxLength: 7.0,
                  minPix: 3)
    range7 = Cuts(minCharge: 5.9e5,
                  maxCharge: 1.0e6,
                  minRms: 0.0,
                  maxRms: 1.1,
                  maxLength: 7.0,
                  minPix: 3)
    ranges = [range0, range1, range2, range3, range4, range5, range6, range7]
    xray_ref = getXrayRefTable()

  result = initTable[string, Cuts]()
  for key, vals in pairs(xray_ref):
    result[vals] = ranges[key]

proc getRegionCut*(region: ChipRegion): CutsRegion =
  const
    xMinChip = 0.0
    xMaxChip = 14.0
    yMinChip = 0.0
    yMaxChip = 14.0

  case region
  of crGold:
    result = CutsRegion(xMin: 4.5,
                        xMax: 9.5,
                        yMin: 4.5,
                        yMax: 9.5,
                        radius: 0.0)
  of crSilver:
    # based on radius of 4.5 from center
    result = CutsRegion(xMin: 0.0,
                        xMax: 0.0,
                        yMin: 0.0,
                        yMax: 0.0,
                        radius: 4.5)
  of crBronze:
    result = CutsRegion(xMin: 0.0,
                        xMax: 0.0,
                        yMin: 0.0,
                        yMax: 0.0,
                        radius: 5.5)
  of crAll:
    result = CutsRegion(xMin: 0.0,
                        xMax: 14.0,
                        yMin: 0.0,
                        yMax: 0.0,
                        radius: 0.0)



################################################################################
##################### HDF5 related helper functions ############################
################################################################################

proc getTrackingEvents*(h5f: var H5FileObj, group: H5Group, num_tracking: int = -1, tracking = true): seq[int] =
  ## given a `group` in a `h5f`, filter out all indices, which are part of
  ## a tracking (`tracking == true`) or not part of a tracking (`tracking == false`)
  ## NOTE: the indices of the whole timestamp array correspond to the event
  ## numbers of a run, since this is a 1:1 mapping
  result = @[]
  # attributes of this group
  var attrs = group.attrs
  try:
    # try except for check of num_trackings
    let ntrackings = attrs["num_trackings", int]
    const date_syntax = getDateSyntax()

    var
      tr_starts = newSeq[DateTime](ntrackings)
      tr_stops  = newSeq[DateTime](ntrackings)
    for i in 0 ..< ntrackings:
      # get start and stop time of each tracking
      tr_starts[i] = attrs[&"tracking_start_{i}", string].parse(date_syntax)
      tr_stops[i]  = attrs[&"tracking_stop_{i}", string].parse(date_syntax)
    # get the timestamp of all events
    let
      tstamp = h5f[(group.name / "timestamp").dset_str][int64]
      # start and stop in seconds
      tr_starts_s = mapIt(tr_starts, int(it.toTime.toSeconds))
      tr_stops_s  = mapIt(tr_stops,  int(it.toTime.toSeconds))
    # filter out all indices of timestamps, which lie inside the tracking

    # first get all indices of all trackings in a seq[seq[int]]
    var allTrackingInds: seq[seq[int]] = @[]
    for k in 0 ..< ntrackings:
      allTrackingInds.add filterIt(toSeq(0 ..< tstamp.len)) do:
        tstamp[it] > tr_starts_s[k] and tstamp[it] < tr_stops_s[k]
    if tracking == true:
      if num_tracking >= 0:
        # simply get the correct indices from allTrackingInds
        result = allTrackingInds[num_tracking]
      else:
        # flatten the allTrackingInds nested seq and return
        result = flatten(allTrackingInds)
    else:
      # all outside trackings are simply the indices, which are not part of a flattened
      # allTrackingInds
      let allTrackingsFlat = flatten(allTrackingInds)
      # and now filter all indices not part of flattened index
      result = filterIt(toSeq(0 ..< tstamp.len)) do:
        it notin allTrackingsFlat
  except KeyError:
    # in this case there is no tracking information. Keep all indices
    echo &"No tracking information in {group.name} found, use all clusters"
    result = @[]


proc filterTrackingEvents*[T: SomeInteger](cluster_events: seq[T], tracking_inds: seq[int]): seq[int] =
  ## filters out all indices (= event numbers) of a reconstructed run for one chip
  ## Need to remove all indices, which are within the tracking indices, but for which
  ## no cluster is found in the datasets, so that we can only read the clusters, which
  ## happened during (or outside) of a tracking
  ## inputs:
  ##   `cluster_events`: all events for one reconstructed chip
  ##   `tracking_inds`: the indices which are part (or not) of a tracking
  # set result to the indices of tracking (non tracking), i.e.
  # all allowed events
  if tracking_inds.len == 0:
    # in case we are handed an empty seq, we simply use all cluster events
    result = toSeq(0 .. cluster_events.high)
  else:
    # create capped sequence of max possible length `cluster_events`
    result = newSeqOfCap[int](cluster_events.len)
    # using allowed events get indices for other events by iterating
    # over all allowed events and removing those, which are not
    # in the events of a chip
    for ind in tracking_inds:
      # remove all events of the allowed events, which are not
      # part of the events for one chip
      if ind in cluster_events:
        # if index in cluster events, add it
        result.add find(cluster_events, ind)

proc filterTrackingEvents*(h5f: var H5FileObj, group: H5Group, tracking_inds: seq[int]): seq[int] =
  ## wrapper around the above proc, which reads the data about which events are allowed
  ## by itself
  ## inputs:
  ##   `h5f`: H5file from which to read the data
  ##   `group`: H5Group object of the specific chip, which contains the clustes
  ##   `tracking_inds`: the indices which are part (or not) of a tracking
  let
    # event numbers of clusters of this chip
    evNumbers = h5f[(group.name / "eventNumber").dset_str][int64]
  result = filterTrackingEvents(evNumbers, tracking_inds)

iterator runs*(h5f: var H5FileObj, data_basename = recoBase()): (string, string) =
  ## simple iterator, which yields the run number and group name of runs in the file.
  ## If reco is true (default) we yield reconstruction groups, else raw groups
  ## Iterator saves the state of `h5f` during the first call to this iterator! If
  ## additional groups are added while iterating, they will be ignored.
  if h5f.visited == false:
    h5f.visit_file

  # get a seq of all groups currently in the H5 file
  # only want to iterate over groups existing at the time, when
  # this proc is being called.
  # If we just use the `keys` iterator for `h5f.groups` we end up
  # skipping runs randomly, since we insert new groups, changing the
  # iterator while iterating. Bad! Solves issue #8
  let groups = toSeq(keys(h5f.groups))
  let runRegex = re(data_basename & r"(\d+)$")
  var run: array[1, string]
  for grp in groups:
    if grp.match(runRegex, run) == true:
      # now read some data. Return value will be added later
      yield (run[0], grp)

iterator dsets*(h5f: var H5FileObj,
                dsetName: string,
                dtype: typedesc,
                chipNumber: int,
                dataBasename = recoBase()): (int, seq[dtype]) =
  ## given a dataset name, its corresponding datatype (needed to define the return type)
  ## and a chip number, read all datasets of all runs stored in the given H5 file.
  ## Choose a base location, by default reconstruction group
  ## NOTE: this cannot yield any datatypes with variable length data!

  if h5f.visited == false:
    h5f.visit_file

  let runChipName = joinPath(dataBasename, r"run_(\d+)/chip_" & $chipNumber)
  let dsetPath = joinPath(runChipName, dsetName)
  let dsetLocationReg = re(dsetPath)
  var runNumber = newSeq[string](1)
  for dset in keys(h5f.datasets):
    if match(dset, dsetLocationReg, runNumber):
      # found a matching dataset, yield the group number as well as the actual
      # data
      var mdset = h5f[dsetPath.dset_str]
      echo mdset.name
      yield (runNumber[0].parseInt, mdset[dtype])

when isMainModule:

  assert combineRawBasenameToT(0, 1) == "/runs/combined/ToT_0_1"
  assert combineRecoBasenameToT(0, 1) == "/reconstruction/combined/ToT_0_1"

  let energies = @[0.1, 0.0, 12.4, 4.4, 2.3, 2.0]
  let inds = [0, 0, 7, 5, 4, 3]
  let refs = ["C-EPIC-0.6kV", "C-EPIC-0.6kV", "Cu-Ni-15kV", "Ti-Ti-9kV", "Ag-Ag-6kV", "Al-Al-4kV"]
  let xray_table = getXrayRefTable()
  let binning = getEnergyBinning()
  forEach e in energies, i in inds, r in refs:
    assert(binning.lowerBound(e) == i)
    assert(toRefDset(e) == r)


  echo "All tests passed!"
