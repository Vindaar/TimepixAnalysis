import strutils, strformat, strscans, terminal, os, parseutils, times, algorithm
import re, tables, memfiles, sequtils, sugar

when compileOption("threads"):
  when not defined(gcDestructors):
    import threadpool_simple
  else:
    import weave
import math
import streams, parsecsv

# cus modules
import helpers/utils
import ../ingrid_types
# zero functional is fine, because it's purely Nim
import zero_functional

const
  # some helper constants
  StartTot* = 20.0
  # constant regex for InGrid type events for the Virtex TOS
  eventRegexVirtex = r".*data\d{4,9}(_1_[0-9]+)?.*\.txt$"
  # constant regex for InGrid type events for the SRS TOS
  eventRegexSrs = r".*run_(\d{6})_data_(\d{6,9})_(\d{6})_(\d{2})-(\d{2})-(\d{2}).txt"

  newVirtexRunRegex = r".*Run_(\d+)_\d{6}-\d{2}-\d{2}.*"
  oldVirtexRunRegex = r".*Run(\d{6})_(\d{2})-(\d{2})-(\d{2}).*"
  srsRunRegex       = r".*Run_(\d{6})_\d{6}_\d{2}-\d{2}-\d{2}.*"

  OldVirtexEventScanf = r"$*data$*_$*_" # needs: dummy, evNumber, chipNumber
  SrsEventScanf = r"$*run_$*_data_$*_${parseSrsTosDate}" # needs: dummy, runNumber, evNumber, date
  OldVirtexEventScanfNoPath = r"data$*_$*_$i." # needs: dummy, evNumber, chipNumber, timestamp
  NewVirtexEventScanfNoPath* = r"data$*.txt$." # needs: evNumber
  FadcEventScanfNoPath = r"data$*.txt-fadc" # needs: evNumber
  SrsEventScanfNoPath = r"run_$*_data_$*_${parseSrsTosDate}" # needs: dummy, runNumber, evNumber, date

  # default chip names, shutter modes and shutter times
  OldShutterMode = "verylong"
  OldShutterTime = "13"
  OldChipName* = "D3 W63"
  OldTosRunDescriptorPrefix* = r".*\/(\d{1,3})-"

proc readToTFileTpx1*(filename: string,
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
  # let nstarts = pulses.filterIt(it == startTot).len

  let lastInd = pulses.len - pulses.reversed.find(startToT) - 2
  if lastInd > 0:
    # if there is only a single StartToT value, lastInd will be -1
    pulses.delete(0, lastInd)
    mean.delete(0, lastInd)
    std.delete(0, lastInd)

  # filter out elements with std == 0.0
  let nonZero = zip(std, pulses, mean) --> filter(it[0] > 0.0)
  # see zips above for indices
  result[1].pulses = nonZero.mapIt(it[1])
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
  result = toTime(parse(date_str, TosDateString))

proc parseRunType*(runType: string): RunTypeKind =
  ## given a string describing a run type, return the correct
  ## `RunTypeKind`
  if runType.normalize in ["calibration", "calib", "c", "rtcalibration"]:
    result = rtCalibration
  elif runType.normalize in ["background", "back", "b", "rtbackground"]:
    result = rtBackground
  elif runType.normalize in ["xrayfinger", "xray", "x", "rtxray"]:
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
  let ti = local(t)
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
  for file in list_of_files:
    let filename = joinPath(path, file)
    let date_str = readDateFromEvent(filename)
    if date_str != "":
      result.add(parseTOSDateString(date_str))
    echoCount(i, msg = " read files and parsed times.")

proc walkRunFolderAndGetTimes*(folder: string): seq[Time] =
  var i: int = 0
  for file in walkFiles(joinPath(folder, "data*.txt-fadc")):
    let path = file
    if "data" in path and "fadc" in path:
      let filename = strip(path, false, true, {'-', 'f', 'a', 'd', 'c'})
      let date_str = readDateFromEvent(filename)
      if date_str != "":
        result.add(parseTOSDateString(date_str))
      echoCount(i, msg = " files walked and parsed times.")

proc writeDateSeqToFile*(date_seq: seq[Time]): void =
  let outfile = "times.txt"
  var f: File
  if open(f, outfile, fmWrite):
    f.write("# date\t unix time\n")
    for time in date_seq:
      let sec = $toUnix(time)
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

proc parseSrsRunInfo*(path: string): Table[string, string] =
  result = initTable[string, string]()
  let runPath = path / "run.txt"
  if fileExists(runPath):
    # all good, can parse it
    # run number, start and end time
    # won't be considered
    var s = openFileStream(runPath)
    defer: s.close()
    var line = ""
    var
      inParameters = false
      inChipIds = false

    while readLine(s, line):
      if line.startsWith("Run parameters:"):
        inParameters = true
        continue
      elif line.startsWith("Chip IDs:"):
        inParameters = false
        inChipIds = true
        continue
      elif line.len == 0 or line[0] != '\t':
        continue
      elif inParameters:
        let val = line.splitWhitespace[^1]
        var key = ""
        if "run mode" in line:
          key = "runMode"
        elif "run time" in line:
          key = "runTime"
        elif "shutter mode" in line:
          key = "shutterMode"
        elif "shutter time" in line:
          key = "shutterTime"
        if key.len > 0:
          result[key] = val
      elif inChipIds:
        # parse the chip ids
        var
          fec = 0
          board = 0
          chip = 0
          colRow = ""
          wafer = ""
        if scanf(line, "$sFEC $i Board $i Chip $i: $w-$w", fec, board, chip, colRow, wafer):
          # subtract 1 from chip number to get 0 indexing
          let
            chipName = &"chip_{chip - 1}"
            nChips = result.getOrDefault("numChips", "0").parseInt
          result[chipName] = &"{colRow} {wafer}"
          # add number of chips
          result["numChips"] = $(nChips + 1)
    # finally correct `runTime` and add `runTimeFrames`
    if result["runMode"] == "1":
      result["runTimeFrames"] = result["runTime"]
      result["runTime"] = "0"
    else:
      result["runTimeFrames"] = "1"

    # if `inChipIds` is still false, there was no chip ID information in the
    # `run.txt` (old format). Add note to run header
    if not inChipIds:
      result[SrsNoChipId] = SrsNoChipIdMsg
  else:
    # doesn't exist, mark run as incomplete and assign
    # dummy values for information stored in run.txt
    result[SrsRunIncomplete] = SrsRunIncompleteMsg
    result["shutterMode"] = "0"
    result["shutterTime"] = "0"

# forward declaration
proc getListOfEventFiles*(folder: string, eventType: EventType,
                          rfKind: RunFolderKind): seq[(int, string)]

proc estimateRunTime(files: seq[(int, string)],
                     runStart: DateTime): (DateTime, int64) =
  ## estimates the total run time of the given run by walking over all files
  ## (to make sure we do not miss a whole day) and returns the stop time of
  ## the run as well as the run time as an integer of seconds
  ## NOTE: the list of files needs to be sorted!
  var
    days = 0
    lastEvNum = -1
    lastHour = -1
    dummy: string
    timestamp: int
  for tup in files:
    let (evNum, evName) = tup
    doAssert lastEvNum < evNum, " List of files MUST be sorted!"
    let (_, tail) = evName.splitPath
    discard scanf(tail, OldVirtexEventScanfNoPath, dummy, dummy, timestamp)
    let tstamp = ($timestamp).align(9, padding = '0')
    let hour = tstamp[0 .. 1].parseInt
    if lastHour == 23 and hour == 0:
      inc days
    lastHour = hour
    lastEvNum = evNum
  # timestamp still points to last event of
  let tstamp = ($timestamp)
    .align(9, padding = '0')
    .parse("HHmmssfff")
  var runStop = runStart
  let daysDur = initDuration(days = days)
  runStop = runStop + daysDur
  runStop.hour = tstamp.hour
  runStop.minute = tstamp.minute
  runStop.second = tstamp.second
  let runTime = (runStop - runStart).inSeconds
  result = (runStop, runTime)

proc getOldRunInformation*(folder: string, runNumber: int, rfKind: RunFolderKind):
  (int, int, int64, int64, int64) =
  ## given a TOS raw data run folder of kind `rfOldTos`, parse information based
  ## on the `*.dat` file contained in it
  case rfKind
  of rfOldTos:
    const oldTosRunDescriptorPrefix = OldTosRunDescriptorPrefix
    let
      (_, tail) = folder.splitPath
      datFile = joinPath(folder, tail & ".dat")
    try:
      let
        lines = readFile(datFile).splitLines
        # now just parse the files correctly. Everything in last column, except
        # runTime
        totalEvents = lines[0].splitWhitespace[^1].parseInt
        numEvents = lines[1].splitWhitespace[^1].parseInt
        startTime = lines[2].splitWhitespace[^1].parseInt
        stopTime = lines[3].splitWhitespace[^1].parseInt
        runTime = lines[4].splitWhitespace[^2].parseInt
      result = (totalEvents, numEvents,
                startTime.int64, stopTime.int64, runTime.int64)
    except IOError:
      # `*.dat` file does not exist for current run. Instead derive some rough
      # guesses
      let files = sortedByIt(getListOfEventFiles(folder, EventType.InGridType, rfOldTos),
                             it[0])
      let numEvents = files.len
      # get upper limit of number of events based on last event number
      let totalEvents = files[^1][0]
      let oldTosRunDescriptor = re(oldTosRunDescriptorPrefix & oldVirtexRunRegex)
      var runStart = fromUnix(0).inZone(utc())
      if folder =~ oldTosRunDescriptor:
        runStart = matches[1].parse("yyMMdd")
        runStart.hour = matches[2].parseInt
        runStart.minute = matches[3].parseInt
        runStart.second = matches[3].parseInt
      let (runStop, runTime) = estimateRunTime(files, runStart)
      result = (totalEvents, numEvents,
                runStart.toTime.toUnix,
                runStop.toTime.toUnix,
                runTime)

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

proc getRunHeader*(ev: Event,
                   runNumber: int,
                   rfKind: RunFolderKind): Table[string, string] =
  ## returns the correct run header based on the `RunFolderKind`
  case rfKind
  of rfNewTos:
    result = ev.evHeader
  of rfOldTos:
    result = ev.evHeader
    # combine evHeader data and data stored (potentially) in
    # the `*.dat` file
    let
      runPath = result["pathName"]
      runInfo = getOldRunInformation(runPath, runNumber, rfOldTos)
      (totalEvents, numEvents, startTime, stopTime, runTime) = runInfo
      start = fromUnix(startTime)
      stop = fromUnix(stopTime)
      mid = ((stop - start).inSeconds div 2 + start.toUnix).fromUnix
    result["numEvents"] = $numEvents
    result["totalEvents"] = $totalEvents
    result["dateTime"] = start.format("yyyy-MM-dd'.'hh:mm:ss")
    result["dateMid"] = mid.format("yyyy-MM-dd'.'hh:mm:ss")
    result["dateStop"] = stop.format("yyyy-MM-dd'.'hh:mm:ss")
    result["runTime"] = $runTime
  of rfSrsTos:
    # here we need to combine evHeader data and data stored in
    # `run.txt`
    # event header as basis
    result = ev.evHeader
    let runPath = result["pathName"]
    let runInfo = parseSrsRunInfo(runPath)
    # echo write (potentially replace) keys read from run.txt in result
    for key, val in pairs(runInfo):
      result[key] = val
  else:
    discard


proc readMemFilesIntoBuffer*(list_of_files: seq[string]): seq[ProtoFile] =
  ## procedure which reads a list of files via memory mapping and as a single
  ## string and returns a seq of `ProtoFile` objects, which can be parsed in parallel
  ## inputs:
  ##    list_of_files: seq[string] = seq of strings containing the filenames to be read
  ## outputs:
  ##    seq[ProtoFile] = seq containing a seq of data for each file
  result = newSeq[ProtoFile](len(list_of_files))
  var badCount = 0
  echo "free memory ", getFreeMem()
  echo "occ memory ", getOccupiedMem()
  var ff: MemFile
  for i, f in list_of_files:
    # add filename to result
    try:
      ff = memfiles.open(f)
      result[i].name = f
      result[i].fileData = newString(ff.size)
      doAssert ff.size > 0, "Size of data in memory mapped file is 0! File " & $f
      copyMem(result[i].fileData[0].addr, ff.mem, ff.size)
      ff.close()
    except OSError:
      # file exists, but is completely empty. Probably HDD ran full!
      # in this case remove the filename from the `dat` seq again
      when not defined(release):
        echo "Warning: Discarding ", f, " since it cannot be opened!"
      inc badCount
      continue
  echo "free memory ", getFreeMem()
  echo "occ memory ", getOccupiedMem()
  if badCount > 0:
    stdout.styledWrite(fgRed, "WARNING: Number of broken files: ", $badCount)
  result.setLen(result.len - badCount)

proc processEventWithScanf*(data: ProtoFile): Event =
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

  let filename = data.name
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
    if likely(line[0] != '[') and
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
        if keyMatch != "numHits":
          raise newException(IOError, "ChipHeader: This shouldn't happen. File is broken!" &
            " filename: " & filepath)
        let nhits = valMatch.parseInt
        pix_to_read = if nhits < 4096: nhits else: 4095
        if pix_to_read == 0:
          var ch_event = ChipEvent(version: Timepix1)
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
        var ch_event = ChipEvent(version: Timepix1)
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
  # first parse the event header (first 18 lines)
  # we can skip line `data[0 .. 1]`, since that:
  # - data[0] == `[General]` header
  var
    fidx, lineCnt = 0
  var line: string
  while fidx < data.fileData.len:
    fidx += parseUntil(data.fileData, line, '\n', fidx)
    inc fidx
    if lineCnt < 1:
      # skip the general header
      inc lineCnt
      continue
    # event header goes until line including 18
    if lineCnt < 19:
      line[3 .. line.high].matchEventHeader(e_header,
                                            keyMatch,
                                            valMatch)
    else:
      # NOTE: match using matching template
      case line[0]
      of '#':
        # Chip header
        try:
          line[2 .. line.high].matchChipHeader(c_header,
                                               pixels,
                                               keyMatch,
                                               valMatch,
                                               pix_counter,
                                               cHeaderCount,
                                               pix_to_read,
                                               chips,
                                               filename)
        except IOError:
          # broken file, `isValid` will be false
          return
      else:
        # in this case we have matched a pixel hit line
        # get number of hits to process
        # match the line with scanf
        try:
          line.matchPixels(c_header,
                           pixels,
                           x, y, ch,
                           pix_counter,
                           pix_to_read,
                           chips,
                           filename)
        except IOError:
          # broken file, `isValid` will be false
          return
    # increase line counter
    inc lineCnt

  # finally add the timestamp from the dateTime to the table as well
  e_header["timestamp"] = $(parseTOSDateString(e_header["dateTime"]).toUnix)

  result.evHeader = e_header
  result.chips = chips
  result.nChips = chips.len
  # classify as a valid event
  result.isValid = true

proc addOldHeaderKeys(e_header, c_header: var Table[string, string],
                      eventNumber, chipNumber, timestamp, pix_counter: int,
                      filepath: string) {.inline.} =
  ## inline adds the keys to tables of event header and chip header
  ## NOTE: timestamp is an int, but it's still in the format, e.g.:
  ## 103647153 == "HHmmssfff", where fff is milliseconds
  ## and we STORE IT AS SUCH A STRING!
  ## Conversion to the proper timestamp can only be done, after the
  ## start time of the whole run is considered! Done after call to `getRunHeader`
  ## in `raw_data_manipulation`
  e_header["eventNumber"] = $eventNumber
  # TODO: in this case simply use the "old default constants"
  e_header["shutterMode"] = OldShutterMode
  e_header["shutterTime"] = OldShutterTime
  # only support 1 chip for old storage format anyways
  e_header["numChips"] = $1
  e_header["timestamp"] = $timestamp
  c_header["numHits"] = $pix_counter
  c_header["chipName"] = OldChipName
  # subtract 1, because in the past TOS was 1 indexed
  c_header["chipNumber"] = $(chipNumber - 1)
  let (head, _) = filepath.splitPath
  e_header["pathName"] = head

proc processOldEventScanf*(data: ProtoFile): OldEvent =
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

  let filepath = data.name
  # check for run folder kind by looking at first line of file
  doAssert data.fileData[0] != '#', "This is not a valid rfOldTos file (old storage format)! " & data.name
  var
    x: int
    y: int
    ch: int

  var
    fidx, lineCnt = 0
  var line: string
  while fidx < data.fileData.len:
    fidx += parseUntil(data.fileData, line, '\n', fidx)
    inc fidx
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
      inc lineCnt

  # in that case we're reading an ``old event files (!)``. Get the event number
  # from the filename
  var
    evNumber: string
    chipNumber: string
    timestamp: int # timestamp as integer, although it actually is a
                   # time in 24h format!
  let fname = filepath.extractFilename
  if scanf(fname, OldVirtexEventScanfNoPath, evNumber, chipNumber, timestamp):
    addOldHeaderKeys(e_header,
                     c_header,
                     evNumber.parseInt,
                     chipNumber.parseInt,
                     timestamp,
                     lineCnt,
                     filepath)

  # now we are reading the last hit, process chip header and pixels
  var ch_event = ChipEvent(version: Timepix1)
  ch_event.chip = (c_header["chipName"], c_header["chipNumber"].parseInt)
  ch_event.pixels = pixels
  # add  the chip event object to the sequence
  chips.add(ch_event)

  result.evHeader = e_header
  result.chips = chips
  result.nChips = chips.len
  # classify as a valid event file
  result.isValid = true

proc parseSrsTosDate(input: string, date: var Time, start: int): int =
  ## proc for `scanf` macro to parse an SRS TOS date string from a filename
  ## example filename:
  ## run_004001_data_079135_181009_23-19-28.txt
  ## we parse from:         ^ here         ^ to here (excl.)
  # length of parsing range is 15
  result = 15 # 6 from year / month / day, 6 from hour / min / sec, 3 '_','-'
  var stop = start + result
  let dateStr = input[start ..< stop]
  const SrsTosDateSyntax = "yyMMdd'_'HH-mm-ss"
  date = dateStr.parse(SrsTosDateSyntax).toTime
  assert input[stop] == '.', "Assertion failed. Token was: " & $input[stop]

proc processSrsEventScanf*(data: ProtoFile): SrsEvent =
  ## Reads an SRS TOS zero suppressed data file using the strscans.scanf
  ## macro
  ## - read the pixel data in an event
  # in this case we have a sequence of strings
  var
    # variables to read data into
    e_header = initTable[string, string]()
    c_header = initTable[string, string]()
    # create a sequence large enough to hold most events, so that only for very large
    # events we need to resize the sequence
    pixels: Pixels = newSeqOfCap[Pix](400)
    # variable to store resulting chip events
    chips: seq[ChipEvent] = newSeqOfCap[ChipEvent](1) # 1 chip at least, pot. more
    # variable to determine, when we are reading the last pixel of one chip
    # important, because we cannot trust numHits in chip header, due to zero
    # suppression! We init by `-1` to deal with old TOS format, in which there
    # is no header
  let filepath = data.name
  # check for run folder kind by looking at first line of file
  doAssert data.fileData[0 .. 2] == "FEC", "This is not a valid rfSrsTos file!" & data.name
  var
    x: int
    y: int
    ch: int
    valMatch: int
    lineCnt = 1
    cnt: int
    line = ""

  proc parseVal(line: string,
                key: string,
                start: int,
                valMatch: var int,
                c_header: var Table[string, string]): int =
    result = start
    result += skipUntil(line, ' ', start = result)
    result += parseInt(line, valMatch, start = result + 1)
    if key == "chipNumber":
      # subtract 1 from chip number, since SRS TOS still starts
      # countint at 1
      valMatch -= 1
    c_header[key] = $valMatch

  template incLine(line: var string,
                   lineCnt: var int,
                   fidx: var int) =
    fidx += parseUntil(data.fileData, line, '\n', fidx)
    inc fidx

  var fidx= 0
  while fidx < data.fileData.len:
    incLine(line, lineCnt, fidx)
    cnt = parseVal(line, "FEC", 0, valMatch, c_header)
    incLine(line, lineCnt, fidx)
    cnt = parseVal(line, "Board", 0, valMatch, c_header)
    incLine(line, lineCnt, fidx)
    cnt = parseVal(line, "chipNumber", 0, valMatch, c_header)
    # start from cnt + 1 to skip space after chip number to get to `,Hits`
    cnt = parseVal(line, "numHits", cnt + 2, valMatch, c_header)
    let nPix = parseInt(c_header["numHits"])
    let pixToRead = if nPix > 4096: 4096 else: nPix
    # now iterate over all hits
    for i in 0 ..< pixToRead:
      incLine(line, lineCnt, fidx)
      if scanf(line, "$i$s$i$s$i", x, y, ch):
        # if the compiler flag (-d:REMOVE_FULL_PIX) is set, we cut all pixels, which have
        # ToT values == 11810, the max ToT value
        when defined(REMOVE_FULL_PIX):
          if ch != 11810:
            pixels.add((x.uint8, y.uint8, ch.uint16))
        else:
          pixels.add((x.uint8, y.uint8, ch.uint16))
    # once we're done with all pixels, add chip header
    # now we are reading the last hit, process chip header and pixels
    var ch_event = ChipEvent(version: Timepix1)
    ch_event.chip = (SrsDefaultChipName, parseInt(c_header["chipNumber"]))
    ch_event.pixels = pixels
    # add  the chip event object to the sequence
    chips.add(ch_event)
    # reset the pixels sequence
    pixels = newSeqOfCap[Pix](400)

  # in that case we're reading an ``old event files (!)``. Get the event number
  # from the filename
  var
    evNumber: string
    runNumber: string
  # need to use `$*` to parse, because if we use $i for integers, we end up
  # parsing the `_` tokens as well, making scanf fail

  # assign chips and nChips to store in event header
  result.chips = chips
  result.nChips = chips.len

  var date: Time
  let fname = filepath.extractFilename
  if scanf(fname, SrsEventScanfNoPath, runNumber, evNumber, date):
    e_header["dateTime"] = $date
    e_header["timestamp"] = $(date.toUnix)
    e_header["eventNumber"] = evNumber.strip(trailing = false, chars = {'0'})
    e_header["runNumber"] = runNumber.strip(trailing = false, chars = {'0'})
    e_header["numChips"] = $result.nChips
    let (head, _) = filepath.splitPath
    e_header["pathName"] = head

    # since we strip the '0' characters, we might strip everything, if the
    # event or run number is number 0. Thus check if nothing is left and
    # if so set to 0 manually.
    if e_header["eventNumber"].len == 0:
      e_header["eventNumber"] = "0"
    if e_header["runNumber"].len == 0:
      e_header["runNumber"] = "0"
  else:
    raise newException(IOError, "SRS filename does not match `scanf` syntax! " &
      "Filename: " & filepath)
  result.evHeader = e_header
  # classify as a valid event
  result.isValid = true

proc processEventWrapper(data: ProtoFile, rfKind: RunFolderKind): Event =
  ## wrapper around both process event procs, which determines which one to call
  ## based on the run folder kind. Need a wrapper, due to usage of spawn in
  ## caling prof `readListOfFiles`.
  case rfKind
  of rfNewTos:
    result = processEventWithScanf(data)
  of rfOldTos:
    result = processOldEventScanf(data)
  of rfSrsTos:
    result = processSrsEventScanf(data)
  of rfUnknown:
    raise newException(IOError, "Unknown run folder kind. Unclear what files " &
      "are events!")

proc pixelsToTOT*(pixels: Pixels): seq[uint16] {.inline.} =
  ## extracts all charge values (ToT values) of a given pixels object (all hits
  ## of a single chip in a given event, filters the full values of 11810
  ## inputs:
  ##    pixels: Pixels: input pixel sequence of tuples containing hits of one event
  ## output:
  ##    seq[int]: all ToT values smaller ToT
  result = filter(map(pixels, (p: tuple[x, y: uint8, ch: uint16]) -> uint16 => p.ch),
                  (ch: uint16) -> bool => ch < 11810'u16)

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

proc getListOfEventFiles*(folder: string, eventType: EventType,
                          rfKind: RunFolderKind): seq[(int, string)] =
  if existsDir(folder) == false:
    return result
  var
    dummy: string
    dummyInt: int
    evNumber: string
  for file in walkDirRec(folder):
    let fname = file.extractFilename
    case event_type:
    of EventType.InGridType:
      case rfKind
      of rfNewTos:
        if scanf(fname, NewVirtexEventScanfNoPath, evNumber, dummy):
          result.add (evNumber.parseInt, file)
      of rfOldTos:
        if scanf(fname, OldVirtexEventScanfNoPath, evNumber, dummy, dummyInt):
          result.add (evNumber.parseInt, file)
      of rfSrsTos:
        var t: Time
        if scanf(fname, SrsEventScanfNoPath, dummy, evNumber, t):
          result.add (evNumber.parseInt, file)
      else:
        raise newException(IOError, "Unknown run folder kind. Unclear what files " &
          "are events!")
    of EventType.FadcType:
      if scanf(fname, FadcEventScanfNoPath, evNumber, dummy, dummy):
        result.add (evNumber.parseInt, file)

proc getSortedListOfFiles*(run_folder: string,
                           sort_type: EventSortType,
                           event_type: EventType,
                           rfKind: RunFolderKind): seq[string] =
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
  # get the list of files from this run folder and sort it
  case sort_type
  of fname:
    # sort by parsed event number from the above proc
    let evNumFiles = sortedByIt(getListOfEventFiles(run_folder, event_type, rfKind),
                        it[0])
    # no need for the event numbers anymore
    result = evNumFiles --> map(it[1])
  of inode:
    let evNumFiles = getListOfEventFiles(run_folder, event_type, rfKind)
    let files = evNumFiles --> map(it[1])
    result = sortByInode(files)

when compileOption("threads"):
  proc unwrapAndRemoveInvalid[T](res: var seq[T], flow: seq[FlowVar[T]]) =
    ## assigns all valid events from `flow` into `res`
    var validCount = 0
    for i, x in mpairs(res):
      when not defined(gcDestructors):
        let ev = ^flow[i]
      else:
        let ev = sync flow[i]
      if ev.isValid:
        x = ev
        inc validCount
    if validCount < res.len:
      res.setLen(validCount)


  proc readListOfFiles*[T](list_of_files: seq[string],
                           rfKind: RunFolderKind = rfUnknown):
                             seq[T] = #{.inline.} =
    ## As this is to be called from a function specifying the datatype, see the calling functions
    ## for descriptions of input and output
    # backward compatible flag to disable reading FADC files straight
    # from memory mapping
    const fadcMemFiles = true

    let nfiles = len(list_of_files)
    echo "Reading files into buffer from " & $0 & " to " & $(nfiles - 1)
    # seq of lines from memmapped files
    result = newSeq[T](nfiles)
    when not defined(gcDestructors):
      var res = newSeq[FlowVar[T]](nfiles)
      var f_count = 0
      let p = newThreadPool()
      when T is Event or not fadcMemFiles:
        let protoFiles = readMemFilesIntoBuffer(list_of_files)
        echo "...done reading"
        for i, s in protoFiles:
          # loop over each file and call work on data function
          if i < len(result):
            when T is Event:
              res[i] = p.spawn processEventWrapper(s, rfKind)
            elif T is FadcFile:
              res[i] = p.spawn readFadcFile(s)
          echoCount(f_count)
      elif T is FadcFile and fadcMemFiles:
        # should be faster than not using memory mapping
        for i, f in list_of_files:
          # loop over each file and call work on data function
          if i < len(result):
            res[i] = p.spawn readFadcFileMem(f)
          echoCount(f_count)
      echo "Now sync and unwrap"
      p.sync()
      # finally await flowvars and assign to result
      unwrapAndRemoveInvalid(result, res)
      echo "Done syncing and unwrapping"
    else:
      init(Weave)
      var resBuf = cast[ptr UncheckedArray[T]](result[0].addr)
      when T is Event:
        let protoFiles = readMemFilesIntoBuffer(list_of_files)
        echo "...done reading"
        let numFiles = protoFiles.len
        var ppBuf = cast[ptr UncheckedArray[ProtoFile]](protoFiles[0].unsafeAddr)
        parallelFor i in 0 ..< numFiles:
          # loop over each file and call work on data function
          captures: {resBuf, ppBuf, rfKind}
          resBuf[i] = processEventWrapper(ppBuf[i], rfKind)
          echoCounted(i)
      else:
        let numFiles = result.len
        let inputBuf = cast[ptr UncheckedArray[string]](list_of_files[0].unsafeAddr)
        parallelFor i in 0 ..< numFiles:
          # loop over each file and call work on data function
          captures: {resBuf, inputBuf, rfKind}
          resBuf[i] = readFadcFileMem(inputBuf[i])
          echoCounted(i)
      exit(Weave)
      echo "Exit Weave!"

  proc readListOfInGridFiles*(list_of_files: seq[string], rfKind: RunFolderKind):
                            seq[Event] =
    ## this procedure receives a list of files, reads them into memory (as a buffer)
    ## and processes the content into a seq of ref Events
    ## inputs:
    ##    list_of_files: seq[string] = a seq of filenames, which are to be read in one go
    ## outputs:
    ##    seq[FlowVar[ref Event]] = a seq of flow vars pointing to events, since we read
    ##                              in parallel
    result = readListOfFiles[Event](list_of_files, rfKind)

#else:
#  proc readListOfFiles*[T](list_of_files: seq[string],
#                           rfKind: RunFolderKind = rfUnknown):
#                             seq[T] = #{.inline.} =
#    {.error: "`readListOfFiles` needs thread support! Compile with `--threads:on`.".}
#
#  proc readListOfInGridFiles*(list_of_files: seq[string], rfKind: RunFolderKind):
#                            seq[Event] =
#    {.error: "`readListOfInGridFiles` needs thread support! Compile with `--threads:on`.".}


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

  let runRegex = re(newVirtexRunRegex)
  # else check for old `Chistoph style` runs
  let oldRunRegex = re(oldVirtexRunRegex)
  let srsRunRegex = re(srsRunRegex)
  # TODO: check whether following works
  const oldTosRunDescriptorPrefix = OldTosRunDescriptorPrefix
  var oldTosRunDescriptor: Regex

  var matches_rf_name: bool = false
  var runNumber: array[1, string]
  result.rfKind = rfUnknown
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
    oldTosRunDescriptor = re(oldTosRunDescriptorPrefix & oldVirtexRunRegex)
    if folder =~ oldTosRunDescriptor:
      result.runNumber = matches[0].parseInt
  elif match(folder, srsRunRegex, runNumber) == true:
    matches_rf_name = true
    # in case of the old tos, extract the run number from the folder name
    result.runNumber = runNumber[0].parseInt
    result.rfKind = rfSrsTos

  var eventRegex: Regex
  case result.rfKind
  of rfOldTos, rfNewTos:
    eventRegex = re(eventRegexVirtex)
  of rfSrsTos:
    eventRegex = re(eventRegexSrs)
  else:
    # in this case `eventRegex` will stay unknown, because we are probably
    # not in a run folder.
    result.is_rf = false
    # set regex to a "default" regex, which matches any TOS kind
    eventRegex = re(r".*data.*\.txt")
    # however, run folder kind is unknown in this case!
    result.rfKind = rfUnknown
  for kind, path in walkDir(folder):
    if kind == pcFile:
      if match(path, eventRegex) == true and matches_rf_name == true:
        result.is_rf = true
        # found an event file -> is run folder -> break
        break
    else:
      # else we deal with a folder. call this function recursively
      let (is_rf, _, _, _) = isTosRunFolder(path)
      # if the underlying folder contains an event file, this folder thus
      # contains a run folder
      if is_rf == true:
        result.contains_rf = true

proc findLastEvent(files: seq[string]): (string, Time) =
  var
    lastFile: string
    timeLast: Time
  for i in countdown(files.high, 0):
    lastFile = files[i]
    try:
      echo "Trying file: ", lastFile, " at index ", i
      # determine last event and get time from it.
      timeLast = getTimeFromEvent(last_file)
      echo "File ", lastFile, " is valid!"
      break
    except ValueError:
      # value error might be raised, if the file is empty and thus parsing the
      # date time string fails. In this case just continue with the next file
      continue
  result = (lastFile, timeLast)

proc getRunTimeInfo*(run_files: seq[string]): RunTimeInfo =
  ## this procdure creates a RunTimeInfo object from a given list of files
  ## in a run folder (not sorted). The list is sorted and then the first and
  ## last event are read, converted to Time objects and the length of the run
  ## is calculated from the difference between both events
  ## sort list of files
  ## NOTE: If for some reason the very last event does not contain information,
  ## because e.g. the file is broken etc., we try to event before that until
  ## we find a suitable file, which will be considered the end of the run
  result = RunTimeInfo()
  let
    sortedFiles = sorted(run_files, system.cmp[string])
    firstFile = sortedFiles[0]
    # and times
    timeFirst = getTimeFromEvent(firstFile)

  let (_, timeLast) = findLastEvent(sortedFiles)

  # calc run length
  let runLength = timeLast - timeFirst

  echo "Time first is $# and time last is $#" % [$timeFirst, $timeLast]
  echo "Time difference in seconds $#" % $((timeLast - timeFirst).inSeconds)
  echo "Time difference start end $#" % $(runLength)

  result.t_start = timeFirst
  result.t_end = timeLast
  result.t_length = runLength

proc getRunInfo*(path: string): RunInfo =
  ## wrapper around the above proc if only the path to the run is known
  let regex = r"^/([\w-_]+/)*data\d{6,9}\.txt$"
  let fadcRegex = r"^/([\w-_]+/)*data\d{6,9}\.txt-fadc$"
  let (_, runNumber, rfKind, _) = isTosRunFolder(path)
  let files = getListOfFiles(path, regex)
  let fadcFiles = getListOfFiles(path, fadcRegex)
  if files.len > 0:
    result.timeInfo = getRunTimeInfo(files)
  result.runNumber = runNumber
  result.rfKind = rfKind
  result.runType = rtNone
  result.path = path
  result.nEvents = files.len
  result.nFadcEvents = fadcFiles.len

proc extractRunNumber*(runFolder: string): int =
  ## given a valid TOS run folder, extract its run Number
  let (_, runNumber, _, _) = isTosRunFolder(runFolder)
  result = runNumber

proc extractRunFolderKind*(runFolder: string): RunFolderKind =
  ## given a valid TOS run folder, extract its run Number
  let (_, _, rfKind, _) = isTosRunFolder(runFolder)
  result = rfKind

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
