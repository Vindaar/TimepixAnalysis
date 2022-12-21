## This program reads a set of log files from CAST (slow control and tracking)
## and appends the data to a HDF5 file, which was created using the `raw_data_manipulation`
## and `reconstruction` programs if desired.
##
## Alternatively, information about the magnet history can be computed.
## CAST log file reader. Reads log files (Slow Control or Tracking Logs)
## and outputs information about the files.
##
## Note: This tool can be compiled with the `-d:pure` flag to avoid pulling
## in the HDF5 dependency.
##
## If this tool is run with a single slow control and / or tracking log
## file (i.e. not a folder containing several log files), the tool only
## prints information from that file. In case of a tracking log this is
## the start and end date of a contained tracking (if any). For the slow
## control log currently no useful output is produced.
##
## If a path to tracking log files is given, and a H5 file is given via
## the `--h5out` option, which  InGrid data for the (some of) the times
## covered by the log files being read. The start and end times of a run
## will be added to the attributes to the TOS runs contained within.
##
## If only a single H5 file is given together with the `--delete` argument
## all tracking information will be deleted from the attributes in the file.


import std / [times, strutils, strformat, os, algorithm,
              sequtils, strutils, parseutils, hashes, options, streams]
import ggplotnim, cligen, flatty, untar
import unchained except Time

when not defined(pure):
  import nimhdf5

import ingrid/tos_helpers

{.deadCodeElim: on.}

type
  TrackingKind* = enum
    rkTracking, # if tracking took place in this log
    rkNoTracking # if no tracking took place in this log

  LogFileKind* = enum
    lkSlowControl, lkTracking
  # object storing all slow control log data
  Temperatures = object
    amb: float
    iron: float
    mrb: float
    mfb: float
    ext: float
    vent: float

  SlowControlLog = object
    filename: string
    # store date of log file as string or Time?
    date: Time
    # time of day, stored as time duration starting from date
    # each log starts at 12am
    timestamps: seq[int] #seq[Duration] Unix timestamp
    # pressures relevant for InGrid
    pmm, p3, p3_ba: seq[float]
    # what's called `I_magnet` which is ``not`` the current, but the
    # magnetic field in tesla!
    B_magnet: seq[float]
    # environmental ambient temperature
    temps: seq[Temperatures]
    humidity: seq[float]
    # pointing positions as angles
    h_angle, v_angle: seq[float]
    # and encoder values
    h_encoder, v_encoder: seq[uint16]
    # MM gas alarm (if on, gas supply closed!)
    mm_gas: seq[bool]
    badLineCount: int # number of broken lines in this file
  # object storing the tracking log data
  TrackingLog* = object
    # date the tracking log covers
    name*: string
    date*: Time
    case kind*: TrackingKind
    of rkNoTracking: discard
    of rkTracking:
      tracking_start*, tracking_stop*: Time
    # and the line based data
    timestamps*: seq[int]
    #speedH: seq[float]
    #speedV: seq[float]
    isMoving*: seq[bool]
    isTracking*: seq[bool]
    magB*: seq[float] # magnetic field
                      # apparently magB is not the real value in the tracking logs anymore!
    badLineCount: int

  ## A single version schema
  Version = distinct string
  VersionSchema = object
    version: Version         ## version number of the given schema
    fTab: Table[string, int] ## field table mapping field name to column index

  ## A file describing multiple version schemas
  VersionSchemaFile = object
    versions: Table[Version, VersionSchema] ## maps version names to a schema

proc hash(v: Version): Hash {.borrow.}
proc `==`(v1, v2: Version): bool {.borrow.}

proc len(v: VersionSchema): int = v.fTab.len

proc `[]=`(v: var VersionSchemaFile, key: Version, schema: VersionSchema) =
  v.versions[key] = schema

proc `[]`(v: VersionSchemaFile, key: Version): VersionSchema =
  v.versions[key]

proc `[]`(v: VersionSchema, key: string): Option[int] =
  ## Either returns the index of the given `key` if it exists or `None`
  if key in v.fTab:
    some(v.fTab[key])
  else:
    none[int]()

proc `$`(v: VersionSchema): string =
  result = "{version: " & v.version.string
  result.add ", fTab: " & $v.ftab & "}"

proc `$`(v: VersionSchemaFile): string =
  result = "{versions: "
  for key in keys(v.versions):
    result.add key.string & ", "
  result.add "}"

proc contains(v: VersionSchemaFile, key: Version): bool = key in v.versions

proc parseField(line: string, tok: var string, idx: var int) =
  inc idx, line.parseUntil(tok, until = "##", start = idx)
  inc idx, 2

proc parseVersionSchemaFile(path: string): VersionSchemaFile =
  for schemaProto in readFile(path).splitLines:
    # first field is version. Parse it separate
    var tok: string
    var idx = 0
    var col = 0
    parseField(schemaProto, tok, idx)
    ## XXX: there is a refc bug that causes the `ver` variable to point to the same
    ## string as `tok`!
    when not defined(gcDestructors):
      let ver = Version(deepCopy tok)
    else:
      let ver = Version(tok)
    var tab = initTable[string, int]()
    # now parse all fields
    while idx < schemaProto.len:
      parseField(schemaProto, tok, idx)
      tab[tok] = col
      inc col
    result[ver] = VersionSchema(version: ver, fTab: tab)

proc getVersionSchema(schemaFile: VersionSchemaFile, fname: string): VersionSchema =
  let (dir, fn, ext) = splitFile(fname)
  if fn.startsWith("SCDV"):
    # should be a schema version file
    var idx = 4
    var verBuf = ""
    inc idx, parseUntil(fn, verBuf, until = '_', idx)
    if idx != 8:
      # failed to parse
      raise newException(IOError, "Could not parse version number from filename: " & $fn)
    elif Version(verBuf) notin schemaFile:
      # bad version
      raise newException(IOError, "Could not find version " & $verBuf & " in table of version numbers.")
    else:
      result = schemaFile[Version(verBuf)]
  else:
    # bad file
    raise newException(IOError, "Input file " & $fn & " does not follow convention.")

proc newSlowControlLog(name: string): SlowControlLog =
  result.filename = name
  result.date = fromUnix(0)

proc hash(x: TrackingLog): Hash =
  ## proc to define a hash value for a tracking log object
  ## needed to store it in a hash Table
  var h: Hash = 0
  # can use hashing of Time
  h = h !& hash(x.name)
  h = h !& hash(x.date.toUnix)
  case x.kind
  of rkTracking:
    h = h !& hash(x.tracking_start.toUnix)
    h = h !& hash(x.tracking_stop.toUnix)
  else: discard
  result = !$h

proc `==`(x, y: TrackingLog): bool =
  ## equality of tracking logs, since variant types do not provide an equality
  ## operator by itself
  # don't compare names, as file names can be different
  result = x.date == y.date
  result = result and (x.kind == y.kind)
  if x.kind == y.kind and x.kind == rkTracking:
    # if tracking true, compare start and end times
    result = result and (x.tracking_start == y.tracking_start)
    result = result and (x.tracking_stop == y.tracking_stop)

proc parseTime(time_str: string): Duration =
  ## proc to parse the time from a string hh:mm:ss returned as a
  ## time duration from midnight
  let t = time_str.split(":")
  # create the time interval based on the date given
  result = initDuration(hours = parseInt(t[0]), minutes = parseInt(t[1]), seconds = parseInt(t[2]))

proc is_magnet_moving(h_me, v_me: (int, int)): bool {.inline.} =
  ## determine whether magnet is moving horizontally and (!) vertically
  ## by checking previous and current motor encoder values
  result = if h_me[0] != h_me[1] and v_me[0] != v_me[1]: true else: false

proc sortTrackingLogs(tr_logs: seq[TrackingLog], order = SortOrder.Ascending): seq[TrackingLog] =
  ## proc to sort a sequence of tracking logs by date
  var
    gt = 0
    lt = 0
  # check the sorting order. Depending on case, we define the variables greater and
  # lesser than to either be 1 or 0
  case order
  of Ascending:
    gt = 1
    lt = 0
  of Descending:
    gt = 0
    lt = 1
  result = sorted(tr_logs) do (r, t: TrackingLog) -> int:
    let c = times.`<`(r.date, t.date)
    result = if c == true: lt else: gt

proc filterTrackingLogs(logs: seq[TrackingLog], startDate, endDate: Time): seq[TrackingLog] =
  echo "Filtering tracking logs to date: ", startDate, " and ", endDate
  result = logs.filterIt(it.kind == rkNoTracking or (it.tracking_start >= startDate and it.tracking_stop < endDate))

proc sortAndFilter(logs: seq[TrackingLog], startTime, endTime: string): seq[TrackingLog] =
  let startDate = if startTime.len > 0: parseTime(startTime, "YYYY/MM/dd", utc())
                  else: fromUnix(0)
  let endDate = if endTime.len > 0: parseTime(endTime, "YYYY/MM/dd", utc())
                else: fromUnix(int32.high)
  result = sortTrackingLogs(logs)
  # filter to stard / end dates
  result = result.filterTrackingLogs(startDate, endDate)

proc print_tracking_logs(logs: seq[TrackingLog], print_type: TrackingKind,
                         startTime, endTime = "") =
  ## proc to pretty print a seq of TrackingLogs using org date format
  ## inputs:
  ##    logs: seq[TrackingLog] = seq of (potentially mixed) tracking logs to be
  ##      printed using org mode date
  ##    print_type: TrackingKind = the kind of logs to be printed. Either tracking
  ##      no tracking
  let s_logs = logs.sortAndFilter(startTime, endTime)
  for log in s_logs:
    case log.kind
    of rkTracking:
      echo "<$#>    <$#>" % [formatAsOrgDate(log.tracking_start), formatAsOrgDate(log.tracking_stop)]
    of rkNoTracking:
      echo "<$#>" % formatAsOrgDate(log.date)

  case print_type
  of rkTracking:
    echo &"There are {s_logs.len} solar trackings found in the log file directory"
    var trackTime: Duration
    for log in s_logs:
      if log.kind == rkTracking:
        doAssert log.tracking_stop > log.tracking_start, "the fuck " & $log
        trackTime += (log.tracking_stop - log.tracking_start)
    echo &"The total time of all trackings: {trackTime.inHours()} h (exact: {tracktime})"
  of rkNoTracking:
    echo &"There are {s_logs.len} runs without solar tracking found in the log file directory"

  # compute total time magnet was on
  when true:
    ## NOTE: this is wrong for modern tracking logs!
    var sumB: Second
    for log in s_logs:
      #echo "log ", log.date, " has Bs: ", log.magB.filterIt(it > 0.0), " in name ", log.name
      #if log.date.utc().year > 2015: quit()
      for i in 1 ..< log.timestamps.len:
        let diff = log.timestamps[i].s - log.timestamps[i-1].s
        #echo "log.times ", log.timestamps[i], " to ", log.timestamps[i-1], " is ", diff
        if log.magB[i] > 8.0:
          sumB = sumB + diff
    echo &"Total time the magnet was on (> 1 T): {sumB.to(Hour)} h"

proc print_slow_control_logs(logs: seq[SlowControlLog], magnetField = 8.0) =
  ## proc to pretty print useful information about the Slow Control data
  ## Mainly related to magnet activity.
  # compute total time magnet was on, here we do it based on difference between last and this
  # timestamp. Then add that diff if magnet is on now.
  var sumB: Second
  var Bvals = newSeq[float]()
  var Tdiffs = newSeq[Second]()
  let sortedLogs = logs.filterIt(it.timestamps.len > 0).sortedByIt(it.date)
  var oldTime = sortedLogs[0].timestamps[0]
  for log in sortedLogs:
    if log.B_magnet.len == 0: continue # no magnet data, nothing to learn
    for i in 1 ..< log.timestamps.len:
      let curTime = log.timestamps[i]
      let diff = (curTime - oldTime).Second
      Tdiffs.add diff
      if log.B_magnet[i] > magnetField:
        sumB = sumB + diff
      if log.B_magnet[i] > 0.0:
        Bvals.add log.B_magnet[i]
      #elif log.B_magnet[i] > 0.0 and log.B_magnet[i] < magnetField:

      oldTime = curTime

  let df = toDf({"B" : Bvals})
  ggplot(df.filter(f{`B` <= magnetField}), aes("B")) + geom_histogram(bins = 100) + ggsave(&"/tmp/B_field_larger_0_smaller_{magnetField}.pdf")
  ggplot(df.filter(f{float -> bool: `B` >= magnetField}), aes("B")) + geom_histogram(bins = 100) + ggsave(&"/tmp/B_field_larger_{magnetField}.pdf")
  ggplot(df, aes("B")) + geom_histogram(bins = 100) + ggsave(&"/tmp/B_field_larger_0.pdf")
  let dfT = toDf({"Tdiff" : Tdiffs.mapIt(it.float)})
    .filter(f{`Tdiff` > 0.1 and `Tdiff` < 500.0})
  echo dfT
  ggplot(dfT, aes("Tdiff")) + geom_histogram(bins = 100) +
    scale_y_continuous() + scale_x_continuous() +
    ggsave("/tmp/T_diffs.pdf")
  echo &"Total time the magnet was on (> {magnetField} T): {sumB.to(Hour)} h"

when not defined(pure):
  proc print_mapped_tracking_logs(logs: seq[TrackingLog], map: Table[TrackingLog, int],
                                  startTime, endTime: string) =
    let logs = sortAndFilter(logs, startTime, endTime)
    # 1. first print the *mapped* tracking logs:
    echo "========== Logs mapped to trackings in the output file: =========="
    for log in logs.filterIt(it in map): ## filter and get the run to output in sorted order!
      let run = map[log]
      doAssert log.kind == rkTracking
      echo "<$#>    <$#>  for run: $#" % [formatAsOrgDate(log.tracking_start),
                                          formatAsOrgDate(log.tracking_stop),
                                          $run]
    # 2. print those runs *not* mapped to anything:
    echo "==================================================================\n"
    echo "========== Logs *not* mapped to a run ============================"
    for log in logs.filterIt(it notin map):
      doAssert log.kind == rkTracking
      echo "<$#>    <$#>" % [formatAsOrgDate(log.tracking_start),
                             formatAsOrgDate(log.tracking_stop)]
    echo "=================================================================="

  proc map_log_to_run(logs: seq[TrackingLog], h5file: string,
                      startTime, endTime: string): Table[TrackingLog, int] =
    ## maps a sequence of tracking logs to run numbers given in a H5 file
    ## inputs:
    ##   logs: seq[TrackingLog] = tracking logs to map
    ##   h5file: string = path to H5 file containing runs
    ## outputs:
    ##   Table[TrackingLog, int] = table mapping each tracking log to the
    ##     run number in which the tracking is contained
    ## throws:
    ##   HDF5LibraryError = if a call to the H5 library fails

    result = initTable[TrackingLog, int]()

    var h5f = H5open(h5file, "r")

    var run_times = initTable[int, (int, int)]()
    for grp in items(h5f, "/reconstruction", depth = 1):
      if "/reconstruction/run" in grp.name:
        # given a run, check its start and end time
        var tstamp = h5f[(grp.name & "/timestamp").dset_str]
        let
          t_start = tstamp[0, int]
          t_stop  = tstamp[tstamp.high, int]
        # for now need mutable attributes object to access
        var attrs = grp.attrs
        let run_num = attrs["runNumber", int]
        # add start and end time to table of run times
        run_times[run_num] = (t_start, t_stop)

    # now iterate over tracking logs and for each log determine in which run
    # it fits
    for log in logs:
      case log.kind
      of rkTracking:
        # now given start time, look for run for which it fits. Filter elements
        # where run covers tracking start - end
        proc filterRuns(runs: seq[int], log: TrackingLog): int =
          result = -1
          for r in runs:
            let
              t_start = int(log.tracking_start.toUnix)
              t_stop  = int(log.tracking_stop.toUnix)
            if run_times[r][0] < t_start and run_times[r][1] > t_stop:
              doAssert result < 0, "There are multiple trackings? Something is wrong. " & $result & " and " & $r
              result = r
        let run_num = filterRuns(toSeq(run_times.keys), log)
        if run_num > 0:
          result[log] = run_num
      else: discard
    discard h5f.close()

    # in addition print the information about all logs we found. Both the
    # logs that _are_ mapped (and which runs they are mapped to), but also those
    # that are _not_ mapped, i.e. runs that we missed!
    print_mapped_tracking_logs(logs, result, startTime, endTime)

  proc deleteTrackingAttributes(h5file: string) =
    ## proc to delete all tracking related attributes in a H5 file
    withH5(h5file, "rw"):
      let baseName = if "/runs" in h5f: rawDataBase()
                     elif "/reconstruction" in h5f: recoBase()
                     else: raise newException(IOError, "Invalid input file " &
          $h5file & ". It contains neither `runs` nor `reconstruction` group!")
      for (run, grpName) in runs(h5f, baseName):
        # get number of tracking related attributes
        # for now need mutable attributes object to access
        let grp = h5f[grpName.grp_str]
        var attrs = grp.attrs
        var num_tr = 0
        var mgrp = grp
        try:
          num_tr = attrs["num_trackings", int]
        except KeyError:
          # no trackings, continue
          continue
        for i in 0 ..< num_tr:
          var deleted = false
          let
            attr_start = "tracking_start_$#" % $i
            attr_stop  = "tracking_stop_$#" % $i
            attr_num   = "num_trackings"
          deleted = mgrp.deleteAttribute(attr_start)
          deleted = mgrp.deleteAttribute(attr_stop)
          deleted = mgrp.deleteAttribute(attr_num)
          if deleted == false:
            echo "Could not delete one of " &
              "$#, $# or $# in group $#" % [$attr_start, $attr_stop, $attr_num, $mgrp.name]

  proc write_tracking_h5(trck_tab: Table[TrackingLog, int], h5file: string) =
    ## proc to write the mapping of tracking logs to run numbers to the appropriate
    ## groups of the H5 file
    ## inputs:
    ##   trck_tab: Table[TrackingLog, int] = table mapping tracking logs to run numbers
    ##   h5file: string = h5 file in which to write the information
    ## throws:
    ##   HDF5LibraryError = if a call to the H5 library fails
    var runTab = initCountTable[int]()

    withH5(h5file, "rw"):
      for log, run_number in trck_tab:
        # increment run table of current run number
        inc(runTab, runNumber)
        let num = runTab[runNumber]

        let baseName = if "/runs" in h5f: rawDataBase()
                       elif "/reconstruction" in h5f: recoBase()
                       else: raise newException(IOError, "Invalid input file " &
            $h5file & ". It contains neither `runs` nor `reconstruction` group!")
        # take run number, build group name and add attributes
        let runName = baseName & $run_number
        var runGrp = h5f[runName.grp_str]
        # check if there is a number of trackings already

        if "num_trackings" notin runGrp.attrs:
          runGrp.attrs["num_trackings"] = 1
        else:
          runGrp.attrs["num_trackings"] = num

        # store trackings zero indexed
        let numIdx = num - 1

        # add tracking. Trackings will be zero indexed
        runGrp.attrs["tracking_start_$#" % $numIdx] = $log.tracking_start
        runGrp.attrs["tracking_stop_$#" % $numIdx] = $log.tracking_stop

proc tryParseTime(s: string, formatStrings: openArray[string]): Time =
  ## Attempts to parse the time string `s` using any of the `formatStrings`
  for fmt in formatStrings:
    try:
      return toTime(parse(s.strip(chars = {'\0'}), fmt))
    except TimeParseError:
      continue
  raise newException(TimeParseError, "Parsing of string: " & $s &
    " failed with all format strings: " & $formatStrings)

proc copyMemStrBuf(s: var string, buf: var string, len: int) =
  s.setLen(len) # set to `len` to fit buffer exactly
  copyMem(s[0].addr, buf[0].addr, len * sizeof(char))

proc splitMinTwo(lineData: var seq[string],
                 idx: var int, # global file index
                 buf: var string, # buffer for the current column field
                 validData: var bool, # whether we have arrived at valid data
                 isHeaderLine: var bool,
                 data: string, # the full file data
                 expected: int,
                 lineCnt: int): int =
  ## line data is the sequence of fields that will be filled
  var bufIdx = 0
  var spaceCount = 0
  var lastWasSpace = false
  var count = 0
  var parseWithoutStore = false
  while idx < data.len:
    case data[idx]
    of '\0':
      discard # ignore null characters
    of '\n', '\r':
      # end of line
      break
    of ' ':
      lastWasSpace = true
      inc spaceCount
    of '\t':
      lastWasSpace = true
      inc spaceCount, 4 # count tab as multiple spaces so that single tab is enough to parse
    of '#':
      isHeaderLine = true
    else:
      if lastWasSpace and spaceCount > 1: # need more than 1 space:
        # add last element
        if count >= lineData.len:
          parseWithoutStore = true
        elif bufIdx > 0: # only count if we actually parsed something (otherwise line starting ' ' will count)
          copyMemStrBuf(lineData[count], buf, bufIdx)
          if not validData and count < 2 and lineData[count].startsWith("DATE"):
            isHeaderLine = true
            validData = true # should be fine from here on?
            parseWithoutStore = true

        if bufIdx > 0:
          inc count
        bufIdx = 0
      buf[bufIdx] = data[idx]
      inc bufIdx
      lastWasSpace = false

    inc idx
  if count >= lineData.len:
    discard
  else:
    copyMemStrBuf(lineData[count], buf, bufIdx)
    inc count
    #lineData[count] = buf[0 ..< bufIdx]
  if not isHeaderLine and count != expected and lineCnt != 0:
    echo "parsed : ", lineData.mapIt(it.strip(chars = {' ', '\0'}))
    echo "in line: ", lineCnt
  elif not isHeaderLine:
    doAssert count == expected, " Count " & $count & " but expected " & $expected
    validData = true
  result = count

proc parse_sc_logfile(content: string,
                      schema: VersionSchema,
                      filename = ""
                     ): SlowControlLog =
  ## proc to read a slow control log file
  ## inputs:
  ##    filename: string = the filename of the log file
  ##    schemaFile: The map of all known schemas. Used to parse the correct columns.
  ## outputs:
  ##    SlowControlLog = an object storing the data from a slow control
  ##    log file.
  # define the indices for each data column. See the ./sc_log_understand.org file
  # look up the correct fields
  let
    date_i = schema["DATE"]
    time_i = schema["TIME"]
    # for these: taking into account offset of 1
    pmm_i = schema["P-MM"]
    p3_i = schema["P-3"]
    p3ba_i = schema["P3_BA"]
    Imag_i = schema["I_magnet"] # is actually the magnetic field in ``T``
    tamb_i = schema["Tenv_Amb"]
    tiron_i = schema["Tenv_Iron"]
    tmrb_i = schema["Tenv_MRB"]
    tmfb_i = schema["Tenv_MFB"]
    text_i = schema["TEnv_Ext"]
    tvent_i = schema["TEnv_Vent"]
    humid_i = schema["Humidity"]
    mm_gas_i = schema["MM GAS"]
    hang_i = schema["Horiz_Angle"]
    vang_i = schema["Verti_Angle"]
    hme_i = schema["Horiz_ME"]
    vme_i = schema["Verti_ME"]

  result = newSlowControlLog(name = filename)
  var parsedDate = false
  var lineCnt = 0
  var badLineCount = 0

  #for line in content.splitLines():
  var lineStart = 0
  var lineBuf = newString(5000)
  var d = newSeqWith(schema.len, newString(100))
  var i = 0
  #var fVal: float
  #var iVal: int16
  var isHeader = false
  var validData = false
  while i < content.len:
    ## accumulate until line break
    case content[i]
    of '\r', '\n':
      # process last line
      lineStart = i+1
      inc i
    else:
      # accumulate
      # parse into seq
      let count = splitMinTwo(d, i, lineBuf,
                              validData, isHeader,
                              content,
                              schema.len, lineCnt)
      # `d` now contains correct data
      # skip the first line (header)
      if isHeader:
        inc lineCnt
        isHeader = false
        continue

      # else add the data to the SlowControlLog object
      if count == 0: # parsing of line unsuccessful
        echo "bad line count+1 ", badLineCount, " in line : ", lineCnt
        inc lineCnt
        inc badLineCount
        continue
      if not parsedDate and not d[0].strip.startsWith("DATE"):
        # set date based on first row
        try:
          result.date = tryParseTime(d[date_i.get], ["MM/dd/yyyy",
                                                     "M/d/yyyy",
                                                     "dd-MMM-yy"]) # DATE always exists, safe to get
          parsedDate = true
          if result.date.utc().year == 1970:
            echo "fuck this file: ", filename
            quit()
        except TimeParseError:
          ## Parsing *can* fail, because some files contain line breaks in the header...
          if lineCnt > 5: # don't output for first few lines, as they often have more header etc
            echo "Failed to parse the following date line: ", d[date_i.get]
            #echo content
            echo "Full line: ", d
            echo "In file ", filename, " Trying again next line"
          inc lineCnt
          continue
          #quit()
      ## Different versions of the files have different fields set. Probably need to
      ## skip files that don't have what we need?

      if count != schema.len:
        echo "INVALID FILE: ", filename, ". ", count, " columns found, but ", schema.len, " expected!"
        echo "Offending line: ", content[lineStart .. i]
        echo "split to: ", d.mapIt(it.strip(chars = {' ', '\0'}))
        inc lineCnt
        inc badLineCount
        continue

      # now parse both date and time & add as a unix timestamp
      let date = tryParseTime(d[date_i.get], ["MM/dd/yyyy",
                                              "M/d/yyyy",
                                              "dd-MMM-yy"]) # DATE always exists, safe to get
      if abs(date - result.date) > initDuration(days = 1): # one day is allowed to account for midnight to next day logs
        # If there's a line with a different date suddenly, skip it. We cannot safely handle such cases!
        echo "Line detected with mismatching date: ", date, " vs ", result.date
        echo "Skipping line. File: ", filename
        inc badLineCount
        inc lineCnt
        continue

      let time = parseTime(d[time_i.get].strip(chars = {'\0'}))
      result.timestamps.add (date + time).toUnix.int

      template addIfAvailable(arg, body: untyped): untyped =
        if arg.isSome:
          let idx {.inject.} = arg.unsafeGet
          body
      template getIf(arg: untyped): untyped =
        if arg.isSome:
          let idx = arg.unsafeGet
          parseFloat d[idx]
        else:
          0.0
      addIfAvailable(pmm_i, result.pmm.add parseFloat(d[idx]))
      addIfAvailable(p3_i, result.p3.add parseFloat(d[idx]))
      addIfAvailable(p3_ba_i, result.p3_ba.add parseFloat(d[idx]))
      addIfAvailable(Imag_i, result.B_magnet.add parseFloat(d[idx]))
      addIfAvailable(humid_i, result.humidity.add parseFloat(d[idx]))
      addIfAvailable(hang_i, result.h_angle.add parseFloat(d[idx]))
      addIfAvailable(vang_i, result.v_angle.add parseFloat(d[idx]))
      addIfAvailable(hme_i, result.h_encoder.add uint16(parseInt(d[idx])))
      addIfAvailable(vme_i, result.v_encoder.add uint16(parseInt(d[idx])))
      addIfAvailable(mm_gas_i):
        let mm_gas_b = d[idx][0] == '0'
        result.mm_gas.add mm_gas_b
      # finally try to add all temperatures
      let temps = Temperatures(amb: getIf(tamb_i),
                               iron: getIf(tiron_i),
                               mrb: getIf(tmrb_i),
                               mfb: getIf(tmfb_i),
                               ext: getIf(text_i),
                               vent: getIf(tvent_i))
      result.temps.add temps

      inc lineCnt
  result.badLineCount = badLineCount

proc read_sc_logfile(filename: string,
                               schemaFile: VersionSchemaFile): SlowControlLog =
  # get schema
  let schema = schemaFile.getVersionSchema(filename)
  result = parse_sc_logfile(readFile(filename), schema, filename)

proc splitWhitespaceBuf(data: var seq[string],
                        buf: var string,
                        line: string): int =
  ## splits the given `line` at spaces and inserts the elements into data fields
  var idx = 0
  var count = 0
  var bufIdx = 0
  var lastWasSpace = false
  while idx < line.len:
    case line[idx]
    of '\0', '\n', '\r':
      # done, break
      break
    of ' ', '\t':
      lastWasSpace = true
    else:
      if lastWasSpace and bufIdx > 0:
        copyMemStrBuf(data[count], buf, bufIdx)
        inc count
        bufIdx = 0
      lastWasSpace = false
      buf[bufIdx] = line[idx]
      inc bufIdx
    inc idx
  result = count

var seenHeaders: set[uint16]
proc parse_tracking_logfile*(content: string, filename: string): TrackingLog =
  ## reads a tracking log file and returns an object, which stores
  ## the relevant data
  const
    tracking_i = 0
    date_i = 6
    time_i = 7
    h_me = 9
    v_me = 10
    magB_i = 22

  var
    lineCnt = 0
    h_me_p = 0
    v_me_p = 0
    tracking_p = false
    tracking_start = fromUnix(0)
    tracking_stop = fromUnix(0)
    # helper bool, needed because some tracking logs start
    # before midnight
    date_set = false
    badLineCount = 0

  var stream = newStringStream(content)
  var line = newString(5000)
  var d = newSeqWith(100, newString(50)) # just enough space for everything
  var buf = newString(100)
  while stream.readLine(line):
    if lineCnt < 2:
      let spl = line.split(",")
      if spl.len.uint16 notin seenHeaders:
        echo "header: ", line
        for i, s in spl:
          echo "index ", i, " = " , s
        seenHeaders.incl spl.len.uint16

      # skip header
      inc lineCnt
      continue
    if line.len == 0:
      break
    let count = splitWhitespaceBuf(d, buf, line)
    #echo d.len
    ## Order of columns is ``kept`` in tracking logs. That means we only need to remove those
    ## files that don't have the information we want, namely the "bare bones" files from the
    ## commissioning
    var date: Time
    if lineCnt > 1 and not date_set and count >= 22:
      date = toTime(parse(d[date_i], "MM/dd/yy"))
      # check whether this date is at the end of day from the log file or already past midnight
      if date.utc().hour > 23:
        continue
      else:
        # set the date of the tracking log
        result.date = date
        date_set = true
    elif count < 22:
      inc badLineCount
      continue
    else:
      # now parse date to make sure we have it
      try:
        date = toTime(parse(d[date_i], "MM/dd/yy"))
      except:
        echo "d was ", d, " with count ", count
        raise

    if abs(date - result.date) > initDuration(days = 7):# one week is allowed, because old tracking logs were sometimes multiple days
      # If there's a line with a different date suddenly, skip it. We cannot safely handle such cases!
      echo "Line detected with mismatching date: ", date, " vs ", result.date
      echo "Skipping line. File: ", filename
      inc badLineCount
      continue
    let
      h_me = int(parseFloat(d[h_me]))
      v_me = int(parseFloat(d[v_me]))
      daytime = parseTime(d[time_i]) # time on the current day since midnight
      tracking = int(parseFloat(d[tracking_i])) == 1
      magB = parseFloat(d[magB_i])
    # determine magnet movement and set old encoder values
    let move = is_magnet_moving((h_me, h_me_p), (v_me, v_me_p))
    h_me_p = h_me
    v_me_p = v_me
    if not tracking_p and tracking:
      tracking_start = date + daytime
      tracking_p = not tracking_p
    elif tracking_p and not tracking:
      tracking_stop = date + daytime
      tracking_p = not tracking_p

    # append seq data
    result.timestamps.add (date + daytime).toUnix.int
    result.isMoving.add move
    result.isTracking.add tracking
    result.magB.add magB
    inc lineCnt

  # now set the tracking variant object depending on whether tracking took place
  # or not
  if tracking_start == tracking_stop:
    result = TrackingLog(kind: rkNoTracking,
                         date: result.date,
                         timestamps: result.timestamps,
                         isMoving: result.isMoving,
                         isTracking: result.isTracking,
                         magB: result.magB)
  elif tracking_start > tracking_stop: # likely file _ends_ while tracking (broken old file)
    doAssert tracking_stop == fromUnix(0)
    result.tracking_start = tracking_start
    result.tracking_stop = result.timestamps[^1].fromUnix
  else:
    result.tracking_start = tracking_start
    result.tracking_stop = tracking_stop
  result.name = filename
  result.badLineCount = badLineCount
  if result.tracking_stop < parseTime("1975/01/01", "YYYY/MM/dd", utc()):
    quit()

proc read_tracking_logfile*(filename: string): TrackingLog =
  result = parse_tracking_logfile(readFile(filename), filename)

proc split_tracking_logs(logs: seq[TrackingLog]): (seq[TrackingLog], seq[TrackingLog]) =
  ## given a sequence of sorted tracking logs, splits them by `kind`, i.e.
  ## tracking or no tracking
  ## inputs:
  ##   logs: seq[TrackingLog] = seq of mixed TrackingLog kinds to be split by kind
  ## outputs:
  ##   (seq[TrackingLog], seq[TrackingLog]) =
  ##     0: TrackingLog.kind == rkTracking
  ##     1: TrackingLog.kind == rkNoTracking
  ##     both seqs are sorted by date
  # init result
  result[0] = @[]
  result[1] = @[]
  for log in logs:
    case log.kind
    of rkTracking: result[0].add log
    of rkNoTracking: result[1].add log

proc read_tracking_log_folder(log_folder: string): seq[TrackingLog] =
  ## reads all log files from `log_folder` and returns a tuple of sorted `TrackingLog`
  ## objects, one set for logs w/ tracking, others without
  ## inputs:
  ##   log_folder: string = folder in which log files are stored
  ## outputs:
  ##   (seq[TrackingLog], seq[TrackingLog]) =
  ##     0: TrackingLog.kind == rkTracking
  ##     1: TrackingLog.kind == rkNoTracking
  ##     both seqs are sorted by date
  var tracking_logs: seq[TrackingLog] = @[]

  # for each file in log_folder we perform a check of which run corresponds
  # to the date of this tracking log
  for log_p in walkDir(log_folder):
    try:
      case log_p.kind
      of pcFile:
        let log = log_p.path
        echo "Reading ", log
        # check whether actual log file (extension fits)
        let (dir, fn, ext) = splitFile(log)
        if ext == ".log":
          tracking_logs.add read_tracking_logfile(log)
        else:
          # skipping files other than log files
          continue
      else: discard
    except:
      continue # drop bad file

  result = sortTrackingLogs(tracking_logs)

proc toDf(sc: SlowControlLog, enforceSameFields = false): DataFrame =
  if enforceSameFields:
    doAssert sc.timestamps.len == sc.B_magnet.len or sc.B_magnet.len == 0, "Enforcement of times & magnet data length failed! " &
      "times.len = " & $sc.timestamps.len & ", magnet.len = " & $sc.B_magnet.len
  let B = if sc.B_magnet.len == 0: sc.timestamps.mapIt(NaN) # fill with NaN to have floats
          else: sc.B_magnet
  ## XXX: even 22/12/2022 cannot use `toDf` here due to overload resolution issue causing it to complain
  ## due to different data types, wtf
  var df = seqsToDf({ "Time / s" : sc.timestamps,
                      "B / T" : B })
  df["Date"] = newTensorWith(df.len, sc.date.toUnix) #constantColumn(sc.date.toUnix, df.len)
  result = df

proc plotTemperatures(scLogs: seq[SlowControlLog]) =
  ## Plots the different temperatures recorded in the CAST hall in the given time range
  let df = seqstoDf({ "Time" : scLogs.mapIt(it.timestamps).flatten,
                  "T_amb" : scLogs.mapIt(it.temps.mapIt(it.amb)).flatten,
                  "T_iron" : scLogs.mapIt(it.temps.mapIt(it.iron)).flatten,
                  "T_mrb" : scLogs.mapIt(it.temps.mapIt(it.mrb)).flatten,
                  "T_mfb" : scLogs.mapIt(it.temps.mapIt(it.mfb)).flatten,
                  "T_ext" : scLogs.mapIt(it.temps.mapIt(it.ext)).flatten,
                  "T_vent" : scLogs.mapIt(it.temps.mapIt(it.vent)).flatten,
                  "Humidity" : scLogs.mapIt(it.humidity).flatten })
    .gather(["T_amb", "T_iron", "T_mrb", "T_mfb", "T_ext", "T_vent"], key = "Temperature", value = "TempVal")
    .arrange("Time")
  ggplot(df, aes("Time", "TempVal", color = "Temperature")) +
    geom_line() +
    scale_x_date(name = "Date", isTimestamp = true,
                 dateSpacing = initDuration(weeks = 12),
                 formatString = "yyyy-MM") +
                        #dateSpacing = initDuration(weeks = 26), dateAlgo = dtaAddDuration) +
    ggtitle("Temperatures in CAST hall & outside") +
    ggsave("/tmp/temperatures_cast.png", width = 1200, height = 800)

proc read_sc_log_folder(log_folder: string,
                        schemaFile: VersionSchemaFile,
                        magnetField = 8.0) =
  ## reads all slow contro log files from `log_folder` and (at the moment)
  ## determines the time the magnet was turned on in total during the time
  ## covered in all the log files in the folder.
  # for each file in log_folder we perform a check of which run corresponds
  # to the date of this tracking log
  var dfDir = newDataFrame()

  echo log_folder
  var scLogs: seq[SlowControlLog]
  for log_p in walkDir(log_folder):
    case log_p.kind
    of pcFile:
      let log = log_p.path
      # check whether actual log file (extension fits)
      let (dir, fn, ext) = splitFile(log)
      echo "File ", log, " and ext ", ext
      if ext == ".daq" and not fn.endsWith(".OLD"): # some weird files with `OLD.daq`
        let sc = read_sc_logfile(log, schemaFile)
        dfDir.add toDf(sc)
        scLogs.add sc
      else:
        # skipping files other than log files
        echo "Skipping file ", log_p
        continue
    else: discard

  ## NOTE: we assume that each row in the slow control file covers 1 minute,
  ## because normally the slow control software updates once a minute and writes
  ## at that interval. A more thorough student should make sure the time duration
  ## between each line is roughly 1 minute each!
  if dfDir.len > 0:
    dfDir = dfDir.arrange("Date", SortOrder.Ascending)
      .filter(f{c"Date" != 0})

    ggplot(dfDir, aes("Time / s", "B / T")) +
      geom_point() +
      ggsave("/tmp/B_against_time.png", width = 1920, height = 1080)

    let firstDate = fromUnix(dfDir["Time / s"][0, int])
    let lastDate = fromUnix(dfDir["Time / s"][dfDir.high, int])
    echo firstDate
    echo lastDate
    let nRowsActive = dfDir.filter(f{float: c"B / T" > magnetField}).len
    echo &"Magnet was turned on for (> {magnetField} T): {nRowsActive.float / 60.0} h between " &
         &"{$firstDate} and {$lastDate} based on number of rows (assuming 1 row == 60 s)."

  print_slow_control_logs(scLogs, magnetField)
  # plot the temperature data
  plotTemperatures(scLogs)

proc process_log_folder(folder: string, logKind: LogFileKind,
                        h5file = "",
                        schemaFile: VersionSchemaFile = VersionSchemaFile(),
                        magnetField = 8.0,
                        startTime, endTime = "",
                        dryRun = false) =
  case logKind
  of lkSlowControl:
    read_sc_log_folder(folder, schemaFile = schemaFile, magnetField = magnetField)
  of lkTracking:
    let
      tracking_logs = read_tracking_log_folder(folder)
      (trk, notrk) = split_tracking_logs(tracking_logs)
    echo "No tracking days : "
    print_tracking_logs(notrk, rkNoTracking, startTime = startTime, endTime = endTime)
    echo "Tracking days :"
    print_tracking_logs(trk, rkTracking, startTime = startTime, endTime = endTime)

    when not defined(pure):
      if h5file.len > 0:
        # given the H5 file, create a referential table connecting
        # the trackings with run numbers
        let trackmap = map_log_to_run(trk, h5file, startTime, endTime)

        # given mapping of tracking logs to run numbers, finally
        # add tracking information to H5 file
        if not dryRun:
          write_tracking_h5(trackmap, h5file)

proc toDf(tr: TrackingLog, enforceSameFields = false): DataFrame =
  if enforceSameFields:
    doAssert tr.timestamps.len == tr.magB.len or tr.magB.len == 0, "Enforcement of times & magnet data length failed! " &
      "timestamps.len = " & $tr.timestamps.len & ", magnet.len = " & $tr.magB.len
  let B = if tr.magB.len == 0: tr.timestamps.mapIt(NaN) # fill with NaN to have floats
          else: tr.magB
  ## XXX: even 22/12/2022 cannot use `toDf` here due to overload resolution issue causing it to complain
  ## due to different data types, wtf
  var df = seqsToDf({ "Time / s" : tr.timestamps,
                      "B / T" : B })
  df["Date"] = newTensorWith(df.len, tr.date.toUnix)
  #df["Date"] = constantColumn(tr.date.toUnix, df.len)
  result = df

proc filterBasedOnLogType[T](df: DataFrame, _: typedesc[T]): DataFrame =
  when T is TrackingLog:
    # extract data up to 2010, because from 2010 on the magnet data is found
    # in the slow control as well
    result = df.filter(f{int: idx("Time / s") < dateTime(year = 2010,
                                                         month = mJan,
                                                         monthday = 1,
                                                         zone = utc()).toTime.toUnix})
  else:
    # filter to larger 2010. Not strictly needed, as the data isn't there in the first
    # place, but well.
    result = df.filter(f{int: idx("Time / s") >= dateTime(year = 2010,
                                                          month = mJan,
                                                          monthday = 1,
                                                          zone = utc()).toTime.toUnix})

proc filterInvalidTimestamps(df: DataFrame): DataFrame =
  let times = df["Time / s", int]
  var last = times[0]
  var idxToKeep = newSeqOfCap[int](df.len)
  for i in 0 ..< times.len:
    let time = times[i]
    if time >= last:
      idxToKeep.add i
      last = time
    # else we drop the line. Only update `last` if it was actually bigger. Otherwise we start
    # counting 2 lines after a time jump again!
  result = df.filterToIdx(idxToKeep)
  # sanity check that indeed the data is now sorted already
  #when true:
  #  doAssert result == result.arrange("Time / s", SortOrder.Ascending)

proc extractCycles[T](df: DataFrame, magnetField: float,
                      _: typedesc[T]): DataFrame =
  ## Extracts the number of magnet cycles found in the data as well as the
  ## magnet on time for the given type of logfiles.
  let df = df.filterBasedOnLogType(T) # filter based on log type
    .arrange("Date", SortOrder.Ascending) # sort by date
    .filterInvalidTimestamps()
  let dates = df["Time / s", int]
  let Bs = df["B / T", float]
  var lastAbove = false
  var magnetStart = 0
  const minTimeCurrent = 10 * 60 # 10 min at least of magnet on
  var cycles = 0
  var magnetCycleStarts = initHashSet[int]() # starting indices to check if we've already increased cycle counter
  var magnetCycleIdxs = newSeq[(int, int)]() # indices of magnet start ⇒ stop
  var magCycleStart = 0
  var sumTime = 0.Second
  var cumTime = newSeq[int]()
  var cycleLength = newSeq[int]()
  for i in 1 ..< Bs.len:
    # Logic dealing with sum of magnet on time
    if lastAbove and Bs[i] > magnetField:
      sumTime += abs(dates[i] - dates[i-1]).Second # safe, we start at idx == 1
    # Logic dealing with magnet cycles
    if not lastAbove and Bs[i] > magnetField:
      # passed threshold
      lastAbove = true
      magnetStart = dates[i]
      magCycleStart = i
      magnetCycleStarts.incl i
    elif lastAbove and Bs[i] < magnetField:
      sumTime += abs(dates[i] - dates[i-1]).Second # safe, we start at idx == 1
      # cycle done
      magnetCycleIdxs.add (magCycleStart, i)
      let length = abs(dates[i] - dates[magCycleStart])
      cycleLength.add length
      cumTime.add sumTime.int # add current time as cumulative
      lastAbove = false
  when false:
    for (start, stop) in magnetCycleIdxs:
      if stop - start < 10:
        echo "Start at: ", start, " stop at ", stop
        echo "Relevant df : ", df[start - 5 .. stop + 5]
  # convert cycles to DF
  result = seqsToDf({ "cumulativeTime / s" : cumTime,
                      "cycleLength / s" : cycleLength,
                      "cycleStart (unix)" : magnetCycleIdxs.mapIt(dates[it[0]]),
                      "cycleStop (unix)" : magnetCycleIdxs.mapIt(dates[it[1]]),
                      "Type" : $T })
    .mutate(f{int -> string: "cycleStart" ~ $idx("cycleStart (unix)").fromUnix },
            f{int -> string: "cycleStop" ~ $idx("cycleStop (unix)").fromUnix })
  echo "Number of magnet cycles ", result.len, " from data: ", $T
  echo "Total time the magnet was on: ", sumTime.to(Hour), " from data: ", $T
  echo "Total time the magnet was on: ", sumTime.to(Day), " from data: ", $T

proc print_total_magnet_information(df: DataFrame) =
  echo "Number of magnet cycles ", df.len
  let totalTime = df["cumulativeTime / s", float][df.high].Second
  echo "Total time the magnet was on: ", totalTime.to(Hour)
  echo "Total time the magnet was on: ", totalTime.to(Day)

proc handleAllLogs(all_logs_path: string, schemaFile: VersionSchemaFile,
                   magnetField = 1.0, colorPlots = true,
                   filterOutliers = false, pretty = false) =
  ## computes everything about the magnet from all log files
  let scLogsTar = newTarFile(all_logs_path / "new_slowcontrol_logs.tar.gz")
  let trLogsTar = newTarFile(all_logs_path / "old_logs.tar.gz")

  var scLogs = newSeq[SlowControlLog]()

  const scBinStore = "./slow_control_logs.bin"
  if existsFile(scBinStore):
    scLogs = fromFlatty(readFile(scBinStore), seq[SlowControlLog])
  else:
    for (fi, content) in walk(scLogsTar):
      if fi.filename.extractFilename.startsWith("SCDV"):
        echo "Trying to parse file: ", fi.filename
        var scLog: SlowControlLog
        try:
          let schema = schemaFile.getVersionSchema(fi.filename)
          scLog = parse_sc_logfile(content, schema, fi.filename)
        except Exception as e:
          echo "Exception: ", e.msg
          continue
        scLogs.add scLog

    writeFile(scBinStore, toFlatty(scLogs))

  var dfs = newSeqOfCap[DataFrame](scLogs.len)
  for i, log in scLogs:
    let df = toDf(log, enforceSameFields = true)
    if df.len > 1:
      dfs.add df

  var df = assignStack(dfs)
  # filter out NaN already for this block of data (files that have missing magnet columns
  # fields are replaced by NaN)
  df = df.filter(f{float -> bool: classify(idx("B / T")) != fcNaN })

  var cycleDf = df.extractCycles(magnetField, SlowControlLog)
  df["From"] = newTensorWith(df.len, "SlowControlLogs")
  ## XXX: once switch back to stashed datamancer code:
  #df["From"] = "SlowControlLogs"
  echo df["From"].kind
  echo df
  #print_slow_control_logs(scLogs, 2.0)

  var trLogs = newSeq[TrackingLog]()
  const trBinStore = "./tracking_logs.bin"
  if existsFile(trBinStore):
    trLogs = fromFlatty(readFile(trBinStore), seq[TrackingLog])
  else:
    for (fi, content) in walk(trLogsTar):
      let (dir, name, ext) = fi.filename.splitFile
      if dir == "SlowControl/tracking-log" and ext == ".log": # skip `/old` subdir and other root dirs
        echo fi.filename
        let trLog = parse_tracking_logfile(content, name & ext)

        trLogs.add trLog
    writeFile(trBinStore, toFlatty(trLogs))

  trLogs = trLogs.sortedByIt(it.date)
  dfs = newSeqOfCap[DataFrame](trLogs.len)
  for log in trLogs:
    let dfLoc = toDf(log)
    if dfLoc.len > 1:
      dfs.add dfLoc

  ## add this data as well
  var dfTr = assignStack(dfs)
  # TrackingLog covers time *before* sloc control. Thus add its cumulative time to cycleDf
  let cycleDfTr = dfTr.extractCycles(magnetField, TrackingLog)
  # get cumulative time so far and add to it
  let cumTimeSc = cycleDfTr["cumulativeTime / s", int].max
  cycleDf = cycleDf
    .mutate(f{int: "cumulativeTime / s" ~ idx("cumulativeTime / s") + cumTimeSc})
  cycleDf.add cycleDfTr
  cycleDf = cycleDf.arrange("cycleStart (unix)", SortOrder.Ascending)

  dfTr["From"] = newTensorWith(dfTr.len, "TrackingLogs")
  #dfTr["From"] = "TrackingLogs"
  df.add dfTr

  # sort the data by date to do computations
  df = df.arrange("Time / s", SortOrder.Ascending)
  let breaks = toSeq(2003 .. 2022)
    .mapIt(dateTime(year = it, month = mJan, monthDay = 1, zone = utc())
      .toTime()
      .toUnixFloat
  )

  if filterOutliers:
    # remove outliers from df
    df = df.filter(f{idx("B / T") > -0.1 and idx("B / T") <= 10.0})


  # Reduce the time series data to be more manageable by the likes of Cairo
  df["Index"] = toSeq(0 ..< df.len)
  df = df.filter(f{int -> bool: `Index` mod 10 == 0}) # keep only every 10th row

  # print the summarized information about cycles and time
  print_total_magnet_information(cycleDf)

  # sort and plot
  var pltSeries: GgPlot
  var pltCum: GgPlot
  if colorPlots:
    pltSeries = ggplot(df, aes("Time / s", "B / T", color = "From")) + margin(right = 6)
    pltCum = ggplot(cycleDf, aes("cycleStart (unix)",
                                 f{float: "cumulativeTime / h" ~ to(Second(idx("cumulativeTime / s")), Hour).float},
                                 color = "Type")) +
      margin(right = 6)
  else:
    # filter duplicate data, as it doesn'd add anything valuable without color
    let yr2010 = dateTime(year = 2010, month = mJan, monthday = 1, zone = utc()).toTime.toUnix
    let dfF = df.filter(f{int: (idx("From", string) == "TrackingLogs" and idx("Time / s") < yr2010) or
                               (idx("From", string) == "SlowControlLogs" and idx("Time / s") >= yr2010)})
    pltSeries = ggplot(dfF, aes("Time / s", "B / T"))
    pltCum = ggplot(cycleDf, aes("cycleStart (unix)",
                                 f{float: "cumulativeTime / Days" ~ to(Second(idx("cumulativeTime / s")), Day).float}))
  let totalTime = (cycleDf["cumulativeTime / s", float][cycleDf.high]).Second.to(Day)
  let yr = 1.Year.to(Second).float
  pltSeries +
    geom_line() +
    xlab("Date", rotate = -45.0) +
    ylim(-2, 12) +
    xlim(breaks[0] - yr, breaks[^1] + yr) +
    ggtitle(&"Magnetic field in the CAST magnet, 2003 to 2021. Total time on (> {magnetField} T): {totalTime}") +
    scale_x_date(name = "Date", isTimestamp = true,
                 breaks = breaks, formatString = "yyyy") +
                        #dateSpacing = initDuration(weeks = 26), dateAlgo = dtaAddDuration) +
    ggsave(&"/tmp/B_time_series_color_{colorPlots}_filtered_{filterOutliers}_magnetic_field_{magnetField}.pdf",
           width = 800, useTex = pretty, standalone = true)
  pltCum +
    geom_line() +
    xlab("Date", rotate = -45.0) +
    xlim(breaks[0], breaks[^1]) +
    ggtitle(&"Cumulative time the CAST magnet was at > {magnetField} T, between 2003 and 2021") +
    scale_x_date(name = "Date", isTimestamp = true,
                 breaks = breaks, formatString = "yyyy") +
                        #dateSpacing = initDuration(weeks = 26), dateAlgo = dtaAddDuration) +
    ggsave(&"/tmp/B_cumulative_time_magnet_{magnetField}_color_{colorPlots}.pdf", width = 800, useTex = pretty, standalone = true)

  #df[0 .. 5000].showBrowser()
  cycleDf.showBrowser()
  when false:
    echo "Writing CSV"
    df.writeCsv("/tmp/all_magnet_timeseries_data.csv")
    cycleDf.writeCsv("/tmp/magnet_cycles_data.csv")

proc sc(path: string, schemas: string, magnetField = 8.0) =
  # parse docopt string and determine
  # correct proc to call
  let schemaFile = parseVersionSchemaFile(schemas)
  # check whether actual log file (extension fits)
  let (dir, fn, ext) = splitFile(path)
  if ext == ".daq":
    # single slow control file
    let sc = read_sc_logfile(path, schemaFile)
    echo "Parsed log file: ", sc
    echo "Lot's of useless output, huh?"
  else:
    doAssert ext.len == 0, "Invalid slow control with extension " & $ext &
      ". Instead we expect .daq"
    process_log_folder(path, lkSlowControl,
                       schemaFile = schemaFile, magnetField = magnetField)

proc tracking(path: string, h5out = "",
              startTime = "", endTime = "",
              dryRun = false) =
  # check whether actual log file (extension fits)
  let (dir, fn, ext) = splitFile(path)
  if ext == ".log":
    # single tracking file
    let t = read_tracking_logfile(path)
    echo "Parsed slow control file: ", t
    echo "Lot's of useless output, huh?"
  else:
    doAssert ext.len == 0, "Invalid tracking with extension " & $ext &
      ". Instead we expect .log"
    process_log_folder(path, lkTracking, h5out,
                       startTime = startTime, endTime = endTime, dryRun = dryRun)

proc h5file(file: string) =
  when not defined(pure):
    deleteTrackingAttributes(file)
  else:
    echo "INFO: The program was compiled with the `pure` flag. The `h5file` subcommand is not available."

proc allLogs(path: string, schemas: string,
             magnetField = 1.0,
             noColorPlots = false,
             filterOutliers = false,
             pretty = false) =
  let schemaFile = parseVersionSchemaFile(schemas)
  handleAllLogs(path, schemaFile, magnetField, not noColorPlots, filterOutliers, pretty)

proc main() =
  dispatchMulti([sc, help={"path" : "A path to a single slow control file or folder of SC log files is needed.",
                           "magnetField" : "A given magnetic field is used to filter data to `> magnetField`",
                           "schemas" : "The path to the `Versions.idx` schema file description is needed."}],
                [tracking, help={ "path" : "A path to a single tracking file or folder of tracking log files is needed.",
                                  "h5out" : "A path to a H5 file is used to add tracking information to the file.",
                                  "dryRun" : """Can be used to simulate writing or just get information about all runs that " &
"fit into the given run file.""",
                                  "startTime" : "Date in YYYY/MM/dd format from which to print tracking information. This is inclusive.",
                                  "endTime" : "Date in YYYY/MM/dd format from which to print tracking information. This is exclusive.",
                }],
                [allLogs, help={"path" : "A path to the location containing all logs in form of `.tar.gz` files is needed.",
                                "schemas" : "The path to the `Versions.idx` schema file description is needed.",
                                "magnetField" : "Optional cut off to determine if the magnet is on.",
                                "noColorPlots" : "Optional flag to disable coloring of magnet plots.",
                                "filterOutliers" : "Optional flag to filter magnet outlier data (< -0.5 T and > 10 T).",
                                "pretty" : "Optional flag to generate pretty plots."}],
                [h5file, help={"file" : "A HDF5 file containing tracking logs, to delete them from."}])

when isMainModule:
  main()
