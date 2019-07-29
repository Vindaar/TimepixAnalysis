import ../ingrid_types
import helpers/utils
import strutils, ospaths, times, strformat, sequtils, tables, re, algorithm, sets
import nimhdf5
import macros
import pure

#####################################################
# Procs describing the data layout in the HDF5 file #
#####################################################

func cdlGroupName*(tfKindStr, year, dset: string): string =
  let dsetName = dset.extractFilename
  var myr = ""
  case year
  of "2014", "2015":
    myr = "apr2014"
  of "2018", "2019":
    myr = "feb2019"
  result = &"calibration-cdl-{myr}-{tfKindStr}/{dsetName}"

func cdlToXrayBinning2014Map(): Table[string, tuple[bins: int, min, max: float]] =
  ## Maps the names of the `XrayReferenceDataSet.h5` (2014) to the
  ## number of bins and min, max values that must be given to the histogram function
  ## to arrive at the result from the `calibration-cdl.h5` (2014) file.
  result = { "lengthdivbyradius" : (bins: 100, min: 0.9950000047683716, max: 1.985000014305115),
             "skewnessl" : (bins: 100, min: -5.050000190734863, max: 4.849999904632568),
             "skewnesst" : (bins: 100, min: -5.050000190734863, max: 4.849999904632568),
             "rmsy" : (bins: 150, min: -0.01666666753590107, max: 4.949999809265137),
             "excentricity" : (bins: 150, min: 0.9700000286102295, max: 9.909999847412109),
             "pixels" : (bins: 250, min: -0.5, max: 497.5),
             "kurtosisl" : (bins: 100, min: -5.050000190734863, max: 4.849999904632568),
             "kurtosist" : (bins: 100, min: -5.050000190734863, max: 4.849999904632568),
             "length" : (bins: 200, min: -0.05000000074505806, max: 19.85000038146973),
             "xrayperevent" : (bins: 6, min: -1.0, max: 4.0),
             "fractionwithin0.5radius" : (bins: 100, min: -0.004999999888241291, max: 0.9850000143051147),
             "radiusdivbyrmsy" : (bins: 100, min: -0.05000000074505806, max: 9.850000381469727),
             "balance" : (bins: 500, min: -0.0003999999898951501, max: 0.3987999856472015),
             "width" : (bins: 100, min: -0.05000000074505806, max: 9.850000381469727),
             "rmsx" : (bins: 150, min: -0.01666666753590107, max: 4.949999809265137),
             "lengthdivbyrmsy" : (bins: 150, min: -0.1000000014901161, max: 29.70000076293945),
             "rotAngle" : (bins: 100, min: -0.0157079640775919, max: 3.094468832015991),
             "energy" : (bins: 100, min: -0.05000000074505806, max: 9.850000381469727),
             "likelihood" : (bins: 200, min: -40.125, max: 9.625),
             "radius" : (bins: 100, min: -0.02500000037252903, max: 4.925000190734863),
             "fractionwithinrmsy" : (bins: 100, min: -0.004999999888241291, max: 0.9850000143051147),
             "charge" : (bins: 200, min: -6250.0, max: 2481250.0) }.toTable

func cdlToXrayBinning2014*(name: string): tuple[bins: int, min, max: float] =
  const map = cdlToXrayBinning2014Map()
  result = map[name]

func cdlToXrayBinning2018Map(): Table[string, tuple[bins: int, min, max: float]] =
  ## Maps the names of the `XrayReferenceDataSet.h5` (2014) to the
  ## number of bins and min, max values that must be given to the histogram function
  ## to arrive at the result from the `calibration-cdl.h5` (2014) file.
  result = { "skewnessLongitudinal" : (bins: 100, min: -5.050000190734863, max: 4.849999904632568),
             "skewnessTransverse" : (bins: 100, min: -5.050000190734863, max: 4.849999904632568),
             "rmsTransverse" : (bins: 150, min: -0.01666666753590107, max: 4.949999809265137),
             "eccentricity" : (bins: 150, min: 0.9700000286102295, max: 9.909999847412109),
             "hits" : (bins: 250, min: -0.5, max: 497.5),
             "kurtosisLongitudinal" : (bins: 100, min: -5.050000190734863, max: 4.849999904632568),
             "kurtosisTransverse" : (bins: 100, min: -5.050000190734863, max: 4.849999904632568),
             "length" : (bins: 200, min: -0.05000000074505806, max: 19.85000038146973),
             "width" : (bins: 100, min: -0.05000000074505806, max: 9.850000381469727),
             "rmsLongitudinal" : (bins: 150, min: -0.01666666753590107, max: 4.949999809265137),
             "lengthDivRmsTrans" : (bins: 150, min: -0.1000000014901161, max: 29.70000076293945),
             "rotationAngle" : (bins: 100, min: -0.0157079640775919, max: 3.094468832015991),
             "energyFromCharge" : (bins: 100, min: -0.05000000074505806, max: 9.850000381469727),
             "likelihood" : (bins: 200, min: -40.125, max: 9.625),
             "fractionInTransverseRms" : (bins: 100, min: -0.004999999888241291, max: 0.9850000143051147),
             "totalCharge" : (bins: 200, min: -6250.0, max: 2481250.0) }.toTable

func cdlToXrayBinning2018*(name: string): tuple[bins: int, min, max: float] =
  const map = cdlToXrayBinning2018Map()
  result = map[name]

func inCdl2018*(name: string): bool =
  const dsetSet = [ "skewnessLongitudinal",
                    "skewnessTransverse",
                    "rmsTransverse",
                    "eccentricity",
                    "hits",
                    "kurtosisLongitudinal",
                    "kurtosisTransverse",
                    "length",
                    "width",
                    "rmsLongitudinal",
                    "lengthDivRmsTrans",
                    "rotationAngle",
                    "energyFromCharge",
                    "likelihood",
                    "fractionInTransverseRms",
                    "totalCharge" ].toHashSet
  result = name in dsetSet

func inCdl2014*(name: string): bool =
  const dsetSet = [ "EnergyFromCharge",
                    "SkewnessLongitudinal",
                    "RotationAngle",
                    "Radius",
                    "Width",
                    "Length",
                    "KurtosisLongitudinal",
                    "Excentricity",
                    "KurtosisTransverse",
                    "RmsLongitudinal",
                    "RmsTransverse",
                    "SkewnessTransverse",
                    "LikelihoodMarlin",
                    "NumberOfPixels",
                    "TotalCharge",
                    "FractionWithinRmsTransverse" ].toHashSet
  result = name in dsetSet

func cdlToXray2014Map(): Table[string, string] =
  ## Maps the datasets from the `calibration-cdl.h5` (2014) file to the
  ## `XrayReferenceDataSet.h5` (2014) file.
  ## The latter is derived from the former via the charge cuts via:
  ## `func getEnergyBinMinMaxVals*(): Table[string, Cuts]`
  ## applied to the datasets.
  ## The result is binned via the binning described in cdlToXrayBinning2014()
  ## The following datasets are unaccounted for in the map, since they must
  ## be calculated from the existing other datasets:
  ## - "xrayperevent"
  ## - "lengthdivbyradius"
  ## - "lengthdivbyrmsy"
  ## - "fractionwithin0.5radius"
  ## - "radiusdivbyrmsy"
  ## - "balance"
  result = { "EnergyFromCharge" : "energy",
             "SkewnessLongitudinal" : "skewnessl",
             "RotationAngle" : "rotAngle",
             "Radius" : "radius",
             "PositionX" : "",
             "PositionY" : "",
             "Width" : "width",
             "Length" : "length",
             "KurtosisLongitudinal" : "kurtosisl",
             "RunType" : "",
             "EnergyFromPixels" : "",
             "Excentricity" : "excentricity",
             "KurtosisTransverse" : "kurtosist",
             "RmsLongitudinal" : "rmsx",
             "RmsTransverse" : "rmsy",
             "SkewnessTransverse" : "skewnesst",
             "EventNumber" : "",
             "LikelihoodMarlin" : "likelihood",
             "NumberOfPixels" : "pixels",
             "TotalCharge" : "charge",
             "FractionWithinRmsTransverse" : "fractionwithinrmsy",
             "Timestamp" : "",
             "RunNumber" : ""}.toTable

# mapping of 2014 calibration-cdl data to XrayReference data
proc cdlToXray2014*(name: string): string =
  const map = cdlToXray2014Map()
  result = map[name]

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

template rawGroupGrpStr*(): grp_str =
  "/runs/".grp_str

template recoGroupGrpStr*(): grp_str =
  "/reconstruction/".grp_str

template likelihoodGroupGrpStr*(): grp_str =
  "/likelihood/".grp_str

template rawDataBase*(): string =
  "/runs/run_"

template rawDataChipBase*(runNumber: int): string =
  "/runs/run_$#/chip_" % $runNumber # & "$#"

template recoDataChipBase*(runNumber: int): string =
  "/reconstruction/run_$#/chip_" % $runNumber # & "$#"

proc recoPath*(runNumber, chipNumber: int): grp_str {.inline.} =
  result = (recoDataChipBase(runNumber) & $chipNumber).grp_str

proc getGroupNameForRun*(runNumber: int): string =
  # generates the group name for a given run number
  result = rawDataBase() & "$#" % $runNumber

template recoBase*(): string =
  "/reconstruction/run_"

proc recoRunGrpStr*(runNumber: int): grp_str {.inline.} =
  result = (recoBase() & $runNumber).grp_str

template likelihoodBase*(): string =
  "/likelihood/run_"

proc getRecoNameForRun*(runNumber: int): string =
  # generates the reconstrution group name for a given run number
  result = recoBase() & "$#" % $runNumber

proc getRawCombineName*(): string =
  # generates the base path for the combine folder
  result = "/runs/combined/"

proc getRecoCombineName*(): string =
  # generates the base path for the combine folder
  result = "/reconstruction/combined/"

proc hasDset*(h5f: var H5FileObj, runNumber, chipNumber: int, dset: string):
                bool =
  ## returns `true` if the given run and chip has the given `dset`
  let path = recoDataChipBase(runNumber) & $chipNumber / dset
  result = if path in h5f: true else: false

proc hasTotalChargeDset*(h5f: var H5FileObj, runNumber, chipNumber: int):
                           bool {.inline.} =
  ## returns `true` if the given run and chip has the
  ## `totalCharge`
  ## dataset
  result = h5f.hasDset(runNumber, chipNumber, "totalCharge")

proc hasRawRun*(h5f: var H5FileObj, runNumber: int): bool =
  ## checks for the existence of the given run in the file
  let path = rawDataBase & $runNumber
  # this will be the implementation once we reran the whole analysis...
  if path in h5f:
    # check if attribute `done` on run
    var grp = h5f[path.grp_str]
    if "rawDataFinished" in grp.attrs:
      let isDone = grp.attrs["rawDataFinished", string]
      if isDone == "true":
        result = true
      else:
        result = false
    else:
      result = false
  else:
    result = false
  when false:
    # old implementation
    result = if path in h5f: true else: false

proc runFinished*(h5f: var H5FileObj, runNumber: int) =
  ## writes the `rawDataFinished` attribute to the run with
  ## `runNumber`
  let path = rawDataBase & $runNumber
  var grp = h5f[path.grp_str]
  grp.attrs["rawDataFinished"] = "true"

proc getCenterChip*(h5f: var H5FileObj, runNumber: int): int =
  ## reads the `centerChip` attribute from the run group corresponding to
  ## `runNumber`
  result = h5f[recoRunGrpStr(runNumber)].attrs["centerChip", int]

macro createCombineTemplates(name, datatype: string): typed =
  ## creates a template, which returns a basename of the type
  ## combineBasename`name`(chip_number, runNumber): string =
  ##   "/runs/combined/`name`_$#_$#" % [$chip_number, $runNumber]
  var source = ""
  case datatype.strVal
  of "runs":
    source &= "template combineRawBasename" & name.strVal & "*(chip_number, runNumber: int): string =\n"
  of "reconstruction":
    source &= "template combineRecoBasename" & name.strVal & "*(chip_number, runNumber: int): string =\n"
  else:
    discard
  case name.strVal
  of "Hits", "ToT":
    source &= "  \"/" & datatype.strVal & "/combined/" & name.strVal
  else:
    source &= "  \"/" & datatype.strVal.toLowerAscii & "/combined/" & name.strVal.toLowerAscii
  source &= "_$#_$#\" % [$chip_number, $runNumber]"
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

template combineRecoBasenameFadc*(): string =
  "/reconstruction/combined/fadc/"

template combineRecoBasenameNoisy*(runNumber: int): string =
  "/reconstruction/combined/fadc/noisy_$#" % [$runNumber]

template combineRecoBasenameMinvals*(runNumber: int): string =
  "/reconstruction/combined/fadc/minvals_$#" % [$runNumber]

proc fadcRunPath*(runNumber: int): string {.inline.} =
  result = getRecoNameForRun(runNumber) / "fadc"

template noiseBasename*(runNumber: int): string =
  fadcRunPath(runNumber) / "noisy"

template minvalsBasename*(runNumber: int): string =
  fadcRunPath(runNumber) / "minvals"

template rawFadcBasename*(runNumber: int): string =
  getGroupNameForRun(runNumber) / "fadc/raw_fadc"

template trigRecBasename*(runNumber: int): string =
  getGroupNameForRun(runNumber) / "fadc/trigger_record"

template fadcDataBasename*(runNumber: int): string =
  fadcRunPath(runNumber) / "fadc_data"

template fadcBaselineBasename*(runNumber: int): string =
  fadcRunPath(runNumber) / "baseline"

template argMinvalBasename*(runNumber: int): string =
  fadcRunPath(runNumber) / "argMinval"

template riseStartBasename*(runNumber: int): string =
  fadcRunPath(runNumber) / "riseStart"

template fallStopBasename*(runNumber: int): string =
  fadcRunPath(runNumber) / "fallStop"

template riseTimeBasename*(runNumber: int): string =
  fadcRunPath(runNumber) / "riseTime"

template fallTimeBasename*(runNumber: int): string =
  fadcRunPath(runNumber) / "fallTime"

template eventNumberBasename*(runNumber: int): string =
  fadcRunPath(runNumber) / "eventNumber"


################################################################################
##################### HDF5 related helper functions ############################
################################################################################

proc getTrackingEvents*(h5f: var H5FileObj, group: H5Group, num_tracking: int = -1, tracking = true): seq[int] =
  ## given a `group` in a `h5f`, filter out all indices, which are part of
  ## a tracking (`tracking == true`) or not part of a tracking (`tracking == false`)
  ## NOTE: the indices of the whole timestamp array correspond to the event
  ## numbers of a run, since this is a 1:1 mapping.
  ## NOTE2: This only works as long as the raw input does ``actually`` contain
  ## all files! While this was true for Virtex 6 TOS files in 2017/18, it does
  ## not at all hold for V6 2014/15 and SRS TOS files!. Thus get the event numbers
  ## and map to them.
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
      evNums = h5f[group.name / "eventNumber", int64].asType(int)
      tstamp = h5f[group.name / "timestamp", int64]
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
        result = allTrackingInds[num_tracking].mapIt(evNums[it])
      else:
        # flatten the allTrackingInds nested seq and return
        result = flatten(allTrackingInds).mapIt(evNums[it])
    else:
      # all outside trackings are simply the indices, which are not part of a flattened
      # allTrackingInds
      let allTrackingsFlat = flatten(allTrackingInds)
      # and now filter all indices not part of flattened index
      result = toSeq(0 ..< tstamp.len)
        .filterIt(it notin allTrackingsFlat)
        .mapIt(evNums[it])
  except KeyError:
    # in this case there is no tracking information. Keep all indices
    echo &"No tracking information in {group.name} found, use all clusters"
    result = @[]


proc filterTrackingEvents*[T: SomeInteger](cluster_events: seq[T], eventsInTracking: seq[int]): seq[int] =
  ## filters out all event numbers of a reconstructed run for one chip
  ## Need to remove all indices, which are within the tracking indices, but for which
  ## no cluster is found in the datasets, so that we can only read the clusters, which
  ## happened during (or outside) of a tracking
  ## inputs:
  ##   `cluster_events`: all events for one reconstructed chip
  ##   `tracking_inds`: the event numbers which are part (or not) of a tracking
  # set result to the indices of tracking (non tracking), i.e.
  # all allowed events
  result = toSeq(0 .. cluster_events.high)
  if eventsInTracking.len == 0:
    # in case we are handed an empty seq, we simply use all cluster events
    discard
  else:
    # now given all indices describing the cluster event numbers, filter out
    # those, which are ``not`` part of `eventsInTracking``
    result = toSeq(0 ..< cluster_events.len)
      .filterIt(cluster_events[it].int in eventsInTracking)

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

proc getRunInfo*(path: string): RunInfo =
  ## wrapper around the above proc if only the path to the run is known
  let regex = r"^/([\w-_]+/)*data\d{6,9}\.txt$"
  let fadcRegex = r"^/([\w-_]+/)*data\d{6,9}\.txt-fadc$"
  let (is_run_folder, runNumber, rfKind, contains_run_folder) = isTosRunFolder(path)
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

proc getFileInfo*(h5f: var H5FileObj, baseGroup = recoGroupGrpStr()): FileInfo =
  ## returns a set of all run numbers in the given file
  # visit file
  #h5f.visitFile()
  var readAux = false
  # get reconstruction group
  let group = h5f[baseGroup]
  if "runType" in group.attrs:
    result.runType = parseEnum[RunTypeKind](group.attrs["runType", string], rtNone)
  if "runFolderKind" in group.attrs:
    result.rfKind = parseEnum[RunFolderKind](group.attrs["runFolderKind", string],
                                             rfUnknown)
  if "centerChip" in group.attrs:
    result.centerChip = group.attrs["centerChip", int]
  if "centerChipName" in group.attrs:
    result.centerChipName = group.attrs["centerChipName", string]

  for runNumber, group in runs(h5f, data_basename = baseGroup.string / "run_"):
    result.runs.add runNumber.parseInt
    if not readAux:
      let grp = h5f[group.grp_str]
      let nChips = grp.attrs["numChips", int]
      result.chips = toSeq(0 ..< nChips)#.mapIt(it)
      readAux = true
  # sort the run numbers
  result.runs.sort
  echo result

proc parseTracking(grp: H5Group, idx: int): RunTimeInfo =
  let
    start = grp.attrs[&"tracking_start_{idx}", string]
    stop = grp.attrs[&"tracking_stop_{idx}", string]
  result = RunTimeInfo(t_start: start.parseTime("yyyy-MM-dd\'T\'HH:mm:sszzz", utc()),
                       t_end: stop.parseTime("yyyy-MM-dd\'T\'HH:mm:sszzz", utc()))
  result.t_length = result.t_end - result.t_start

proc getExtendedRunInfo*(h5f: var H5FileObj, runNumber: int,
                         runType: RunTypeKind,
                         rfKind: RunFolderKind = rfNewTos): ExtendedRunInfo =
  ## reads the extended run info from a H5 file for `runNumber`
  result.runNumber = runNumber
  let
    grp = h5f[(recoBase() & $runNumber).grp_str]
    tstamp = h5f[grp.name / "timestamp", int64]
    evDuration = h5f[grp.name / "eventDuration", float64]
    nEvents = h5f[(grp.name / "eventNumber").dset_str].shape[0]
  var nFadcEvents = 0
  if minvalsBasename(runNumber) in h5f:
    nFadcEvents = h5f[minvalsBasename(runNumber).dset_str].shape[0]

  result.activeTime = initDuration(seconds = evDuration.foldl(a + b).round.int)
  result.nEvents = nEvents
  result.nFadcEvents = nFadcEvents
  var tInfo = RunTimeInfo(t_start: tstamp[0].fromUnix,
                          t_end: tstamp[^1].fromUnix)
  tInfo.t_length = tInfo.t_end - tInfo.t_start
  result.timeInfo = tInfo

  # calculate background time from trackings
  var trackingDuration = initDuration()

  var numTracking = 0
  if "num_trackings" in grp.attrs:
    numTracking = grp.attrs["num_trackings", int]
  for i in 0 ..< numTracking:
    let tracking = parseTracking(grp, i)
    result.trackings.add tracking
    trackingDuration = trackingDuration + tracking.t_length

  result.trackingDuration = trackingDuration
  result.nonTrackingDuration = result.timeInfo.t_length - trackingDuration

  result.rfKind = rfKind
  result.runType = runType


when isMainModule:
  assert combineRawBasenameToT(0, 1) == "/runs/combined/ToT_0_1"
  assert combineRecoBasenameToT(0, 1) == "/reconstruction/combined/ToT_0_1"
