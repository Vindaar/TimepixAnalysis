import std / [re, strutils, times]
from std / algorithm import sortedByIt, SortOrder
from std / json import `$`, pretty

from std / os import `/`
import pkg / [nimhdf5, parsetoml]

import ./databaseDefinitions

################################################################################
## Helpers related to reading information from reconstruction HDF5 files
################################################################################

proc getConstraints*(runPeriod: H5Group, chipName: string): seq[Constraint] =
  ## If there are constraints in this `runPeriod` we will return the specific
  ## constraint values for the `chipName`.
  ##
  ## If there are none, the result will be an empty seq.
  ##
  ## Constraints are stored as a (string, string) dataset of the sort:
  ## - RunPeriod
  ##   - Constraints   <-- group for constraints
  ##     - Chip XYZ    <-- dataset of strings for each chip with the values
  ##     - Chip ABC    <--                     """
  if ConstraintsGroup notin runPeriod: return # no constraints at all in the file
  elif ConstraintsGroup / chipName notin runPeriod: return
  result = runPeriod[(ConstraintsGroup / chipName).dset_str][Constraint]

################################################################################
## General InGrid Database helpers
################################################################################

template withDatabase*(actions: untyped): untyped =
  ## read only template to open database as `hf5`, work with it
  ## and close it properly
  ## NOTE: The database is only opened, if no other variable
  ## of the name `h5f` is declared in the calling scope.
  ## This allows to have nested calls of `withDebug` without
  ## any issues (otherwise we get will end up trying to close
  ## already closed objects)
  var openedHere = false
  when not declaredInScope(h5f):
    echo "Opening db at ", dbPath
    var h5f {.inject.} = H5open(dbPath, "r")
    openedHere = true
  actions
  if openedHere:
    let err = h5f.close()
    if err < 0:
      echo "Could not properly close database! err = ", err

proc parseChipName*(chipName: string): ChipName =
  ## parses the given `chipName` as a string to a `ChipName` object
  var mChip: array[3, string]
  # parse the chip name
  if match(chipName, ChipNameReg, mChip) == true:
    result.col   = mChip[0][0]
    result.row   = mChip[1].parseInt
    result.wafer = mChip[2].parseInt
  else:
    raise newException(ValueError, "Bad chip name: $#" % chipName)

proc formatChipName(name: string): string = $(name.parseChipName)

proc inRunPeriod(chip: string, grp: H5Group): bool =
  if RunPeriodChipsAttr in grp.attrs:
    let chips = grp.attrs[RunPeriodChipsAttr, seq[string]]
    result = chip.formatChipName in chips
  else:
    result = false # means no chips in this period so far

proc inRunPeriod(run: int, grp: H5Group): bool =
  ## Note: Either the `run` is contained in the valid runs for the
  ## current chips *or* any run is valid indicated by the `runs` just
  ## being a dataset with only `-1`.
  let runs = grp[RunPeriodRunDset.dset_str][int]
  result = run in runs or runs == @[-1]

proc findRunPeriodFor*(h5f: H5File, chipName: string, run: int, chipGroup: H5Group = nil): string =
  ## Determines the correct run period for the given chip and run number,
  ## taking into account potential constraints required for the chip.
  ##
  ## A run period is considered to match, if it contains `chipName` and `run` number and if:
  ## - it either defines no constraints at all
  ## - all its constraints are satisfied by the `chipGroup`
  ##
  ## All matching run periods are finally compared based on which period matches the
  ## most constraints.
  ## Note that this implies if you have a run period with 2 constraints and then add
  ## a new run period with only 1 constraint and you wish to use the latter, you'll
  ## either need to keep the other 2 constraints as well _or_ delete the old run period.
  ## Keep in mind that this will only be an issue if the 2 additional constraints are
  ## actually matching the new data for the new constraint.
  ##
  ## If `chipGroup` is `nil` we always ignore any possible constraints for it.
  # 3. find run period based on chip + run + constraint
  type T = tuple[rp: string, cs: int]
  var matches = newSeq[T]()
  for grp in items(h5f, start_path = "/", depth = 1):
    if chipName.inRunPeriod(grp) and
       run.inRunPeriod(grp):
      # check constraints if any
      if chipGroup != nil:
        let rpConstraints = getConstraints(grp, chipName) # get all constraints of this RP
        var allConstraintsMatch = true
        for (c, v) in rpConstraints:         # iter through them
          if c notin chipGroup.attrs:        # does not exist for chip in run
            allConstraintsMatch = false      # -> no match
          else:                              # exists
            chipGroup.attrs.withAttr(c):     # read attribute, depending on type
              if $attr != v:                 # compare as strings. All constraints are strings
                allConstraintsMatch = false  # -> differ, no match
        if allConstraintsMatch:              # all match, keep this RP
          matches.add (rp: grp.name, cs: rpConstraints.len)
      else: # if no chip group, we ignore any constraints
        matches.add (rp: grp.name, cs: 0)
  if matches.len == 0:
    raise newException(ValueError, "Cannot find any run period matching run " & $run)

  matches = matches.sortedByIt(it[1]) # sorts in ascending order
  result = matches[^1].rp             # pick last, i.e. most constraints

proc findRunPeriodFor*(chipName: string, run: int, chipGroup: H5Group = nil): string =
  ## returns the ``first`` run period that matches the condition
  ## `contains chipName and run in validRuns`
  ##
  ## If the given `chipGroup` is `nil`, we do *NOT* check for any constraints
  ## and instead return the default run period based on chip name and run.
  ## This is mainly for the `databaseTool`.
  withDatabase:
    result = h5f.findRunPeriodFor(chipName, run, chipGroup = chipGroup)

proc chipNameToGroup*(chipName: string, period: string): string =
  ## given a `chipName` will return the correct name of the corresponding
  ## chip's group within the given run `period` of the database.
  ## done by parsing the given string to a `ChipName` and using `ChipNames`
  ## `$` proc and prepending the `period`.
  result = "/" & period & "/" & chipName.formatChipName

proc parseTomlTime*(t: TomlValueRef): DateTime =
  case t.kind
  of TomlValueKind.Date:
    let d = t.dateVal
    result = initDateTime(monthday = d.day, month = Month(d.month), year = d.year,
                          hour = 0, minute = 0, second = 0, zone = utc())
  of TomlValueKind.Time:
    raise newException(ValueError, "Invalid time found in toml file! Time needs " &
      "to include a full date, not only a time!")
  of TomlValueKind.Datetime:
    let dt = t.dateTimeVal
    let d = dt.date
    let t = dt.time
    result = initDateTime(monthday = d.day, month = Month(d.month), year = d.year,
                          hour = t.hour, minute = t.minute, second = t.second,
                          zone = utc())
  else:
    doAssert false, "Invalid TOML kind for time"

proc checkAndGetInt*(it: TomlValueRef): int =
  if not (it.kind == TomlValueKind.Int):
    raise newException(ValueError, "Runs in `runsAvailable` have to be integers!")
  result = it.getInt
