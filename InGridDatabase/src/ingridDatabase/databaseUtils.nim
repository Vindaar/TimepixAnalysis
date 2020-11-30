import nimhdf5, re, strutils, times
import parsetoml

import databaseDefinitions

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

proc inRunPeriod(chip: string, grp: H5Group): bool =
  let chips = grp.attrs[RunPeriodChipsAttr, seq[string]]
  result = chip in chips

proc inRunPeriod(run: int, grp: H5Group): bool =
  let runs = grp[RunPeriodRunsAvailable.dset_str][int]
  result = run in runs

proc findRunPeriodFor*(h5f: H5FileObj, chipName: string, run: int): string =
  ## returns the ``first`` run period that matches the condition
  ## `contains chipName and run in validRuns`
  for grp in items(h5f, start_path = "/", depth = 1):
    if chipName.inRunPeriod(grp) and
       run.inRunPeriod(grp):
      return grp.name
  if result.len == 0:
    raise newException(ValueError, "Cannot find any run period matching run " & $run)

proc chipNameToGroup*(chipName: string, period: string): string =
  ## given a `chipName` will return the correct name of the corresponding
  ## chip's group within the given run `period` of the database.
  ## done by parsing the given string to a `ChipName` and using `ChipNames`
  ## `$` proc and prepending the `period`.
  result = "/" & period & "/" & $(parseChipName(chipName))

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
  if result < 0:
    raise newException(ValueError, "Runs in `runsAvailable` have to positive numbers!")
