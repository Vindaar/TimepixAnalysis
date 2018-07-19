import nimhdf5, re, strutils

import databaseDefinitions

template withDatabase*(actions: untyped): untyped =
  ## read only template to open database as `hf5`, work with it
  ## and close it properly
  ## NOTE: The database is only opened, if no other variable
  ## of the name `h5f` is declared in the calling scope.
  ## This allows to have nested calls of `withDebug` without
  ## any issues (otherwise we get will end up trying to close
  ## already closed objects)
  when not declaredInScope(h5f):
    echo "Opening db at ", dbPath
    var h5f {.inject.} = H5File(dbPath, "r")
  actions
  when not declaredInScope(h5f):
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

proc chipNameToGroup*(chipName: string): string =
  ## given a `chipName` will return the correct name of the corresponding
  ## chip's group
  ## done by parsing the given string to a `ChipName` and using `ChipNames`
  ## `$` proc.
  result = $(parseChipName(chipName))

