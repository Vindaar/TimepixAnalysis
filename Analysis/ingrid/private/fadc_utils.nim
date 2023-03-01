from std / algorithm import lowerBound
import ingrid / ingrid_types

func toFadcSetting*(run: int): FadcSetting =
  ## Convert the given run number to the `FadcSetting` by the ordinal value
  ## associated with the lowest insertion index in the array of run numbers
  ## describing the setting.
  let idx = @[75.5, 100.5, 120.5, 238.5].lowerBound(run.float) - 1
  result = FadcSetting( idx )
