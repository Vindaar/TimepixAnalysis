from std / algorithm import lowerBound
import ingrid / ingrid_types

proc toFadcSetting*(run: int): FadcSetting =
  ## Convert the given run number to the `FadcSetting` by the ordinal value
  ## associated with the lowest insertion index in the array of run numbers
  ## describing the setting.
  result = FadcSetting( @[76, 101, 121, 239].lowerBound(run) )
