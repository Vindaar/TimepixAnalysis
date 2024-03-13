import ../ingrid_types
import pkg / unchained

func getCapacitance*(timepix: TimepixVersion): FemtoFarad =
  ## The values are from the manual of the Timepix1 (and 2) as well as
  ## from the Timepix3 manual.
  case timepix
  of Timepix1: result = 8.fF
  of Timepix3: result = 3.fF
