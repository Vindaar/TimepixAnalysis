import options, strutils

proc getBinRangeForDset*(dset: string): Option[(float, float)] =
  if dset == "hits":
    result = some((0.0, 500.0))
  elif dset == "energyFromPixel":
    result = some((0.0, 10000.0))
  elif dset == "sumTot":
    # TODO: check
    result = some((0.0, 20000.0))
  elif dset == "kurtosisLongitudinal":
    result = some((-5.0, 10.0))
  elif dset == "kurtosisTransverse":
    result = some((-2.0, 8.0))
  elif dset == "eccentricity":
    result = some((1.0, 3.5))
  elif dset == "ToT":
    result = some((0.0, 250.0))
  elif dset == "length_rmsTransverse":
    result = some((2.0, 8.0))
  elif dset == "energyCut":
    result = some((0.0, 10000.0))
  elif dset == "length":
    result = some((0.0, 14.0))
  elif "minvals" in dset:
    result = some((-0.6, 0.0))
  elif "riseTime" in dset:
    result = some((50.0, 120.0))
  elif "fallTime" in dset:
    result = some((100.0, 700.0))
  else:
    result = none((float, float))

proc getNumBinsForDset*(dset: string): Option[int] =
  if dset == "hits":
    result = some(500)
  elif dset == "energyFromPixel":
    result = some(100)
  elif dset == "sumTot":
    # TODO: check
    result = some(100)
  elif dset == "kurtosisLongitudinal":
    result = some(100)
  elif dset == "kurtosisTransverse":
    result = some(100)
  elif dset == "eccentricity":
    result = some(100)
  elif dset == "ToT":
    result = some(250)
  elif dset == "length_rmsTransverse":
    result = some(100)
  elif dset == "energyCut":
    result = some(100)
  elif dset == "length":
    result = some(300)
  elif "minvals" in dset:
    result = some(100)
  elif "riseTime" in dset:
    result = some(100)
  elif "fallTime" in dset:
    result = some(50)
  else:
    result = none(int)
