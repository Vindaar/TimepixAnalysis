import std / os
import pkg / nimhdf5

import ../ingrid_types
from ./geometry import inRegion

## XXX: replace at some point by handling this manually!
converter toGenericCut*(x: tuple[dset: string, lower, upper: float]): GenericCut =
  result = GenericCut(dset: x.dset, min: x.lower, max: x.upper,
                      inverted: false)

converter toGenericCut*(x: varargs[tuple[dset: string, lower, upper: float]]): seq[GenericCut] =
  for el in x:
    result.add el.toGenericCut()

proc cutOnProperties*(h5f: H5File,
                      group: H5Group,
                      region: ChipRegion,
                      cuts: varargs[GenericCut]): seq[int] =
  ## applies the cuts from `cuts` and returns a sequence of indices, which pass
  ## the cut.
  ## Any datasets given will be converted to float after reading.
  ## For usage with FADC data, be careful to extract the event numbers using the
  ## indices first, before applying the indices on InGrid data.
  # first get data
  var dsets = newSeqOfCap[seq[float]](cuts.len)
  for c in cuts:
    let dset = h5f[(group.name / c.dset).dset_str]
    # use `convertType` proc from nimhdf5
    let convert = dset.convertType(float)
    dsets.add dset.convert
  # if chip region not all, get `posX` and `posY`
  var
    posX: seq[float]
    posY: seq[float]
  case region
  of crAll:
    discard
  else:
    try:
      # TODO: this is a workaround for now. `centerX` is the name used for
      # TimepixAnalysis, but the `calibration-cdl.h5` data from Marlin uses
      # PositionX. I don't want to add a `year` field or something to this proc,
      # so for now we just depend on an exception.
      posX = h5f.readAs(group.name / "centerX", float)
      posY = h5f.readAs(group.name / "centerY", float)
    except KeyError:
      posX = h5f.readAs(group.name / "PositionX", float)
      posY = h5f.readAs(group.name / "PositionY", float)

  # take any sequence for number of events, since all need to be of the same
  # type regardless
  let nEvents = if dsets.len > 0: dsets[0].len
                elif posX.len > 0: posX.len
                else: 0
  doAssert nEvents > 0, "crAll cannot be combined with no `cuts`"
  for i in 0 ..< nEvents:
    # cut on region if applicable
    if region != crAll and not inRegion(posX[i], posY[i], region):
      continue
    var skipThis = false
    for j in 0 ..< cuts.len:
      let
        el = cuts[j]
      let
        d = dsets[j][i]
      if not el.inverted:
        if d < el.min or d > el.max:
          skipThis = true
          break
      else: # if inverted, remove everything *inside* the cut
        if d > el.min and d < el.max:
          skipThis = true
          break
    if not skipThis:
      # else add this index
      result.add i

proc cutOnProperties*(h5f: H5File,
                      group: H5Group,
                      cuts: varargs[tuple[dset: string,
                                          lower, upper: float]],): seq[int] {.inline.} =
  ## wrapper around the above for the case of the whole chip as region
  result = h5f.cutOnProperties(group, crAll, cuts.toGenericCut())

proc cutOnProperties*(h5f: H5File,
                      group: H5Group,
                      cuts: seq[tuple[dset: string,
                                      lower, upper: float]]): seq[int] {.inline.} =
  result = h5f.cutOnProperties(group, crAll, cuts.toGenericCut())
