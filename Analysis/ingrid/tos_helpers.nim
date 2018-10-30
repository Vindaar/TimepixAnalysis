import os, ospaths
import strutils, strformat, strscans
import parseutils
import times
import algorithm
import re
import tables
import memfiles
import sequtils, future
import threadpool
import math
import streams, parsecsv
import sets
import typetraits

# cus modules
import helpers/utils
import ingrid_types
import seqmath

# other modules
import loopfusion
import zero_functional

import macros

import private / [pure, geometry]
export pure, geometry

when not defined(pure):
  import nimhdf5
  import arraymancer
  import private / [hdf5_utils, arraymancer_utils, python_utils]
  export hdf5_utils, arraymancer_utils, python_utils


{.deadCodeElim: on.}

const
  # some helper constants
  StartTot* = 20.0
  # constant regex for InGrid type events for the Virtex TOS
  eventRegexVirtex = r".*data\d{4,6}(_1_[0-9]+)?.*\.txt$"
  # constant regex for InGrid type events for the SRS TOS
  eventRegexSrs = r".*run_(\d{6})_data_(\d{6})_(\d{6})_(\d{2})-(\d{2})-(\d{2}).txt"

  newVirtexRunRegex = r".*Run_(\d+)_\d{6}-\d{2}-\d{2}.*"
  oldVirtexRunRegex = r".*Run\d{6}_\d{2}-\d{2}-\d{2}.*"
  srsRunRegex       = r".*Run_(\d{6})_\d{6}_\d{2}-\d{2}-\d{2}.*"

  # default chip names, shutter modes and shutter times
  OldShutterMode = "verylong"
  OldShutterTime = "13"
  OldChipName = "D3 W63"
  OldTosRunDescriptorPrefix* = r".*\/(\d{1,3})-"

#####################################################
#### Procs specifically realted to hardware #########
#####################################################

proc getSeptemHChip*(chipNumber: int): string =
  ## returns the name of a given SeptemH chip
  const names = ["E6 W69",
                 "K6 W69",
                 "H9 W69",
                 "H10 W69",
                 "G10 W69",
                 "D9 W69",
                 "L8 W69"]
  result = names[chipNumber]

################################################################################
##################### procs related to X-ray reference datasets ################
################################################################################

template cdlPrefix*(): string =
  ## return the prefix of the group names in the `calibration-cdl.h5` file
  ## part of the names of the reference names below
  "calibration-cdl-apr2014-"

proc getChristophCutVals*(): Table[string, float] =
  ## returns the cut values used by Christoph in his PhD thesis
  ## to compare with these values
  result = { "C-EPIC-0.6kV" :  11.7,
             "Cu-EPIC-0.9kV" : 10.7,
             "Cu-EPIC-2kV" :   9.7,
             "Al-Al-4kV" :     9.1,
             "Ag-Ag-6kV" :     8.1,
             "Ti-Ti-9kV" :     7.7,
             "Mn-Cr-12kV" :    7.6,
             "Cu-Ni-15kV" :    7.4 }.toTable()

proc getXrayRefTable*(): Table[int, string] =
  ## returns a table mapping the different energy bins to the correct
  ## datasets in the X-ray reference file
  # NOTE: we could also simply store this in a seq...
  result = { 0: "C-EPIC-0.6kV",
             1: "Cu-EPIC-0.9kV",
             2: "Cu-EPIC-2kV",
             3: "Al-Al-4kV",
             4: "Ag-Ag-6kV",
             5: "Ti-Ti-9kV",
             6: "Mn-Cr-12kV",
             7: "Cu-Ni-15kV" }.toTable()

proc getEnergyBinning(): seq[float] =
  ## returns the binning of the energy (upper range for each bin)
  ## as a sequence of floats
  result = @[0.4, 0.7, 1.2, 2.1, 3.2, 4.9, 6.9, Inf]

proc toRefDset*(energy: float): string =
  ## returns the correct X-ray reference table for a given
  ## `energy`
  # define xray table as global to only initialize it once
  const
    xray_table = getXrayRefTable()
    binning = getEnergyBinning()
  let ind = binning.lowerBound(energy)
  result = xray_table[ind]

proc findReplace(cImpl: NimNode, assgn: NimNode): NimNode =
  ## iterates through `cImpl`, looks for the assignment in `assgn`
  ## and if found replaces the value within. If not found appends
  ## the given statement to the NimNode
  result = copy(cImpl)
  for i in 1 ..< cImpl.len:
    let curNode = cImpl[i]
    case curNode.kind
    of nnkExprColonExpr:
      # check if field is `assgn` field
      if eqIdent(curNode[0].strVal, assgn[0].strVal):
        # replace this nodes value
        result[i][1] = assgn[1]
        # return now
        return result
    else:
      let msg = "No, shouldn't be here: " & curNode.repr
      error(msg)
  # if we end up here we didn't find the identifier
  let newExpr = nnkExprColonExpr.newTree(
    assgn[0],
    assgn[1]
  )
  result.add newExpr

macro replace(c: typed, x: untyped): untyped =
  ## mini dsl to replace specific fields of a given object
  ## with values given in body
  ##
  ## Example:
  ## .. code-block:
  ##   type
  ##     MyObj = object
  ##       a: int
  ##       b: seq[int]
  ##       c: bool
  ##       d: float
  ##   let obj = MyObj(a: 5,
  ##                   b: @[1, 2, 3]
  ##                   c: false)
  ##   let nObj = replace(obj):
  ##     c = true
  ##     d = 5.5
  ## # will be rewritten to
  ##   let nObj = MyObj(a: 5,
  ##                   b: @[1, 2, 3]
  ##                   c: true,
  ##                   d: 5.5)
  var cImpl = c.getImpl
  # now just replace the correct cImpl fields
  # by the given fields of the user
  for ch in x:
    cImpl = findReplace(cImpl, ch)
  result = cImpl
  echo result.repr

func getXraySpectrumCutVals*(): Table[string, Cuts] =
  ## returns a table of Cuts (kind ckXray) objects, one for each energy bin
  let baseCut = Cuts(kind: ckXray,
                     minPix: 3,
                     cutTo: crSilver,
                     maxLength: Inf,
                     minRms: 0.1,
                     maxRms: 1.1,
                     maxEccentricity: Inf)
  let range0 = replace(baseCut):
    maxLength = 6.0
  let range1 = replace(baseCut):
    maxEccentricity = 2.0
  let range2 = replace(baseCut):
    maxEccentricity = 2.0
  let range3 = replace(baseCut):
    maxEccentricity = 2.0
  let range4 = replace(baseCut):
    maxEccentricity = 1.4
    maxRms = 1.0
    maxLength = 6.0
  let range5 = replace(baseCut):
    maxEccentricity = 1.3
    maxRms = 1.0
  let range6 = replace(baseCut):
    maxEccentricity = 1.3
    maxRms = 1.0
  let range7 = replace(baseCut):
    maxEccentricity = 1.3
    maxRms = 1.0
  let
    ranges = [range0, range1, range2, range3, range4, range5, range6, range7]
    xray_ref = getXrayRefTable()

  result = initTable[string, Cuts]()
  for key, vals in pairs(xray_ref):
    result[vals] = ranges[key]


func getEnergyBinMinMaxVals*(): Table[string, Cuts] =
  ## returns a table of Cuts (kind ckReference) objects, one for each energy bin
  let baseCut = Cuts(kind: ckReference,
                     minRms: 0.1,
                     maxRms: 1.1,
                     maxLength: 7.0,
                     minPix: 3,
                     minCharge: -Inf,
                     maxCharge: Inf)
  let range0 = replace(baseCut):
    minCharge = 0.0
    maxCharge = 5e4
    minRms = -Inf
    maxRms = Inf
    maxLength = 6.0
  let range1 = replace(baseCut):
    minCharge = 3.0e4
    maxCharge = 8.0e4
    maxLength = 6.0
  let range2 = replace(baseCut):
    minCharge = 7.0e4
    maxCharge = 1.3e5
  let range3 = replace(baseCut):
    minCharge = 5.9e4
    maxCharge = 2.1e5
  let range4 = replace(baseCut):
    minCharge = 2.0e5
    maxCharge = 4.0e5
  let range5 = replace(baseCut):
    minCharge = 2.9e5
    maxCharge = 5.5e5
  let range6 = replace(baseCut):
    minCharge = 3.5e5
    maxCharge = 6.0e5
  let range7 = replace(baseCut):
    minCharge = 5.9e5
    maxCharge = 1e6
  let
    ranges = [range0, range1, range2, range3, range4, range5, range6, range7]
    xray_ref = getXrayRefTable()

  result = initTable[string, Cuts]()
  for key, vals in pairs(xray_ref):
    result[vals] = ranges[key]

func getRegionCut*(region: ChipRegion): CutsRegion =
  const
    xMinChip = 0.0
    xMaxChip = 14.0
    yMinChip = 0.0
    yMaxChip = 14.0

  case region
  of crGold:
    result = CutsRegion(xMin: 4.5,
                        xMax: 9.5,
                        yMin: 4.5,
                        yMax: 9.5,
                        radius: 0.0)
  of crSilver:
    # based on radius of 4.5 from center
    result = CutsRegion(xMin: 0.0,
                        xMax: 0.0,
                        yMin: 0.0,
                        yMax: 0.0,
                        radius: 4.5)
  of crBronze:
    result = CutsRegion(xMin: 0.0,
                        xMax: 0.0,
                        yMin: 0.0,
                        yMax: 0.0,
                        radius: 5.5)
  of crAll:
    result = CutsRegion(xMin: 0.0,
                        xMax: 14.0,
                        yMin: 0.0,
                        yMax: 0.0,
                        radius: 0.0)

func inRegion*(centerX, centerY: float, region: ChipRegion): bool {.inline.} =
  ## returns the result of a cut on a certain chip `region`. Inputs the
  ## `centerX` and `centerY` position of a cluster and returns true if
  ## the cluster is within the region
  const centerChip = 7.0
  case region
  of crGold:
    # make sure this is only initialized once somehow...
    let regCut = getRegionCut(region)
    result = if centerX >= regCut.xMin and
                centerX <= regCut.xMax and
                centerY >= regCut.yMin and
                centerY <= regCut.yMax:
               true
             else:
               false
  of crAll:
    # simply always return good
    result = true
  else:
    # make sure this is only initialized once somehow...
    let regCut = getRegionCut(region)
    # silver and bronze region only different by radius
    let
      xdiff = (centerX - centerChip)
      ydiff = (centerY - centerChip)
      radius = distance(xdiff, ydiff)
    # TODO: gold cut is NOT part of the silver region (see C. Krieger PhD p. 133)
    result = if radius <= regCut.radius: true else : false

when isMainModule:

  assert combineRawBasenameToT(0, 1) == "/runs/combined/ToT_0_1"
  assert combineRecoBasenameToT(0, 1) == "/reconstruction/combined/ToT_0_1"

  let energies = @[0.1, 0.0, 12.4, 4.4, 2.3, 2.0]
  let inds = [0, 0, 7, 5, 4, 3]
  let refs = ["C-EPIC-0.6kV", "C-EPIC-0.6kV", "Cu-Ni-15kV", "Ti-Ti-9kV", "Ag-Ag-6kV", "Al-Al-4kV"]
  let xray_table = getXrayRefTable()
  let binning = getEnergyBinning()
  forEach e in energies, i in inds, r in refs:
    assert(binning.lowerBound(e) == i)
    assert(toRefDset(e) == r)


  echo "All tests passed!"
