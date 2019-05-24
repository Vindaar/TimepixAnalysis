import helpers / utils

import ../ingrid_types

import tables, algorithm
export tables, algorithm

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

when isMainModule:
  let energies = @[0.1, 0.0, 12.4, 4.4, 2.3, 2.0]
  let inds = [0, 0, 7, 5, 4, 3]
  let refs = ["C-EPIC-0.6kV", "C-EPIC-0.6kV", "Cu-Ni-15kV", "Ti-Ti-9kV", "Ag-Ag-6kV", "Al-Al-4kV"]
  let xray_table = getXrayRefTable()
  let binning = getEnergyBinning()
  forEach e in energies, i in inds, r in refs:
    assert(binning.lowerBound(e) == i)
    assert(toRefDset(e) == r)


  echo "All tests passed!"
