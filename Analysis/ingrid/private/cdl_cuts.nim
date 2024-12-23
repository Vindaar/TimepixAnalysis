import helpers / utils

import ../ingrid_types

import tables, algorithm
export tables, algorithm

################################################################################
##################### procs related to X-ray reference datasets ################
################################################################################

template cdlPrefix*(year: string): string =
  ## return the prefix of the group names in the `calibration-cdl.h5` file
  ## part of the names of the reference names below
  if "2014" in year:
    "calibration-cdl-apr2014-"
  else:
    "calibration-cdl-feb2019-"

proc getChristophCutVals*(): CutValueInterpolator =
  ## returns the cut values used by Christoph in his PhD thesis
  ## to compare with these values
  var cutTab = { "C-EPIC-0.6kV" :  11.7,
                 "Cu-EPIC-0.9kV" : 10.7,
                 "Cu-EPIC-2kV" :   9.7,
                 "Al-Al-4kV" :     9.1,
                 "Ag-Ag-6kV" :     8.1,
                 "Ti-Ti-9kV" :     7.7,
                 "Mn-Cr-12kV" :    7.6,
                 "Cu-Ni-15kV" :    7.4 }.toOrderedTable()
  result = initCutValueInterpolator(cutTab)

proc getXrayRefTable*(): OrderedTable[int, string] =
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
             7: "Cu-Ni-15kV" }.toOrderedTable()

proc getXrayTfKindTable*(): OrderedTable[int, TargetFilterKind] =
  ## returns a table mapping the different energy bins to the correct
  ## datasets in the X-ray reference file
  # NOTE: we could also simply store this in a seq...
  result = { 0: tfCEpic0_6,
             1: tfCuEpic0_9,
             2: tfCuEpic2,
             3: tfAlAl4,
             4: tfAgAg6,
             5: tfTiTi9,
             6: tfMnCr12,
             7: tfCuNi15 }.toOrderedTable()

proc getInverseXrayRefTable*(): OrderedTable[string, int] =
  ## returns the table mapping the different energy bins to the correct
  ## datasets in the X-ray reference file
  # NOTE: we could also simply store this in a seq...
  result = { "C-EPIC-0.6kV"  : 0,
             "Cu-EPIC-0.9kV" : 1,
             "Cu-EPIC-2kV"   : 2,
             "Al-Al-4kV"     : 3,
             "Ag-Ag-6kV"     : 4,
             "Ti-Ti-9kV"     : 5,
             "Mn-Cr-12kV"    : 6,
             "Cu-Ni-15kV"    : 7 }.toOrderedTable()

proc getEnergyBinning*(): seq[float] =
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

proc toRefTfKind*(energy: float): TargetFilterKind =
  ## returns the correct target filter kind associated with the given energy
  # define xray table as global to only initialize it once
  const
    xrayTab = getXrayTfKindTable()
    binning = getEnergyBinning()
  let ind = binning.lowerBound(energy)
  result = xrayTab[ind]

func getXrayFluorescenceTable*(): OrderedTable[string, string] =
  ## returns a table mapping the target filter combinations to the
  ## names of the fluorescence lines
  result = { "C-EPIC-0.6kV" :  "C Kα",
             "Cu-EPIC-0.9kV" : "O Kα",
             "Cu-EPIC-2kV" :   "Cu Lα" ,
             "Al-Al-4kV" :     "Al Kα",
             "Ag-Ag-6kV" :     "Ag Lα",
             "Ti-Ti-9kV" :     "Ti Kα",
             "Mn-Cr-12kV" :    "Mn Kα",
             "Cu-Ni-15kV" :    "Cu Kα"}.toOrderedTable()

func getXrayFluorescenceLines*(): seq[float] =
  ## returns a sequence of the energey of the fluorescence lines in each target
  ## filter combination. Energy in keV
  result = @[8.04, 5.89, 4.51, 2.98, 1.49, 0.930, 0.525, 0.277].reversed

proc toXrayLineEnergy*(target: TargetFilterKind): float =
  ## Returns the energy of the given target/filter kind associated fluorescence line
  const invXrayTab = getInverseXrayRefTable()
  const energies = getXrayFluorescenceLines()
  result = energies[invXrayTab[$target]]

from std / strutils import parseEnum
proc toXrayLineEnergy*(target: string): float =
  ## Returns the energy of the given target/filter kind associated fluorescence line
  result = toXrayLineEnergy(parseEnum[TargetFilterKind](target))

func getXrayCleaningCuts*(): OrderedTable[string, Cuts] =
  ## returns a table of Cuts (kind ckXray) objects, one for each energy bin. These
  ## cuts are used to remove likely background candidates contaminating the CDL
  ## dataset. They are applied to the raw CDL data before computing the logL
  ## distributions from them.
  ##
  ## `NOTE:` These cuts are the ones used to select the clusters that are used
  ## to `generate the log likelihood distributions` from the CDL data file (typically
  ## called `calibration-cdl.h5` or similar). These are used `in addition` to the
  ## `getEnergyBinMinMaxVals201*`!
  ## This is why thes are `kind: cfXray`, because they refer to likely Xrays to be
  ## used to determine what the LogL distribution is supposed to look like.
  let baseCut = Cuts(kind: ckXray,
                     minPix: 3,
                     cutTo: crSilver,
                     maxLength: Inf,
                     minRms: 0.1,
                     maxRms: 1.3, # 1.1, # old valu
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
    maxRms = 1.2 # 1.0 # old value
    maxLength = 6.0
  let range5 = replace(baseCut):
    maxEccentricity = 1.3
    maxRms = 1.2 # 1.0 # old value
  let range6 = replace(baseCut):
    maxEccentricity = 1.3
    maxRms = 1.2 # 1.0 # old value
  let range7 = replace(baseCut):
    maxEccentricity = 1.3
    maxRms = 1.2 # 1.0 # old value
  let
    ranges = [range0, range1, range2, range3, range4, range5, range6, range7]
    xray_ref = getXrayRefTable()

  result = initOrderedTable[string, Cuts]()
  for key, vals in pairs(xray_ref):
    result[vals] = ranges[key]

func getXrayCleaningCutArray*(): array[TargetFilterKind, Cuts] =
  const xrayTab = getXrayCleaningCuts()
  for tf in TargetFilterKind:
    result[tf] = xrayTab[$tf]

func getEnergyBinMinMaxVals2018*(): OrderedTable[string, Cuts] =
  ## returns a table of Cuts (kind ckReference) objects, one for each energy bin for the
  ## CDL data from February 2019.
  ## The charge cut values are derived from the fits to the main peaks in the spectra.
  ## It's
  ## minCharge = gmu - 3 * gs
  ## maxCharge = gmu + 3 * gs
  ## to cover 99.7 % of the spectrum
  ## for an overview of all spectras that were used to derive these values, see
  ## master thesis of Hendrik Schmick
  ## (that might be important in case the cdl_spectrum_creation code might change, making
  ## these used values void for some reason!).
  ##
  ## `NOTE:` These cuts are the ones used to select the clusters that are used
  ## to `generate the X-ray reference spectra` from the CDL data file (typically
  ## called `calibration-cdl.h5` or similar).
  ## In addition they are applied to the raw CDL data file to filter out clusters to
  ## be used to generate the LogL distributions (this is in combination with `getXraySpectrumCuts`.
  ## This is why thes are `kind: cfReference`!
  let baseCut = Cuts(kind: ckReference,
                     minRms: 0.1,
                     maxRms: 1.2, # 1.1 # old value
                     maxLength: 7.0,
                     minPix: 3,
                     minCharge: -Inf,
                     maxCharge: Inf)
  # the charges are assigned _after_ by reading from the correct group
  # of the CDL data with the charge low / high attributes!
  let range0 = replace(baseCut):
    minRms = -Inf
    maxRms = Inf
    maxLength = 6.0
  let range1 = replace(baseCut):
    maxLength = 6.0
  let range2 = baseCut
  let range3 = baseCut
  let range4 = baseCut
  let range5 = baseCut
  let range6 = baseCut
  let range7 = baseCut
  let
    ranges = [range0, range1, range2, range3, range4, range5, range6, range7]
    xray_ref = getXrayRefTable()

  result = initOrderedTable[string, Cuts]()
  for key, vals in pairs(xray_ref):
    result[vals] = ranges[key]

func getEnergyBinMinMaxVals2014*(): OrderedTable[string, Cuts] =
  ## returns a table of Cuts (kind ckReference) objects, one for each energy bin
  ##
  ## `NOTE:` See the docstring of `getEnergyBinMinMaxVals2018` for more information on the
  ## meaning of these cuts.
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

  result = initOrderedTable[string, Cuts]()
  for key, vals in pairs(xray_ref):
    result[vals] = ranges[key]

func getRegionCut*(region: ChipRegion): CutsRegion =
  const
    xMinChip = 0.0
    xMaxChip = TimepixSize
    yMinChip = 0.0
    yMaxChip = TimepixSize

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
                        xMax: TimepixSize,
                        yMin: 0.0,
                        yMax: 0.0,
                        radius: 0.0)

when isMainModule:
  import sequtils
  let energies = @[0.1, 0.0, 12.4, 4.4, 2.3, 2.0]
  let inds = [0, 0, 7, 5, 4, 3]
  let refs = ["C-EPIC-0.6kV", "C-EPIC-0.6kV", "Cu-Ni-15kV", "Ti-Ti-9kV", "Ag-Ag-6kV", "Al-Al-4kV"]
  let xray_table = getXrayRefTable()
  let binning = getEnergyBinning()
  for i in 0 ..< energies.len:
    assert(binning.lowerBound(energies[i]) == inds[i])
    assert(toRefDset(energies[i]) == refs[i])


  echo "All tests passed!"
