import std / [strutils, parsecsv, strformat, streams, os, sequtils]

import pkg / nimhdf5
import ../ingrid_types
import hdf5_utils, cdl_cuts, cut_utils

proc removeSuff*(s, suff: string): string =
  result = s
  result.removeSuffix(suff)

# helper converters to get targets, filters and kV values from a `TargetFilterKind`
proc toTarget*(tfKind: TargetFilterKind): TargetKind = ($tfKind).split("-")[0].parseEnum[:TargetKind]()
proc toFilter*(tfKind: TargetFilterKind): FilterKind = ($tfKind).split("-")[1].parseEnum[:FilterKind]()
proc toHV*(tfKind: TargetFilterKind, withSuffix = true): string =
  if withSuffix: ($tfKind).split("-")[2]
  else: ($tfKind).split("-")[2].removeSuff("kV")

proc toCutStr*(run: CdlRun): string =
  let hv = block:
    if run.hv > 1.0 and run.hv < 10.0:
      &"{run.hv:1}"
    elif run.hv > 10.0:
      &"{run.hv:2}"
    else:
      &"{run.hv:1.1f}"
  result = &"{run.target}-{run.filter}-{hv}kV"

proc toTfKind*(run: CdlRun): TargetFilterKind =
  result = parseEnum[TargetFilterKind](&"{toCutStr(run)}")

proc readRuns*(fname: string): seq[CdlRun] =
  var s = newFileStream(fname, fmRead)
  var parser: CsvParser
  if not s.isNil:
    parser.open(s, fname, separator = '|')
    parser.readHeaderRow()
    discard parser.readRow()
    while parser.readRow:
      let row = parser.row
      let run = CdlRun(number: row[1].strip.parseInt,
                       runType: parseEnum[RunTypeKind](row[2].strip, rtNone),
                       hasFadc: row[3].strip.parseBool,
                       target: parseEnum[TargetKind](row[4].strip, tEmpty),
                       filter: parseEnum[FilterKind](row[5].strip, fEmpty),
                       hv: if row[6].strip.len > 0: row[6].strip.parseFloat else: 0.0)
      result.add run

iterator tfRuns*(h5f: H5File, tfKind: TargetFilterKind,
                 filename: string): (int, H5Group) =
  ## Yields the center chip group of all runs from `filename`,
  ## which match `tfKind`
  let runs = readRuns(filename)
  for r in runs:
    case r.runType
    of rtXrayFinger:
      let tfk = r.toTfKind
      if tfk == tfKind:
        let runGrp = h5f[recoRunGrpStr(r.number)]
        let centerChip = runGrp.attrs["centerChip", int]
        let chpGrp = h5f[(runGrp.name / "chip_" & $centerChip).grp_str]
        yield (r.number, chpGrp)
    else:
      # nothing to yield if not an "XrayFinger" (read CDL) run
      discard

proc getCdlCutIdxs*(h5f: H5File, runNumber, chip: int, tfKind: TargetFilterKind,
                    eMin = 0.0, eMax = 0.0, energyDset = igEnergyFromCharge): seq[int] =
  ## Returns a sequence of all indices which match the X-ray cleaning cuts from the
  ## CDL data of the given run and target/filter.
  let cutTab = getXrayCleaningCuts()
  let grp = h5f[(recoDataChipBase(runNumber) & $chip).grp_str]
  let cut = cutTab[$tfKind]
  if eMin > 0.0 or eMax > 0.0:
    result = cutOnProperties(h5f,
                             grp,
                             cut.cutTo,
                             ("rmsTransverse", cut.minRms, cut.maxRms),
                             ("length", 0.0, cut.maxLength),
                             ("hits", cut.minPix, Inf),
                             ("eccentricity", 0.0, cut.maxEccentricity),
                             (energyDset.toDset(), eMin, eMax))
  else:
    result = cutOnProperties(h5f,
                             grp,
                             cut.cutTo,
                             ("rmsTransverse", cut.minRms, cut.maxRms),
                             ("length", 0.0, cut.maxLength),
                             ("hits", cut.minPix, Inf),
                             ("eccentricity", 0.0, cut.maxEccentricity))

proc readCutCDL*[T](h5f: H5File, runNumber, chip: int, dset: string,
                    tfKind: TargetFilterKind, _: typedesc[T]): seq[T] =
  ## Reads the desired dataset `dset` from `runNumber` by applying the X-ray cleaning cuts
  ## to the data. I.e. get the data valid for the CDL fits.
  let passIdx = h5f.getCdlCutIdxs(runNumber, chip, tfKind)
  let data = h5f.readAs(recoDataChipBase(runNumber) & $chip / dset, T)
  result = passIdx.mapIt(data[it])

proc readCutCDL*[T](h5f: H5File, runNumber, chip: int, dset: string,
                    tfKind: TargetFilterKind, passIdx: seq[int], _: typedesc[T]): seq[T] =
  ## Reads the desired dataset `dset` from `runNumber` by applying the X-ray cleaning cuts
  ## to the data. I.e. get the data valid for the CDL fits.
  let data = h5f.readAs(recoDataChipBase(runNumber) & $chip / dset, T)
  result = passIdx.mapIt(data[it])
