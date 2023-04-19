import std / [os, sequtils, tables, strformat, strutils]
import nimhdf5, seqmath, helpers/utils, arraymancer

import databaseUtils
import databaseDefinitions

import ingrid / ingrid_types
import ingrid / private / [hdf5_utils, pure, timepix_utils]


# procs to get data from the file
proc getTotCalib*(h5f: H5File, chipName: string, runPeriod: string): Tot =
  let dsetName = joinPath(chipNameToGroup(chipName, runPeriod), TotPrefix)
  var dset = h5f[dsetName.dset_str]
  let data = dset[float64].reshape2D(dset.shape).transpose
  result.pulses = data[0]
  result.mean = data[1]
  result.std = data[2]

proc getTotCalib*(chipName: string, runPeriod: string): Tot =
  ## reads the TOT calibration data from the database for `chipName`
  withDatabase:
    result = h5f.getTotCalib(chipName, runPeriod)

proc getTotCalib*(chipName: string, run: int): Tot =
  var runPeriod: string
  withDatabase:
    runPeriod = h5f.findRunPeriodFor(chipName, run)
    result = h5f.getTotCalib(chipName, runPeriod)

proc getScurve*[T: SomeInteger](h5f: var H5FileObj,
                                chipName: string,
                                runPeriod: string,
                                voltage: T):
                                  SCurve =
  ## reads the given voltages' SCurve for `chipName`
  result.name = SCurvePrefix & $voltage
  let dsetName = joinPath(joinPath(chipNameToGroup(chipName, runPeriod),
                                   SCurveFolder),
                          result.name)
  result.voltage = voltage
  var dset = h5f[dsetName.dset_str]
  let data = dset[float64].reshape2D(dset.shape).transpose
  result.thl = data[0].asType(int)
  result.hits = data[1]

proc getScurve*[T: SomeInteger](h5f: var H5FileObj,
                                chipName: string,
                                run: int,
                                voltage: T):
                                  SCurve =
  let runPeriod = h5f.findRunPeriodFor(chipName, run)
  result = h5f.getScurve(chipName, runPeriod, voltage)

proc getScurve*(chipName: string,
                run: int,
                voltage: int): SCurve =
  ## reads the given voltages' SCurve for `chipName`
  ## overload of above for the case of non opened H5 file
  withDatabase:
    result = h5f.getScurve(chipName, run, voltage)

proc getScurve*(chipName: string,
                runPeriod: string,
                voltage: int): SCurve =
  ## reads the given voltages' SCurve for `chipName`
  ## overload of above for the case of non opened H5 file
  withDatabase:
    result = h5f.getScurve(chipName, runPeriod, voltage)

proc getScurveSeq*(chipName: string, run: int | string): SCurveSeq =
  ## read all SCurves of `chipName` and return an `SCurveSeq`
  withDatabase:
    h5f.visitFile()
    when typeof(run) is string:
      let runPeriod = run
    else:
      let runPeriod = h5f.findRunPeriodFor(chipName, run)
    var groupName = joinPath(chipNameToGroup(chipName, runPeriod),
                             SCurveFolder)
    var grp = h5f[grp_str(groupName)]
    for dset in items(grp):
      var mdset = dset
      let voltage = mdset.attrs["voltage", int64].int
      let curve = h5f.getSCurve(chipName, runPeriod, voltage)
      result.files.add dset.name
      result.curves.add curve

proc getThreshold*(chipName: string, runPeriod: string): Threshold =
  ## reads the threshold matrix of `chipName`
  withDatabase:
    let dsetName = joinPath(chipNameToGroup(chipName, runPeriod), ThresholdPrefix)
    var dset = h5f[dsetName.dset_str]
    result = dset[int64].toTensor.reshape(256, 256).asType(int)

proc getThreshold*(chipName: string, run: int): Threshold =
  ## reads the threshold matrix of `chipName`
  var runPeriod: string
  withDatabase:
    runPeriod = h5f.findRunPeriodFor(chipName, run)
  result = getThreshold(chipName, runPeriod)

proc getTotCalibParameters*(chipName: string, runPeriod: string):
                          (float, float, float, float)
    {.raises: [KeyError, IOError, Exception].} =
  ## returns the factors of the TOT calibration result:
  ## `a`, `b`, `c`, `t`
  ## may raise a `KeyError` if the calibration wasn't performed
  ## for the given chip
  withDatabase:
    let groupName = chipNameToGroup(chipName, runPeriod).grp_str
    var grp = h5f[groupName]
    let
      a = grp.attrs["a", float64]
      b = grp.attrs["b", float64]
      c = grp.attrs["c", float64]
      t = grp.attrs["t", float64]
    result = (a, b, c, t)

proc getTimepixVersion*(chipName, runPeriod: string): TimepixVersion =
  withDatabase:
    let groupName = chipNameToGroup(chipName, runPeriod).grp_str
    let grp = h5f[groupName]
    if "timepixVersion" in grp.attrs:
      result = parseEnum[TimepixVersion](grp.attrs["timepixVersion", string])
    else:
      result = Timepix1

proc getChipNumber*(chipName, runPeriod: string): int =
  withDatabase:
    let groupName = chipNameToGroup(chipName, runPeriod).grp_str
    let grp = h5f[groupName]
    if "chipNumber" notin grp.attrs:
      raise newException(Exception, "Chip number not found in database for chip of name " &
        $chipName & " in run period " & $runPeriod & ".")
    result = grp.attrs["chipNumber", string].parseInt

proc getTotCalibParameters*(chipName: string, run: int):
                          (float, float, float, float)
    {.raises: [KeyError, IOError, Exception].} =
  ## returns the factors of the TOT calibration result:
  ## `a`, `b`, `c`, `t`
  ## may raise a `KeyError` if the calibration wasn't performed
  ## for the given chip
  var runPeriod: string
  withDatabase:
    runPeriod = h5f.findRunPeriodFor(chipName, run)
  result = getTotCalibParameters(chipName, runPeriod)

proc getCalibVsGasGainFactors*(chipName: string, run: int, suffix = ""): tuple[b, m: float] =
  ## returns the fit parameters (no errors) for the given chip
  ## of the calibration of Fe charge spectrum vs gas gain
  withDatabase:
    h5f.visitFile()
    let runPeriod = h5f.findRunPeriodFor(chipName, run)
    let grpName = chipNameToGroup(chipName, runPeriod)
    let dsetName = ChargeCalibGasGain & suffix
    if hasKey(h5f.datasets, grpName / dsetName):
      var dset = h5f[(grpName / dsetName).dset_str]
      let
        b = dset.attrs["b", float64]
        m = dset.attrs["m", float64]
      result = (b: b, m: m)
    else:
      discard h5f.close()
      raise newException(Exception, "Charge calibration vs gas gain dataset " &
                         &"does not exist for chip {parseChipName(chipName)}")

proc inDatabase*(chipName: string, run: int): bool =
  ## used to check whether a chip is contained in the InGridDatabase
  if chipName == SrsDefaultChipName:
    return false
  withDatabase:
    h5f.visitFile()
    let runPeriod = h5f.findRunPeriodFor(chipName, run)
    result = chipNameToGroup(chipName, runPeriod) in h5f

proc readToTFile*(filename: string,
                  startRead = 0.0,
                  totPrefix = "TOTCalib"): (int, Tot) =
  if filename.endsWith(".txt"):
    result = readToTFileTpx1(filename, startRead, totPrefix)
  elif filename.endsWith(".h5"):
    result = readToTFileTpx3(filename, startRead, totPrefix)
  else:
    raise newException(IOError, "Given filename " & $filename & " does not seem to be a " &
      "ToT calibration file.")

from pkg / unchained import FemtoFarad, fF
proc initCalibInfo*(h5f: H5File,
                    runNumber: int,
                    chipName: string, chipNumber: int,
                    capacitance: FemtoFarad,
                    basePath = recoBase()): CalibInfo =
  # get factors for charge calibration
  let (a, b, c, t) = getTotCalibParameters(chipName, runNumber)
  echo "Getting parameters for ", runNumber, " of hip ", chipName, " a ", a, " b ", b, " c ", c, " t ", t

  # get factors for charge / gas gain fit
  let (bL, mL) = getCalibVsGasGainFactors(chipName, runNumber, suffix = $gcIndividualFits)
  # now compute gas gain to use here by computing mean of all gas gain slices in this run (most sensible)
  let group = basePath & $runNumber
  let gain = h5f[group / &"chip_{chipNumber}/gasGainSlices", GasGainIntervalResult].mapIt(it.G).mean
  result = CalibInfo(a: a, b: b, c: c, t: t, mL: mL, bL: bL, capacitance: capacitance, gain: gain)

proc initCalibInfo*(h5f: H5File): CalibInfo =
  let fileInfo = h5f.getFileInfo()
  let run = fileInfo.runs[0]
  let capacitance = getCapacitance(fileInfo.timepix)
  result = h5f.initCalibInfo(run, fileInfo.centerChipName, fileInfo.centerChip, capacitance)
