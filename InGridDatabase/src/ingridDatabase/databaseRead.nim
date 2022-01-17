import ospaths, sequtils, tables, strformat
import nimhdf5, seqmath, helpers/utils, arraymancer

import databaseUtils
import databaseDefinitions

import ingrid/ingrid_types


# procs to get data from the file
proc getTotCalib*(h5f: H5File, chipName: string, runPeriod: string): Tot =
  let dsetName = joinPath(chipNameToGroup(chipName, runPeriod), TotPrefix)
  var dset = h5f[dsetName.dset_str]
  let data = dset[float64].reshape2D(dset.shape).transpose
  result.pulses = data[0].asType(int)
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

proc getTotCalibParameters*(chipName: string, run: int):
                          (float, float, float, float)
    {.raises: [KeyError, IOError, Exception].} =
  ## returns the factors of the TOT calibration result:
  ## `a`, `b`, `c`, `t`
  ## may raise a `KeyError` if the calibration wasn't performed
  ## for the given chip
  withDatabase:
    let runPeriod = h5f.findRunPeriodFor(chipName, run)
    let groupName = chipNameToGroup(chipName, runPeriod).grp_str
    var grp = h5f[groupName]
    let
      a = grp.attrs["a", float64]
      b = grp.attrs["b", float64]
      c = grp.attrs["c", float64]
      t = grp.attrs["t", float64]
    result = (a, b, c, t)

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
