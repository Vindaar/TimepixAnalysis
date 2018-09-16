import ospaths, sequtils, tables, strformat
import nimhdf5, seqmath, helpers/utils, arraymancer

import databaseUtils
import databaseDefinitions

import ingrid/ingrid_types


# procs to get data from the file
proc getTotCalib*(chipName: string): Tot =
  ## reads the TOT calibration data from the database for `chipName`
  let dsetName = joinPath(chipNameToGroup(chipName), TotPrefix)
  withDatabase:
    var dset = h5f[dsetName.dset_str]
    let data = dset[float64].reshape2D(dset.shape).transpose
    result.pulses = data[0].asType(int)
    result.mean = data[1]
    result.std = data[2]

proc getScurve*[T: SomeInteger](h5f: var H5FileObj,
                                chipName: string,
                                voltage: T):
                                  SCurve =
  ## reads the given voltages' SCurve for `chipName`
  result.name = SCurvePrefix & $voltage
  let dsetName = joinPath(joinPath(chipNameToGroup(chipName),
                                   SCurveFolder),
                          result.name)
  result.voltage = voltage
  var dset = h5f[dsetName.dset_str]
  let data = dset[float64].reshape2D(dset.shape).transpose
  result.thl = data[0].asType(int)
  result.hits = data[1]

proc getScurve*(chipName: string, voltage: int): SCurve =
  ## reads the given voltages' SCurve for `chipName`
  ## overload of above for the case of non opened H5 file
  withDatabase:
    result = h5f.getScurve(chipName, voltage)

proc getScurveSeq*(chipName: string): SCurveSeq =
  ## read all SCurves of `chipName` and return an `SCurveSeq`
  # init result seqs (for some reason necessary?!)
  result.files = @[]
  result.curves = @[]
  var groupName = joinPath(chipNameToGroup(chipName),
                           SCurveFolder)
  echo "done ", groupName
  withDatabase:
    h5f.visitFile()
    groupName = joinPath(chipNameToGroup(chipName),
                         SCurveFolder)
    var grp = h5f[groupName.grp_str]
    for dset in grp:
      var mdset = dset
      let voltage = mdset.attrs["voltage", int64].int
      let curve = h5f.getSCurve(chipName, voltage)
      result.files.add dset.name
      result.curves.add curve

proc getThreshold*(chipName: string): Threshold =
  ## reads the threshold matrix of `chipName`
  withDatabase:
    let dsetName = joinPath(chipNameToGroup(chipName), ThresholdPrefix)
    var dset = h5f[dsetName.dset_str]
    result = dset[int64].toTensor.reshape(256, 256).asType(int)

proc getTotCalibParameters*(chipName: string):
                          (float, float, float, float)
    {.raises: [KeyError, IOError, Exception].} =
  ## returns the factors of the TOT calibration result:
  ## `a`, `b`, `c`, `t`
  ## may raise a `KeyError` if the calibration wasn't performed
  ## for the given chip
  withDatabase:
    let groupName = chipNameToGroup(chipName).grp_str
    var grp = h5f[groupName]
    let
      a = grp.attrs["a", float64]
      b = grp.attrs["b", float64]
      c = grp.attrs["c", float64]
      t = grp.attrs["t", float64]
    result = (a, b, c, t)

proc getCalibVsGasGainFactors*(chipName: string): (float, float) =
  ## returns the fit parameters (no errors) for the given chip
  ## of the calibration of Fe charge spectrum vs gas gain
  withDatabase:
    h5f.visitFile()
    let grpName = chipNameToGroup(chipName)
    if hasKey(h5f.datasets, grpName / ChargeCalibGasGain):
      var dset = h5f[(grpName / ChargeCalibGasGain).dset_str]
      let
        b = dset.attrs["b", float64]
        m = dset.attrs["m", float64]
      result = (b, m)
    else:
      discard h5f.close()
      raise newException(Exception, "Charge calibration vs gas gain dataset " &
                         &"does not exist for chip {parseChipName(chipName)}")
