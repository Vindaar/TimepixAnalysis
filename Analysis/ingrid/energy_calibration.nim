import sequtils, strutils
import ospaths
import future
import seqmath
import nimhdf5
import tables

import tos_helpers

proc cutFeSpectrum(data: array[4, seq[float64]], event_num, hits: seq[int64]): seq[int64] =
  ## proc which receives the data for the cut, performs the cut and returns the
  ## event numbers of the passing elements
  ## inputs:
  ##    data: array[4, seq[float64]] = array containing 4 sequences
  ##      - pos_x, pos_y, eccentricity, rms_transverse which we need for cuts
  ##    event_num: seq[int] = sequence containing event numbers of data stored in
  ##        other seqs
  ##    hits: seq[int] = sequence containing the hits of the corresponding event
  ##        which are the data in the final spectrum
  ## outputs:
  ##    seq[int] = hits of the events passing the cuts
  # constants whihc define the cuts
  result = @[]
  const
    cut_x = 7.0
    cut_y = 7.0
    cut_r = 4.5
    cut_ecc_high = 1.3
    cut_rms_trans_high = 1.2

  let
    pos_x = data[0]
    pos_y = data[1]
    ecc = data[2]
    rms_trans = data[3]

  for i in 0 .. pos_x.high:
    let dist = distance( (pos_x[i] - cut_x), (pos_y[i] - cut_y) )
    if dist > cut_r:
      continue
    if ecc[i] > cut_ecc_high:
      continue
    if rms_trans[i] > cut_rms_trans_high:
      continue

    # else we keep these events, hence add event number to output
    result.add hits[i]

proc createFeSpectrum*(h5f: var H5FileObj, run_number: int) =
  ## proc which reads necessary reconstructed data from the given H5 file,
  ## performs cuts (as Christoph) and writes resulting spectrum to H5 file
  ## NOTE: currently does not perform a check whether the run is actually a
  ## calibration run or not...
  ## throws:
  ##    HDF5LibraryError = in case a call to the H5 library fails, this may happen
  # spectrum will be calculated for center chip, everything else makes no
  # sense, since we don't have X-rays there
  # obviously means that calibration will only be good for center chip,
  # but that is fine, since we only care about details on that one
  const chip = 3

  # what we need:
  # pos_x, pos_y, eccentricity < 1.3, transverse RMS < 1.2
  var reco_group = recoDataChipBase(run_number) & $chip
  # get the group from file
  var group = h5f[reco_group.grp_str]
  # get the chip number from the attributes of the group
  let chip_number = group.attrs["chipNumber", int]
  # sanity check:
  assert chip_number == chip
  var
    pos_x_dset = h5f[(group.name / "centerX").dset_str]
    pos_y_dset = h5f[(group.name / "centerY").dset_str]
    ecc_dset   = h5f[(group.name / "eccentricity").dset_str]
    rms_trans_dset = h5f[(group.name / "rmsTransverse").dset_str]
    event_num_dset = h5f[(group.name / "eventNumber").dset_str]
    hits_dset = h5f[(group.name / "hits").dset_str]
  let
    pos_x  = pos_x_dset[float64]
    pos_y  = pos_y_dset[float64]
    ecc    = ecc_dset[float64]
    rms_trans = rms_trans_dset[float64]
    event_num = event_num_dset[int64]
    hits = hits_dset[int64]

  # given this data, filter all events which don't conform
  let hits_spectrum = cutFeSpectrum([pos_x, pos_y, ecc, rms_trans], event_num, hits)
  # with the events to use for the spectrum
  echo "Elements passing cut : ", hits_spectrum.len

  # given hits, write spectrum to file
  var spectrum_dset = h5f.create_dataset(group.name & "/FeSpectrum", hits_spectrum.len, dtype = int)
  spectrum_dset[spectrum_dset.all] = hits_spectrum

proc applyEnergyCalibration*(h5f: var H5FileObj, run_number: int, calib_factor: float) =
  ## proc which applies an energy calibration based on the number of hit pixels in an event
  ## using a conversion factor of unit eV / hit pixel to the run given by run_number contained
  ## in file h5f
  ## throws:
  ##     HDF5LibraryError = in case a call to the H5 library fails, this might be raised

  # what we need:
  # the hits of the clusters is all we need
  var chip_base = recoDataChipBase(run_number)
  # get the group from file
  for grp in keys(h5f.groups):
    if chip_base in grp:
      # now can start reading, get the group containing the data for this chip
      var group = h5f[grp.grp_str]
      # get the chip number from the attributes of the group
      let chip_number = group.attrs["chipNumber", int]
      # get dataset of hits
      var hits_dset = h5f[(grp / "hits").dset_str]
      let hits = hits_dset[int64]

      # now calculate energy for all hits
      let energy = mapIt(hits, float(it) * calib_factor)
      # create dataset for energy
      var energy_dset = h5f.create_dataset(grp / "energyFromPixel", energy.len, dtype = float)
      energy_dset[energy_dset.all] = energy
      # attach used conversion factor to dataset
      energy_dset.attrs["conversionFactorUsed"] = calib_factor




  
