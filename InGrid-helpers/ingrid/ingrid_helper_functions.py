# a combination of some useful functions dealing with InGrid related helper functions
# mostly related to plotting data, which was produced by Nim

import numpy as np
import h5py
import os

def readFileColumn(filename):
    lines = open(filename, "r").readlines()
    data = []
    for line in lines:
        if "#" not in line and "nan" not in line:
            data.append(float(line))
    return np.asarray(data)
    
def readFadcSpectrumFile(filename):
    return readFileColumn(filename)

def readHitsFile(filename):
    return readFileColumn(filename)

def readTotsFile(filename):
    return readFileColumn(filename)


def recoDataChipBase(run_number):
    return "/reconstruction/run_{}/chip_".format(run_number)

def readH5Data(h5file, group_name, dset_names):
    # given a H5 file, a run number one, a chip number
    # and one ore more dataset names, read the data and return
    # a list of the `dset_names` arrays

    print("Opening h5file {}".format(h5file))

    h5f = h5py.File(h5file, "r")
    print(h5f)
    print("Reading group {}".format(group_name))
    group = h5f[group_name]
    print(group)
    result = []
    for dset in dset_names:
        print("Opening dset {}".format(dset))
        dset_h5 = group[dset]
        result.append(dset_h5[:])
        
    h5f.close()
    return result


def writeFitParametersH5(h5file, fit_results, group_name, dset_name):
    # writes the fit results to the H5 file given
    # currently only writing the energy calibration
    # as the results for the spectrum are not actully used in the
    # reconstruction or analysis
    h5f = h5py.File(h5file, "w")
    group = h5f[group_name]
    dset_h5 = group[dset_name]

    popt, pcov, popt_E, pcov_E = fit_results

    dset_h5.attrs["eV_per_pix"] = popt_E[0]
    dset_h5.attrs["d_eV_per_pix"] = np.sqrt(pcov_E[0][0])

    h5f.close()
    
