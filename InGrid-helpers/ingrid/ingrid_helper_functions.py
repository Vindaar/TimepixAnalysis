# a combination of some useful functions dealing with InGrid related helper functions
# mostly related to plotting data, which was produced by Nim

import numpy as np
import h5py
import os
import matplotlib.pyplot as plt

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

def recoBase():
    return "/reconstruction/"

def recoDataChipBase(run_number):
    return "/reconstruction/run_{}/chip_".format(run_number)

def binData(data, cuts):
    if cuts is not None:
      for cut_str in cuts:
          # careful, this is somewhat dangerous, because we just evaluate any string
          # this could do something evil!
          cut = eval(cut_str)
          data = data[np.where(cut)[0]]
          print data

    # assume binsize of 1 for now
    data = np.concatenate(data).flatten()
    binning = np.linspace(-0.5, np.max(data) + 0.5, np.max(data) + 2)
    
    hist, bin_edges = np.histogram(data, binning)
    bins = np.arange(np.max(data) + 1)

    # return data as tuple (bin content / binning)
    # we remove the last element due to the way we create the bins
    return (hist, bin_edges[:-1])

def plotData(hist, binning, outfile, title, xlabel, ylabel, save_plot = True):
    fig, ax = plt.subplots(1, 1)
    #ax.hist(data, bins = 199)
    if binning is None:
        ax.hist(hist, 500)
    else:
        print(np.shape(hist), np.shape(binning))
        ax.bar(binning, hist, 1., align='edge', linewidth=0.2)
    if title is not None:
        ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    #ax.set_xlim(0, np.percentile(data, 98))
    outname = outfile + ".pdf"

    if save_plot == True:
        plt.savefig(outname)
        plt.show()
        return None
    else:
        return fig, ax

def readH5DataSingle(h5f, group_name, dset_names):
    # proc which reads datasets from the a single group
    # with name `dset_name` in the already opened
    # file h5file
    print("Reading group {}".format(group_name))
    print("Reading dset {}".format(dset_names))
    group = h5f[group_name]
    print(group)
    result = []
    for dset in dset_names:
        print("Opening dset {}".format(dset))
        dset_h5 = group[dset]
        result.append(dset_h5[:])
    return result

def readH5Data(h5file, group_name, chip, dset_names):
    # given a H5 file, a run number one, a chip number
    # and one ore more dataset names, read the data and return
    # a list of the `dset_names` arrays

    print("Opening h5file {}".format(h5file))
    h5f = h5py.File(h5file, "r")
    print(h5f)
    
    # first check whether we read one run or all
    result = []
    if group_name == recoBase():
        basegroup = h5f[group_name]
        # in this case have to iterate over all runs
        for grp in basegroup.keys():
            print "Reading group ", grp
            if "run_" in grp:
                result.append(readH5DataSingle(h5f, recoBase() + grp + "/chip_{}".format(chip), dset_names))
                print result
        result = np.asarray(result).flatten()
    else:
        result = readH5DataSingle(h5f, group_name, dset_names)
        
    h5f.close()
    return result


def writeFitParametersH5(h5file, fit_results, group_name, dset_name):
    # writes the fit results to the H5 file given
    # currently only writing the energy calibration
    # as the results for the spectrum are not actully used in the
    # reconstruction or analysis
    h5f = h5py.File(h5file, "r+")
    group = h5f[group_name]
    dset_h5 = group[dset_name]

    popt, pcov, popt_E, pcov_E = fit_results

    a_inv = 1.0 / popt_E[0] * 1000
    da_inv = a_inv * np.sqrt(pcov_E[0][0]) / popt_E[0]

    dset_h5.attrs["eV_per_pix"] = a_inv
    dset_h5.attrs["d_eV_per_pix"] = da_inv

    h5f.close()

    # return the conversion factor, so that this script can be called
    # from external programs and the factor used
    return a_inv
