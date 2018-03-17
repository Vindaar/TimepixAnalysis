# a combination of some useful functions dealing with InGrid related helper functions
# mostly related to plotting data, which was produced by Nim

import numpy as np
import h5py
import os
import matplotlib.pyplot as plt
from datetime import datetime

def getListOfDsets():
    dsets = ["hits", "sumToT", "centerX", "centerY",
             "energyFromPixel", "rmsLongitudinal", "rmsTransverse", 
             "skewnessLongitudinal", "skewnessTransverse", "kurtosisLongitudinal",
             "kurtosisTransverse", "eccentricity", "rotationAngle",
             "length", "width", "fractionInTransveseRms"]
    return dsets

def getBinRangeForDset(dset):
    if dset == "hits":
        return (0, 500)
    elif dset == "energyFromPixel":
        return (0, 10000)
    elif dset == "sumToT":
        # TODO: check
        return (0, 20000)
    elif dset == "kurtosisLongitudinal":
        return (-5, 10)
    elif dset == "kurtosisTransverse":
        return (-2, 8)
    elif dset == "eccentricity":
        return (0, 7.5)
    elif dset == "ToT":
        return (0, 250)
    elif dset == "length_rmsTransverse":
        return (2, 8)
    elif "minvals" in dset:
        return (-0.6, 0.0)
    elif "riseTime" in dset:
        return (50, 120)
    elif "fallTime" in dset:
        return (100, 700)
    else:
        return None
    

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

def rawBase():
    return "/runs/"

def recoBase():
    return "/reconstruction/"

def recoFadcBase(run_number):
    return "/reconstruction/run_{}/".format(run_number)

def recoDataChipBase(run_number):
    return "/reconstruction/run_{}/chip_".format(run_number)

def rawDataChipBase(run_number):
    return "/runs/run_{}/chip_".format(run_number)

def binData(data, cuts):
    if cuts is not None:
      for cut_str in cuts:
          # careful, this is somewhat dangerous, because we just evaluate any string
          # this could do something evil!
          cut = eval(cut_str)
          data = data[np.where(cut)[0]]
          #print data

    # assume binsize of 1 for now
    data = np.concatenate(data).flatten()
    binning = np.linspace(-0.5, np.max(data) + 0.5, np.max(data) + 2)
    
    hist, bin_edges = np.histogram(data, binning)
    bins = np.arange(np.max(data) + 1)

    # return data as tuple (bin content / binning)
    # we remove the last element due to the way we create the bins
    return (hist, bin_edges[:-1])

def plotData(hist, binning, range, outfile, title, xlabel, ylabel, save_plot = True, fitting_only = False):
    fig, ax = plt.subplots(1, 1)
    #ax.hist(data, bins = 199)
    if binning is None:
        try:
            if len(hist) > 1 and len(hist) < 10:
                # in this case we plot 2 files instead of 1
                colors = ["red", "blue", "green", "purple", "sienna"]
                for i in xrange(len(hist)):
                    hist_i = hist[i] #np.concatenate(hist[i]).flatten()
                    ax.hist(hist_i,
                            bins = 249,
                            range = range,
                            normed = True,
                            linewidth = 0.0,
                            alpha = 0.5,
                            color = colors[i])
            else:
                hist = hist #np.concatenate(hist).flatten()
                ax.hist(hist, bins = 500, range = range, linewidth = 0.0)
        except ValueError:
            print("something broken on outfile {}".format(outfile))
            #print hist
            raise
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
        if fitting_only == False:
            plt.show()
        return None
    else:
        return fig, ax

def plotScatter(x_data, y_data, x_range, y_range, outfile, title, xlabel, ylabel, save_plot = True):
    # function to plot a scatter plot of x and y data.
    fig, ax = plt.subplots(1, 1)

    ax.plot(x_data, y_data)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    outname = outfile + ".pdf"

    if save_plot == True:
        plt.savefig(outname)

    plt.show()

def plotVsTime(times, y_data, x_range, y_range, outfile, title, xlabel, ylabel, save_plot = True):
    # function to plot values against time. Times needs to be an array, which stores the time
    # information as a unix timestamp
    fig, ax = plt.subplots(1, 1)

    print "Converting..."
    print np.size(times)
    print np.size(y_data)
    # TODO: the following is a hack to get rid of broken
    # timestamp in some events contained in calibratino_w_timestamps.h5 and
    # background_w_timestamps.h5
    # NOTE: Only works for data from 2017 obviously...
    datetimes = []
    yvals = y_data
    for i, el in enumerate(times):
        #if el < 1514764800:
        datetimes.append(el)
        #yvals.append(y_data[i])

    print "Plotting..."
    ax.plot(datetimes, yvals, linestyle = "", marker = ".")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    outname = outfile + ".pdf"

    if save_plot == True:
        plt.savefig(outname)

    plt.show()    
    

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
        print result
    return result

def iterH5RunGroups(h5f, reco = True):
    # given a H5 group, yields (a tuple; not right now) of group name
    # (and group object; no) for each run in the file
    # if reco is False, we iterate over raw data groups
    # instead
    # NOTE: needs an opened H5 file as input!

    group_name = ""
    if reco == True:
        group_name = recoBase()
    else:
        group_name = rawBase()
    basegroup = h5f[group_name]
    # in this case have to iterate over all runs
    for grp in basegroup.keys():
        if "run_" in grp:
            print "Reading group ", grp
            yield grp

def iterH5DatasetAndChps(h5f, dset, run_number = None):
    # this function yields a combined dataset (unless run_number is specified)
    # from a H5 file for all chips one after another
    nchips = 7
    for chip in xrange(nchips):
        data = readH5Data(h5f, recoBase(), chip, [dset])
        yield (chip, data)

def iterTwoH5DatasetAndChps(h5_files, dset, run_number = None):
    # this function yields a combined dataset (unless run_number is specified)
    # from a list of H5 files
    nchips = 7

    for chip in xrange(nchips):
        data = []
        for h5f in h5_files:
            data.append(readH5Data(h5f, recoBase(), chip, [dset]))
        yield (chip, data)

def readFadcInGridDset(h5f, group_name, dset_name):
    # fn which distinguishes datasets based on whether
    # it's for FADC or InGrid
    if "fadc" in dset_name:
        dset = dset_name.split("_")[1]
        result = readH5DataSingle(h5f, group_name, [dset])
        
    else:
        result = readH5DataSingle(h5f, group_name, [dset_name])

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
        for name in iterH5RunGroups(h5f):
            if "fadc" in dset_names[0]:
                grp_name = recoBase() + name + "/fadc"
                result.append(readH5DataSingle(h5f, grp_name, [dset_names[0].lstrip("fadc_")]))
            else:
                grp_name = recoBase() + name + "/chip_{}".format(chip)
                try:
                    single_data = readH5DataSingle(h5f, grp_name, dset_names)
                except KeyError as e:
                    print("Could not find dataset or group, error was {}".format(str(e)))
                    continue
                result.extend(single_data)
        #print result[0][0]
        #print("Sum of all data is ", np.sum(np.asarray(result).flatten()))
        result = np.concatenate(result).flatten()
    else:
        if type(group_name) is list:
            for group in group_name:
                result.append(readFadcInGridDset(h5f, group, dset_names[0]))
            result = np.concatenate(result).flatten()
        else:
            result = readFadcInGridDset(h5f, group_name, dset_names[0])
            result = np.asarray(result)
        
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
