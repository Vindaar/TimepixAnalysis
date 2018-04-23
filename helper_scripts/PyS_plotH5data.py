#!/usr/bin/env python

import os
import matplotlib
matplotlib.use("TKagg")
import matplotlib.pyplot as plt
import pylab
import datetime
import dateutil.parser
import argparse
from ingrid.ingrid_helper_functions import *
from scipy.optimize import curve_fit

import numpy as np

from fit_fe_spectrum import *


def fancy_plotting():
    # set up some LaTeX plotting parameters
    # still need to change parameters
    # next line is for standard article document
    # fig_width_pt = 478.00812#246.0  # Get this from LaTeX using \showthe\columnwidth
    # next line is for thesis after Brock thesis guide
    fig_width_pt = 451.58598
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = 1.3*fig_width_pt*inches_per_pt  # width in inches
    fig_height = 1.3*fig_width*golden_mean      # height in inches
    fig_size =  [fig_width,fig_height]
    params = {'backend': 'ps',
              'axes.labelsize':      10,#10,
              'axes.titlesize':      10,
              'font.size':           10,
              'legend.fontsize':     8,#10,
              'xtick.labelsize':     8,#8,
              'ytick.labelsize':     8,#8,
              'text.usetex':         True,
              'text.latex.preamble': [r'\usepackage{siunitx} \usepackage{mhchem}'],
              'font.family':         'serif',
              'font.serif':          'cm',
              'figure.figsize':      fig_size}
    pylab.rcParams.update(params)

def plotsForAllDsets(h5file, all_runs = False, h5file2 = None):
    # function which creates histograms for all datasets
    # in the /reconstruction/run_$#/chip_$# folders

    nchips = 7
    dsets = getListOfDsets()
    # open h5file
    if all_runs == False:
        h5f = h5py.File(h5file, "r")
        # basically we will iterate over all groups, and all
        # datasets, which are defined in ingrid_helper_functions
        for grp in iterH5RunGroups(h5f):
            run_number = grp.split("_")[-1]
            for chip in xrange(nchips):
                grp_name = recoBase() + grp + "/chip_{}".format(chip)
                for dset in dsets:
                    print("Reading dataset {}".format(dset))
                    range = getBinRangeForDset(dset)
                    try:
                        data = readH5DataSingle(h5f, grp_name, [dset])
                                 "plots/{0}_run_{1}_chip_{2}".format(dset, run_number, chip),
                                 "{0} of run {1} for chip {2}".format(dset, run_number, chip),
                                 dset,
                                 "# hits",
                                 fitting_only = True)
                    except KeyError:
                        print("Key {} could not be found in file {}".format(dset, h5file))
    if all_runs == True:
        # give recoBase as group name for the argument
        if h5file2 is None:
            for dset in dsets:
                for chip, data in iterH5DatasetAndChps(h5file, dset):
                    range = getBinRangeForDset(dset)
                    nbins = getNumBinsForDset(dset)
                    plotData(data, nbins, range,#binning,
                             "plots/{0}_all_runs_chip_{1}".format(dset, chip),
                             "{0} for chip {1} for all runs".format(dset, chip),
                             dset,
                             "# hits",
                             fitting_only = True)
        else:
            for dset in dsets:
                for chip, data in iterTwoH5DatasetAndChps([h5file, h5file2], dset):
                    range = getBinRangeForDset(dset)
                    plotData(data, None, range,#binning,
                             "plots/{0}_all_runs_chip_{1}".format(dset, chip),
                             "{0} for chip {1} for all runs".format(dset, chip),
                             dset,
                             "# hits",
                             fitting_only = True)


def createGroupName(run_number, all_plots, dset_name, chip, all_runs):
    # this function parses the run_number (either number, list, or None)
    # as well as the all_plots flag
    group_name = []
    if type(run_number) is not list:
        print "broken"
        if "fadc" not in dset_name:
            group_name = recoDataChipBase(run_number) + str(chip)
        else:
            group_name = recoFadcBase(run_number) + "fadc/"
    else:
        for run in run_number:
            if "fadc" not in dset_name:
                group_name.append(recoDataChipBase(run) + str(chip))
            else:
                group_name.append(recoFadcBase(run) + "fadc/")
                        
    if all_runs == True:
        group_name = recoBase()
    return group_name

def parseDsetName(dset_name):
    # function to parse the dataset name for potential ratio of two
    # datasets
    ratio = False
    if "/" in dset_name:
        ratio = True
        dset_name = dset_name.split("/")
    return dset_name, ratio

def filterFullFrames(data, h5file, group_name, chip):
    # given some dataset and the file + group + chip it was read from
    # reads the number of hits in the cluster and filters all events
    # have > 4095 pixels
    hits = np.asarray(readH5Data(h5file, group_name, chip, ["hits"]))
    # flatten concat and flatten, in case we read several runs
    non_full = np.where(hits < 4095)[0]
    # concat and flatten potential list of datasets
    if data.ndim > 1:
        data = np.concatenate(data).flatten()
    # return indices of non full arrays
    return data[non_full]

def readAndPlotRatioDsets(h5file, h5file2, group_name, chip, dset_name, run_number, ratio, ignore_full_frames):
    dataf1 = []
    dataf2 = []
    for dset in dset_name:
        d = readH5Data(h5file, group_name, chip, [dset])
        dataf1.append(d)
        if h5file2 != None:
            d2 = readH5Data(h5file2, group_name, chip, [dset])
            dataf2.append(d2)
    dataf1 = np.asarray(dataf1)
    dataf2 = np.asarray(dataf2)
    if ignore_full_frames == True:
        dataf1 = [filterFullFrames(d, h5file, group_name, chip) for d in dataf1]
        if h5file2 != None:
            dataf2 = [filterFullFrames(d, h5file2, group_name, chip) for d in dataf2]

    # of the remaining datasets, calculate the ratio of the two
    dataf1 = dataf1[0] / dataf1[1]
    if h5file2 != None:
        dataf2 = dataf2[0] / dataf2[1]    

    if dataf2.size > 0:
        data = [dataf1, dataf2]
    else:
        data = [dataf1]
    dset_name = dset_name[0] + "_" + dset_name[1]
    range = getBinRangeForDset(dset_name)
    print range
    plotData(data, None, range,#binning,
             "{0}_{1}".format(dset_name, run_number),
             "{0} of run {1}".format(dset_name, run_number),
             dset_name,
             "# hits")
        

def main(args):

    parser = argparse.ArgumentParser(description = 'H5 Data plotter')
    parser.add_argument('file',
                        help = "The H5 file from which to read data")
    parser.add_argument('--cuts',
                        default = None,
                        dest = "cuts",
                        help = "The cuts to be applied on the data before plotting (unstable and potentially dangerous!)")
    parser.add_argument('--file2',
                        default = None,
                        dest = "file2",
                        help = "A potential second file to plot agains")
    parser.add_argument('--dset_name',
                        default = "FeSpectrum",
                        dest = "dset_name",
                        help = """The data to be plotted. 
                        Note: if you want to plot an FADC
                        dataset (minvals etc.), prepend it with `fadc_` such as:
                        --dset_name fadc_minvals
                        Note2: if you wish to plot a ratio of two datasets, simply 
                        provide the two datasets to divide separated by a '/' without
                        spaces, e.g.:
                        --dset_name length/rmsTransverse""")
    parser.add_argument('--run_number',
                        nargs = "+",
                        default = None,
                        type = int,
                        dest = "run_number",
                        help = """The run number(s) to read from the file. If multiple are provided
                        and --all_runs is not set, we compare the two runs for the given dataset.
                        Currently not supported in combination with --all_plots.""")
    parser.add_argument('--all_plots',
                        action = 'store_true',
                        help = "Toggle this, if you want to create all plots for properties")
    parser.add_argument('--ignore_full_frames',
                        action = 'store_true',
                        help = """Toggle this, if you want to ignore clusters with > 4096 hits.
                        Note: name is taken from event display; reason for 'full_frames' instead
                        of 'clusters'""")
    parser.add_argument('--no_fit',
                        action = 'store_true',
                        help = """Toggle this, if you do not want to call the fitting procedure, despite
                        using the FeSpectrum dataset""")
    parser.add_argument('--all_runs',
                        action = 'store_true',
                        help = "Toggle this, if you want to run over all runs in the file")
    parser.add_argument('--chip',
                        default = None,
                        dest = "chip",
                        help = "The chip for which to read data from `run_number`")
    parser.add_argument('--outfolder',
                        default = ".",
                        dest = "outfolder",
                        help = "The folder in which to save the plots")
    parser.add_argument('--fitting_only',
                        action = 'store_true',
                        help = """Set this if you don't want any plots showed (only saved). 
                        Useful to call script externally, without getting blocking behavior.""")
    parser.add_argument('--fancy',
                        default = False,
                        action = 'store_true',
                        help = "Flag to activate fancy plotting via LaTeX output")

    args_dict = vars(parser.parse_args())
    h5file = os.path.abspath(args_dict["file"])
    print(args_dict)

    if args_dict["fancy"] == True:
        fancy_plotting()

    # in order to define cuts, we can do it the `ROOT` way and parse strings, which
    # are evaluated and then used. E.g.
    # ["data < 300", "data > 100"] as argument
    # which is then on the fly quoted
    # this does work, but of course demands to use the correct names the variables
    # have when the stuff is evaluated etc... Well, plus its obvsiously unsafe, because
    # the code might do evil things!
    cuts = args_dict["cuts"]
    run_number = args_dict["run_number"]
    outfolder = args_dict["outfolder"]
    chip = args_dict["chip"]
    fitting_only = args_dict["fitting_only"]
    all_runs = args_dict["all_runs"]
    all_plots = args_dict["all_plots"]
    ignore_full_frames = args_dict["ignore_full_frames"]
    h5file2 = args_dict["file2"]
    dset_name = args_dict["dset_name"]


    if (chip is None or (run_number is None and all_runs == False)) and all_plots == False and "fadc" not in dset_name:
        import sys
        sys.exit("Please provide a run and chip number!")
    
    group_name = createGroupName(run_number, all_runs, dset_name, chip, all_runs)

    # given the dataset name, parse the user input and determine if we
    # wish to calculate a ratio of two sets
    dset_name, ratio = parseDsetName(dset_name)

    if all_plots == True:
        # in this case almost all other flags and arguments will be ignored
        plotsForAllDsets(h5file, all_runs, h5file2)
    else:
        data = []
        if ratio == False:
            data.append(readH5Data(h5file, group_name, chip, [dset_name]))
            if ignore_full_frames == True:
                data[0] = filterFullFrames(data[0], h5file, group_name, chip)
            if h5file2 != None:            
                data.append(readH5Data(h5file2, group_name, chip, [dset_name]))
                if ignore_full_frames == True:
                    data[1] = filterFullFrames(data[1], h5file2, group_name, chip)
            if args_dict["no_fit"] == True or dset_name != "FeSpectrum":
                #hist, binning = binData(data, cuts)
                range = getBinRangeForDset(dset_name)
                plotData(data, None, range,#binning,
                         "{0}_{1}".format(dset_name, run_number),
                         "{0} of run {1}".format(dset_name, run_number),
                         dset_name,
                         "# hits")
            else:
                fit_results = fitAndPlotFeSpectrum(data, cuts, outfolder, run_number, fitting_only)

                if all_runs == False:
                    # now finally write the fit results back to the H5 file
                    calib_factor = writeFitParametersH5(h5file, fit_results, group_name, dset_name)

                    if fitting_only == True:
                        # now print the calibration factor again so that it can be read as the last line
                        # from a calling process
                        print("\n{}".format(calib_factor))
        else:
            # in this case we plot the ratio of two datasets
            readAndPlotRatioDsets(h5file, h5file2, group_name, chip, dset_name, run_number, ratio, ignore_full_frames)
                
if __name__=="__main__":
    import sys
    main(sys.argv[1:])
