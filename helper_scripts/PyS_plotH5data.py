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

def plotsForAllDsets(h5file):
    # function which creates histograms for all datasets
    # in the /reconstruction/run_$#/chip_$# folders
    pass

def main(args):

    parser = argparse.ArgumentParser(description = 'H5 Data plotter')
    parser.add_argument('file',
                        help = "The H5 file from which to read data")
    parser.add_argument('--cuts',
                        default = None,
                        dest = "cuts",
                        help = "The cuts to be applied on the data before plotting (unstable and potentially dangerous!)")
    parser.add_argument('--dset_name',
                        default = "FeSpectrum",
                        dest = "dset_name",
                        help = "The data to be plotted.")
    parser.add_argument('--run_number',
                        default = None,
                        dest = "run_number",
                        help = "The run number to read from the file")
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
    
    cuts = args_dict["cuts"]
    run_number = args_dict["run_number"]
    outfolder = args_dict["outfolder"]
    chip = args_dict["chip"]
    fitting_only = args_dict["fitting_only"]
    all_runs = args_dict["all_runs"]
    if chip is None or (run_number is None and all_runs == False):
        import sys
        sys.exit("Please provide a run and chip number!")
    
    # combine = args_dict["combine"]
    # filter_s = args_dict["filter_s"]
    to_plot = {}

    group_name = recoDataChipBase(run_number) + str(chip)
    dset_name = args_dict["dset_name"]
                        
    if all_runs == True:
        group_name = recoBase()
    data = readH5Data(h5file, group_name, chip, [dset_name])
    #print("Read the following data {}".format(len(data[0])))
        

    # in order to define cuts, we can do it the `ROOT` way and parse strings, which
    # are evaluated and then used. E.g.
    # ["data < 300", "data > 100"] as argument
    # which is then on the fly quoted
    # this does work, but of course demands to use the correct names the variables
    # have when the stuff is evaluated etc... Well, plus its obvsiously unsafe, because
    # the code might do evil things!


    if dset_name != "FeSpectrum":
        #hist, binning = binData(data, cuts)
        plotData(data, None,#binning,
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


if __name__=="__main__":
    import sys
    main(sys.argv[1:])
