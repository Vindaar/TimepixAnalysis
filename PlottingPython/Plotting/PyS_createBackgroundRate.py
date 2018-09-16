#!/usr/bin/env python

import os
import matplotlib
matplotlib.use("TKagg")
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import pylab
import datetime
import argparse
from ingrid.ingrid_helper_functions import *
import numpy as np

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
              'axes.labelsize':      40,#10,
              'axes.titlesize':      40,
              'font.size':           40,
              'legend.fontsize':     32,#10,
              'xtick.labelsize':     32,#8,
              'ytick.labelsize':     32,#8,
              'text.usetex':         True,
              'text.latex.preamble': [r'\usepackage{siunitx} \usepackage{mhchem}'],
              'font.family':         'serif',
              'font.serif':          'cm',
              'figure.figsize':      fig_size}
    pylab.rcParams.update(params)

def readXrayData(h5file, chip):
    # read data
    group_name = likelihoodBase()
    chipNumber = chip
    energy = readH5Data(h5file, group_name, chipNumber, ["energyFromCharge"])
    return energy

def main(args):
    parser = argparse.ArgumentParser(description = 'H5 Data plotter')
    parser.add_argument('file',
                        help = "The H5 file from which to read data")
    parser.add_argument('--region',
                        default = "gold",
                        dest = "region",
                        help = "The chip region considered")
    parser.add_argument('--chip',
                        default = 3,
                        dest = "chip",
                        help = "The chip to plot data for")
    parser.add_argument('--log',
                        default = False,
                        action = 'store_true',
                        help = "Flag to plot data in semi log y")
    parser.add_argument('--fancy',
                        default = False,
                        action = 'store_true',
                        help = "Flag to activate fancy plotting via LaTeX output")

    args_dict = vars(parser.parse_args())
    h5file = os.path.abspath(args_dict["file"])
    print(args_dict)

    fancy = args_dict["fancy"]
    region = args_dict["region"]
    chip = args_dict["chip"]
    logY = args_dict["log"]
    if fancy == True:
        fancy_plotting()

    energy = readXrayData(h5file, chip) # / 1000.0 division needed for E from P since in eV
    print(np.shape(energy))
    #print energy
    hist, bin_edges = np.histogram(energy, bins=25, range=(0.2, 10))
    hist_err = np.sqrt(hist)

    # scale hist
    if logY == False:
        factor = 1e5
    else:
        factor = 1.0
    time_back = 0
    shutter_open = 0
    year = ""
    if "2014" in h5file:
        year = "2014"
        # Christoph 2014 / 15
        time_back = 4000 * 3600
        shutter_open = 0.97
    elif "2017" in h5file:
        # 2017
        year = "2017"
        time_back = 1123 * 3600
        shutter_open = 0.88
    else:
        import sys
        sys.exit("File needs to state if 2014 or 2017 data!")
    area = (0.95 - 0.45)**2
    bin_width = 0.392
    scale = factor / (time_back * shutter_open * area * bin_width)
    print("Scale is ", scale)
    print("Hist is ", hist)
    hist = hist * scale
    print("Hist is now ", hist)
    hist_err = hist_err * scale

    bins = [(bin_edges[i+1] + bin_edges[i]) / 2.0 for i in range(len(bin_edges) - 1)]
    plt.errorbar(bins, hist, yerr = hist_err, xerr = bin_width / 2.0,
                 linestyle = '',
                 marker = '.',
                 markersize = 15,
                 #color = (147/255.0, 88/255.0, 254/255.0),#'blue',
                 #color = (186/255.0, 31/255.0, 86/255.0),#'blue',a
                 #color = (220/255.0, 37/255.0, 102/255.0),
                 color = (99/255.0, 59/255.0, 171/255.0),
                 #color = (253 / 255.0, 151 / 255.0, 31 / 255.0),
                 linewidth = 2)
    if fancy == True:
        plt.xlabel('Energy / $\\si{\\keV}$')
        if logY == False:
            plt.ylabel('Rate / $\\SI{1e-5}{\\keV \\per \\cm^2 \\per \\s}$')
        else:
            plt.ylabel('Rate / $\\si{\\keV \\per \\cm^2 \\per \\s}$')
    plt.xticks(np.arange(0, 11, 1))
    if region != "all":
        plt.title("Background rate of {} in {} region for chip {}".format(year, region, chip))
    else:
        plt.title("Background rate of {} over whole chip {}".format(year, chip))
    plt.xlim(0, 10)
    if logY == True:
        plt.semilogy()
    plt.grid()
    plt.show()


if __name__=="__main__":
    import sys
    main(sys.argv[1:])
