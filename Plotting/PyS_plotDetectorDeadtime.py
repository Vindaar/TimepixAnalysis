#!/usr/bin/env python

# this script is used to plot the dead time of the detector
# for each run, extracted via fadc_analysis --noise_analysis

import os
import matplotlib
matplotlib.use("TKagg")
import matplotlib.pyplot as plt
import pylab
import dateutil.parser
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
              'axes.labelsize':      30,#10,
              'axes.titlesize':      30,
              'font.size':           30,
              'legend.fontsize':     24,#10,
              'xtick.labelsize':     24,#8,
              'ytick.labelsize':     24,#8,
              'text.usetex':         True,
              'text.latex.preamble': [r'\usepackage{siunitx} \usepackage{mhchem}'],
              'font.family':         'serif',
              'font.serif':          'cm',
              'figure.figsize':      fig_size}
    pylab.rcParams.update(params)

def plotTime(t50, t100):
    """
    Proc to plot the dead times for the two different FADC integration times
    t50, t100: dict w/ key == Date string, val: live time
    """
    dates50_s = t50.keys()
    dates100_s = t100.keys()
    dates50 = []
    dates100 = []
    dates50 = map(dateutil.parser.parse, t50.keys())
    dates100 = map(dateutil.parser.parse, t100.keys())
    plt.plot(dates50, t50.values(),
             linestyle = '', marker = '.',
             markersize = 15, color = (147/255.0, 88/255.0, 254/255.0),
             label = "$\\SI{50}{\\nano\\second}$ integration time")
    plt.plot(dates100, t100.values(), linestyle = '', marker = '.',
             markersize = 15, color = (220/255.0, 37/255.0, 102/255.0),
             label = "$\\SI{100}{\\nano\\second}$ integration time")
    plt.legend(loc = 3)
    plt.ylabel("\% shutter was open")
    plt.title("\% of time the shutter was open for different FADC integration times during tracking")
    plt.show()
    
def main(args):

    parser = argparse.ArgumentParser(description = 'InGrid noise analysis plotter')
    parser.add_argument('file',
                        help = "The data file from which to take dead times")
    parser.add_argument('--fancy',
                        default = False,
                        action = 'store_true',
                        help = "Flag to activate fancy plotting via LaTeX output")

    args_dict = vars(parser.parse_args())
    infile = os.path.abspath(args_dict["file"])
    print(args_dict)

    if args_dict["fancy"] == True:
        fancy_plotting()

    flines = open(infile, "r").readlines()
    t50 = {}
    t100 = {}
    for line in flines:
        if "#" not in line:
            line = line.split()
            if line[0] == "fk50":
                t50[line[2]] = line[-1]
            elif line[0] == "fk100":
                t100[line[2]] = line[-1]

    plotTime(t50, t100)

    
if __name__=="__main__":
    import sys
    main(sys.argv[1:])
