#!/usr/bin/env python

# This script is used to create a scatter plot of
# a certain dataset in a H5 file, given an X and Y
# dataset

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


def main(args):

    parser = argparse.ArgumentParser(description = 'H5 Data plotter')
    parser.add_argument('file',
                        help = "The H5 file from which to read data")
    parser.add_argument('--x_dset',
                        default = None,
                        dest = "x_dset",
                        help = "The dataset to use for the X axis")
    parser.add_argument('--y_dset',
                        default = None,
                        dest = "y_dset",
                        help = "The dataset to use for the X axis")
    parser.add_argument('--dates',
                        default = False,
                        action = 'store_true',
                        help = "Flag to activate plotting against time instead")
    parser.add_argument('--fancy',
                        default = False,
                        action = 'store_true',
                        help = "Flag to activate fancy plotting via LaTeX output")

    args_dict = vars(parser.parse_args())
    h5file = os.path.abspath(args_dict["file"])
    h5f = h5py.File(h5file, "r")    
    print(args_dict)

    x_dset = args_dict["x_dset"]
    y_dset = args_dict["y_dset"]
    if x_dset == None or y_dset == None:
        import sys
        sys.exit("Please provide both an X and Y dataset")

    if args_dict["fancy"] == True:
        fancy_plotting()

    dates = args_dict["dates"]


    # TODO: change group name appropriately
    group_name = rawBase() + "run_102"

    print group_name
    print x_dset
    print y_dset
    
    data_x = readH5DataSingle(h5f, group_name, [x_dset])
    # add chip to group name
    # TODO: remove this... or fit into proper program
    grp_chip = group_name + "/chip_3"
    data_y = readH5DataSingle(h5f, grp_chip, [y_dset])    

    if dates == True:
        plotVsTime(data_x[0], data_y[0], None, None, "hits_vs_time", None, "Time", y_dset, True)
    else:
        plotScatter(data_x, data_y, None, None, "hits_vs_time", None, "Time", y_dset, True)
    
    h5f.close()

    
if __name__=="__main__":
    import sys
    main(sys.argv[1:])
