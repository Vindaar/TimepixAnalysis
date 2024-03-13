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

def plotHitFile(hits, filename, title):

    #hits = hits[hits < np.percentile(hits, 99)]

    fname = os.path.basename(filename)

    binning = np.linspace(-0.5, np.max(hits) + 0.5, np.max(hits) + 2)
    #binning = np.linspace(-0.5, np.max(hits) + 0.5, (np.max(hits) + 2)/2)
    hist, bin_edges = np.histogram(hits, binning)
    bins = np.arange(np.max(hits) + 1)

    #print len(binning), len(hist)
    binning = binning[:-1]
    
    fig, ax = plt.subplots(1, 1)
    ax.bar(binning, hist, 1.)
    if title is not None:
        ax.set_title(fname + title)
    else:
        ax.set_title(fname)
    #ax.hist(hits, bins = 100)
    ax.set_xlabel("Pixels hit per event")
    ax.set_ylabel("\# events")
    ax.set_xlim(0, np.percentile(hits, 98))    
    outname = filename.replace(".txt", ".pdf")
    plt.savefig(outname)
    plt.show()

def plotToTFile(tots, filename, title):

    fname = os.path.basename(filename)
    
    tots = tots[np.where(tots <= 200)[0]]

    binning = np.linspace(-0.5, np.max(tots) + 0.5, np.max(tots) + 2)
    hist, bin_edges = np.histogram(tots, binning)
    bins = np.arange(np.max(tots) + 1)

    fig, ax = plt.subplots(1, 1)
    #ax.hist(tots, bins = 199)
    print(np.shape(hist), np.shape(binning))
    print(binning)
    print bin_edges
    ax.bar(bin_edges[:-1], hist, 1., align='edge')
    if title is not None:
        ax.set_title(fname + title)
    else:
        ax.set_title(fname)
    ax.set_xlabel("ToT values per pixel")
    ax.set_ylabel("\# pixel")
    #ax.set_xlim(0, np.percentile(tots, 98))
    outname = filename.replace(".txt", ".pdf")
    plt.savefig(outname)
    plt.show()

def plotRotAngleFile(angles, filename, title):

    fname = os.path.basename(filename)

    angles = angles[np.where(angles != np.nan)[0]]

    print max(angles)
    #angles = angles[angles < np.percentile(angles, 99)]

    binning = np.linspace(-0.5, max(angles) + 0.5, max(angles) + 2)
    #binning = np.linspace(-0.5, np.max(angles) + 0.5, (np.max(angles) + 2)/2)
    print angles, binning
    print max(angles), min(angles)
    # hist, bin_edges = np.histogram(angles, binning)

    # #print len(binning), len(hist)
    # binning = binning[:-1]
    
    fig, ax = plt.subplots(1, 1)
    #ax.bar(binning, hist, 1.)
    if title is not None:
        ax.set_title(fname + title)
    else:
        ax.set_title(fname)
    ax.hist(angles, bins = 100)
            
    ax.set_xlabel("Rot angle")
    ax.set_ylabel("\# pixel")
    #ax.set_xlim(0, np.percentile(tots, 98))
    outname = filename.replace(".txt", ".pdf")
    plt.savefig(outname)
    plt.show()

def plotFadcSpectrumFile(hits, filename, title):

    fname = os.path.basename(filename)

    #hits = hits[hits < np.percentile(hits, 99)]

    # nbins = 100
    # start = -0.5
    # end = 0
    # width = end - start 
    # binning = np.linspace(start, end, nbins)
    # #binning = np.linspace(-0.5, np.max(hits) + 0.5, (np.max(hits) + 2)/2)
    # hist, bin_edges = np.histogram(hits, binning)
    # bins = np.arange(np.max(hits) + 1)

    # print len(binning), len(hist)
    # binning = binning[:-1]
    
    fig, ax = plt.subplots(1, 1)
    #ax.bar(binning, hist, nbins / width)
    ax.hist(hits, bins = 50)
    ax.set_xlabel("Voltage of dip / $\\si{\\milli\\volt}$")
    if title is not None:
        ax.set_title(fname + title)
    else:
        ax.set_title(fname)
    ax.set_ylabel("\# events")
    ax.set_xlim(-0.5, 0)
    outname = filename.replace(".txt", ".pdf")
    plt.savefig(outname)
    plt.show()

    
def main(args):

    parser = argparse.ArgumentParser(description = 'Nim Data plotter')
    parser.add_argument('folder',
                        help = "The folder containing files with column data (Hits, ToTs, Dips / Peaks, RotAngles)")
    parser.add_argument('--filter',
                        default = '.',
                        dest = "filter_s",
                        help = "Filter to be applied to files, which are in given folder, e.g. 'Hit'")
    parser.add_argument('--title',
                        default = None,
                        dest = "title",
                        help = "Title to be added to plot")
    parser.add_argument('--combine',
                        default = False,
                        action = 'store_true',
                        help = "Flag to combine files matching filter_s in folder into 1 plot")
    parser.add_argument('--fancy',
                        default = False,
                        action = 'store_true',
                        help = "Flag to activate fancy plotting via LaTeX output")

    args_dict = vars(parser.parse_args())
    folder = os.path.abspath(args_dict["folder"])
    files = os.listdir(folder)
    print(args_dict)

    if args_dict["fancy"] == True:
        fancy_plotting()
    
    title = args_dict["title"]
    combine = args_dict["combine"]
    filter_s = args_dict["filter_s"]


    hits = []
    tots = []
    angles = []
    dips = []
    to_plot = {}
    
    for f in files:
        if ".txt" in f:
            # for each file we need to read the second column containing the dates
            print("Starting to read file %s" % f)
            filename = os.path.join(folder, f)
            if filter_s in f and ("hits" in f or "tots" in f or "rot" in f or "fadc" in f):
                data = readFileColumn(filename)
            if "hits" in f and filter_s in f:
                if combine == False:
                    plotHitFile(data, filename, title)
                else:
                    hits.extend(data)
                    to_plot["hits"] = hits
            elif "tot" in f and filter_s in f:
                if combine == False:
                    plotToTFile(data, filename, title)
                else:
                    tots.extend(data)
                    to_plot["tots"] = tots
            elif "rot" in f and filter_s in f:
                if combine == False:
                    plotRotAngleFile(data, filename, title)
                else:
                    angles.extend(data)
                    to_plot["angles"] = angles
            elif "fadc" in f and filter_s in f:
                print("Actually reading file %s" % f)
                if combine == False:
                    plotFADCspectrum(data, filename, title)
                else:
                    dips.extend(data)
                    to_plot["dips"] = dips

    if combine == True:
        for k, v in to_plot.iteritems():
            vals = np.asarray(v)
            print("Vals is %s and its length %s" % (vals, np.shape(vals)))
            if k == "hits":
                plotHitFile(vals, folder, title)
            elif k == "tots":
                plotToTFile(vals, folder, title)
            elif k == "angles":
                plotRotAngleFile(vals, folder, title)
            elif k == "dips":
                plotFadcSpectrumFile(vals, folder, title)
                
        


if __name__=="__main__":
    import sys
    main(sys.argv[1:])
