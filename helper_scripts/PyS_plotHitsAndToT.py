#!/usr/bin/env python

import os
import matplotlib
matplotlib.use("TKagg")
import matplotlib.pyplot as plt
import pylab
import datetime
import dateutil.parser

import numpy as np

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

def plotHitFile(filename, title):

    lines = open(filename, "r").readlines()
    hits = []
    for line in lines:
        if "#" not in line:
            hits.append(int(line))

    hits = np.asarray(hits)
    #hits = hits[hits < np.percentile(hits, 99)]

    binning = np.linspace(-0.5, np.max(hits) + 0.5, np.max(hits) + 2)
    #binning = np.linspace(-0.5, np.max(hits) + 0.5, (np.max(hits) + 2)/2)
    hist, bin_edges = np.histogram(hits, binning)
    bins = np.arange(np.max(hits) + 1)

    #print len(binning), len(hist)
    binning = binning[:-1]
    
    fig, ax = plt.subplots(1, 1)
    ax.bar(binning, hist, 1.)
    if title is not None:
        ax.set_title(title)
    #ax.hist(hits, bins = 100)
    ax.set_xlabel("Pixels hit per event")
    ax.set_ylabel("\# events")
    ax.set_xlim(0, np.percentile(hits, 98))    
    outname = filename.replace(".txt", ".pdf")
    plt.savefig(outname)
    plt.show()

def plotToTFile(filename):

    lines = open(filename, "r").readlines()
    tots = []
    for line in lines:
        if "#" not in line:
            tots.append(int(line))

    binning = np.linspace(-0.5, np.max(tots) + 0.5, np.max(tots) + 2)
    hist, bin_edges = np.histogram(tots, binning)
    bins = np.arange(np.max(tots) + 1)            

    tots = np.asarray(tots)
    fig, ax = plt.subplots(1, 1)
    ax.hist(tots, bins = np.max(tots))
    ax.set_xlabel("ToT values per pixel")
    ax.set_ylabel("\# pixel")
    ax.set_xlim(0, np.percentile(tots, 98))
    outname = filename.replace(".txt", ".pdf")
    plt.savefig(outname)
    plt.show()

def main(args):

    if len(args) > 0:
        folder = args[0]
    files = os.listdir(folder)
    if len(args) > 1:
        chip_select = args[1]
    else:
        chip_select = "."
    if len(args) > 2:
        title = args[2]
    else:
        title = None


    for f in files:
        if ".txt" in f:
            # for each file we need to read the second column containing the dates
            print("Starting to read file %s" % f)
            if "hits" in f and chip_select in f:
                plotHitFile(os.path.join(folder, f), title)
            elif "tot" in f and chip_select in f:
                plotToTFile(os.path.join(folder, f))


if __name__=="__main__":
    import sys
    main(sys.argv[1:])
