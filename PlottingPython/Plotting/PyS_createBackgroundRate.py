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

color2017 = (99/255.0, 59/255.0, 171/255.0) #(220/255.0, 37/255.0, 102/255.0)#(147/255.0, 88/255.0, 254/255.0) #(99/255.0, 59/255.0, 171/255.0)
color2014 = (253/255.0, 151/255.0, 31/255.0) #(147/255.0, 88/255.0, 254/255.0) #(220/255.0, 37/255.0, 102/255.0)
colorMarlin = (143/255.0, 192/255.0, 41/255.0)

colorYellow = (255.0/255.0, 193.0/255.0, 7.0/255.0)
colorDarkBlue = (2.0/255.0, 136.0/255.0, 209.0/255.0)
colorBlue = (3.0/255.0, 169.0/255.0, 244.0/255.0)
colorLightBlue = (179.0/255.0, 229.0/255.0, 252.0/255.0)
colorGrey = (189.0/255.0, 189.0/255.0, 189.0/255.0)
colorLightGrey = (117.0/255.0, 117.0/255.0, 117.0/255.0)
colorDarkGrey = (33.0/255.0, 33.0/255.0, 33.0/255.0)

colorOrange = (255.0/255.0, 80.0/255.0, 47.0/255.0)
colorWine = (98.0/255.0, 55.0/255.0, 78.0/255.0)
colorCyan = (0.0/255.0, 255.0/255.0, 255.0/255.0)


def fancy_plotting():
    # set up some LaTeX plotting parameters
    # still need to change parameters
    # next line is for standard article document
    # fig_width_pt = 478.00812#246.0  # Get this from LaTeX using \showthe\columnwidth
    # next line is for thesis after Brock thesis guide
    fig_width_pt = 451.58598
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    #golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    golden_mean = 10.0 / 20.0 #(np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = 1.3*fig_width_pt*inches_per_pt  # width in inches
    fig_height = 1.3*fig_width*golden_mean      # height in inches
    fig_size =  [fig_width,fig_height]
    #params = {'backend': 'ps',
    #          'axes.labelsize':      40,#10,
    #          'axes.titlesize':      40,
    #          'font.size':           40,
    #          'legend.fontsize':     32,#10,
    #          'xtick.labelsize':     32,#8,
    #          'ytick.labelsize':     32,#8,
    #          'text.usetex':         True,
    #          'text.latex.preamble': [r'\usepackage{siunitx} \usepackage{mhchem}'],
    #          'font.family':         'serif',
    #          'font.serif':          'cm',
    #          'figure.figsize':      fig_size}
    params = {'backend': 'ps',
              'axes.labelsize':      20,#10,
              'axes.titlesize':      20,
              'font.size':           20,
              'legend.fontsize':     16,#10,
              'xtick.labelsize':     16,#8,
              'ytick.labelsize':     16,#8,
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
    if not isinstance(chipNumber, int):
        raise TypeError("Chip number MUST be an integer!")
    energy = readH5Data(h5file, group_name, chipNumber, ["energyFromCharge"])
    return energy

def prepareChristophFrameworkPlot(logY):
    ## creates the plot for the background rate of Christoph's analysis framework
    lines = open("../../resources/background-rate-gold.2014+2015.dat", 'r').readlines()[1:]
    energy = []
    rate = []
    ratePlus = []
    rateMinus = []
    for l in lines:
        spl = l.split()
        energy.append(float(spl[0]))
        rate.append(float(spl[1]))
        ratePlus.append(float(spl[2]))
        rateMinus.append(float(spl[3]))
    energy = np.asarray(energy)
    rate = np.asarray(rate)
    ratePlus = np.asarray(ratePlus)
    rateMinus = np.asarray(rateMinus)

    # scale hist
    if logY == False:
        factor = 1e5
    else:
        factor = 1.0
    rate = rate * factor
    ratePlus = ratePlus * factor
    rateMinus = rateMinus * factor
    dRate = (rateMinus, ratePlus)

    year = "Marlin2014"

    plt.errorbar(energy, rate, yerr = dRate,
                 linestyle = '',
                 marker = '.',
                 markersize = 8.,
                 color = colorMarlin,
                 label = year,
                 linewidth = 1.0
    )
    return year, 0

def preparePlot(h5file, chip, logY, CK_binning):
    year = ""
    h5f = h5py.File(h5file, "r")
    lhGrp = h5f["/likelihood"]
    time_back = lhGrp.attrs["totalDuration"]
    shutter_open = 1.0
    if "2014" in h5file:
        year = "2014/15"
        # Christoph 2014 / 15
        #time_back = 4000 * 3600
        #shutter_open = 0.97
        color = color2014
    elif "2017" in h5file:
        # 2017
        year = "2017/18"
        #time_back = 1123 * 3600
        #shutter_open = 0.88
        color = color2017
    elif "2018" in h5file and not "fadc" in h5file and not "scinti" in h5file:
        # 2017
        year = "2018"
        #time_back = 1123 * 3600
        #shutter_open = 0.88
        color = color2017
    elif "2018" in h5file and not "scinti" in h5file:
        # Case of FADC veto
        year = "2018+FADC veto"
        color = (1.0, 0.0, 0.0)#colorBlue
    elif "2018" in h5file and not "fadc" in h5file:
        # Case of scintillators veto
        year = "2018+scinti veto"
        color = colorYellow#DarkGrey
    elif "2018" in h5file:
        # Case of FADC + scintillators veto
        year = "2018+scinti+FADC veto"
        color = colorCyan
    else:
        import sys
        sys.exit("File needs to state if 2014, 2017 or 2018 data!")

    if year == "2014/15":
        chip = 0

    # close H5 file before we read from it again in `readH5Data`
    h5f.close()

    energy = readXrayData(h5file, chip) # / 1000.0 division needed for E from P since in eV
    print(np.shape(energy))
    #print energy
    hist = None
    bin_edges = None
    if CK_binning:
        hist, bin_edges = np.histogram(energy, bins=50, range=(0.2, 10.2))
    else:
        hist, bin_edges = np.histogram(energy, bins=25, range=(0.2, 10))
    print("Bin edges are ", bin_edges)
    hist_err = np.sqrt(hist)

    # scale hist
    if logY == False:
        factor = 1e5
    else:
        factor = 1.0

    area = (0.95 - 0.45)**2
    bin_width = None
    if CK_binning:
        bin_width = 0.2
    else:
        bin_width = 0.392
    scale = factor / (time_back * shutter_open * area * bin_width)
    print("Total duration is {} h".format(time_back / 3600.0))
    print("Scale is ", scale)
    print("Hist is ", hist)
    hist = hist * scale
    print("Hist is now ", hist)
    hist_err = hist_err * scale

    bins = [(bin_edges[i+1] + bin_edges[i]) / 2.0 for i in range(len(bin_edges) - 1)]
    plt.errorbar(bins, hist, yerr = hist_err, xerr = bin_width / 2.0,
                 linestyle = '',
                 marker = '.',
                 markersize = 8,
                 color = color,
                 #color = (147/255.0, 88/255.0, 254/255.0),#'blue',
                 #color = (186/255.0, 31/255.0, 86/255.0),#'blue',a
                 #color = (220/255.0, 37/255.0, 102/255.0),
                 #color = (99/255.0, 59/255.0, 171/255.0),
                 #color = (253 / 255.0, 151 / 255.0, 31 / 255.0),
                 linewidth = 1,
                 label = year)

    return year, chip

def addRegion(fname):
    if "all" in fname:
        return "all_"
    elif "bronze" in fname:
        return "bronze_"
    elif "silver" in fname:
        return "silver_"
    elif "gold" in fname:
        return "gold_"
    else:
        return ""

def main(args):
    parser = argparse.ArgumentParser(description = 'H5 Data plotter')
    parser.add_argument('file',
                        help = "The H5 file from which to read data")
    parser.add_argument('--files',
                        nargs='*',
                        dest = "fileList",
                        help = "More files to be plotted against `file`")
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
    parser.add_argument('--ck_binning',
                        default = False,
                        action = 'store_true',
                        help = "Flag to plot data same binning as Marlin2014 data")
    parser.add_argument('--title',
                        default = None,
                        dest = "title",
                        help = "The plot title to override the auto generated title")
    parser.add_argument('--preliminary',
                        default = False,
                        action = 'store_true',
                        help = "Flag to activate `preliminary` in center of plot")
    parser.add_argument('--show2014',
                        default = False,
                        action = 'store_true',
                        help = "Flag to activate plotting of 2014/15 data w/ Marlin")
    parser.add_argument('--fancy',
                        default = False,
                        action = 'store_true',
                        help = "Flag to activate fancy plotting via LaTeX output")

    args_dict = vars(parser.parse_args())
    h5file = os.path.abspath(args_dict["file"])
    print(args_dict)

    fancy = args_dict["fancy"]
    customTitle = args_dict["title"]
    preliminary = args_dict["preliminary"]
    region = args_dict["region"]
    chip = args_dict["chip"]
    logY = args_dict["log"]
    show2014 = args_dict["show2014"]
    CK_binning = args_dict["ck_binning"]

    if fancy == True:
        fancy_plotting()

    year, chip = preparePlot(h5file, chip, logY, CK_binning)

    # now handle all additional files
    fileList = args_dict["fileList"]
    if fileList is not None:
        for f in fileList:
            year2, chip2 = preparePlot(f, chip, logY, CK_binning)
            year = str(year) + " and " + str(year2)
            chipStr = str(chip) + " and " + str(chip2)

    if show2014:
        year3, chip3 = prepareChristophFrameworkPlot(logY)
        year = str(year) + " and " + str(year3)
        chipStr = str(chip) + " and " + str(chip3)

    if fancy:
        plt.xlabel('Energy / $\\si{\\keV}$')
        if logY == False:
            plt.ylabel('Rate / $\\SI{1e-5}{\\keV \\per \\cm^2 \\per \\s}$')
        else:
            plt.ylabel('Rate / $\\si{\\keV \\per \\cm^2 \\per \\s}$')
    plt.xticks(np.arange(0, 11, 1))
    if customTitle != None:
        plt.title(customTitle)
    elif region != "all":
        plt.title("Background rate of {} in {} region for chip {}".format(year, region, chipStr))
    else:
        plt.title("Background rate of {} over whole chip {}".format(year, chipStr))
    plt.xlim(0, 10)
    if logY == True:
        plt.semilogy()
    plt.grid()
    plt.legend()

    axes = plt.gca()
    axes.grid(color = (1.0, 1.0, 1.0))

    if preliminary:
        ax = plt.gca()
        fontsize = 30
        if fancy:
            fontsize = 75
        plt.text(0.5, 0.5, "preliminary", fontsize = fontsize, rotation = 30,
                 horizontalalignment = "center",
                 verticalalignment = "center",
                 transform = ax.transAxes)

    logStr = "log_" if logY == True else ""
    fname = "background_rate_" + addRegion(h5file) + logStr
    fname = fname + year
    fname = fname.replace("/", "_").replace(" ", "_")
    print("Saving plot as ", fname + ".pdf")
    plt.savefig(fname + ".pdf")
    plt.show()

if __name__=="__main__":
    import sys
    main(sys.argv[1:])
