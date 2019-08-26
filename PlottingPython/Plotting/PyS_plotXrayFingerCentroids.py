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

def readXrayData(h5file):
    # read data
    #group_name = likelihoodBase()
    group_name = "/reconstruction/"
    chipNumber = 3
    center_x = readH5Data(h5file, group_name, chipNumber, ["centerX"])
    center_y = readH5Data(h5file, group_name, chipNumber, ["centerY"])

    print(len(center_x))
    print(len(center_y))
    centers = np.zeros((256, 256))
    x_inds = []
    y_inds = []
    for x, y in zip(center_x, center_y):
        scale = 256.0 / 14.0
        x_ind = 256 - x * scale
        y_ind = 256 - y * scale
        x_inds.append(x_ind)
        y_inds.append(y_ind)
        centers[int(y_ind)][int(x_ind)] += 1
    centroid_x = np.mean(x_inds)
    centroid_y = np.mean(y_inds)
    #print(centers)
    print("Centroid is at {} / {}".format(centroid_x, centroid_y))

    return centroid_x, centroid_y, centers

def addFakeData(centers, aroundVal, sigma = 0.3, posx = 27, posy = 27):
    """
    This func can be used to add some fake data to `centers` with a gaussian distr
    around `aroundVal`.
    This is used in order to highlight deficiencies in certain color maps, which vary
    too much in hue and are not perceptually uniform.
    """
    print("Around val: ", aroundVal)
    size = 150
    noise = np.sort(np.random.normal(aroundVal, sigma, 150))
    half = int(np.shape(noise)[0] / 2.0)
    noise[:half] = np.roll(noise, half)[:half]
    noise[half:] = np.flip(np.roll(noise, half))[half:]
    xvals = np.sort(np.random.normal(posx, 4, 150))
    yvals = np.random.normal(posy, 4, 150)
    print("vals ", xvals)
    print("vals ", yvals)
    res = centers
    count = 0
    for x, y in zip(xvals, yvals):
        if noise[count] > 0:
            res[int(np.round(x)), int(np.round(y))] = noise[count]
        count += 1
    return res

def main(args):
    parser = argparse.ArgumentParser(description = 'H5 Data plotter')
    parser.add_argument('file',
                        help = "The H5 file from which to read data")
    parser.add_argument('--file2',
                        default = None,
                        dest = "file2",
                        help = "A potential second file to plot against")
    parser.add_argument('--cmap',
                        default = 'viridis',
                        dest = "cmap",
                        help = "The colormap to use, default viridis")
    parser.add_argument('--stack',
                        action = "store_true",
                        help = "If toggled, the occupancies are stacked")
    parser.add_argument('--centroid',
                        action = "store_true",
                        help = "If toggled, the centroid of all events is shown")
    parser.add_argument('--fakeData',
                        action = "store_true",
                        dest = "fakeData",
                        help = "If toggled, will add some fake noise to the plot")
    parser.add_argument('--noLabels',
                        action = "store_true",
                        dest = "noLabels",
                        help = "If toggled, won't show any labels, axes or titles")
    parser.add_argument('--maxval',
                        default = 6.0,
                        dest = 'maxval',
                        help = "Use to set the maximum colormap value")
    parser.add_argument('--fancy',
                        default = False,
                        action = 'store_true',
                        help = "Flag to activate fancy plotting via LaTeX output")

    args_dict = vars(parser.parse_args())
    h5file = os.path.abspath(args_dict["file"])
    h5file2 = args_dict["file2"]
    cmap = args_dict["cmap"]
    toAddFakeData = args_dict["fakeData"]
    noLabels = args_dict["noLabels"]
    stack = args_dict["stack"]
    maxval = args_dict["maxval"]
    show_centroid = args_dict["centroid"]
    print(args_dict)

    if args_dict["fancy"] == True:
        fancy_plotting()

    centroid_x, centroid_y, centers = readXrayData(h5file)
    if h5file2 != None:
        centroid_x2, centroid_y2, centers2 = readXrayData(h5file2)
        if stack == True:
            centers += centers2

    circle = Circle((centroid_x, centroid_y), 5.0, color = 'red')
    print("Max values of data: ", np.max(centers))
    ax = None
    if h5file2 != None and stack == True or h5file2 == None:
        fig, ax = plt.subplots()
        if toAddFakeData:
            nonEmpty = np.nonzero(centers)
            centers = addFakeData(centers, 1.4, posy = 87, posx = 125)
            centers = addFakeData(centers, 0.2, posx = 221, posy = 221)
        ax = ax.imshow(centers, cmap = cmap, vmax = maxval)
        if show_centroid:
            ax.add_artist(circle)
        if h5file2 != None and not noLabels:
            ax.set_title("Stacked occupancies of 2017/2018")
        if h5file2 != None:
            circle2 = Circle((centroid_x2, centroid_y2), 5.0, color = 'red')
            if show_centroid:
                ax.add_artist(circle2)
        print("Adding both")

        if noLabels:
            plt.gca().get_xaxis().set_visible(False)
            plt.gca().get_yaxis().set_visible(False)

    elif h5file2 != None and stack == False:
        fig, ax = plt.subplots(1, 2)
        circle2 = Circle((centroid_x2, centroid_y2), 5.0, color = 'red')

        ax[0].imshow(centers, cmap = cmap, vmax = 4)
        ax[1].imshow(centers2, cmap = viridis, vmax = 4)
        if show_centroid:
            ax[0].add_artist(circle)
            ax[1].add_artist(circle2)
        if not noLabels:
            ax[0].set_title(os.path.basename(h5file))
            ax[1].set_title(os.path.basename(h5file2))
        else:
            plt.gca().get_xaxis().set_visible(False)
            plt.gca().get_yaxis().set_visible(False)

    if h5file2 == None and not noLabels:
        plt.title(os.path.basename(h5file))


    plt.savefig("xray_finger_clusters_" + str(cmap) + "_max_" + str(maxval) + ".pdf")
    fig.colorbar(ax)

    plt.show()





if __name__=="__main__":
    import sys
    main(sys.argv[1:])
