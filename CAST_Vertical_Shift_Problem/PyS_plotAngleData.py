#!/usr/bin/env python

import argparse
import os
import numpy as np
import scipy.signal
import math
import matplotlib.pyplot as plt
import pylab
from time import mktime
from collections import OrderedDict
from scipy.optimize import curve_fit
from scipy.stats import linregress
from datetime import datetime

def congrid(a, newdims, method='linear', centre=False, minusone=False):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).
    
    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if not a.dtype in [np.float64, np.float32]:
        a = np.cast[float](a)

    m1 = np.cast[int](minusone)
    ofs = np.cast[int](centre) * 0.5
    old = np.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print "[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions."
        return None
    newdims = np.asarray( newdims, dtype=float )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = np.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = np.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = np.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [np.arange(i, dtype = np.float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + range( ndims - 1 )
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = np.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = np.mgrid[nslices]

        newcoords_dims = range(np.ndim(newcoords))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (np.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print "Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported."
        return None

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

def convert_datetime_str_to_datetime(datetime_str):
    # convert datetime string to datetime object
    syntax = "%H:%M:%S"
    date = datetime.strptime(datetime_str, syntax)
    return date

def read_prox_sens_from_sclog(fname):
    # reads proximity sensor values from a slow control file
    # returns a tuple of timestamps, proxJ, proxA

    f = open(fname, "r").readlines()
    t = []
    pJ = []
    pA = []
    # indices for prox sensors and time
    t_i  = 1
    pJ_i = 60
    pA_i = 61
    
    for i, line in enumerate(f):
        line = line.split()
        if i > 0 and len(line) > 0:
            t.append(line[t_i])
            pJ.append(line[pJ_i])
            pA.append(line[pA_i])

    t = [convert_datetime_str_to_datetime(x) for x in t]
    pJ = np.asarray(pJ, dtype=float)
    pA = np.asarray(pA, dtype=float)
            
    return t, pJ, pA

def append_line(lists, vals):
    # appends a single line to the lists
    for lst, val in zip(lists, vals):
        lst.append(val)

def read_angle_file(fname, args_dict):
    f = open(fname, "r").readlines()
    
    cut_tracking = True
    # indices of each column
    ind_dict = OrderedDict([("h_ae_ind", 1),
                            ("h_me_ind", 3),
                            ("v_ae_ind", 2),
                            ("v_me_ind", 4),
                            ("t_l_ind" , 0)])

    # empty lists for data
    h_ae = []
    v_ae = []
    h_me = []
    v_me = []
    t_l = []

    only_up = args_dict["only_up"]
    only_down = args_dict["only_down"]
    print("Only up is {} and only down is {}".format(only_up, only_down))
    
    t0 = convert_datetime_str_to_datetime(f[1].split()[-1])
    v_prev = 0
    for i, line in enumerate(f):
        line = line.split()
        if i > 0 and len(line) > 0:
            t = convert_datetime_str_to_datetime(line[-1])
            t_diff = (t - t0).total_seconds()
            if t_diff < 3 and t_diff > 0:
                # parse line
                val_lst = [float(line[ind]) for ind in ind_dict.values()]
                lst_lst = [h_ae, h_me, v_ae, v_me, t_l]
                if only_up == True and float(line[ind_dict["v_me_ind"]]) - v_prev > 0:
                    # if less than 3 seconds between time, we are moving, since we
                    # are updating every second
                    if cut_tracking == True and t.hour <= 9:
                        append_line(lst_lst, val_lst)
                        # h_ae.append(float(line[h_ae_ind]))
                        # v_ae.append(float(line[v_ae_ind]))
                        # h_me.append(float(line[h_me_ind]))
                        # v_me.append(float(line[v_me_ind]))
                        # t_l.append(float(line[t_l_ind]))
                elif only_down == True and float(line[ind_dict["v_me_ind"]]) - v_prev < 0:
                    if cut_tracking == True and t.hour <= 9:
                        append_line(lst_lst, val_lst)
                elif only_up == False and only_down == False:
                    if cut_tracking == True and t.hour <= 9:
                        append_line(lst_lst, val_lst)
                elif cut_tracking == False:
                    append_line(lst_lst, val_lst)
                    
                v_prev = float(line[ind_dict["v_me_ind"]])
            t0 = t
            
    # create numpy arrays
    print "Elements : ",np.size(v_ae)
    h_ae = np.asarray(h_ae)
    v_ae = np.asarray(v_ae)
    h_me = np.asarray(h_me)
    v_me = np.asarray(v_me)
    t_l = np.asarray([datetime.fromtimestamp(x) for x in t_l])
    
    return h_ae, v_ae, h_me, v_me, t_l

def plot_h_vs_v(h_ae, v_ae, h_me, v_me, show = False):
    # rescale m encoder values
    h_me = h_me - (max(h_me) - min(h_me)) / 2.0
    v_me = v_me - (max(v_me) - min(v_me)) / 2.0

    h_me = h_me / max(h_ae)
    v_me = v_me / max(v_ae)

    print np.max(h_me), np.max(h_ae)
    print np.min(h_me), np.min(h_ae)

    if show == True:
        plt.plot(h_ae, v_ae)
        plt.plot(h_me, v_me)        
        plt.show()

def linear(x, a, m):
    return a + m * x

def fit_linear(v_ae, v_me):    
    slope, intercept, r_v, p_v, std_err = linregress(v_me, v_ae)
    print "Slope is : ", slope
    print "Intersection is : ", intercept

    #x_lin = np.linspace(3300, 53000, 1000)
    y_lin = linear(v_me, intercept, slope)
    #return x_lin, y_lin
    return v_me, y_lin

def plot_vae_vs_vme(v_ae, v_me, show = False):

    print "Plotting vertical ME vs vertial AE"
    x_lin, y_lin = fit_linear(v_ae, v_me)
    if show == True:
        plt.plot(v_me, v_ae)
        plt.plot(x_lin, y_lin, color = 'red')
        plt.show()

def plot_fit_diff(v_ae, v_me, title = "", show = True, from_prox = False):
    x_lin, y_lin = fit_linear(v_ae, v_me)

    y_diff = v_ae - y_lin

    y_diff = 10000 * np.tan(y_diff * np.pi / 180.0)

    #y_diff /= np.max(y_diff)
    y_diff /= np.percentile(y_diff, 95)

    print title
    l_name = os.path.basename(title).lstrip("Angles_data_.")
    if from_prox == False:
        plt.plot(x_lin, y_diff, label = l_name)
        plt.title(title)
    else:
        plt.plot(x_lin, y_diff, label = "V_ae")
        plt.title(l_name)
    plt.xlabel("Vertical ME units")
    plt.ylabel("Diff Vertical AE vals - Vertical AE fit")
    if show == True:    
        plt.show()
    #else:
    #    plt.clf()

def get_n_n1_data(v_ae, v_me, t_l):
    n_el = np.size(v_ae) - 1
    v_ae_diff = np.zeros(n_el)
    v_me_diff = np.zeros(n_el)
    t_l_diff = np.zeros(n_el)
    for i in xrange(n_el):
        v_ae_diff[i] = v_ae[i+1] - v_ae[i]
        v_me_diff[i] = v_me[i+1] - v_me[i]
        t_l_diff[i] = (t_l[i+1] - t_l[i]).total_seconds()

    print "Encoder deviation larger -30: ", np.where(v_me_diff < -30)[0]

    # v_me_diff /= np.max(v_me_diff)
    # v_ae_diff /= np.max(v_ae_diff)
    v_me_diff /= np.percentile(v_me_diff, 95)
    v_ae_diff /= np.percentile(v_ae_diff, 95)

    return v_ae_diff, v_me_diff, t_l_diff


def plot_n_n1_diff(v_ae, v_me, t_l, show = False):

    # get the data 
    v_ae_diff, v_me_diff, t_l_diff = get_n_n1_data(v_ae, v_me, t_l)

    if show == True:
        plt.plot(np.arange(np.size(v_ae_diff)), v_ae_diff)
        plt.plot(np.arange(np.size(v_me_diff)), v_me_diff, color='red')        
        plt.show()

    #v_me_diff = v_me_diff[1000:-500]
    #t_l_diff = t_l_diff[1000:-500]

    #v_me_diff /= np.max(v_me_diff)
    #t_l_diff /= np.max(t_l_diff)
    
    if show == True:
        plt.plot(np.arange(np.size(v_me_diff)), v_me_diff, color='red')
        plt.plot(np.arange(np.size(v_me_diff)), t_l_diff, color='blue')
        plt.show()

def processProxVme(fname, prox_fname, args_dict):
    only_up = args_dict["only_up"]
    only_down = args_dict["only_down"]
    
    
    # processes proximity versus V_me plot
    t, pJ, pA = read_prox_sens_from_sclog(prox_fname)

    # now also need angle encoder data
    h_ae, v_ae, h_me, v_me, t_l = read_angle_file(fname, args_dict)

    # since sc logs are 1 minute apart each, whereas angle```s file
    # happens every ~second once tracking starts, need to bin the
    # angle data accordingly. Questionable whether this is good enough
    # anyways? Once per minute is not a lot

    # get difference between 1904 timestamp and 1970 timestamp
    # for correction of t_l times
    dt_base = datetime(1904, 1, 1, 0, 0, 0)
    dt = datetime(1970, 1, 1, 0, 0, 0)
    diff_time = (dt - dt_base)

    # apply correctionto get datetime object of correct dates
    t_l = np.asarray([(x - diff_time) for x in t_l])


    t_l_date = t_l[0].date()
    # use date of t_l and replace its date 
    t = [(x.replace(year = t_l_date.year, month = t_l_date.month, day = t_l_date.day)) for x in t]
    t = np.asarray(t)

    # cut away everything of pJ, pA and t that is not in v_ae, v_me and t_l
    # get indices
    # this only works properly is only_up and only_down are both False
    if only_up == True or only_down == True:
        # in this case cut to interesting parts..
        # given rather cut on v_me
        v_me_HIGH = 51000
        v_me_LOW = 27000
        v_me_inds = np.where(np.logical_and(v_me <= v_me_HIGH, v_me >= v_me_LOW))[0]
        print v_me_inds
        v_me = v_me[v_me_inds]
        v_ae = v_ae[v_me_inds]
        # now use this to cut on time
        t_l = t_l[v_me_inds]

    # if only up or down set, we have cut to specific ME values, can now perform normal
    # time cut
    t_inds = np.where(np.logical_and([x >= t_l[0] for x in t], [x <= t_l[-1] for x in t]))
    # cut
    t = t[t_inds]
    pJ = pJ[t_inds]
    pA = pA[t_inds]

    print len(pJ), len(t_l)
    # pJ_norm = (pJ - np.mean(pJ)) / np.max(pJ - np.mean(pJ))
    # pA_norm = (pA - np.mean(pA)) / np.max(pA - np.mean(pA))
    pJ_norm = (pJ - np.mean(pJ)) / np.percentile(pJ - np.mean(pJ), 95)
    pA_norm = (pA - np.mean(pA)) / np.percentile(pA - np.mean(pA), 95)        

    print t
    plt.plot(t, pJ_norm, label = "prox Jura")
    plt.plot(t, pA_norm, label = "prox Airport")
    # get n_n1 data as comparison
    v_ae_diff, v_me_diff, t_l_diff = get_n_n1_data(v_ae, v_me, t_l)
    print t_l.size, v_me_diff.size
    plt.plot(t_l[1:], v_me_diff, label = "V_me")
    plt.plot(t_l[1:], v_ae_diff, label = "V_ae")
    plt.legend(loc = 4)
    plt.show()

    if only_up == False and only_down == False:
        import sys
        sys.exit("Warning: plotting the whole range (up and down) against the vertical encoders" +
                 "does not work, since slow control sometimes updates every second, which breaks" +
                 "possiblity to interpolate properly (without big effort).")

    # now bin the data from angle encoders
    # congrid fails us, because the target space is non uniform in distance :/
    # sometimes the slow control does indeed update every second instead of
    # every minute... (makes sense considering that pressures are sometimes updated
    # that often)
    v_me_b = congrid(v_me, t.shape, 'linear')
    v_ae_b = congrid(v_ae, t.shape, 'linear')
    # plt.plot(pJ - np.mean(pJ), v_me_b)
    # plt.plot(pA - np.mean(pA), v_me_b)
    #plt.plot(v_me_b, pJ - np.mean(pJ))    
    #plt.plot(v_me_b, pA - np.mean(pA))

    # linear space suffers from the same problem as congrid
    v_me_lin = np.linspace(v_me[0], v_me[-1], t.size)
    
    # calculation by linear interpolation based on real data
    # allows us to deal with cases where slow control updates per second
    # however, problem if we use neither only_up nor only_down!
    t_l_s = [mktime(x.timetuple()) for x in t_l]
    t_s = [mktime(x.timetuple()) for x in t]    
    f_vme = scipy.interpolate.interp1d(t_l_s, v_me)
    v_me_f = f_vme(t_s)
    
    # now bin the data from angle encoders
    plt.plot(v_me_f, pJ_norm, label = "prox Jura")
    plt.plot(v_me_f, pA_norm, label = "prox Airport")
    plt.title("Proximity sensors vs vertical motor encoder values")
    title_add = ""
    if only_up == True:
        title_add = " moving upwards"
    elif only_down == True:
        title_add = " moving downwards"
    
    plot_fit_diff(v_ae, v_me, fname + title_add, show = False, from_prox = True)
    plt.legend(loc = 4)
    plt.show()
    

def processLogFile(fname, args_dict, single = True):
    h_ae, v_ae, h_me, v_me, t_l = read_angle_file(fname, args_dict)
    # plot_h_vs_v(h_ae, v_ae, h_me, v_me)
    plot_vae_vs_vme(v_ae, v_me, show = False)

    title = fname
    plot_fit_diff(v_ae, v_me, title, show = single)

    if args_dict["n_n1"] == True:
        plot_n_n1_diff(v_ae, v_me, t_l, show = True)
    else:
        plot_n_n1_diff(v_ae, v_me, t_l, show = False)

def main(args):

    parser = argparse.ArgumentParser(description = 'H5 Data plotter')
    parser.add_argument('file',
                        help = "The H5 file from which to read data")
    parser.add_argument('--n_n1',
                        default = False,
                        action = 'store_true',
                        help = "Flag to activate plot of N - N1 diff")
    parser.add_argument('--prox_ME',
                        default = None,
                        dest = 'prox_log',
                        help = "Flag to only create plot of proximity sensor against vertical motor encoder")
    parser.add_argument('--only_up',
                        default = False,
                        action = 'store_true',
                        help = "Flag to only read tracking upwards")
    parser.add_argument('--only_down',
                        default = False,
                        action = 'store_true',
                        help = "Flag to only read tracking downwards")    
    parser.add_argument('--fancy',
                        default = False,
                        action = 'store_true',
                        help = "Flag to activate fancy plotting via LaTeX output")
    args_dict = vars(parser.parse_args())
    if args_dict["fancy"] == True:
        fancy_plotting()

    fname = args_dict["file"]
    v_ae_dict = {}
    v_me_dict = {}
    if os.path.isfile(fname):
        prox_fname = args_dict["prox_log"]
        if prox_fname == None:
            processLogFile(fname, args_dict)
        else:
            processProxVme(fname, prox_fname, args_dict)
    else:
        for logf in os.listdir(fname):
            processLogFile(os.path.join(fname, logf), args_dict, False)
        plt.legend(loc = 1)
        plt.title("Combined V_ae_corr vs V_me")
        plt.show()

    

if __name__=="__main__":
    import sys
    main(sys.argv[1:])
