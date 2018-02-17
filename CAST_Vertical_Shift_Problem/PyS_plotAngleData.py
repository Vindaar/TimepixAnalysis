#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
from datetime import datetime

def convert_datetime_str_to_datetime(datetime_str):
    # convert datetime string to datetime object

    syntax = "%H:%M:%S"
    date = datetime.strptime(datetime_str, syntax)
    return date

def read_angle_file(fname):
    f = open(fname, "r").readlines()

    cut_tracking = True
    
    # indices of each column
    h_ae_ind = 1
    v_ae_ind = 2
    h_me_ind = 3
    v_me_ind = 4
    t_l_ind = 0

    # empty lists for data
    h_ae = []
    v_ae = []
    h_me = []
    v_me = []
    t_l = []
    
    t0 = convert_datetime_str_to_datetime(f[1].split()[-1])
    for i, line in enumerate(f):
        line = line.split()
        if i > 0 and len(line) > 0:
            t = convert_datetime_str_to_datetime(line[-1])
            t_diff = (t - t0).total_seconds()
            if t_diff < 3 and t_diff > 0:
                # if less than 3 seconds between time, we are moving, since we
                # are updating every second
                if cut_tracking == True and t.hour <= 9:
                    h_ae.append(float(line[h_ae_ind]))
                    v_ae.append(float(line[v_ae_ind]))
                    h_me.append(float(line[h_me_ind]))
                    v_me.append(float(line[v_me_ind]))
                    t_l.append(float(line[t_l_ind]))
                elif cut_tracking == False:
                    h_ae.append(float(line[h_ae_ind]))
                    v_ae.append(float(line[v_ae_ind]))
                    h_me.append(float(line[h_me_ind]))
                    v_me.append(float(line[v_me_ind]))
                    t_l.append(float(line[t_l_ind]))                    
            t0 = t
            
    # create numpy arrays
    print "Elements : ",np.size(v_ae)
    h_ae = np.asarray(h_ae)
    v_ae = np.asarray(v_ae)
    h_me = np.asarray(h_me)
    v_me = np.asarray(v_me)
    t_l = np.asarray(t_l)
    
    return h_ae, v_ae, h_me, v_me, t_l

def plot_h_vs_v(h_ae, v_ae, h_me, v_me):
    # rescale m encoder values
    h_me = h_me - (max(h_me) - min(h_me)) / 2.0
    v_me = v_me - (max(v_me) - min(v_me)) / 2.0

    h_me = h_me / max(h_ae)
    v_me = v_me / max(v_ae)

    print np.max(h_me), np.max(h_ae)
    print np.min(h_me), np.min(h_ae)
            
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

def plot_vae_vs_vme(v_ae, v_me):

    print "Plotting vertical ME vs vertial AE"
    plt.plot(v_me, v_ae)
    x_lin, y_lin = fit_linear(v_ae, v_me)
    plt.plot(x_lin, y_lin, color = 'red')
    
    plt.show()

def plot_fit_diff(v_ae, v_me):
    x_lin, y_lin = fit_linear(v_ae, v_me)

    y_diff = v_ae - y_lin

    y_diff = 10000 * np.tan(y_diff * np.pi / 180.0)
    
    plt.plot(x_lin, y_diff)
    plt.show()

def plot_n_n1_diff(v_ae, v_me, t_l):

    n_el = np.size(v_ae) - 1
    v_ae_diff = np.zeros(n_el)
    v_me_diff = np.zeros(n_el)
    t_l_diff = np.zeros(n_el)
    for i in xrange(n_el):
        v_ae_diff[i] = v_ae[i+1] - v_ae[i]
        v_me_diff[i] = v_me[i+1] - v_me[i]
        t_l_diff[i] = t_l[i+1] - t_l[i]

    print "Encoder deviation larger -30: ", np.where(v_me_diff < -30)[0]

    #v_me_diff /= np.max(v_me_diff)
    #v_ae_diff /= np.max(v_ae_diff)

    plt.plot(np.arange(np.size(v_ae_diff)), v_ae_diff)
    plt.plot(np.arange(np.size(v_me_diff)), v_me_diff, color='red')
    plt.show()


    
    v_me_diff = v_me_diff[1000:-500]
    t_l_diff = t_l_diff[1000:-500]

    v_me_diff /= np.max(v_me_diff)
    #t_l_diff /= np.max(t_l_diff)
    
    plt.plot(np.arange(np.size(v_me_diff)), v_me_diff, color='red')
    plt.plot(np.arange(np.size(v_me_diff)), t_l_diff, color='blue')
    plt.show()
    
        

def main(args):

    fname = args[0]

    h_ae, v_ae, h_me, v_me, t_l = read_angle_file(fname)
    print v_ae
    print v_me

    # plot_h_vs_v(h_ae, v_ae, h_me, v_me)

    plot_vae_vs_vme(v_ae, v_me)

    plot_fit_diff(v_ae, v_me)

    plot_n_n1_diff(v_ae, v_me, t_l)
    

if __name__=="__main__":
    import sys
    main(sys.argv[1:])
