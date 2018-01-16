from ingrid.ingrid_helper_functions import *
from scipy.optimize import curve_fit
import numpy as np

# this whole fit is probably going to be pretty slow.
# convert to Nim functions, which we compile to lib and call
# in python via ctypes
# see: https://akehrer.github.io/posts/connecting-nim-to-python/
# given the good starting parameters and that we only have to fit
# /once/ per calibration run, this is not a priorty for now, although
# it would be fun :)
def gauss(x, mean, sigma, norm = False):
    # based on the ROOT implementation of TMath::Gaus:
    # https://root.cern.ch/root/html524/src/TMath.cxx.html#dKZ4iB
    if sigma == 0:
        return 1.e30
    arg = (x-mean)/sigma;
    res = np.exp(-0.5 * arg * arg)
    if norm == False:
        return res
    else:
        return res / (2.50662827463100024 * sigma) # sqrt(2*Pi)=2.5066282746310002

def expGauss(x_ar, p):
    # exponential * times (?!) from Christoph's expogaus.c
    if len(p) != 5:
        return np.nan
    p_val = 2 * (p[1] * np.power(p[4], 2) - p[3])
    q_val = 2 * np.power(p[4], 2) * p[0] + np.power(p[3], 2) - np.log(p[2]) * 2 * np.power(p[4], 2)

    threshold = - p_val / 2.0 - np.sqrt( np.power(p_val, 2) / 4.0 - q_val )

    x_out = np.array(x_ar, copy = True)
    for i, x in enumerate(x_ar):
        if x < threshold:
            x_out[i] = np.exp(p[0] + p[1] * x)
        else:
            x_out[i] = p[2] * gauss(x, p[3], p[4])

    return x_out

def feSpectrumFunc(x, *p_ar):#p_ar):
    # NOTE: should we init with 0? or 1?

    parAll = np.zeros(20, dtype = np.float)
    if len(p_ar) != 15:
        print("Lenght of params is ", len(p_ar))
        import sys
        sys.exit("Bad parameter set! Sorry for killing the progam :)")

    # while this is the actual fitting function, it is in fact only
    # a wrapper around 4 Gaussians * exponential decay
    # It's a mixture of the K_alpha, K_beta lines as well as
    # the escape peaks
    # see fitFeSpectrum() (where bounds are set) for discussion of different
    # parameters
    parAll[0] = p_ar[0]
    parAll[1] = p_ar[1]
    parAll[2] = p_ar[2]
    parAll[3] = p_ar[3]
    parAll[4] = p_ar[4]

    parAll[5] = p_ar[5]
    parAll[6] = p_ar[6]
    parAll[7] = p_ar[14]*p_ar[2]
    parAll[8] = 3.5/2.9*p_ar[3]
    parAll[9] = p_ar[4]

    parAll[10] = p_ar[7]
    parAll[11] = p_ar[8]
    parAll[12] = p_ar[9]
    parAll[13] = p_ar[10]
    parAll[14] = p_ar[11]

    parAll[15] = p_ar[12]
    parAll[16] = p_ar[13]
    parAll[17] = p_ar[14]*p_ar[9]
    parAll[18] = 6.35/5.75*p_ar[10]
    parAll[19] = p_ar[11]

    value = 0
    value += expGauss(x, parAll[0:5])
    value += expGauss(x, parAll[5:10])
    value += expGauss(x, parAll[10:15])
    value += expGauss(x, parAll[15:20])

    return value

def fitFeSpectrum(hist, binning, cuts):
    # given our histogram and binning data
    # for fit := (y / x) data
    # fit a double gaussian to the data

    # define center, std and amplitude of K_alpha line
    # as well as escape peak
    mu_kalpha = np.argmax(hist)
    sigma_kalpha = mu_kalpha / 10.0
    n_kalpha = hist[mu_kalpha]
    print mu_kalpha
    print sigma_kalpha
    print n_kalpha
    mu_kalpha_esc = mu_kalpha * 2.9/5.75
    sigma_kalpha_esc = mu_kalpha_esc / 10.0
    n_kalpha_esc = n_kalpha / 10.0

    params = np.zeros(15, dtype = np.float)
    params[2] = n_kalpha_esc
    params[3] = mu_kalpha_esc
    params[4] = sigma_kalpha_esc

    params[9] = n_kalpha
    params[10] = mu_kalpha
    params[11] = sigma_kalpha

    params[14] = 17.0 / 150.0

    l_bounds = []
    u_bounds = []
    for i in xrange(15):
        l_bounds.append(-np.inf)
        u_bounds.append(np.inf)
    # set bound on paramerters
    # constrain amplitude of K_beta to some positive value    
    l_bounds[2] = 0
    u_bounds[2] = 10000
    # constrain amplitude of K_alpha to some positive value
    l_bounds[9] = 0
    u_bounds[9] = 10000

    # location of K_alpha escape peak, little more than half of K_alpha location
    l_bounds[3] = mu_kalpha_esc * 0.8
    u_bounds[3] = mu_kalpha_esc * 1.2

    # bounds for center of K_alpha peak (should be at around 220 electrons, hits)
    l_bounds[10] = mu_kalpha*0.8
    u_bounds[10] = mu_kalpha*1.2

    # some useful bounds for K_alpha escape peak width
    l_bounds[4] = sigma_kalpha_esc*0.5
    u_bounds[4] = sigma_kalpha_esc*1.5

    # some useful bounds for K_alpha width
    l_bounds[11] = sigma_kalpha*0.5
    u_bounds[11] = sigma_kalpha*1.5
    # param 14: "N_{K_{#beta}}/N_{K_{#alpha}}"
    # known ratio of two K_alpha and K_beta, should be in some range
    l_bounds[14] = 0.01
    u_bounds[14] = 0.3

    # this leaves parameter 7 and 8, as well as 12 and 13 without bounds
    # these describe the exponential factors contributing, since they will
    # be small anyways...
    bounds = (l_bounds, u_bounds)
    print(len(bounds))
    print bounds
    print(len(params))

    # only fit in range up to 350 hits. Can take index 350 on both, since we
    # created the histogram for a binning with width == 1 pixel per hit
    data_tofit = hist[0:350]
    bins_tofit = binning[0:350]
    #print data_tofit
    #print bins_tofit
    lb, ub = [np.asarray(b, dtype=float) for b in bounds]
    print ub
    
    result = curve_fit(feSpectrumFunc, bins_tofit, data_tofit, p0=params, bounds = bounds)#, full_output=True)
    popt = result[0]
    pcov = result[1]
    # Calculate the reduced Chi^2:
    # n_dof: # degrees of freedom
    #n_dof  = (np.size(infodict[0]['fvec']) - np.size(popt))
    #chi_sq = np.sum(infodict[0]['fvec']**2) / n_dof
    print '--------------------------------------------------'
    print 'Parameters of calibration fit: '
    for i, p in enumerate(popt):
        print 'p_{} ='.format(i), popt[i], '+-', np.sqrt(pcov[i][i])
    #print 'Chi^2 / dof =', chi_sq

    print "yay :)"

    return popt, pcov

def linear_func(x, a):
    return a * x

def fitAndPlotFeSpectrum(data, cuts, outfolder, run_number, fitting_only = False):

    # bin the data
    hist, binning = binData(data, cuts)

    popt, pcov = fitFeSpectrum(hist, binning, cuts)

    # now get some nice x / y points from the fit parameters
    # given these values, plot
    x_pl = np.linspace(0, 350, 3500)
    y_pl = feSpectrumFunc(x_pl, *popt)

    fig, ax = plotData(hist, binning, "", "Fe spectrum", "\\# pixels hit", "\\# events", False)

    # TODO: add fit parameter results to plots as legend!

    ax.set_xlim(0, 350)
    ax.plot(x_pl, y_pl, linewidth = 2, color = "red")
    plt.savefig(os.path.join(outfolder, "fe_spectrum_{}.pdf".format(run_number)))
    if fitting_only == False:
        plt.show()

    # now we can fit the energy calibration function
    energies = [2.925, 5.755]
    pixels_peaks = [popt[3], popt[10]]
    pixels_err = [np.sqrt(pcov[3][3]), np.sqrt(pcov[10][10])]
    result = curve_fit(linear_func, energies, pixels_peaks, sigma = pixels_err, full_output=True)
    popt_E = result[0]
    pcov_E = result[1]
    infodict = result[2:]
    # n_dof: # degrees of freedom
    n_dof  = (np.size(infodict[0]['fvec']) - np.size(popt_E))
    chi_sq = np.sum(infodict[0]['fvec']**2) / n_dof
    print '--------------------------------------------------'
    print 'Parameters of calibration fit: '
    for i, p in enumerate(popt_E):
        print 'p_{} ='.format(i), popt_E[i], '+-', np.sqrt(pcov_E[i][i])
    a_inv = 1.0 / popt_E[0] * 1000
    da_inv = a_inv * np.sqrt(pcov_E[0][0]) / popt_E[0]
    print("a^-1 = {0} +- {1}".format(a_inv, da_inv))
    print 'Chi^2 / dof =', chi_sq
                                       
    E_calc = np.linspace(2, 7, 1000)
    H_calc = linear_func(E_calc, popt_E[0])
    plt.errorbar(energies, pixels_peaks, yerr = pixels_err, marker = ".", markersize = 8, linestyle = "", color = "red")
    plt.plot(E_calc, H_calc, marker = "", linestyle = "-")
    plt.xlabel("Energy / keV")
    plt.ylabel("# pixels hit")
    plt.title("Energy calibration function based on # pix in Fe55 spectrum")
    plt.grid()
    plt.savefig(os.path.join(outfolder, "fe_energy_calibration_{}.pdf".format(run_number)))
    if fitting_only == False:
        plt.show()

    # return fit results so that we can write them to the H5 file
    return (popt, pcov, popt_E, pcov_E)
