from ingrid.ingrid_helper_functions import *
from scipy.optimize import curve_fit
import numpy as np
import matplotlib as mpl

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

def feSpectrumFuncCharge(x, *p_ar):
    ## the fitting function used for the Fe spectrum with charge bins, see `XrayCalib.c`
    # fix N_Kbeta/N_Kalpha to theoretical value from XDB
    p6 = 17.0 / 150.0
    t1 = p_ar[0] * np.exp(-(x - p_ar[1]) * (x - p_ar[1]) / (2 * p_ar[2] * p_ar[2]))
    t2exp = -(x - (p_ar[1] * 3.53 / 2.94)) * (x - (p_ar[1] * 3.53/2.94)) / (2 * p_ar[2] * p_ar[2])
    t2 = p_ar[0] * p6 * np.exp(t2exp)
    t3 = p_ar[3] * np.exp(-(x - p_ar[4]) * (x - p_ar[4]) / (2 * p_ar[5] * p_ar[5]))
    t4exp = -(x - (p_ar[4] * 6.49 / 5.90)) * (x - (p_ar[4] * 6.49 / 5.90)) / (2 * p_ar[5] * p_ar[5])
    t4 = p_ar[3] * p6 * np.exp(t4exp)
    return (t1 + t2 + t3 + t4)

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

def getLines(hist):
    # define center, std and amplitude of K_alpha line
    # as well as escape peak
    mu_kalpha = np.argmax(hist)
    sigma_kalpha = mu_kalpha / 10.0
    n_kalpha = hist[mu_kalpha]

    mu_kalpha_esc = mu_kalpha * 2.9/5.75
    sigma_kalpha_esc = mu_kalpha_esc / 10.0
    n_kalpha_esc = n_kalpha / 10.0
    return mu_kalpha, sigma_kalpha, n_kalpha, mu_kalpha_esc, sigma_kalpha_esc, n_kalpha_esc

def getBoundsList(nParams):
    l_bounds = []
    u_bounds = []
    for i in range(nParams):
        l_bounds.append(-np.inf)
        u_bounds.append(np.inf)
    return l_bounds, u_bounds

def fitFeSpectrumToCharge(hist, binning, cuts):
    # perform the fit of the Fe spectrum to the charge data instead
    # of the number of hit pixels data
    mu_kalpha, sigma_kalpha, n_kalpha, mu_kalpha_esc, sigma_kalpha_esc, n_kalpha_esc = getLines(hist)
    params = np.zeros(6, dtype = np.float)
    # parameter names
    # 0 : "N^{esc}_{K_{#alpha}}"
    # 1 : "#mu^{esc}_{K_{#alpha}} [10^{3} e]"
    # 2 : "#sigma^{esc}_{K_{#alpha}} [10^{3} e]"
    # 3 : "N_{K_{#alpha}}"
    # 4 : "#mu_{K_{#alpha}} [10^{3} e]"
    # 5 : "#sigma_{K_{#alpha}} [10^{3} e]"
    # 6 : "N_{K_{#beta}}/N_{K_{#alpha}}"
    params[0] = n_kalpha_esc
    params[1] = mu_kalpha_esc
    params[2] = sigma_kalpha_esc

    params[3] = n_kalpha
    params[4] = mu_kalpha
    params[5] = sigma_kalpha
    # get lists defining the bounds of the parameters
    l_bounds, u_bounds = getBoundsList(6)
    # set the bounds
    l_bounds[0] = 0
    u_bounds[0] = 10000
    l_bounds[3] = 0
    u_bounds[3] = 10000

    l_bounds[1] = mu_kalpha_esc*0.8
    u_bounds[1] = mu_kalpha_esc*1.2
    l_bounds[4] = mu_kalpha*0.8
    u_bounds[4] = mu_kalpha*1.2

    l_bounds[2] = sigma_kalpha_esc*0.5
    u_bounds[2] = sigma_kalpha_esc*1.5
    l_bounds[5] = sigma_kalpha*0.5
    u_bounds[5] = sigma_kalpha*1.5

    # combine bounds
    bounds = (l_bounds, u_bounds)
    #print(len(bounds))
    print("Bounds for charge fit: ", bounds)
    print("N params for charge fit: ", len(params))
    # only fit in range up to 350 hits. Can take index 350 on both, since we
    # created the histogram for a binning with width == 1 pixel per hit
    inds = np.where(np.logical_and(binning > 200, binning < 2000))[0]
    data_tofit = hist[inds]
    #print("Data to fit: ", data_tofit)
    bins_tofit = binning[inds]
    print("Binning to fit: ", bins_tofit)
    print("l bounds : ", l_bounds)
    #lb, ub = [np.asarray(b, dtype=float) for b in bounds]
    print("Bounds are : ", bounds)
    result = curve_fit(feSpectrumFuncCharge, bins_tofit, data_tofit, p0=params, bounds = bounds)#, full_output=True)
    popt = result[0]
    pcov = result[1]
    # Calculate the reduced Chi^2:
    # n_dof: # degrees of freedom
    #n_dof  = (np.size(infodict[0]['fvec']) - np.size(popt))
    #chi_sq = np.sum(infodict[0]['fvec']**2) / n_dof
    print('--------------------------------------------------')
    print('Parameters of calibration fit: ')
    for i, p in enumerate(popt):
        print('p_{} ='.format(i), popt[i], '+-', np.sqrt(pcov[i][i]))
    #print('Chi^2 / dof =', chi_sq)

    print("yay :)")

    return popt, pcov

def fitFeSpectrum(hist, binning, cuts):
    # given our histogram and binning data
    # for fit := (y / x) data
    # fit a double gaussian to the data
    mu_kalpha, sigma_kalpha, n_kalpha, mu_kalpha_esc, sigma_kalpha_esc, n_kalpha_esc = getLines(hist)

    params = np.zeros(15, dtype = np.float)
    params[2] = n_kalpha_esc
    params[3] = mu_kalpha_esc
    params[4] = sigma_kalpha_esc

    params[9] = n_kalpha
    params[10] = mu_kalpha
    params[11] = sigma_kalpha

    params[14] = 17.0 / 150.0

    l_bounds, u_bounds = getBoundsList(15)
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
    print("Bounds for pixel fit: ", bounds)
    print("N params for pixel fit: ", len(params))

    # only fit in range up to 350 hits. Can take index 350 on both, since we
    # created the histogram for a binning with width == 1 pixel per hit
    data_tofit = hist[0:350]
    bins_tofit = binning[0:350]

    result = curve_fit(feSpectrumFunc, bins_tofit, data_tofit, p0=params, bounds = bounds)#, full_output=True)
    popt = result[0]
    pcov = result[1]
    # Calculate the reduced Chi^2:
    # n_dof: # degrees of freedom
    #n_dof  = (np.size(infodict[0]['fvec']) - np.size(popt))
    #chi_sq = np.sum(infodict[0]['fvec']**2) / n_dof
    print('--------------------------------------------------')
    print('Parameters of calibration fit: ')
    for i, p in enumerate(popt):
        print('p_{} ='.format(i), popt[i], '+-', np.sqrt(pcov[i][i]))
    #print('Chi^2 / dof =', chi_sq)

    print("yay :)")

    return popt, pcov

def linear_func(x, a):
    return a * x

def fitAndPlotFeSpectrumCharge(data, cuts, outfolder, run_number, fitting_only = False):
    # bin the data
    # before we bin data, divide it by 1000
    data = np.asarray(data) / 1000
    hist, binning = binData(data, cuts)
    popt, pcov = fitFeSpectrumToCharge(hist, binning, cuts)

    # now get some nice x / y points from the fit parameters
    # given these values, plot
    x_pl = np.linspace(0, max(binning), 3500)
    y_pl = feSpectrumFuncCharge(x_pl, *popt)
    # create both plots
    fig, ax = plotData(hist, binning, None, "", "Fe spectrum", "charge / e-", "\\# events", False)
    # TODO: add fit parameter results to plots as legend!
    #ax.set_xlim(0, )
    ax.plot(x_pl, y_pl, linewidth = 2, color = "red")
    plt.show()

    return popt

def fitAndPlotFeSpectrumPixels(data, cuts, outfolder, run_number, fitting_only = False):
    # bin the data
    hist, binning = binData(data, cuts)
    popt, pcov = fitFeSpectrumToPixels(hist, binning, cuts)
    print("POPT ", popt)
    return popt

def fitAndPlotFeSpectrum(data, cuts, outfolder, run_number, fitting_only = False):
    # bin the data
    hist, binning = binData(data, cuts)
    popt, pcov = fitFeSpectrum(hist, binning, cuts)
    k_alpha = popt[10]

    # now get some nice x / y points from the fit parameters
    # given these values, plot
    x_pl = np.linspace(0, 350, 3500)
    y_pl = feSpectrumFunc(x_pl, *popt)
    # before we plot the Fe spectrum, perform the calculations and fit of
    # the fit to the spectrum peaks
        # now we can fit the energy calibration function
    energies = [2.925, 5.755]
    pixels_peaks = [popt[3], k_alpha]
    pixels_err = [np.sqrt(pcov[3][3]), np.sqrt(pcov[10][10])]
    result = curve_fit(linear_func, energies, pixels_peaks, sigma = pixels_err, full_output=True)
    popt_E = result[0]
    pcov_E = result[1]
    infodict = result[2:]
    # n_dof: # degrees of freedom
    n_dof  = (np.size(infodict[0]['fvec']) - np.size(popt_E))
    chi_sq = np.sum(infodict[0]['fvec']**2) / n_dof
    print('--------------------------------------------------')
    print('Parameters of linear fit to peaks: ')
    for i, p in enumerate(popt_E):
        print('p_{} ='.format(i), popt_E[i], '+-', np.sqrt(pcov_E[i][i]))

    a_inv = 1.0 / popt_E[0] * 1000
    da_inv = a_inv * np.sqrt(pcov_E[0][0]) / popt_E[0]
    print("a^-1 = {0} +- {1}".format(a_inv, da_inv))
    print('Chi^2 / dof =', chi_sq)

    E_calc = np.linspace(2, 7, 1000)
    H_calc = linear_func(E_calc, popt_E[0])

    # create both plots
    fig, ax = plotData(hist, binning, None, "", "Fe spectrum", "\\# pixels hit", "\\# events", False)

    # TODO: add fit parameter results to plots as legend!

    ax.set_xlim(0, 350)
    ax.plot(x_pl, y_pl, linewidth = 2, color = "red")

    if mpl.rcParams["text.usetex"] == True:
        text  = "$\mu = \SI{" + "{0:.1f}".format(k_alpha) + "}{pix}$"
        text2 = "$\sim\SI{" + "{0:.1f}".format(a_inv) + "}{\electronvolt \per pix}$"
    else:
        text  = "mu = " + "{0:.1f}".format(k_alpha) + "pix"
        text2 = "{0:.1f}".format(a_inv) + "ev / pix"
    ax.text(120, 150, text, fontsize = 20)
    ax.text(120, 135, text2, fontsize = 20)


    plt.savefig(os.path.join(outfolder, "fe_spectrum_{}.pdf".format(run_number)))
    if fitting_only == False:
        plt.show()
    else:
        # in that case clear the current figure to not end up with
        # both plots in one
        plt.clf()

    plt.errorbar(energies, pixels_peaks, yerr = pixels_err, marker = ".", markersize = 8, linestyle = "", color = "red")
    plt.plot(E_calc, H_calc, marker = "", linestyle = "-")
    plt.xlabel("Energy / keV")
    plt.ylabel("# pixels hit")
    plt.title("Energy calibration function based on \# pix in Fe55 spectrum")
    plt.grid()
    plt.savefig(os.path.join(outfolder, "fe_energy_calibration_{}.pdf".format(run_number)))
    if fitting_only == False:
        plt.show()
    else:
        plt.close()

    # return fit results so that we can write them to the H5 file
    return (popt, pcov, popt_E, pcov_E)










def feSpectrumFuncPixels(x, *p_ar):
    ## to be implemented
    pass

def fitFeSpectrumToPixels(hist, binning, cuts):
    # perform the fit of the Fe spectrum to the charge data instead
    # of the number of hit pixels data
    mu_kalpha, sigma_kalpha, n_kalpha, mu_kalpha_esc, sigma_kalpha_esc, n_kalpha_esc = getLines(hist)
    params = np.zeros(8, dtype = np.float)
    # parameter names are the following:
    # 0 : "N^{esc}_{K_{#alpha}}"
    # 1 : "#mu^{esc}_{K_{#alpha}}"
    # 2 : "#sigma^{esc}_{K_{#alpha}}"
    # 3 : "N_{K_{#beta}}/N_{K_{#alpha}}");
    # 4 : "a_{K_{#alpha}}"
    # 5 : "b_{K_{#alpha}}"
    # 6 : "N_{K_{#alpha}}"
    # 7 : "#mu_{K_{#alpha}}"
    # 8 : "#sigma_{K_{#alpha}}"
    # assign lines as parameters
    params[0] = n_kalpha_esc
    params[1] = mu_kalpha_esc
    params[2] = sigma_kalpha_esc
    params[5] = n_kalpha
    params[6] = mu_kalpha
    params[7] = sigma_kalpha

    # get lists defining the bounds of the parameters
    l_bounds, u_bounds = getBoundsList(8)
    # set parameter bounds
    l_bounds[0] = 0
    u_bounds[0] = 10000
    l_bounds[5] = 0
    u_bounds[5] = 10000

    l_bounds[1] = mu_kalpha_esc*0.8
    u_bounds[1] = mu_kalpha_esc*1.2
    l_bounds[6] = mu_kalpha*0.8
    u_bounds[6] = mu_kalpha*1.2

    l_bounds[2] = sigma_kalpha_esc*0.5
    u_bounds[2] = sigma_kalpha_esc*1.5
    l_bounds[7] = sigma_kalpha*0.5
    u_bounds[7] = sigma_kalpha*1.5

    bounds = (l_bounds, u_bounds)
    print(len(bounds))
    print("Bounds for pixel fit: ", bounds)
    print("N params for pixel fit: ", len(params))

    # only fit in range up to 350 hits. Can take index 350 on both, since we
    # created the histogram for a binning with width == 1 pixel per hit
    print(hist)
    data_tofit = hist[80:320]
    print(binning)
    bins_tofit = binning[80:320]
    #print(data_tofit)
    #print(bins_tofit)
    lb, ub = [np.asarray(b, dtype=float) for b in bounds]
    print(ub)

    # NOTE NOTE NOTE the fit function used is WRONG!!! still points to old XrayCalibPixel.C instead
    # of new one from 2014_15 folder!
    result = curve_fit(feSpectrumFuncPixels, bins_tofit, data_tofit, p0=params, bounds = bounds)#, full_output=True)
    popt = result[0]
    pcov = result[1]
    # Calculate the reduced Chi^2:
    # n_dof: # degrees of freedom
    #n_dof  = (np.size(infodict[0]['fvec']) - np.size(popt))
    #chi_sq = np.sum(infodict[0]['fvec']**2) / n_dof
    print('--------------------------------------------------')
    print('Parameters of calibration fit: ')
    for i, p in enumerate(popt):
        print('p_{} ='.format(i), popt[i], '+-', np.sqrt(pcov[i][i]))
    #print('Chi^2 / dof =', chi_sq)

    print("yay :)")

    return popt, pcov

def fitAndPlotFeSpectrumPixels(data, cuts, outfolder, run_number, fitting_only = False):
    # bin the data
    hist, binning = binData(data, cuts)
    popt, pcov = fitFeSpectrumToPixels(hist, binning, cuts)
    print("POPT ", popt)
    return popt
