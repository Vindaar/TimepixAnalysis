from ingrid.ingrid_helper_functions import *
from scipy.optimize import curve_fit
import numpy as np
import matplotlib as mpl
import sys
sys.path.append("./")
from procsForPython import polyaPython
from scipy.special import gamma
import math

def polyaInPython(x, p0, p1, p2):
  ## Polya function to fit to TOT histogram / charge in electrons of a
  ## run. This is the actual implementation of the polya distribution.
  ## Parameters:
  ## N     = p[0]    scaling factor
  ## G     = p[1]    gas gain
  ## theta = p[2]    parameter, which describes distribution (?! I guess it makes sens
  ##                 since we take its power and it enters gamma)
  thetaDash = p2 + 1
  coeff1 = (p0 / p1) * np.power((thetaDash), thetaDash) / gamma(thetaDash)
  coeff2 = np.power((x / p1), p2) * np.exp(-thetaDash * x / p1)
  result = coeff1 * coeff2
  return result

def fitPolyaFunc(chToFit, countsToFit, p, bounds):
    result = curve_fit(polyaPython, chToFit, countsToFit,
                       p0=p, bounds = bounds)#, full_output=True)
    return result
