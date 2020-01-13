#!/usr/bin/env python3
import numpy as np
from scipy.optimize import curve_fit

########################################################################
# This packages contains functions that are used for fitting data
########################################################################

def poly_fct(x, c0, c1, c2):
    """
    polynominal fitting
    """
    y = c0 + c1*x + c2*np.power(x,2) # definition of function
    return(y)

def best_vals_of_poly_fct(l_x, l_y):
    """
    looking for the best values of polynominal fitting parameters
    """
    # fitting part
    init_vals = [0.5, 0.5, 0.5] # for c0, c1, c2
    best_vals, covar = curve_fit(poly_fct, l_x, l_y, p0=init_vals)     
    print("best_vals: {}".format(best_vals))
    return(best_vals)



def lin_fct(x, c1, c2):
    """
    linear fitting
    """
    y = c1 + c2*x
    return(y)

def best_vals_of_lin_fct(l_x, l_y):
    """
    looking for the best values of linear fitting parameters
    """
    init_vals = [0.5, 0.5]
    best_vals, covar = curve_fit(lin_fct, l_x, l_y, p0=init_vals)
    print("best_vals: {}".format(best_vals))
    return(best_vals)
