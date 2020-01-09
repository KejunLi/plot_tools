#!/usr/bin/env python3
import numpy as np

########################################################################
# This packages contains functions that are used for fitting data
########################################################################


def poly_fit(x, c1, c2, c3):
    """
    polynominal fitting
    """
    y = c1 + c2*x + c3*np.power(x,2)
    return(y)

