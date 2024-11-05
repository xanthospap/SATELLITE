## translated from:
## SATELLITE/satellite/satellite/norm_flux_error_script.py

import numpy as np

def fluxError(a1,sa1,a2,sa2):
    if (a1 != 0e0) and (a2 != 0e0):
        par1 = sa1 / a1
        par2 = sa2 / a2
        error = (a1 * 100e0 / a2) * np.sqrt(par1 * par1 + par2 * par2)
    else:
        error = 0
    return error

def ratioError(a1,sa1,a2,sa2):
    if (a1 != 0e0) and (a2 != 0e0):
        par1 = sa1 / a1
        par2 = sa2 / a2
    else:
        par1 = par2 = 0e0
    return np.sqrt(par1*par1+par2*par2)

