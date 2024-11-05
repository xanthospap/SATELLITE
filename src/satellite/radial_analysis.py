## translated from:
## SATELLITE/satellite/satellite/radial_analysis_script.py

import numpy as np
import scipy.ndimage as nd
import satellite.fitsutils as fs

angle_step = 
slit_width = 
slit_height = 

def radialAnalysis(inlist: list):
    for angle in range(0, 360, angle_step):
        for entry in inlist:
# load FITS instance
            fits = fs.loadFitsImageData(entry['fn'])
# reshape so that (row, col) pixel is at the centre
            fits = fs.setCenterPixel(fits, row, col, False, np.inf)
# rotate the matrix by angle
            fits = nd.rotate(fits, angle, axis=(1,0), reshape=False, order=5, cval=np.inf)
# get vertical slit (1D-array)
            n, m = fits.shape
            slit = fs.getVerticalSlit(n//2, m//2, slit_width, m//2)
# remove np.inf values from array and get mean value
            mean = np.mean(slit[~np.isinf(slit)])
