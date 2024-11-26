import numpy as np
import scipy.ndimage as nd
import satellite.fitsutils as fs

def getFitsSlit(fits_fn: str, slit: dict):
    mat = fs.loadFitsImageData(fits_fn)
    return fs.getVerticalSlit(mat, slit['y'], slit['x'], slit['w'], slit['h'])

def specific_slit_analysis(fitsd: dict, slits: list):
    # for every slit
    for slit in slits:
        # for every FITS get the slit specified
        for idx, fits in enumerate(fitsd):
            ar = getFitsSlit(fits['fn'], slit)
            print(ar)
