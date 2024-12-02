import numpy as np
import scipy.ndimage as nd
import copy

import satellite.fitsutils as fs
import satellite.intensity as si

def getFitsSlit(fits_fn: str, slit: dict):
    mat = fs.loadFitsImageData(fits_fn)
    return fs.getVerticalSlit(mat, slit['y'], slit['x'], slit['w'], slit['h'])

def indexOf(atom: str, spectrum: str, atomic_number: int, fitsd: list):
    for idx, lst in enumerate(fitsd):
        if lst['element'] == atom and lst['spectrum'] == spectrum and lst['atomic'] == atomic_number:
            return idx
    return -1

def specific_slit_analysis(fitsd: list, slits: list):
# Reference element for intensity: Hi_4861
    ref_index = indexOf('H', 'i', 4861, fitsd)
    if ref_index < 0:
        print('[ERROR] Fits list does not contain reference element for intensities!', file=sys.stderr)
        raise RuntimeError
# copy of dictionary
    cpd = copy.deepcopy(fitsd)
# for every slit
    for slit in slits:
# for every FITS get the slit specified
        for idx, fits in enumerate(fitsd):
            ar = getFitsSlit(fits['fn'], slit)
            # print(ar)
# sum all elements of slit
            sm = np.sum(np.sum(ar))
            # print("sum={:.2f}".format(sm))
            # print(fits)
            cpd[idx]['slit_sum'] = sm
            # print(cpd[idx])
# compile the intensities data file (for PyNeb)
        si.makeIntensitiesDataFile(cpd, ref_index, 'slit_sum', 'test.dat')
