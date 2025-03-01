import sys
import pyneb as pn
import numpy as np

# from satellite import cfgio
from satellite import astroflux
import satellite.roman as sr
import satellite.cfgio as sc

def makeIntensitiesDataFile(fitsd: list, reference_element: dict, value_keys: list, fn: str, factor=100e0):
    ref_index = sc.indexOf(
        reference_element['element'], reference_element['spectrum'], reference_element['atomic'], fitsd)
    ref_sval = fitsd[ref_index][value_keys[0]]
    ref_eval = fitsd[ref_index][value_keys[1]]
    with open(fn, 'w') as fout:
        print('LINE test err', file=fout)
        for obj in fitsd:
            pnlabel = sc.objectIntensityPyNebCode(
                obj['element'], obj['spectrum'], obj['atomic'])
            print('{:} {:+9e} {:+9e}'.format(pnlabel, factor*obj[value_keys[0]]/ref_sval, astroflux.fluxError(
                obj[value_keys[0]], obj[value_keys[1]], ref_sval, ref_eval)), file=fout)
    return fn


def computeIntensities(fitsd: dict, pnObs, pnErrObs, pnRC, reference_element: dict, logger):
    ref_pnstr = sc.objectIntensityPyNebCode(
        reference_element['element'], reference_element['spectrum'], reference_element['atomic'], logger)
    iref = float(pnObs.getIntens()[ref_pnstr])
    eref = float(pnObs.getError()[ref_pnstr])
    intensities_list = []
    for idx, fits in enumerate(fitsd):
        scor = pnRC.getCorr(fits['atomic'], reference_element['atomic'])
        ecor = pnRC.getErrCorr(fits['atomic'], np.std(
            pnErrObs.extinction.E_BV), reference_element['atomic'])
        pnstr = sc.objectIntensityPyNebCode(
            fits['element'], fits['spectrum'], fits['atomic'], logger)
        iele = float(pnObs.getIntens()[pnstr])
        eele = float(pnObs.getError()[pnstr])
        err = np.sqrt(eele**2 + eref**2 + float(ecor/scor)**2)
        intensities_list.append({'element_pn': pnstr, 'element': '{:}{:}_{:}'.format(
            fits['element'], sr.roman2int(fits['spectrum']), fits['atomic']), 'intensity': iele, 'intensity_err': err})
    return intensities_list
