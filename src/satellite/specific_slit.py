import numpy as np
import scipy.ndimage as nd
import copy

import satellite.fitsutils as fs
import satellite.intensity as si
import satellite.cfgio as sc

import pyneb as pn

def getFitsSlit(fits_fn: str, slit: dict):
    mat = fs.loadFitsImageData(fits_fn)
    return fs.getVerticalSlit(mat, slit['y'], slit['x'], slit['w'], slit['h'])

ext_law = 'S79 H83 CCM89'
monte_carlo_fake_obs = 3
reference_element = {'element': 'H', 'spectrum': 'i', 'atomic': 4861}

def specific_slit_analysis(fitsd: list, slits: list):
# copy of dictionary
    cpd = copy.deepcopy(fitsd)
# for every slit
    for slit in slits:
# for every FITS get the slit specified
        for idx, fits in enumerate(fitsd):
            ar = getFitsSlit(fits['fns'], slit)
# sum all elements of slit
            sm = np.sum(np.sum(ar))
            cpd[idx]['sslit_sum'] = sm
            ar = getFitsSlit(fits['fne'], slit)
# sum all elements of slit
            sm = np.sum(np.sum(ar))
            cpd[idx]['eslit_sum'] = sm
            # print(cpd[idx])
# compile the intensities data file (for PyNeb)
        si.makeIntensitiesDataFile(cpd, reference_element, ['sslit_sum', 'eslit_sum'], 'test.dat')
# PyNeb stuff
        sobs = pn.Observation()
        sobs.readData('test.dat', fileFormat='lines_in_rows_err_cols', errIsRelative=False)
        sobs.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)
        sobs.extinction.law = ext_law
        sobs.correctData()

        eobs = pn.Observation()
        eobs.readData('test.dat', fileFormat='lines_in_rows_err_cols', errIsRelative=False)
        eobs.addMonteCarloObs(N=monte_carlo_fake_obs)
        eobs.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)
        eobs.extinction.law = ext_law
        eobs.correctData()

        RC = pn.RedCorr(E_BV=sobs.extinction.E_BV[0], law=ext_law)
        ref_pnstr = si.objectIntensityPyNebCode(reference_element['element'], reference_element['spectrum'], reference_element['atomic'])
        iref = float(sobs.getIntens()[ref_pnstr])
        eref = float(sobs.getError()[ref_pnstr])

        for idx, fits in enumerate(fitsd):
            scor  = RC.getCorr(fits['atomic'], reference_element['atomic'])
            ecor  = RC.getErrCorr(fits['atomic'], np.std(eobs.extinction.E_BV), reference_element['atomic'])
            pnstr = si.objectIntensityPyNebCode(fits['element'], fits['spectrum'], fits['atomic'])
            iele  = float(sobs.getIntens()[pnstr])
            eele  = float(sobs.getError()[pnstr])
            err   = np.sqrt(eele**2 + eref**2 + float(ecor/scor)**2)

