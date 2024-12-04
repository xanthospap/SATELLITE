import numpy as np
import scipy.ndimage as nd
import copy
import sys

import satellite.fitsutils as fs
import satellite.intensity as si
import satellite.cfgio as sc
import satellite.roman as sr

import pyneb as pn

def getFitsSlit(fits_fn: str, slit: dict):
    mat = fs.loadFitsImageData(fits_fn)
    return fs.getVerticalSlit(mat, slit['y']-1, slit['x']-1, slit['w'], slit['h'])

def computeRatio(ratio: str, intensity_list: list):
    def getIntensity(element):
        for obj in intensity_list:
            if obj['element'] == element:
                return obj['intensity'], obj['intensity_err']
        raise RuntimeError('[ERROR] Cannot find element {:} for ratio {:} in intensities list!'.format(element, ratio))
        return None,None
    ar, par = sc.resolveRatioStr(ratio)
    var = 0e0; vpar = 0e0; par1 = 0e0; par2 = 0e0;
    for idx in range(1, len(ar), 2):
        # val, err = getIntensity(sc.satellite_str2pyneb_str(ar[idx]))
        val, err = getIntensity(ar[idx])
        var  += ar[idx-1] * val
        par1 += ar[idx-1] * (err / val * np.log(10))
    for idx in range(1, len(par), 2):
        # val, err = getIntensity(sc.satellite_str2pyneb_str(par[idx]))
        val, err = getIntensity(par[idx])
        vpar  += par[idx-1] * val
        par2  += par[idx-1] * (err / val * np.log(10))
    return var / vpar, np.sqrt(par1**2 + par2**2), ratio

monte_carlo_fake_obs = 3
reference_element = {'element': 'H', 'spectrum': 'i', 'atomic': 4861}

def specific_slit_analysis(fitsd: list, slits: list, ratios: list, ext_law: str, intensities_out: str, ratios_out):

    global_intensities = {}
    def add_global_intensities(new_list, new_index):
        if new_index == 0:
            for element in new_list:
                global_intensities[element['element']] = [element['intensity'], element['intensity_err']] + [np.nan]*(len(slits)-1)*2
        else:
            for element in new_list:
                if element['element'] not in global_intensities:
                    raise RuntimeError("[ERROR] Element {:} not found in stored list at slit nr {:}".format(element['element'], new_index))
                global_intensities[element['element']][new_index*2  ] = element['intensity']
                global_intensities[element['element']][new_index*2+1] = element['intensity_err']

    global_ratios = {}
    def add_global_ratio(val, err, ratio, new_index):
        if new_index == 0:
            global_ratios[ratio] = [val, err] + [np.nan]*(len(slits)-1)*2
        else:
            if ratio not in global_ratios:
                raise RuntimeError("[ERROR] Element {:} not found in stored list at slit nr {:}".format(ratio, new_index))
            global_ratios[ratio][new_index*2  ] = val
            global_ratios[ratio][new_index*2+1] = err

# for every slit
    for slit_idx, slit in enumerate(slits):
# copy of dictionary
        cpd = copy.deepcopy(fitsd)
# for every FITS get the slit specified
        for idx, fits in enumerate(fitsd):
            ar = getFitsSlit(fits['fns'], slit)
            np.nan_to_num(ar, False)
# sum all elements of slit
            sm = np.sum(np.sum(ar))
            cpd[idx]['sslit_sum'] = sm
            ar = getFitsSlit(fits['fne'], slit)
            np.nan_to_num(ar, False)
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
        sobs.correctData(normWave=4861.)

        eobs = pn.Observation()
        eobs.readData('test.dat', fileFormat='lines_in_rows_err_cols', errIsRelative=False)
        eobs.addMonteCarloObs(N=monte_carlo_fake_obs)
        eobs.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)
        eobs.extinction.law = ext_law
        eobs.correctData(normWave=4861.)

        RC = pn.RedCorr(E_BV=sobs.extinction.E_BV[0], law=ext_law)
        ref_pnstr = sc.objectIntensityPyNebCode(reference_element['element'], reference_element['spectrum'], reference_element['atomic'])
        iref = float(sobs.getIntens()[ref_pnstr])
        eref = float(sobs.getError()[ref_pnstr])

        intensities_list = []
        for idx, fits in enumerate(fitsd):
            scor  = RC.getCorr(fits['atomic'], reference_element['atomic'])
            ecor  = RC.getErrCorr(fits['atomic'], np.std(eobs.extinction.E_BV), reference_element['atomic'])
            pnstr = sc.objectIntensityPyNebCode(fits['element'], fits['spectrum'], fits['atomic'])
            iele  = float(sobs.getIntens()[pnstr])
            eele  = float(sobs.getError()[pnstr])
            err   = np.sqrt(eele**2 + eref**2 + float(ecor/scor)**2)
            intensities_list.append({'element_pn': pnstr, 'element': '{:}{:}_{:}'.format(fits['element'], sr.roman2int(fits['spectrum']), fits['atomic']), 'intensity': iele, 'intensity_err': err})
        add_global_intensities(intensities_list, slit_idx)
        
        for ratio in ratios:
            try:
                val, err, rstr = computeRatio(ratio, intensities_list)
                add_global_ratio(val, err, rstr, slit_idx)
            except:
                print('[WRNNG] Skipping ratio {:}'.format(ratio))
                
    with open(intensities_out, 'w') as fout:
        for key, lst in global_intensities.items():
            print("{:<10} {:}".format(key, ' '.join(['{:10.4f}'.format(x) for x in lst])), file=fout)
    
    with open(ratios_out, 'w') as fout:
        for key, lst in global_ratios.items():
            print("{:<45} {:}".format(key, ' '.join(['{:10.4f}'.format(x) for x in lst])), file=fout)
