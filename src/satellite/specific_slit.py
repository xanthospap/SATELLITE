import numpy as np
import scipy.ndimage as nd
import copy
import sys
import re
import logging
import warnings

import satellite.fitsutils as fs
import satellite.intensity as si
import satellite.cfgio as sc
import satellite.roman as sr
import satellite.satlogger as sl
import satellite.tene as st
import satellite.ionic_abundancies as sa
import satellite.abundance as sb
import satellite.icf as sf

import pyneb as pn


def getFitsSlit(fits_fn: str, slit: dict, logger=None):
    mat = fs.loadFitsImageData(fits_fn)
    return fs.getVerticalSlit(mat, slit['y']-1, slit['x']-1, slit['w'], slit['h'], logger)


def computeRatio(ratio: str, intensity_list: list, logger=None):
    def getIntensity(element):
        for obj in intensity_list:
            if obj['element'] == element:
                return obj['intensity'], obj['intensity_err']
        if logger:
            logger.error("Cannot find element {:} for ratio {:} in intensities list!".format(
                element, ratio))
        raise RuntimeError(
            '[ERROR] Cannot find element {:} for ratio {:} in intensities list!'.format(element, ratio))
        return None, None
    ar, par = sc.resolveRatioStr(ratio)
    var = 0e0
    vpar = 0e0
    par1 = 0e0
    par2 = 0e0
    for idx in range(1, len(ar), 2):
        # val, err = getIntensity(sc.satellite_str2pyneb_str(ar[idx]))
        val, err = getIntensity(ar[idx])
        var += ar[idx-1] * val
        par1 += ar[idx-1] * (err / val * np.log(10))
    for idx in range(1, len(par), 2):
        # val, err = getIntensity(sc.satellite_str2pyneb_str(par[idx]))
        val, err = getIntensity(par[idx])
        vpar += par[idx-1] * val
        par2 += par[idx-1] * (err / val * np.log(10))
    return var / vpar, np.sqrt(par1**2 + par2**2), ratio


def extract_ion(label):
    """ Example:
        t = "[OI] 5577/6300+"
        ion = extract_ion(t)
        print(ion)  # Output: "[OI]"
    """
    match = re.match(r"\[(.*?)\]", label)  # Find text inside square brackets
    return match.group(0) if match else None  # Return full [Ion] if found


def err2scalar(err_array):
    """ Define the fucntion to produce a single (i.e. scalar) error value, 
        when we have an array of error values (i.e. from Monte-Carlo 
        simulations.
    """
    if np.isnan(err_array).any():
        print('WARNING  array contains nan! {:}'.format(err_array))
    return np.std(err_array)


monte_carlo_fake_obs = 3
reference_element = {'element': 'H', 'spectrum': 'i', 'atomic': 4861}


def specific_slit_analysis(fitsd: list, slits: list, ratios: list, density_diagnostics: list, tempterature_diagnostics: list, ext_law: str, intensities_out: str, ratios_out: str, logger):

    # attach PyNeb's logger to our, if we have one!
    if logger:
        pyneb_logger = logging.getLogger("pyneb")
        # Remove existing PyNeb handlers (Fixes duplicate log outputs)
        for handler in pyneb_logger.handlers[:]:
            pyneb_logger.removeHandler(handler)
        pyneb_logger.setLevel(logger.level)
        for handler in logger.handlers:
            pyneb_logger.addHandler(handler)
        pyneb_logger.propagate = False
        # Redirect PyNeb warnings (stderr) to our logger
        # Redirects warnings.warn() to logging.WARNING
        logging.captureWarnings(True)
        warnings.simplefilter("always")  # Ensure all warnings are captured
        sys.stderr = sl.PyNebLogRedirector(logger, logging.WARNING)

    pn.log_.open_file('pyneblog.log')

    global_intensities = {}
    global_ratios = {}
    global_tene = {}

    def add_global_intensities(new_list, new_index):
        if new_index == 0:
            for element in new_list:
                global_intensities[element['element']] = [
                    element['intensity'], element['intensity_err']] + [np.nan]*(len(slits)-1)*2
        else:
            for element in new_list:
                if element['element'] not in global_intensities:
                    logger.error("Element {:} not found in stored list at slit nr {:}".format(
                        element['element'], new_index))
                    raise RuntimeError("[ERROR] Element {:} not found in stored list at slit nr {:}".format(
                        element['element'], new_index))
                global_intensities[element['element']
                                   ][new_index*2] = element['intensity']
                global_intensities[element['element']
                                   ][new_index*2+1] = element['intensity_err']

    def add_global_ratio(val, err, ratio, new_index):
        if new_index == 0:
            global_ratios[ratio] = [val, err] + [np.nan]*(len(slits)-1)*2
        else:
            if ratio not in global_ratios:
                logger.error(
                    "Element {:} not found in stored list at slit nr {:}".format(ratio, new_index))
                raise RuntimeError(
                    "[ERROR] Element {:} not found in stored list at slit nr {:}".format(ratio, new_index))
            global_ratios[ratio][new_index*2] = val
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

        ## <-- End Looping FITS --> ##

# compile the intensities data file (for PyNeb)
        si.makeIntensitiesDataFile(cpd, reference_element, [
                                   'sslit_sum', 'eslit_sum'], 'test.dat')

# PyNeb stuff
        sobs = pn.Observation()
        sobs.readData(
            'test.dat', fileFormat='lines_in_rows_err_cols', errIsRelative=False)
        sobs.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)
        sobs.extinction.law = ext_law
        sobs.correctData(normWave=4861.)

        eobs = pn.Observation()
        eobs.readData(
            'test.dat', fileFormat='lines_in_rows_err_cols', errIsRelative=False)
        eobs.addMonteCarloObs(N=monte_carlo_fake_obs)
        eobs.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)
        eobs.extinction.law = ext_law
        eobs.correctData(normWave=4861.)

        RC = pn.RedCorr(E_BV=sobs.extinction.E_BV[0], law=ext_law)

# Compute intensity for each FITS/atom
        add_global_intensities(si.computeIntensities(
            fitsd, sobs, eobs, RC, reference_element, logger), slit_idx)

# Comptue intensity rations
        for ratio in ratios:
            try:
                val, err, rstr = computeRatio(ratio, intensities_list)
                add_global_ratio(val, err, rstr, slit_idx)
            except:
                logger.info("Skipping ratio {:}".format(ratio))

        tene_dict = st.computeTeNePairs(
            density_diagnostics, tempterature_diagnostics, sobs, eobs, logger)
        global_tene[slit_idx] = tene_dict

        ionic_abundancies = sa.computeIonicAbundancies(
            cpd, tene_dict, sobs, eobs, logger)

# ICFs
        elemspec_abundancies = sb.computeAbundancies(
            cpd, ionic_abundancies, logger)
        elem_abundancies = sf.computeIcfs(elemspec_abundancies, logger)
        print(elem_abundancies)

    ## <-- End Looping Slits --> ##
    with open(intensities_out, 'w') as fout:
        for key, lst in global_intensities.items():
            print("{:<10} {:}".format(key, ' '.join(
                ['{:10.4f}'.format(x) for x in lst])), file=fout)

    with open(ratios_out, 'w') as fout:
        for key, lst in global_ratios.items():
            print("{:<45} {:}".format(key, ' '.join(
                ['{:10.4f}'.format(x) for x in lst])), file=fout)

    with open('tene', 'w') as fout:
        all_rows = []
        for idx, elist in global_tene.items():
            for entry in elist:
                # each elist is in the following form:
                # {'tene_pair': (t,d), 'sT': st, 'sN': sn, 'eT': et, 'eN': en}
                if entry['tene_pair'] not in all_rows:
                    all_rows.append(entry['tene_pair'])
        print("{:45s}{:}".format('#', ''.join(["Slit {:5d}{:16s}{:26s}".format(
            d, ' ', ' ') for d in range(len(slits))])), file=fout)
        print("{:45s}{:11s}   {:11s} {:11s}   {:11s}".format(
            "Te/Ne Pair", "Tempature", "", "Density", ""), file=fout)
        print('-'*(45+26*2*len(slits)), file=fout)
        for pair in all_rows:
            print("{:45s}".format('/'.join(pair)), end='', file=fout)
# value of TeNe pair (pair) for all slits ...
            for sj in range(len(slits)):
                elist = global_tene[sj]
                edict = next(
                    (item for item in elist if item["tene_pair"] == pair), None)
# TODO handle case where edict is None!!
                print("{:.5e} \u00B1 {:.5e} {:.5e} \u00B1 {:.5e} ".format(edict['sT'], err2scalar(
                    edict['eT']), edict['sN'], err2scalar(edict['eN'])), end='', file=fout)
# pair done !
            print('', file=fout)

    pn.log_.close_file()
