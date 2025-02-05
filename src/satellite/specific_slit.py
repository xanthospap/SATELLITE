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

import pyneb as pn


def getFitsSlit(fits_fn: str, slit: dict, logger=None):
    mat = fs.loadFitsImageData(fits_fn)
    return fs.getVerticalSlit(mat, slit['y']-1, slit['x']-1, slit['w'], slit['h'], logger)

def computeRatio(ratio: str, intensity_list: list, logger=None):
    def getIntensity(element):
        for obj in intensity_list:
            if obj['element'] == element:
                return obj['intensity'], obj['intensity_err']
        if logger: logger.error("Cannot find element {:} for ratio {:} in intensities list!".format(element, ratio))
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


def get_atom_model(element, ion, logger=None):
    """
    Determines whether to use pn.RecAtom() (for recombination lines)
    or pn.Atom() (for collisional excitation lines).
    
    Parameters:
        element (str): The element symbol (e.g., 'He', 'O')
        ion (int): The ionization state (e.g., 1 for He I, 2 for He II)

    Returns:
        A PyNeb Atom or RecAtom object

    Example:
        atoms = [
            ("H", 1),  # H I (Recombination)
            ("He", 1), # He I (Recombination)
            ("He", 2), # He II (Recombination)
            ("O", 1),  # O I (Collisional)
            ("O", 2),  # O II (Collisional)
            ("N", 2),  # N II (Collisional)
        ]
    """
    rec_atom_elements = {'H', 'He'}  # Hydrogen & Helium use RecAtom
    if element in rec_atom_elements:
        return pn.RecAtom(element, ion)  # Use RecAtom for H and He
    else:
        return pn.Atom(element, ion)  # Use Atom for all other elements

def extract_wavelength(line_id):
    """
    Convert a line identifier (e.g., 'O2_3729A') to the corresponding 'L(wavelength)' format.

    Parameters:
        line_id (str): Line identifier (e.g., 'O2_3729A')

    Returns:
        str: The to_eval string (e.g., 'L(3729)') or None if invalid
    """
    try:
        # Extract the wavelength (everything after the last underscore)
        wavelength = line_id.split('_')[-1].replace('A', '')  # Remove 'A' if present
        # Ensure the extracted value is a valid number
        if wavelength.isdigit():
            return f"L({wavelength})"
    except Exception as e:
        print(f"Error processing {line_id}: {e}")
    return None  # Return None if parsing fails

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
        for handler in logger.handlers: pyneb_logger.addHandler(handler)
        pyneb_logger.propagate = False
        # Redirect PyNeb warnings (stderr) to our logger
        logging.captureWarnings(True)  # Redirects warnings.warn() to logging.WARNING
        warnings.simplefilter("always")  # Ensure all warnings are captured
        sys.stderr = sl.PyNebLogRedirector(logger, logging.WARNING)

    pn.log_.open_file('pyneblog.log')

    global_intensities = {}

    def add_global_intensities(new_list, new_index):
        if new_index == 0:
            for element in new_list:
                global_intensities[element['element']] = [
                    element['intensity'], element['intensity_err']] + [np.nan]*(len(slits)-1)*2
        else:
            for element in new_list:
                if element['element'] not in global_intensities:
                    logger.error("Element {:} not found in stored list at slit nr {:}".format(element['element'], new_index))
                    raise RuntimeError("[ERROR] Element {:} not found in stored list at slit nr {:}".format(
                        element['element'], new_index))
                global_intensities[element['element']
                                   ][new_index*2] = element['intensity']
                global_intensities[element['element']
                                   ][new_index*2+1] = element['intensity_err']

    global_ratios = {}

    def add_global_ratio(val, err, ratio, new_index):
        if new_index == 0:
            global_ratios[ratio] = [val, err] + [np.nan]*(len(slits)-1)*2
        else:
            if ratio not in global_ratios:
                logger.error("Element {:} not found in stored list at slit nr {:}".format(ratio, new_index))
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
        ref_pnstr = sc.objectIntensityPyNebCode(
            reference_element['element'], reference_element['spectrum'], reference_element['atomic'], logger)
        iref = float(sobs.getIntens()[ref_pnstr])
        eref = float(sobs.getError()[ref_pnstr])

        intensities_list = []
        for idx, fits in enumerate(fitsd):
            scor = RC.getCorr(fits['atomic'], reference_element['atomic'])
            ecor = RC.getErrCorr(fits['atomic'], np.std(
                eobs.extinction.E_BV), reference_element['atomic'])
            pnstr = sc.objectIntensityPyNebCode(
                fits['element'], fits['spectrum'], fits['atomic'], logger)
            iele = float(sobs.getIntens()[pnstr])
            eele = float(sobs.getError()[pnstr])
            err = np.sqrt(eele**2 + eref**2 + float(ecor/scor)**2)
            intensities_list.append({'element_pn': pnstr, 'element': '{:}{:}_{:}'.format(
                fits['element'], sr.roman2int(fits['spectrum']), fits['atomic']), 'intensity': iele, 'intensity_err': err})
        add_global_intensities(intensities_list, slit_idx)

        for ratio in ratios:
            try:
                val, err, rstr = computeRatio(ratio, intensities_list)
                add_global_ratio(val, err, rstr, slit_idx)
            except:
                logger.info("Skipping ratio {:}".format(ratio))

# TeNe_specific_slits_script
        diags = pn.Diagnostics()
        # Register all diagnostics with PyNeb
        diags.addDiag(density_diagnostics+tempterature_diagnostics)
        # Retrieve the list of valid diagnostics recognized by PyNeb
        all_diags = diags.getDiagLabels()
        # Filter out only valid diagnostics
        valid_temp = [t for t in tempterature_diagnostics if t in all_diags]
        valid_dens = [d for d in density_diagnostics if d in all_diags]
        # Automatically generate valid (temperature, density) pairs
        valid_pairs = [(t, d) for t in valid_temp for d in valid_dens]
        # Iterate valid pairs and add results to dictionary
        tene_slit_dict = []
        for t, d in valid_pairs:
            try:
                st, sn = diags.getCrossTemDen(t, d, obs=sobs)
                et, en = diags.getCrossTemDen(t, d, obs=eobs)
                tene_slit_dict.append({'tene_pair': (t,d), 'sT': st, 'sN': sn, 'eT': et, 'eN': en})
                print("PyNeb: Te {:} = {:.1f} Ne {:} = {:.2f}".format(
                    extract_ion(t), st, extract_ion(d), sn))
            except:
                logger.info(f"Skipping Tem/Den pair {t} (Temp) ↔ {d} (Density)")

# Ionic Abundancies
        # print(cpd)
        def ref_tene_pair():
            for entry in tene_slit_dict:
                # print(entry['tene_pair'])
                if entry['tene_pair'] == ('[SIII] 6312/9069', '[ClIII] 5538/5518'):
                    return entry
            return None
        reftene = ref_tene_pair()
        for entry in cpd:
            pn_atom = get_atom_model(entry['element'], sr.roman2int(entry['spectrum']))
            pn_element = sc.objectIntensityPyNebCode(entry['element'], entry['spectrum'], entry['atomic'], logger)
            logger.debug("int_ratio={:}, tem={:}, den={:}, to_eval={:}, Hbeta={:}".format(sobs.getIntens(0)[pn_element][0], reftene['sT'], reftene['sN'], extract_wavelength(pn_element), 100.))
            sabd = pn_atom.getIonAbundance(int_ratio=sobs.getIntens(0)[pn_element], tem=reftene['sT'], den=reftene['sN'], to_eval=extract_wavelength(pn_element), Hbeta=100.)[0]
            eabd = pn_atom.getIonAbundance(int_ratio=eobs.getIntens(0)[pn_element], tem=reftene['eT'], den=reftene['eN'], to_eval=extract_wavelength(pn_element), Hbeta=100.)[0]
            print('{:}+{:}({:})/H+ = {:.2e} \u00B1 {:.4e}'.format(entry['element'], sr.roman2int(entry['spectrum']), extract_wavelength(pn_element), sabd, eabd))


    ## <-- End Looping Slits --> ##
    with open(intensities_out, 'w') as fout:
        for key, lst in global_intensities.items():
            print("{:<10} {:}".format(key, ' '.join(
                ['{:10.4f}'.format(x) for x in lst])), file=fout)

    with open(ratios_out, 'w') as fout:
        for key, lst in global_ratios.items():
            print("{:<45} {:}".format(key, ' '.join(
                ['{:10.4f}'.format(x) for x in lst])), file=fout)

    pn.log_.close_file()
