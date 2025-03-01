import numpy as np
import pyneb as pn

import satellite.cfgio as sc
import satellite.roman as sr


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
        wavelength = line_id.split(
            '_')[-1].replace('A', '')  # Remove 'A' if present
        # Ensure the extracted value is a valid number
        if wavelength.isdigit():
            return f"L({wavelength})"
    except Exception as e:
        print(f"Error processing {line_id}: {e}")
    return None  # Return None if parsing fails


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


def computeIonicAbundancies(fitsd, tene_dict, pnObs, pnErrObs, logger):
    ionic_abundancies_dict = []

    def ref_tene_pair(tene):
        for entry in tene_dict:
            if entry['tene_pair'] == tene:
                return entry
        return None

    for entry in fitsd:
        reftene = ref_tene_pair(sc.configRefTenNe2PyNebPair(entry['ref_tene']))
        if reftene is None:
            print("ERROR Failed finding reference Te/Ne pair for ionic abundancies!")
            return
        pn_atom = get_atom_model(
            entry['element'], sr.roman2int(entry['spectrum']))
        pn_element = sc.objectIntensityPyNebCode(
            entry['element'], entry['spectrum'], entry['atomic'], logger)
        # logger.debug("int_ratio={:}, tem={:}, den={:}, to_eval={:}, Hbeta={:}".format(sobs.getIntens(0)[pn_element][0], reftene['sT'], reftene['sN'], extract_wavelength(pn_element), 100.))
        sabd = pn_atom.getIonAbundance(int_ratio=pnObs.getIntens(
            0)[pn_element], tem=reftene['sT'], den=reftene['sN'], to_eval=extract_wavelength(pn_element), Hbeta=100.)[0]
        eabd = pn_atom.getIonAbundance(int_ratio=pnErrObs.getIntens(
            0)[pn_element], tem=reftene['eT'], den=reftene['eN'], to_eval=extract_wavelength(pn_element), Hbeta=100.)[0]
        # print('{:}+{:}({:})/H+ = {:.2e} \u00B1 {:.4e}'.format(entry['element'], sr.roman2int(entry['spectrum']), extract_wavelength(pn_element), sabd, eabd))
        ionic_abundancies_dict.append({'element': entry['element'], 'spectrum': entry['spectrum'],
                                       'atomic': entry['atomic'], 'pn_element': pn_element, 'abundance': sabd, 'abundance_error': eabd})
    return ionic_abundancies_dict
