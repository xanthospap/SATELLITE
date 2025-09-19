import re

import numpy as np
import pyneb as pn

import satellite.cfgio as sc
import satellite.roman as sr
import satellite.nomenclature as sn


def extract_wavelength(line_id):
    """
    Convert a line identifier (e.g., 'O2_3729A') to the corresponding 'L(wavelength)' format.

    Parameters
    ----------
    line_id: str
        Line identifier (e.g., 'O2_3729A')

    Return
    ------
    str:
        The to_eval string (e.g., 'L(3729)') or None if invalid
    """
    try:
        # Extract the wavelength (everything after the last underscore)
        wavelength = line_id.split("_")[-1].replace("A", "")  # Remove 'A' if present
        # Ensure the extracted value is a valid number
        if wavelength.isdigit():
            return f"L({wavelength})"
    except Exception as e:
        print(f"Error processing {line_id}: {e}")
    return None  # Return None if parsing fails


def get_atom_model_obsolete(element: str, ion: int, logger=None):
    """
    Determines whether to use pn.RecAtom() (for recombination lines)
    or pn.Atom() (for collisional excitation lines).

    Parameters
    ----------
        element: str
            The element symbol (e.g., 'He', 'O')
        ion: int
            The ionization state (e.g., 1 for He I, 2 for He II)

    Returns
    -------
        A PyNeb Atom or RecAtom object

    Example
    -------
        >>> atoms = [
            ("H", 1),  # H I (Recombination)
            ("He", 1), # He I (Recombination)
            ("He", 2), # He II (Recombination)
            ("O", 1),  # O I (Collisional)
            ("O", 2),  # O II (Collisional)
            ("N", 2),  # N II (Collisional)
        ]
    """
    rec_atom_elements = {"H", "He"}
    if element in rec_atom_elements:
        logger.debug(
            f"Creating recombination atom from entry {element} {ion} for ionic abundancies. pn.RecAtom({element}, {ion})"
        )
        return pn.RecAtom(element, ion)
    else:
        logger.debug(
            f"Creating (default) atom from entry {element} {ion} for ionic abundancies."
        )
        return pn.Atom(element, ion)


def get_atom_model(dct_entry, logger=None):
    """
    Determines whether to use pn.RecAtom() (for recombination lines)
    or pn.Atom() (for collisional excitation lines).
    """
    ref_type = "CEL" if "ref_type" not in dct_entry else dct_entry["ref_type"]
    if ref_type.upper() not in ["RL", "CEL"]:
        if logger:
            logger.error(f"Invalid ref_type for entry: {dct_entry}")
        raise RuntimeError(f"ERROR Invalid ref_type for entry: {dct_entry}")
    element = dct_entry["element"]
    ion = sr.roman2int(dct_entry["spectrum"])
    if ref_type.upper() == "CEL":
        if "ref_type" not in dct_entry:
            logger.debug(
                f'Assuming ref_type=CEL for {element}{dct_entry["spectrum"]} beacuse ref_type entry is missing.'
            )
        logger.debug(
            f'Creating (default) atom from entry {element} {ion} (or {element}{dct_entry["spectrum"]}) for ionic abundancies.'
        )
        return pn.Atom(element, ion)
    else:
        logger.debug(
            f'Creating recombination atom from entry {element} {ion} (or {element}{dct_entry["spectrum"]}) for ionic abundancies. pn.RecAtom({element}, {ion})'
        )
        return pn.RecAtom(element, ion)


def refTenNe2PyNebPair(ref_tene):
    """
    Example:
        ref_tene="[SIII] 6312/9069 [ClIII] 5538/5518"
    Return:
        "[SIII] 6312/9069", "[ClIII] 5538/5518"
    """
    match = re.match(
        r"\[([A-Za-z0-9]*)\]\s*([0-9a-zA-Z\+]*[/0-9a-zA-Z\+]*)\s*\[([A-Za-z0-9]*)\]\s*([0-9a-zA-Z\+]*[/0-9a-zA-Z\+]*)",
        ref_tene,
    )
    return "[{:}] {:}".format(match.group(1), match.group(2)), "[{:}] {:}".format(
        match.group(3), match.group(4)
    )


def computeIonicAbundancies(fitsd, tene_dict, pnObs, pnErrObs, logger):
    ionic_abundancies_dict = []

    def ref_tene_pair(tene):
        for entry in tene_dict:
            if entry["tene_pair"] == tene:
                return entry
        return None

    for entry in fitsd:
        reftene = ref_tene_pair(refTenNe2PyNebPair(entry["ref_tene"]))
        if reftene is None:
            logger.error(
                "ERROR Failed finding reference Te/Ne pair for ionic abundancies!"
            )
            return
        # pn_atom = get_atom_model_obsolete(
        #    entry["element"], sr.roman2int(entry["spectrum"]), logger
        # )
        pn_atom = get_atom_model(entry, logger)
        pn_element = sn.objectIntensityPyNebCode(
            entry["element"], entry["spectrum"], entry["atomic"], logger
        )
        sabd = pn_atom.getIonAbundance(
            int_ratio=pnObs.getIntens(0)[pn_element],
            tem=reftene["sT"],
            den=reftene["sN"],
            to_eval=extract_wavelength(pn_element),
            Hbeta=100.0,
        )[0]
        # Uncertainty, loop over MC simulated inttensity ratios ...
        int_mc = pnErrObs.getIntens()[pn_element]
        tar = []
        for idx in range(len(int_mc)):
            tar.append(
                pn_atom.getIonAbundance(
                    int_ratio=int_mc[idx],
                    tem=reftene["eT"][idx],
                    den=reftene["eN"][idx],
                    to_eval=extract_wavelength(pn_element),
                    Hbeta=100.0,
                )
            )
        eabd_array = np.array(tar)
        assert not (np.isnan(eabd_array).any() or np.isnan(np.std(eabd_array)))
        ionic_abundancies_dict.append(
            {
                "element": entry["element"],
                "spectrum": entry["spectrum"],
                "atomic": entry["atomic"],
                "pn_element": pn_element,
                "abundance": sabd,
                # "abundance_error": eabd_array.std(),
                "abundance_error": np.std(eabd_array),
            }
        )
    return ionic_abundancies_dict


def printIonicAbundancies(dict_of_abundancies, fn, logger):
    element_format = "pn_element"

    # given an (inner) dictionary of 'intensities' return the specific
    # dictionary of a line, if it exists, else None
    def inDictOf(abundanciesDict, line):
        for entry in abundanciesDict:
            if entry[element_format] == line:
                return entry
        return None

    # extract unique line and store columns
    unique_lines = []
    columns = []
    for k, v in dict_of_abundancies.items():
        columns += [k]
        unique_lines = list(set(unique_lines + [entry["pn_element"] for entry in v]))

    # sort rows & columns
    unique_lines = sorted(unique_lines)
    columns = sorted(columns)

    # if columns are numeric values and start at 0, then add an 1 offset so they start from 1
    offset = ""
    try:
        [int(c) for c in columns]
        if columns[0] == 0:
            offset = 1
    except:
        pass

    with open(fn, "w") as fout:
        # write first line, i.e. column keys
        print("{:15s}".format("Slit Nr."), file=fout, end="")
        for col in columns:
            print("{:->6d}{:25s} ".format(col + offset, "-" * 25), file=fout, end="")
        print("", file=fout)

        # iterate for every line in unique_lines
        for line in unique_lines:
            # do not print abundance for H1r_4861A and H1r_6563A
            if line not in ["H1r_4861A", "H1r_6563A"]:
                print("{:15s}".format(line), file=fout, end="")
                for col in columns:
                    entry = inDictOf(dict_of_abundancies[col], line)
                    if entry is not None:
                        print(
                            "{:15.9e} {:15.9e} ".format(
                                entry["abundance"], entry["abundance_error"]
                            ),
                            file=fout,
                            end="",
                        )
                    else:
                        print("{:31s} ".format(" "), file=fout, end="")
                print("", file=fout)

    return fn
