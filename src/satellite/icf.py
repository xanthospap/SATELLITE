import re
import pyneb as pn
import numpy as np
import satellite.roman as sr
from typing import Union
from math import isnan
from os import remove


def ctrue(*args):
    return True


def omega(elemspec_abundancies: dict) -> float:
    """Compute Omega value, relevant for ICF checking:
    Omega = O2+ / (O+ + O2+)

    Parameters
    ----------
    elemspec_abundancies: dict
        A dictionary that should include the required values, as:
        O+  : 'O2'
        O2+ : 'O3'

    Return
    ------
    float:
        The value of Omega
    """
    return elemspec_abundancies["O3"] / (
        elemspec_abundancies["O2"] + elemspec_abundancies["O3"]
    )


def U(elemspec_abundancies: dict) -> float:
    """Compute U value, relevant for ICF checking:
    U = He2+ / (He+ + He2+)

    Parameters
    ----------
    elemspec_abundancies: dict
        A dictionary that should include the required values, as:
        He+  : 'He2'
        He2+ : 'He3'

    Return
    ------
    float:
        The value of U
    """
    return elemspec_abundancies["He3"] / (
        elemspec_abundancies["He2"] + elemspec_abundancies["He3"]
    )


def omega1435(ea: dict) -> float:
    return omega(ea) <= 0.5


def omega1436(ea: dict) -> float:
    return omega(ea) > 0.5


def omega1414(ea: dict) -> float:
    return omega(ea) <= 0.95


def omega1414b(ea: dict) -> float:
    return omega(ea) > 0.95


def omega1429b(ea: dict) -> float:
    return omega(ea) > 0.02 and omega(ea) < 0.95


def omega1432(ea: dict) -> float:
    return omega(ea) > 0.95


def u1417c(ea: dict) -> float:
    return U(ea) < 0.015


""" Spectrum notation to abundance/Ionic notation """
elemspec2ionstr_dict = {
    "HeI": "He+",
    "HeII": "He2+",
    "HI": "H+",
    "NI": "N0",
    "NII": "N+",
    "OI": "O0",
    "OII": "O+",
    "OIII": "O2+",
    "SII": "S+",
    "SIII": "S2+",
    "ArIII": "Ar2+",
    "ArIV": "Ar3+",
    "ClII": "Cl+",
    "ClIII": "Cl2+",
    "CI": "C0",
    "CII": "C+",
    "NeIII": "Ne2+",
    "NiII": "Ni+",
    "CaII": "Ca+",
    "FeII": "Fe+",
    "FeIII": "Fe2+",
}


def elemspec2ionstr(element: str, spectrum: Union[int, str]) -> str:
    """This function actually applies the dictionary elemspec2ionstr_dict.
    Given an element and a spectrum, we construct a (possible) key entry
    of the dictionary, and return the corresponding value.

    Parameters
    -----------
    element: str
        The element, e.g. 'He', 'H', 'O', etc
    spectrum:
        The spectrum, which can be either an integer value (i.e. 2) or
        an (upper- or lower-case) roman letter; e.g. ('i', 'II', etc)

    Returns
    -------
    str:
        The corresponding entry of the elemspec2ionstr_dict dictionary.

    Examples
    --------
    >>> elemspec2ionstr('He', 1)     -> 'He+'
    >>> elemspec2ionstr('Fe', 'iii') -> 'Fe2+'
    """
    try:
        int(spectrum)
        spectrum = sr.int2roman(int(spectrum))
    except:
        pass
    return elemspec2ionstr_dict["{:}{:}".format(element, spectrum.upper())]


def ionstr2pyneb(ion):
    """Given an ion string (e.g. 'Fe2+') the function will return the
    corresponding pyneb literal.

    Parameters
    ----------
    ion: str
        The ion string (e.g 'Fe2+', 'He+') etc.

    Returns
    -------
    str:
        The corresponding pyneb entry/string

    Examples:
    """
    g = re.match(r"([A-Za-z]+)(0)", ion)
    try:
        return g[1]
    except:
        pass
    g = re.match(r"([A-Za-z]+)([0-9]?)\+", ion)
    try:
        return g[1] + (str(int(g[2]) + 1) if g[2] != "" else "2")
    except:
        raise RuntimeError(
            "ERROR Failed transforming Ion {:} to PyNeb string".format(ion)
        )


""" A dictionary holding all ICFs/DIMS that may be computed, organized by element.
"""
icf_dict = {
    "Ar": {
        "kb": [
            {"name": "KB94_A32", "cond": ctrue},
            {"name": "KB94_A30.10", "cond": ctrue},
        ],
        "dims": [
            {"name": "DIMS14_35", "cond": omega1435},
            {"name": "DIMS14_36", "cond": omega1436},
        ],
    },
    "N": {
        "kb": [{"name": "KB94_A1.10", "cond": ctrue}],
        "dims": [
            {"name": "DIMS14_14", "cond": omega1414},
            {"name": "DIMS14_14b", "cond": omega1414b},
        ],
    },
    "O": {
        "kb": [
            {"name": "KB94_A10", "cond": ctrue},
            {"name": "KB94_A8", "cond": ctrue},
            {"name": "KB94_A6", "cond": ctrue},
        ],
        "dims": [{"name": "DIMS14_12", "cond": ctrue}],
    },
    "S": {
        "kb": [
            {"name": "KB94_A36.10", "cond": ctrue},
            {"name": "KB94_A38.10", "cond": ctrue},
        ],
        "dims": [
            {"name": "DIMS14_23", "cond": ctrue},
            {"name": "DIMS14_26", "cond": ctrue},
        ],
    },
    "Cl": {
        "kb": [],
        "dims": [
            {"name": "DIMS14_29b", "cond": omega1429b},
            {"name": "DIMS14_32", "cond": omega1432},
        ],
    },
    "Ne": {
        "kb": [
            {"name": "KB94_A28.10", "cond": ctrue},
            {"name": "KB94_A27", "cond": ctrue},
        ],
        "dims": [
            {"name": "DIMS14_17a", "cond": ctrue},
            {"name": "DIMS14_17b", "cond": ctrue},
            {"name": "DIMS14_17c", "cond": u1417c},
        ],
    },
    "C": {
        "kb": [
            {"name": "KB94_A12", "cond": ctrue},
            {"name": "KB94_A13.10", "cond": ctrue},
            {"name": "KB94_A16", "cond": ctrue},
            {"name": "KB94_A19", "cond": ctrue},
            {"name": "KB94_A21", "cond": ctrue},
            {"name": "KB94_A26", "cond": ctrue},
        ],
        "dims": [{"name": "DIMS14_39", "cond": ctrue}],
    },
    "Fe": {
        "kb": [
            {"name": "RR05_2", "cond": ctrue},
            {"name": "RR05_3", "cond": ctrue},
            {"name": "RR05_4", "cond": ctrue},
        ],
        "dims": [],
    },
}


def icfIsOfElement(icf, element):
    if element not in icf_dict:
        return False
    for law_type in icf_dict[element]:
        for entry in icf_dict[element][law_type]:
            if entry["name"] == icf:
                return True
    return False


def getElementFromIcfName(icf_name):
    for elem, law_sets in icf_dict.items():
        for laws in law_sets.values():
            for entry in laws:
                if entry["name"] == icf_name:
                    return elem
    return None  # not found


def IcfNameList() -> list:
    """Return a list with all values of 'name' keys, from the
    icf_dict dictionary.
    """
    return [
        entry["name"]
        for element in icf_dict.values()
        for section in element.values()
        for entry in section
    ]


def renameIons(abundancies: dict) -> dict:
    """Returns a copy of the input dictionary where the the keys are changed
    accroding to the transformation:
    key -> elemspec2ionstr() -> ionstr2pyneb() -> new_key

    Note
    ----
    This transformation makes sure that the following key renamings are
    applied:
        HeI  -> becomes -> He2
        HeII -> becomes -> He3

    Parameters
    ----------
    abundancies: dict
        A dictionary where keys are elements, (e.g. 'HeI', 'OII',
    """
    acopy = {}

    def es_split(elemsp):
        g = re.match(r"([A-Za-z]+)([0-9]?)", elemsp)
        return g[1], g[2]

    for k, v in abundancies.items():
        acopy[ionstr2pyneb(elemspec2ionstr(*es_split(k)))] = v
    return acopy


def computeIcfs(abundancies, logger=None):
    icf = pn.ICF()
    allIcf = icf.getElemAbundance(
        renameIons({a[0]: a[1][0] for a in abundancies.items()}), IcfNameList()
    )
    return allIcf


def computeIcfsWithErrors(abundancies, logger=None):
    """Assume that we are dealing with a formula like:
    X = ICF * (A+B+C), then:
    (sigma_X/X)^2 = (sigma_ICF/ICF)^2 + [sigma_A/A * A/(A+B+C)]^2 +
        [sigma_B/B * B/(A+B+C)]^2 +
        [sigma_C/C * C/(A+B+C)]^2
    """
    ne = computeIcfs(abundancies, logger)

    icf = pn.ICF()
    renamed_abunds = renameIons({k: v[0] for k, v in abundancies.items()})
    renamed_errors = renameIons({k: v[1] for k, v in abundancies.items()})

    results = {}

    for icf_name, total_abund in ne.items():
        if not np.isfinite(total_abund):
            continue

        element = getElementFromIcfName(icf_name)

        # Find all contributing ions for this element
        contributing_ions = [ion for ion in renamed_abunds if ion.startswith(element)]
        if not contributing_ions:
            continue

        ionic_sum = sum(renamed_abunds[ion] for ion in contributing_ions)
        if ionic_sum == 0:
            continue

        # Step 1: ionic part of uncertainty (quadrature of fractional contributions)
        ionic_var = 0.0
        for ion in contributing_ions:
            x = renamed_abunds[ion]
            dx = renamed_errors[ion]
            if x == 0:
                continue
            weight = x / ionic_sum
            ionic_var += (dx / x) ** 2 * weight**2

        # Step 2: get ICF and its uncertainty
        try:
            icf_val = icf.getICF(element, renamed_abunds, icf_name)
            icf_err = icf.getICFError(element, renamed_abunds, icf_name)
        except Exception:
            icf_val = 1.0
            icf_err = 0.0

        # Step 3: combine relative variances
        rel_var = ionic_var + (icf_err / icf_val) ** 2
        abs_uncertainty = total_abund * np.sqrt(rel_var)

        results[icf_name] = (total_abund, abs_uncertainty)

    return results


def ionicAbundance2elementAbundance(abundancies, logger=None):
    """
    Example
    -------
    >>> abundancies = {'O1': (5.2507134485026286e-06, 5.2507134485026286e-06), 'O2': (9.374517342869009e-05, 6.784694600365557e-05), 'O3': (0.0005404778598225776, 0.00038217561305225653), 'He1': (0.10051405354336723, 0.07111350159031496), 'He2': (0.0066002184982703, 0.0066002184982703), 'Ar3': (2.1068655865313865e-06, 2.1068655865313865e-06), 'H1': (1.001583596941619, 0.7082274385540935), 'N1': (1.3020037638389353e-06, 1.3020037638389353e-06), 'N2': (1.5042597487660062e-05, 9.072192621667611e-06), 'Cl3': (1.1600426563474183e-07, 8.202740287769977e-08), 'S2': (5.923084867006183e-07, 4.1906409773640354e-07), 'S3': (6.652299608059526e-06, 4.703886163496009e-06)}
    >>>
    """

    def stripElement(elemspec):
        g = re.match(r"([A-Za-z]+)[0-9]+", elemspec)
        return g[1]

    elemdct = {}
    # for entry, vals in abundancies.items():
    #    element = stripElement(entry)
    #    if element in elemdct:
    #        elemdct[element] += vals[0]
    #    else:
    #        elemdct[element] = vals[0]
    for ion_label, (val, err) in abundancies.items():
        element = stripElement(ion_label)
        if element not in elemdct:
            elemdct[element] = {"abundance": 0e0, "unc_sq_sum": 0e0}

        elemdct[element]["abundance"] += val
        elemdct[element]["unc_sq_sum"] += err**2  # Add variance

    # Finalize: compute sqrt of summed variances
    for element in elemdct:
        elemdct[element]["uncertainty"] = np.sqrt(elemdct[element]["unc_sq_sum"])
        del elemdct[element]["unc_sq_sum"]  # clean up

    return elemdct


def printIcfs(icfs, elem_abundancies, fn, logger):
    # extract unique elements and store columns
    unique_elements = []
    columns = []
    for k, v in elem_abundancies.items():
        columns += [k]
        unique_elements = list(set(unique_elements + [k for k in v.keys()]))

    # sort rows & columns
    unique_elements = sorted(unique_elements)
    columns = sorted(columns)

    def inDictOfElements(slit_abundancies, element):
        return slit_abundancies[element] if element in slit_abundancies else None

    def elementIcf(slit_icfs, element):
        element_icfs = {}
        for k, v in slit_icfs.items():
            if not (np.isnan(v[0]) or isnan(v[0])):
                if icfIsOfElement(k, element):
                    element_icfs[k] = v
        return element_icfs

    # first pass: write everything but the first line. count max width per slit
    max_col_widths = [0] * len(columns)
    lines = [[] for _ in unique_elements]

    # iterate for every element in unique_elements
    for line_nr, element in enumerate(unique_elements):
        # do not print abundancies for element H
        if element != "H":
            lines[line_nr].append("{:5s}".format(element))
            for j, col in enumerate(columns):
                # first write element total abundance
                entry = inDictOfElements(elem_abundancies[col], element)
                if entry is not None:
                    colstr = "{:15.9e}/{:15.9} ".format(
                        entry["abundance"], entry["uncertainty"]
                    )
                else:
                    colstr = "{:31s} ".format(" ")
                # write any ICFs/DIMS
                element_icfs = elementIcf(icfs[col], element)
                for k, v in element_icfs.items():
                    msg = "{:}:{:15.9e}/{:15.9e} ".format(k, v[0], v[1])
                    colstr += f"{msg}"
                max_col_widths[j] = max(max_col_widths[j], len(colstr) - 1)
                lines[line_nr].append(colstr)

    # remove empty sublists
    lines = [line for line in lines if line != []]

    # second pass, write header
    with open(fn, "w") as fout:
        print("Slit ", file=fout, end="")
        for j, c in enumerate(max_col_widths):
            print(
                "-----{:2d}{:s} ".format(columns[j], "-" * (c - 7)), file=fout, end=""
            )
        print("", file=fout)
        # print(lines)
        for j, line in enumerate(lines):
            print(f"{line[0]:<5s} ", file=fout, end="")
            for lc in zip(line[1:], max_col_widths):
                print(f"{lc[0]:<{lc[1]}},", file=fout, end="")
            print("", file=fout)
        print("", file=fout)
    return fn
