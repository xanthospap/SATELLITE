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
    try:
        return elemspec_abundancies["O3"] / (
            elemspec_abundancies["O2"] + elemspec_abundancies["O3"]
        )
    except:
        return elemspec_abundancies["O3"][0] / (
            elemspec_abundancies["O2"][0] + elemspec_abundancies["O3"][0]
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
    try:
        return elemspec_abundancies["He3"] / (
            elemspec_abundancies["He2"] + elemspec_abundancies["He3"]
        )
    except:
        return elemspec_abundancies["He3"][0] / (
            elemspec_abundancies["He2"][0] + elemspec_abundancies["He3"][0]
        )


def omega1435(omegav, uv) -> float:
    return omegav <= 0.5


def omega1436(omegav, uv) -> float:
    return omegav > 0.5


def omega1414(omegav, uv) -> float:
    return omegav <= 0.95


def omega1414b(omegav, uv) -> float:
    return omegav > 0.95


def omega1429b(omegav, uv) -> float:
    return omegav > 0.02 and omegav < 0.95


def omega1432(omegav, uv) -> float:
    return omegav > 0.95


def u1417c(omegav, uv) -> float:
    return uv < 0.015


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

    Example call:
    >>> IcfNameList()
    >>> ['KB94_A32', 'KB94_A30.10', 'DIMS14_35', 'DIMS14_36', 'KB94_A1.10', 'DIMS14_14', 'DIMS14_14b', 'KB94_A10', 'KB94_A8', 'KB94_A6', 'DIMS14_12', 'KB94_A36.10', 'KB94_A38.10', 'DIMS14_23', 'DIMS14_26', 'DIMS14_29b', 'DIMS14_32', 'KB94_A28.10', 'KB94_A27', 'DIMS14_17a', 'DIMS14_17b', 'DIMS14_17c', 'KB94_A12', 'KB94_A13.10', 'KB94_A16', 'KB94_A19', 'KB94_A21', 'KB94_A26', 'DIMS14_39', 'RR05_2', 'RR05_3', 'RR05_4']
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
        A dictionary where keys are elements, (e.g. 'HeI', 'OII', ...)

    Example
    --------
    >>> renameIons({'O1': np.float64(2.4326264721336407e-06), 'O2': np.float64(2.2954914792701307e-05), 'O3': np.float64(0.0002963497082918657), 'N1': np.float64(2.742370113432223e-07), 'N2': np.float64(2.3649163749694575e-06), 'He1': np.float64(0.11989790014907455), 'He2': np.float64(0.012375089973927573), 'Ar3': np.float64(7.696229948812361e-07), 'Cl3': np.float64(4.697062601788208e-08), 'H1': np.float64(1.0080578213512084), 'S2': np.float64(8.724656708931813e-08), 'S3': np.float64(1.4236060230320906e-06)})
    >>> {'O': np.float64(2.4326264721336407e-06), 'O2': np.float64(2.2954914792701307e-05), 'O3': np.float64(0.0002963497082918657), 'N': np.float64(2.742370113432223e-07), 'N2': np.float64(2.3649163749694575e-06), 'He2': np.float64(0.11989790014907455), 'He3': np.float64(0.012375089973927573), 'Ar3': np.float64(7.696229948812361e-07), 'Cl3': np.float64(4.697062601788208e-08), 'H2': np.float64(1.0080578213512084), 'S2': np.float64(8.724656708931813e-08), 'S3': np.float64(1.4236060230320906e-06)}
    """
    acopy = {}

    def es_split(elemsp):
        g = re.match(r"([A-Za-z]+)([0-9]?)", elemsp)
        return g[1], g[2]

    for k, v in abundancies.items():
        acopy[ionstr2pyneb(elemspec2ionstr(*es_split(k)))] = v

    return acopy


def computeIcfs(abundancies, logger=None):
    """
    Example
    -------
    >>> computeIcfs({'He1': (np.float64(0.11989790014907455), np.float64(0.0024290548109774427)), 'He2': (np.float64(0.012375089973927573), np.float64(0.00030644222010208356)), 'H1': (np.float64(1.0080578213512084), np.float64(0.010954363783206234)), 'S2': (np.float64(8.724656708931813e-08), np.float64(3.6368619630816025e-09)), 'S3': (np.float64(1.4236060230320906e-06), np.float64(6.474589949531195e-08)), 'N1': (np.float64(2.742370113432223e-07), np.float64(3.088339543010921e-08)), 'N2': (np.float64(2.3649163749694575e-06), np.float64(7.846853280518262e-08)), 'Ar3': (np.float64(7.696229948812361e-07), np.float64(5.118305350315916e-08)), 'Cl3': (np.float64(4.697062601788208e-08), np.float64(1.7198167692195738e-09)), 'O1': (np.float64(2.4326264721336407e-06), np.float64(1.703462882653106e-07)), 'O2': (np.float64(2.2954914792701307e-05), np.float64(2.122843805043223e-06)), 'O3': (np.float64(0.0002963497082918657), np.float64(1.3243705440243252e-05))})
    >>> {'KB94_A32': np.float64(4.416853398761392e-07), 'KB94_A10': np.float64(0.00020091915455175606), 'KB94_A8': np.float64(nan), 'KB94_A6': np.float64(nan), 'DIMS14_12': np.float64(0.0001992916808640836), 'KB94_A36.10': np.float64(6.402147642804062e-07), 'KB94_A38.10': np.float64(2.754948762720838e-07), 'DIMS14_23': np.float64(2.763227013462463e-07), 'DIMS14_26': np.float64(7.934917520559565e-07), 'DIMS14_29b': np.float64(nan), 'DIMS14_32': np.float64(nan), 'KB94_A28.10': np.float64(nan), 'KB94_A27': nan, 'DIMS14_17a': np.float64(nan), 'DIMS14_17b': np.float64(nan), 'DIMS14_17c': np.float64(nan), 'KB94_A12': np.float64(nan), 'KB94_A13.10': np.float64(nan), 'KB94_A16': np.float64(nan), 'KB94_A19': np.float64(nan), 'KB94_A21': np.float64(nan), 'KB94_A26': np.float64(nan), 'DIMS14_39': np.float64(nan), 'RR05_2': np.float64(nan), 'RR05_3': np.float64(nan), 'RR05_4': np.float64(nan), 'DIMS14_35': np.float64(3.2846467418925716e-07), 'DIMS14_36': np.float64(5.511584261396814e-07), 'KB94_A1.10': np.float64(5.058821131007465e-06), 'DIMS14_14': np.float64(4.000525659076789e-06), 'DIMS14_14b': np.float64(2.054823990511968e-05), 'KB94_A30.10': np.float64(nan)}
    """
    icf = pn.ICF()
    allIcf = icf.getElemAbundance(
        renameIons({a[0]: a[1][0] for a in abundancies.items()}), IcfNameList()
    )

    return allIcf


def computeManualIcfs(abundancies, logger=None):
    # compute all possible ICFs
    ne = computeIcfs(abundancies, logger)

    # get total abundancies
    total_abundancies = ionicAbundance2elementAbundance(abundancies)

    # PyNeb
    icf = pn.ICF()
    renamed_abunds = renameIons({k: v[0] for k, v in abundancies.items()})
    renamed_errors = renameIons({k: v[1] for k, v in abundancies.items()})

    # get omega and u
    eo = omega(renamed_abunds)
    eu = U(renamed_abunds)

    results = {}

    for elem, allicfs in icf_dict.items():
        if elem in total_abundancies:
            # total abundance (and error) for element
            ta_val, ta_err = (
                total_abundancies[elem]["abundance"],
                total_abundancies[elem]["uncertainty"],
            )
            # walk through ICF's for this element ...
            for tp, lst in allicfs.items():
                for icfdetails in lst:
                    if icfdetails["cond"](eo, eu):
                        logger.debug(
                            f'Element {elem} computing icf {icfdetails["name"]}'
                        )
                        try:
                            icf_val = icf.getElemAbundance(
                                renamed_abunds, [icfdetails["name"]]
                            )[icfdetails["name"]]
                            print(f"\t>>> value: {icf_val}")
                            if not (np.isnan(icf_val) or isnan(icf_val)):
                                results[icfdetails["name"]] = (icf_val, 0e0)
                        except:
                            pass
                    else:
                        logger.debug(
                            f'Skipping computation of {icfdetails["name"]} for element {elem} due to omega, U values ({eo:.4f}, {eu:.4f})'
                        )
    return results


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

    koko = {k: v[0] for k, v in abundancies.items()}
    print(f"{koko}")
    print(f"{renamed_abunds}")

    results = {}

    # get omega and u
    eo = omega(renamed_abunds)
    eu = U(renamed_abunds)

    for icf_name, total_abund in ne.items():
        if not np.isfinite(total_abund):
            continue

        element = getElementFromIcfName(icf_name)

        # Find all contributing ions for this element
        contributing_ions = [ion for ion in renamed_abunds if ion.startswith(element)]
        # if not contributing_ions:
        #     continue
        logger.debug(f"Contributing ions for {element}: {contributing_ions}")

        ionic_sum = sum(renamed_abunds[ion] for ion in contributing_ions)
        # if ionic_sum == 0:
        #    continue

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

    results["omega"] = eo
    results["U"] = eu
    return results


def ionicAbundance2elementAbundance(abundancies, logger=None):
    """
    Example
    -------
    >>> ionicAbundance2elementAbundance({'He1': (np.float64(0.07289117410488913), np.float64(0.0014473831132186096)), 'He2': (np.float64(0.025174496324910436), np.float64(0.00038447148833430333)), 'N1': (np.float64(2.4613015709276348e-08), np.float64(5.312011544474936e-09)), 'N2': (np.float64(5.097540067692014e-07), np.float64(5.666130804139922e-08)), 'S2': (np.float64(3.059625788523793e-08), np.float64(3.678995609088767e-09)), 'S3': (np.float64(8.109451468778343e-07), np.float64(1.889202176560218e-08)), 'H1': (np.float64(1.010407099346665), np.float64(0.0059810498635109646)), 'Cl3': (np.float64(2.166739769269426e-08), np.float64(7.071592860594983e-10)), 'Ar3': (np.float64(5.035899412128091e-07), np.float64(1.5991809652913825e-08)), 'O1': (np.float64(2.5751403644020295e-08), np.float64(1.0104436788415689e-08)), 'O2': (np.float64(1.1789314340737157e-05), np.float64(4.635885046610191e-06)), 'O3': (np.float64(0.0002431990283552331), np.float64(5.314201887002895e-06))})
    >>> {'He': {'abundance': np.float64(0.09806567042979956), 'uncertainty': np.float64(0.001497576776586893)}, 'N': {'abundance': np.float64(5.343670224784777e-07), 'uncertainty': np.float64(5.690976450145411e-08)}, 'S': {'abundance': np.float64(8.415414047630721e-07), 'uncertainty': np.float64(1.9246908715003583e-08)}, 'H': {'abundance': np.float64(1.010407099346665), 'uncertainty': np.float64(0.0059810498635109646)}, 'Cl': {'abundance': np.float64(2.166739769269426e-08), 'uncertainty': np.float64(7.071592860594983e-10)}, 'Ar': {'abundance': np.float64(5.035899412128091e-07), 'uncertainty': np.float64(1.5991809652913825e-08)}, 'O': {'abundance': np.float64(0.0002550140940996143), 'uncertainty': np.float64(7.052111312284563e-06)}}
    """

    def stripElement(elemspec):
        g = re.match(r"([A-Za-z]+)[0-9]+", elemspec)
        return g[1]

    elemdct = {}
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
            if k != "omega" and k != "U":
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

        # print omega and u values
        print(f"Omega ", end="", file=fout)
        for j, c in enumerate(max_col_widths):
            print(f"{icfs[j]['omega']:30.9e}", end="", file=fout)
        print("", file=fout)
        print(f"U     ", end="", file=fout)
        for j, c in enumerate(max_col_widths):
            print(f"{icfs[j]['U']:30.9e}", end="", file=fout)
        print("", file=fout)

        # print(lines)
        for j, line in enumerate(lines):
            print(f"{line[0]:<5s} ", file=fout, end="")
            for lc in zip(line[1:], max_col_widths):
                print(f"{lc[0]:<{lc[1]}},", file=fout, end="")
            print("", file=fout)
        print("", file=fout)
    return fn
