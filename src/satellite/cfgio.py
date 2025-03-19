import yaml
import os
import re
import sys

import pyneb as pn
from satellite import roman
import satellite.nomenclature as sn


def parseConfigInout(fn: str):
    """
    Parse an input SATELLITE configuration file (YAML format) and return
    the corresponding options as a dictionary.
    """
    with open(fn, "r") as fin:
        ymldata = fin.read()
    return yaml.safe_load(ymldata)


def configElements(dct: dict):
    """
    Once we have parsed the YAML config file into a dictionary, we can use
    the following function to get a list of all the elements that are included.

    Results should be something like:
    ['H', 'He', 'N', 'O']
    """
    return [k["element"] for k in dct["data_list"]["element_list"]]


def configFitsFileList(dct: dict):
    """
    Once we have parsed the YAML config file into a dictionary, we can use
    the following function to get a list of all the FITS images and related
    info that should be present.

    Results should be something like:
    [
      {'element': 'H', 'spectrum' : 'i',  'atomic': 6563, 'fns': '/home/.../Hi_6563s.fit', 'fne': '/home/.../Hi_6563e.fit'}
      {'element': 'H', 'spectrum' : 'i',  'atomic': 4861, 'fns': '/home/.../Hi_4861s.fit', 'fne': '/home/.../Hi_4861e.fit'}
      {'element': 'H', 'spectrum' : 'i',  'atomic': 4340, 'fns': '/home/.../Hi_4340s.fit', 'fne': '/home/.../Hi_4340e.fit'}
      {'element': 'H', 'spectrum' : 'i',  'atomic': 4101, 'fns': '/home/.../Hi_4101s.fit', 'fne': '/home/.../Hi_4101e.fit'}
      {'element': 'He', 'spectrum': 'i',  'atomic': 5876, 'fns': '/home/.../Hei_5876s.fit','fne': '/home/.../Hei_5876e.fit'}
    ]
    """
    fns = []
    did = dct["data_list"]["data_img_id"]
    eid = dct["data_list"]["error_img_id"]
    for element in dct["data_list"]["element_list"]:
        for spec in element["spectrums"]:
            for atomic in spec["atomic_numbers"]:
                sfn = (
                    element["atom"]
                    + spec["spectrum"]
                    + "_"
                    + str(atomic["atomic_number"])
                    + did
                    + "."
                    + dct["data_list"]["suffix"]
                )
                efn = (
                    element["atom"]
                    + spec["spectrum"]
                    + "_"
                    + str(atomic["atomic_number"])
                    + eid
                    + "."
                    + dct["data_list"]["suffix"]
                )
                fns.append(
                    {
                        "element": element["atom"],
                        "spectrum": spec["spectrum"].lower(),
                        "atomic": atomic["atomic_number"],
                        "ref_tene": atomic["ref_TeNe"],
                        "fns": os.path.join(dct["data_list"]["prefix"], sfn),
                        "fne": os.path.join(dct["data_list"]["prefix"], efn),
                    }
                )
    return fns


def elementInputDict(fitsd: list, reference_element: dict, logger=None):
    for obj in fitsd:
        matched = True
        for key in ["element", "spectrum", "atomic"]:
            if obj[key] != reference_element[key]:
                matched = False
        if matched == True:
            return obj
    if logger is not None:
        logger.error(
            "Failed matching element {:}/{:}/{:} to available input file".format(
                reference_element["element"],
                reference_element["spectrum"],
                reference_element["atomic"],
            )
        )
    raise RuntimeError(
        "[ERROR] Failed matching element {:}/{:}/{:} to available input file".format(
            reference_element["element"],
            reference_element["spectrum"],
            reference_element["atomic"],
        )
    )


def configSpecificSlitAnalysis(dct: dict):
    d = dct["analysis"]["specific_slit_analysis"]
    # no specific-slit analysis
    if d["skip"] in ["1", "True", "true"]:
        return None
    # return a list of dictionaries, one for each slit
    slits = []
    for slit in d["slits"]:
        ps = [int(x) for x in slit.split(",")]
        slits.append({"PA": ps[0], "w": ps[1], "h": ps[2], "x": ps[3], "y": ps[4]})
    return slits


def checkInputFits(fitsd: list, logger=None):
    fitsd_out = []
    missing_files = []
    for idx, fits in enumerate(fitsd):
        obs_or_error_missing = False
        for ftype in ["fns", "fne"]:
            file_is_missing = False
            fitsfn = fits[ftype]
            if not os.path.isfile(fitsfn):
                if logger:
                    logger.warning("Missing Fits file {:}".format(fitsfn))
                bn = os.path.basename(fitsfn)
                file_is_missing = True
                # find spectrum in filename
                sstr = fits["spectrum"]
                # replace roman spectrum with int and see if file exists
                if bn.find(sstr) >= 0:
                    gfn = bn.replace(sstr, str(roman.roman2int(sstr)), 1)
                    guess = os.path.join(os.path.dirname(fitsfn), gfn)
                    if os.path.isfile(guess):
                        # change the name in the return list
                        if logger:
                            logger.debug(
                                "Fits filename {:} is missing; using {:}".format(
                                    os.path.basename(fitsfn), gfn
                                )
                            )
                        fitsd[idx][ftype] = guess
                        file_is_missing = False
                if file_is_missing and bn.find(sstr) >= 0:
                    gfn = bn.replace(sstr, sstr.upper(), 1)
                    guess = os.path.join(os.path.dirname(fitsfn), gfn)
                    if os.path.isfile(guess):
                        if logger:
                            logger.debug(
                                "Fits filename {:} is missing; using {:}".format(
                                    os.path.basename(fitsfn), gfn
                                )
                            )
                        fitsd[idx][ftype] = guess
                        file_is_missing = False
            if file_is_missing:
                missing_files.append(fitsfn)
                obs_or_error_missing = True
        if not obs_or_error_missing:
            fitsd_out.append(fits)
    return fitsd_out, missing_files


def indexOf(atom: str, spectrum: str, atomic_number: int, fitsd: list):
    for idx, lst in enumerate(fitsd):
        if (
            lst["element"] == atom
            and lst["spectrum"] == spectrum
            and lst["atomic"] == atomic_number
        ):
            return idx
    return -1


def configElementRatiosList(dct: dict):
    return dct["analysis"]["specific_slit_analysis"]["log_ratios"]


def configDensityDiagnostics(dct: dict):
    return dct["analysis"]["specific_slit_analysis"]["density_diagnostics"]


def configTemperatureDiagnostics(dct: dict):
    return dct["analysis"]["specific_slit_analysis"]["temperature_diagnostics"]


def configRefTenNe2PyNebPair(ref_tene):
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
