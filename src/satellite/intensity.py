import sys
import pyneb as pn
import numpy as np
from typing import List, Tuple, Dict, Any
from typing import Optional

from satellite import astroflux
import satellite.roman as sr
import satellite.cfgio as sc
import satellite.nomenclature as sn


def match_line(
    idx: Dict[str, List[Dict[str, Any]]],
    ion: str,
    target_wave_A: float,
    tol_A: float,
    logger=None,
) -> Dict[str, Any]:
    """
    Find the closest measured line for a given ion and target wavelength.
    Uses precomputed wavelength_A field; works with fractional labels.
    """
    rows = idx.get(ion, [])
    if not rows:
        raise KeyError(f"No measured lines for ion {ion}")

    best = min(rows, key=lambda r: abs(float(r["wavelength_A"]) - target_wave_A))
    diff = abs(float(best["wavelength_A"]) - target_wave_A)
    if diff > tol_A:
        raise KeyError(
            f"No {ion} line within {tol_A}A of {target_wave_A} (closest diff={diff:.3f}A)"
        )
    return best


def build_intensity_index(
    intensities_payload: Dict[str, Any],
) -> Dict[str, List[Dict[str, Any]]]:
    """
    ---- build an index: ion -> list of line rows (sorted by wavelength)

    Example in:
    -----------------------------------------------------------------------------
    intensities_payload = [{'element_pn': 'Fe2_5333.3A', 'wavelength_A': 5333.3,
    'intensity': 0.2754380852937433, 'intensity_err': 0.2655770619186089}, {...}]


    idx={'H1r': [{'element_pn': 'H1r_4861A', 'wavelength_A': 4861.0, 'intensity': 100.0, 'intensity_err': 0.0},
    {'element_pn': 'H1r_6563A', 'wavelength_A': 6563.0, 'intensity': 285.0, 'intensity_err': 1.9945059249406321e-16}],

    'He1r': [{'element_pn': 'He1r_4922A', 'wavelength_A': 4922.0, 'intensity': 1.1164356287394488, 'intensity_err': 0.24307189639613375},
    {'element_pn': 'He1r_5016A', 'wavelength_A': 5016.0, 'intensity': 0.00029613419757631315, 'intensity_err': 0.027502310616225586},
    {'element_pn': 'He1r_5876A', 'wavelength_A': 5876.0, 'intensity': 10.860156981173906, 'intensity_err': 0.032491327185563815},
    {'element_pn': 'He1r_6678A', 'wavelength_A': 6678.0, 'intensity': 3.0653452035057693, 'intensity_err': 0.05025473346631531},
    {'element_pn': 'He1r_7065A', 'wavelength_A': 7065.0, 'intensity': 1.7520544059314276, 'intensity_err': 0.05204130629965182},
    {'element_pn': 'He1r_7281A', 'wavelength_A': 7281.0, 'intensity': 0.7246900823217778, 'intensity_err': 0.1587545571323489}],

    'He2r': [{'element_pn': 'He2r_5411A', 'wavelength_A': 5411.0, 'intensity': 0.2553722409536337, 'intensity_err': 0.43714336293822537}],

    'N1': [{'element_pn': 'N1_5197.9A', 'wavelength_A': 5197.9, 'intensity': 0.7512108730829774, 'intensity_err': 0.32354189843424747}],

    'N2': [{'element_pn': 'N2_5754.6A', 'wavelength_A': 5754.6, 'intensity': 0.9708804818911101, 'intensity_err': 0.12154216685657077},
    {'element_pn': 'N2_6548A', 'wavelength_A': 6548.0, 'intensity': 33.75846323002547, 'intensity_err': 0.03207146748490434},
    {'element_pn': 'N2_6583.5A', 'wavelength_A': 6583.5, 'intensity': 105.06208897767499, 'intensity_err': 0.02077326030239454}],

    'O1': [{'element_pn': 'O1_6300.3A', 'wavelength_A': 6300.3, 'intensity': 3.017453953521789, 'intensity_err': 0.06428001838601238},
    {'element_pn': 'O1_6363.8A', 'wavelength_A': 6363.8, 'intensity': 0.954158013758621, 'intensity_err': 0.11602159109986818},
    {'element_pn': 'O1_6391.7A', 'wavelength_A': 6391.7, 'intensity': 0.41173526112965403, 'intensity_err': 0.08581529065355169}],

    'O2': [{'element_pn': 'O2_7320A', 'wavelength_A': 7320.0, 'intensity': 2.2260970489807304, 'intensity_err': 1.8386362370122522},
    {'element_pn': 'O2_7329.7A', 'wavelength_A': 7329.7, 'intensity': 1.8713260957442526, 'intensity_err': 0.04999735850071672}],

    'O3': [{'element_pn': 'O3_4958.9A', 'wavelength_A': 4958.9, 'intensity': 17.235278845414193, 'intensity_err': 0.03773048312571993},
    {'element_pn': 'O3_5006.8A', 'wavelength_A': 5006.8, 'intensity': 48.9291228149472, 'intensity_err': 0.036392822826505536}],

    'S2': [{'element_pn': 'S2_6716.4A', 'wavelength_A': 6716.4, 'intensity': 24.699821355804577, 'intensity_err': 0.026962961056976917},
    {'element_pn': 'S2_6730.8A', 'wavelength_A': 6730.8, 'intensity': 28.416030881947997, 'intensity_err': 0.020001730247727608}],

    'S3': [{'element_pn': 'S3_6312.1A', 'wavelength_A': 6312.1, 'intensity': 0.7694181747701979, 'intensity_err': 0.16082353231862348},
    {'element_pn': 'S3_9068.6A', 'wavelength_A': 9068.6, 'intensity': 19.2744009654915, 'intensity_err': 0.05108465107917644}],

    'Ar3': [{'element_pn': 'Ar3_7135.8A', 'wavelength_A': 7135.8, 'intensity': 7.243554122355608, 'intensity_err': 0.028024407544087658},
    {'element_pn': 'Ar3_7751.1A', 'wavelength_A': 7751.1, 'intensity': 1.707043702456526, 'intensity_err': 0.0707610172499901}],
    """
    idx: Dict[str, List[Dict[str, Any]]] = {}
    for row in intensities_payload["intensities"]:
        lab = row["element_pn"]
        if "_" not in lab:
            continue
        ion = lab.split("_", 1)[0]  # e.g. "N2", "S2", "Cl3"
        idx.setdefault(ion, []).append(row)

    for ion in idx:
        idx[ion].sort(key=lambda r: float(r["wavelength_A"]))
    return idx


def makeIntensitiesDataFile(
    fitsd: list,
    reference_element: dict,
    value_keys: list,
    fn: str,
    factor=100e0,
    logger=None,
) -> str:
    ref_index = sc.indexOf(
        reference_element["element"],
        reference_element["spectrum"],
        reference_element["atomic"],
        fitsd,
    )

    ref_sval = fitsd[ref_index][value_keys[0]]
    ref_eval = fitsd[ref_index][value_keys[1]]
    with open(fn, "w") as fout:
        print("LINE test err", file=fout)
        for idx, obj in enumerate(fitsd):
            val = factor * obj[value_keys[0]] / ref_sval
            err = astroflux.fluxError(
                obj[value_keys[0]], obj[value_keys[1]], ref_sval, ref_eval
            )
            print(f"{obj['pnstr']} {val:+9e} {err:+9e}", file=fout)
    return fn


def computeIntensities(
    fitsd: dict, pnObs, pnErrObs, pnRC, reference_element: dict, logger
) -> list:
    ref_index = sc.indexOf(
        reference_element["element"],
        reference_element["spectrum"],
        reference_element["atomic"],
        fitsd,
    )
    ref_pnstr = fitsd[ref_index]["pnstr"]
    reference_element = fitsd[ref_index]

    iref = float(pnObs.getIntens()[ref_pnstr])
    eref = float(pnObs.getError()[ref_pnstr])
    intensities_list = []
    for idx, fits in enumerate(fitsd):
        scor = pnRC.getCorr(fits["pn_line"], reference_element["pn_line"])
        ecor = pnRC.getErrCorr(
            fits["pn_line"],
            np.std(pnErrObs.extinction.E_BV),
            reference_element["pn_line"],
        )
        pnstr = fits["pnstr"]
        iele = float(pnObs.getIntens()[pnstr])
        eele = float(pnObs.getError()[pnstr])
        err = np.sqrt(eele**2 + eref**2 + float(ecor / scor) ** 2)
        intensities_list.append(
            {
                "element_pn": pnstr,
                "element": "{:}{:}_{:}".format(
                    fits["element"], sr.roman2int(fits["spectrum"]), fits["atomic"]
                ),
                "intensity": iele,
                "intensity_err": err,
            }
        )
    return intensities_list


def printIntensities(dict_of_intensities, fn, logger):
    """Example: global_intensities[1] = {'intensities': [...], 'E_BV': rc.E_BV, 'cHbeta': rc.cHbeta}"""
    element_format = "element_pn"

    # given an (inner) dictionary of 'intensities' return the specific
    # dictionary of a line, if it exists, else None
    def inDictOf(intensitiesDict, line):
        for entry in intensitiesDict:
            if entry[element_format] == line:
                return entry
        return None

    # extract unique lines and store columns
    unique_lines = []
    columns = []
    for k, v in dict_of_intensities.items():
        aplist = []
        columns += [k]
        for entry in v["intensities"]:
            aplist += [entry["element_pn"]]
        unique_lines = list(set(unique_lines + aplist))

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
            print("{:15s}".format(line), file=fout, end="")
            for col in columns:
                entry = inDictOf(dict_of_intensities[col]["intensities"], line)
                if entry is not None:
                    print(
                        "{:15.9e} {:15.9e} ".format(
                            entry["intensity"], entry["intensity_err"]
                        ),
                        file=fout,
                        end="",
                    )
                else:
                    print("{:31s} ".format(" "), file=fout, end="")
            print("", file=fout)

        # cHbeta line
        print("{:15s}".format("cHbeta"), file=fout, end="")
        for col in columns:
            print(
                "{:15.9e} {:15.9e} ".format(
                    dict_of_intensities[col]["cHbeta"],
                    dict_of_intensities[col]["cHbetaError"],
                ),
                file=fout,
                end="",
            )
        # Total F(Hb) line
        print("\n{:15s}".format("F(Hb)"), file=fout, end="")
        for col in columns:
            print(
                "{:15.9e} {:15.9e} ".format(
                    dict_of_intensities[col]["FHb"],
                    dict_of_intensities[col]["FHb_error"],
                ),
                file=fout,
                end="",
            )
    return fn
