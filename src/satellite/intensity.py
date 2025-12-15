import sys
import pyneb as pn
import numpy as np
from satellite import astroflux
import satellite.roman as sr
import satellite.cfgio as sc
import satellite.nomenclature as sn


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
