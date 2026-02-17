import re
from typing import Any, Dict, Tuple, Optional, List
import numpy as np
import pyneb as pn

import satellite.cfgio as sc
import satellite.roman as sr
import satellite.tene as tn
import satellite.nomenclature as sn
from satellite import intensity as si


def extract_wavelength(line_id: str) -> Optional[str]:
    """
    'O2_3729A' or 'O2_3729.7A' -> 'L(3729)' or 'L(3729.7)'
    """
    try:
        w = line_id.split("_")[-1]
        w = w.replace("A", "").strip()
        float(w)  # validate
        return f"L({w})"
    except Exception:
        return None


_LINE_PN_RE = re.compile(r"^(?P<ion>[A-Za-z0-9]+)_(?P<w>\d+(?:\.\d+)?)A$")


def ion_key_from_entry(entry: Dict[str, Any]) -> str:
    """
    Returns the ion key used by your intensities index, e.g.
      N + II -> 'N2'
      H + I (RL) -> 'H1r'
    """
    element = entry["element"]
    ion_int = sr.roman2int(entry["spectrum"])  # you already import roman
    ref_type = entry.get("ref_type", "CEL").upper()
    suffix = "r" if (ref_type == "RL" and element in {"H", "He"}) else ""
    return f"{element}{ion_int}{suffix}"


def matched_line_intensity_for_entry(
    idx: Dict[str, list],
    entry: Dict[str, Any],
    tol_A: float,
    logger=None,
) -> Tuple[str, float, float]:
    """
    Returns (pn_element_label, intensity, sigma_intensity) for this fit entry.
    """
    ion = ion_key_from_entry(entry)
    # target_A = float(entry["atomic"])  # your line wavelength (can be float)
    # (prefer the canonical PyNeb line wavelength)
    target_A = float(entry.get("pn_line", entry["atomic"]))

    row = si.match_line(idx, ion, target_A, tol_A, logger)  # <-- your existing function

    pn_element = row["element_pn"]  # e.g. 'N2_5754.6A'
    I_c = float(row["intensity"])
    rel = float(row["intensity_err"])  # relative
    sigma_I = abs(I_c) * rel  # absolute sigma

    if logger:
        logger.debug(
            f"Abundance line match: {entry['element']}{entry['spectrum']}@{target_A}A "
            f"-> {pn_element} (meas {row['wavelength_A']}A)"
        )
    return pn_element, I_c, sigma_I


def _parse_pn_line_label(pn_element: str) -> Tuple[str, float]:
    m = _LINE_PN_RE.fullmatch(pn_element.strip())
    if not m:
        raise ValueError(f"Bad pn line label: {pn_element!r}")
    return m.group("ion"), float(m.group("w"))


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

        # return pn.Atom(element, ion)
        # Validate that PyNeb actually has usable atomic data
        atom = pn.Atom(element, ion)
        try:
            w = getattr(atom, "wave_Ang", None)
            if w is None or not hasattr(w, "shape"):
                if logger:
                    logger.warning(
                        f"Skipping {element}{ion}: no usable PyNeb atomic data."
                    )
                return None
        except Exception:
            if logger:
                logger.warning(f"Skipping {element}{ion}: invalid PyNeb atomic model.")
            return None

        return atom

    else:
        logger.debug(
            f'Creating recombination atom from entry {element} {ion} (or {element}{dct_entry["spectrum"]}) for ionic abundancies. pn.RecAtom({element}, {ion})'
        )
        return pn.RecAtom(element, ion)


_PAIR_RE = re.compile(
    r"""^
    (?:                                   # option A: [Ion] expr [Ion] expr
        \[([A-Za-z0-9]+)\]\s*([0-9A-Za-z+]+(?:/[0-9A-Za-z+]+)?)\s+
        \[([A-Za-z0-9]+)\]\s*([0-9A-Za-z+]+(?:/[0-9A-Za-z+]+)?)
    |
        (\S+)\s+(\S+)                      # option B: token1 token2
    )
    $""",
    re.VERBOSE,
)


def refTenNe2PyNebPair(ref_tene):
    """
    Example:
        ref_tene="[SIII] 6312/9069 [ClIII] 5538/5518"
        ref_tene="string1 strings3"
    Return:
        "[SIII] 6312/9069", "[ClIII] 5538/5518"
        "string1", "string2"
    """
    m = _PAIR_RE.match(ref_tene.strip())
    if not m:
        raise ValueError(f"Unrecognized ref_tene format: {ref_tene!r}")

    if m.group(1) is not None:
        # matched the [Ion] expr [Ion] expr form
        return f"[{m.group(1)}] {m.group(2)}", f"[{m.group(3)}] {m.group(4)}"
    else:
        # matched the generic "token token" form
        return m.group(5), m.group(6)


def computeIonicAbundancies(
    fitsd,
    tene_dict,
    intensities_payload,
    min_percentage,
    logger,
    *,
    tol_A: float = 1.0,
    mc_N: int = 1000,
    seed: int = 0,
):
    ionic_abundancies_dict = []

    idx = si.build_intensity_index(intensities_payload)
    rng = np.random.default_rng(seed)

    def ref_tene_pair(tene):
        for entry in tene_dict:
            if entry["tene_pair"] == tene:
                return entry
        return None

    def inspect_mc_err_local(arr: List[float]) -> Tuple[np.ndarray, bool]:
        a = np.asarray(arr, dtype=float).ravel()
        ok = np.isfinite(a)
        frac = ok.mean() if a.size else 0.0
        if frac < min_percentage:
            return a[ok], False
        return a[ok], True

    for entry in fitsd:
        reftene = ref_tene_pair(refTenNe2PyNebPair(entry["ref_tene"]))
        if reftene is None:
            logger.error(
                f'ERROR Failed finding reference Te/Ne pair {entry["ref_tene"]} for ionic abundancies!'
            )
            return []

        pn_atom = get_atom_model(entry, logger)
        if pn_atom is None:
            logger.warning(
                f"Skipping abundance for {entry['element']}{entry['spectrum']}_{entry['atomic']} "
                f"(no usable atomic model)"
            )
            continue

        try:
            pn_element, I_c, sigma_I = matched_line_intensity_for_entry(
                idx, entry, tol_A, logger
            )
        except:
            logger.warning(
                f"Cannot match line for {entry['element']}{entry['spectrum']}_{entry['atomic']} aka {entry['pnstr']}, skipping."
            )
            continue

        # to_eval = extract_wavelength(pn_element)
        to_eval = extract_wavelength(pn_element)
        if to_eval is None:
            logger.warning(f"Cannot build to_eval for {pn_element}, skipping.")
            continue

        # central abundance
        sabd = pn_atom.getIonAbundance(
            int_ratio=I_c,
            tem=reftene["sT"],
            den=reftene["sN"],
            to_eval=to_eval,
            Hbeta=100.0,
        )

        if not np.isfinite(sabd):
            logger.warning(
                f"Central ionic abundance for {pn_element} is non-finite: {sabd!r} "
                f"(T={reftene['sT']}, Ne={reftene['sN']}, int_ratio={I_c}, to_eval={to_eval!r})"
            )
            continue

        # ---------- MC uncertainty (replaces pnErrObs.getIntens()[pn_element] + reftene['eT']['eN'] arrays) ----------
        # We only have scalar Te/Ne errors now (std), so we sample Te/Ne around the central values.
        eT = float(reftene.get("eT", np.nan))
        eN = float(reftene.get("eN", np.nan))

        I_mc = rng.normal(I_c, sigma_I, size=mc_N)
        I_mc = np.clip(I_mc, 0.0, None)

        if np.isfinite(eT) and eT > 0:
            T_mc = rng.normal(float(reftene["sT"]), eT, size=mc_N)
            T_mc = np.clip(T_mc, 1.0, None)
        else:
            T_mc = np.full(mc_N, float(reftene["sT"]))

        if np.isfinite(eN) and eN > 0:
            N_mc = rng.normal(float(reftene["sN"]), eN, size=mc_N)
            N_mc = np.clip(N_mc, 1e-6, None)
        else:
            N_mc = np.full(mc_N, float(reftene["sN"]))

        tar = []
        for k in range(mc_N):
            tar.append(
                pn_atom.getIonAbundance(
                    int_ratio=float(I_mc[k]),
                    tem=float(T_mc[k]),
                    den=float(N_mc[k]),
                    to_eval=to_eval,
                    Hbeta=100.0,
                )
            )

        eabd_array, ok = inspect_mc_err_local(tar)
        if not ok:
            logger.error(
                f"Too many non-finite MC getIonAbundance for {pn_element}. "
                f"Limit was {100.0*(1.0-min_percentage):.1f}% bad draws"
            )
            raise RuntimeError(
                f"Too many non-finite MC getIonAbundance for {pn_element}."
            )

        ionic_abundancies_dict.append(
            {
                "element": entry["element"],
                "spectrum": entry["spectrum"],
                "atomic": entry["atomic"],
                "pn_element": pn_element,
                "abundance": float(sabd),
                "abundance_error": float(np.std(eabd_array)),
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
            try:
                cpo = str(col + offset)
            except:
                cpo = str(col) + str(offset)
            print("{:->6s}{:25s} ".format(cpo, "-" * 25), file=fout, end="")
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
