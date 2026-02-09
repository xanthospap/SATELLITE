import copy
import re
import math
import pyneb as pn
import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional, Any

from satellite import satdebug as db
from satellite import intensity
from satellite import roman


# ---- parse the built-in diagnostic expressions: e.g.
# "(L(6548)+L(6584))/L(5755)"  -> list of waves in numerator and denominator
_L_RE = re.compile(r"L\(\s*(\d+(?:\.\d+)?)\s*\)")


def _split_top_level_div(expr: str) -> Tuple[str, str]:
    depth = 0
    for i, ch in enumerate(expr):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == "/" and depth == 0:
            return expr[:i], expr[i + 1 :]
    raise ValueError(f"No top-level '/' in expression: {expr!r}")


def _waves_from_expr(expr: str) -> Tuple[List[float], List[float]]:
    num_expr, den_expr = _split_top_level_div(expr)
    num = [float(x) for x in _L_RE.findall(num_expr)]
    den = [float(x) for x in _L_RE.findall(den_expr)]
    if not num or not den:
        raise ValueError(f"Could not extract waves from: {expr!r}")
    return num, den


def sum_intensities(
    idx: Dict[str, List[Dict[str, Any]]],
    ion: str,
    waves_A: List[float],
    tol_A: float,
    logger,
) -> Tuple[float, float, List[Tuple[float, str, float]]]:
    """
    Sum intensities for a list of target wavelengths (Å).
    Error combined in quadrature (assumes independent).
    Returns (sum_I, sum_err, matches) where matches records what got matched.
    """
    I = 0.0
    e2 = 0.0
    matches: List[Tuple[float, str, float]] = []
    for w in waves_A:
        row = intensity.match_line(idx, ion, w, tol_A, logger)
        Ii = float(row["intensity"])
        rel = float(row["intensity_err"])  # relative (fractional) error
        si = abs(Ii) * rel  # absolute sigma
        I += Ii
        e2 += si * si
        matches.append((w, row["element_pn"], float(row["wavelength_A"])))
    return I, math.sqrt(e2), matches


def ratio_and_err(
    num: float, den: float, num_err: float, den_err: float
) -> Tuple[float, float]:
    """
    R = num/den, error via standard propagation (independent).
    """
    if num <= 0.0 or den <= 0.0:
        return float("nan"), float("nan")
    R = num / den
    rel2 = 0.0
    if num_err > 0 and num > 0:
        rel2 += (num_err / num) ** 2
    if den_err > 0 and den > 0:
        rel2 += (den_err / den) ** 2
    return R, abs(R) * math.sqrt(rel2)


def observed_ratio_for_diag(
    diag_label: str,
    idx: Dict[str, List[Dict[str, Any]]],
    tol_A: float,
    logger=None,
) -> Tuple[float, float, Dict[str, Any]]:
    """
    Example in:
    ------------
    diag_label=[NII] 5755/6584+

    Uses pn.diags_dict to expand '+' (because it's encoded in the expression),
    then matches measured lines by nearest wavelength.
    """
    if diag_label not in pn.diags_dict:
        raise KeyError(
            f"{diag_label!r} not found in pn.diags_dict. "
            f"If it's custom, you must addDiag(label, tuple) first."
        )

    ion, expr, _err_expr = pn.diags_dict[diag_label]
    # e.g. ion=N2, expr=L(5755)/(L(6548)+L(6584)),
    # _err_expr=RMS([E(6548)*L(6548)/(L(6548)+L(6584)), E(6584)*L(6584)/(L(6584)+L(6548)), E(5755)])

    num_w, den_w = _waves_from_expr(expr)
    # e.g. num_w=[5755.0], den_w=[6548.0, 6584.0]
    logger.info(
        f"Computing diagnostic from {expr}; nominator waves:{num_w}, denominator waves:{den_w}"
    )

    # sum intensities for the different lines in nominator and denominator
    num_I, num_e, num_matches = sum_intensities(idx, ion, num_w, tol_A, logger)
    den_I, den_e, den_matches = sum_intensities(idx, ion, den_w, tol_A, logger)

    R, Rerr = ratio_and_err(num_I, den_I, num_e, den_e)
    meta = {
        "ion": ion,
        "expr": expr,
        "num_waves": num_w,
        "den_waves": den_w,
        "num_matches": num_matches,
        "den_matches": den_matches,
        "num_sum": num_I,
        "den_sum": den_I,
    }
    logger.info(f"Macthed waves in nominator   (summed): {[n[1] for n in num_matches]}")
    logger.info(f"Macthed waves in denominator (summed): {[n[1] for n in den_matches]}")

    return R, Rerr, meta


# does exactly the same as observed_ratio_for_diag but for custom diagnostics
def observed_ratio_for_custom_expr(
    *,
    label: str,
    ion: str,
    expr: str,
    idx,
    tol_A: float,
    logger=None,
):
    num_w, den_w = _waves_from_expr(expr)

    if logger:
        logger.info(
            f"Computing CUSTOM diagnostic {label} from {expr}; "
            f"num={num_w}, den={den_w}"
        )

    num_I, num_e, num_matches = sum_intensities(idx, ion, num_w, tol_A, logger)
    den_I, den_e, den_matches = sum_intensities(idx, ion, den_w, tol_A, logger)

    R, Rerr = ratio_and_err(num_I, den_I, num_e, den_e)
    meta = {
        "ion": ion,
        "expr": expr,
        "num_waves": num_w,
        "den_waves": den_w,
        "num_matches": num_matches,
        "den_matches": den_matches,
        "num_sum": num_I,
        "den_sum": den_I,
    }
    return R, Rerr, meta


def computeTeNePairs(
    density_diagnostics: List[str],
    temperature_diagnostics: List[str],
    intensities_payload: Dict[str, Any],
    tol_A: float,
    logger=None,
    *,
    mc_N: int = 1000,
    seed: int = 0,
    min_percentage: float = 0.7,  # like our old inspectMcErr threshold
) -> Dict[str, Any]:
    """
    - add all diags to pn.Diagnostics
    - compute observed ratios from your intensity list (no Observation)
    - call getCrossTemDen(value_tem=..., value_den=...) for each pair
    """

    def _normalize_diag_entry(d):
        # returns (label, ion, expr) where ion/expr may be None for built-ins
        if isinstance(d, str):
            return d, None, None
        if isinstance(d, dict):
            return d["label"], d["ion"], d["expr"]
        raise TypeError(f"Unsupported diagnostic entry: {d!r}")

    def _inspect_mc_err(arr: np.ndarray) -> Tuple[float, bool]:
        a = np.asarray(arr, dtype=float).ravel()
        ok = np.isfinite(a)
        frac = ok.mean() if a.size else 0.0
        if frac < min_percentage:
            if logger:
                logger.warning(
                    f"MC coverage too low: {frac:.3f} < {min_percentage:.3f}"
                )
            return float("nan"), False
        return float(np.nanstd(a)), True

    # for MC simulations
    rng = np.random.default_rng(seed)

    # build index
    idx = intensity.build_intensity_index(intensities_payload)

    temp_entries = [_normalize_diag_entry(d) for d in temperature_diagnostics]
    dens_entries = [_normalize_diag_entry(d) for d in density_diagnostics]

    # keep the original ordering, but unique by label
    def _unique_by_label(entries):
        seen = set()
        out = []
        for label, ion, expr in entries:
            if label not in seen:
                out.append((label, ion, expr))
                seen.add(label)
        return out

    all_entries = _unique_by_label(temp_entries + dens_entries)

    # collect/combine user-defined diagnostics
    # diags = pn.Diagnostics()
    # all_diags = list(dict.fromkeys(temperature_diagnostics + density_diagnostics))
    # for d in all_diags:
    #    diags.addDiag(d)  # built-in labels like [SII], [NII], etc.
    diags = pn.Diagnostics()
    for label, ion, expr in all_entries:
        if ion is None:
            diags.addDiag(label)  # built-in
        else:
            # Option A (often works): inject into pn.diags_dict, then add by label
            pn.diags_dict[label] = (ion, expr, None)
            diags.addDiag(label)

    # compute ratios once
    # ratio_info: Dict[str, Dict[str, Any]] = {}
    # for d in all_diags:
    #    R, Rerr, meta = observed_ratio_for_diag(d, idx, tol_A, logger)
    #    ratio_info[d] = {"R": float(R), "Rerr": float(Rerr), **meta}
    ratio_info = {}
    for label, ion, expr in all_entries:
        if ion is None:
            R, Rerr, meta = observed_ratio_for_diag(label, idx, tol_A, logger)
        else:
            R, Rerr, meta = observed_ratio_for_custom_expr(
                label=label, ion=ion, expr=expr, idx=idx, tol_A=tol_A, logger=logger
            )
        ratio_info[label] = {"R": float(R), "Rerr": float(Rerr), **meta}

    # to be returned ...
    tene_slit_dict: List[Dict[str, Any]] = []

    temp_labels = [x[0] for x in temp_entries]
    dens_labels = [x[0] for x in dens_entries]
    # solve Te/Ne for every T x n combination
    for t in temp_labels:
        for n in dens_labels:
            Rt = ratio_info[t]["R"]
            Rn = ratio_info[n]["R"]
            Rt_err = ratio_info[t]["Rerr"]
            Rn_err = ratio_info[n]["Rerr"]

            # skip invalid ratios
            if not (np.isfinite(Rt) and np.isfinite(Rn)):
                if logger:
                    logger.warning(
                        f"Skipping {(t, n)} due to invalid ratio(s): Rt={Rt}, Rn={Rn}"
                    )
                continue

            # central Te/Ne (sT, sN)
            st, sn = diags.getCrossTemDen(t, n, value_tem=Rt, value_den=Rn)
            st = float(np.atleast_1d(st)[0])
            sn = float(np.atleast_1d(sn)[0])

            # MC on ratios -> MC Te/Ne
            Rt_mc = rng.normal(Rt, Rt_err, size=mc_N)
            Rn_mc = rng.normal(Rn, Rn_err, size=mc_N)
            # ratios should be positive
            Rt_mc = np.clip(Rt_mc, 1e-30, None)
            Rn_mc = np.clip(Rn_mc, 1e-30, None)
            Te_mc, Ne_mc = diags.getCrossTemDen(t, n, value_tem=Rt_mc, value_den=Rn_mc)

            et, okT = _inspect_mc_err(Te_mc)
            en, okN = _inspect_mc_err(Ne_mc)
            if not (okT and okN):
                logger.warning(f"MC diagnostics failed for TeNe pair {t}-{n}")

            tene_slit_dict.append(
                {
                    "tene_pair": (t, n),
                    "sT": st,
                    "sN": sn,
                    "eT": et,
                    "eN": en,
                }
            )

    return tene_slit_dict, (okT and okN)


def printDiagnostics(dict_of_diagnostics, fn, logger):
    def pair2str(diag):
        return "{:}/{:}".format(diag[0], diag[1])

    def inDictOf(diags_list, diag):
        for d in diags_list:
            if d["tene_pair"] == diag:
                return d
        return None

    # extract unique diagnostics and store columns
    unique_diags = []
    columns = []
    for k, v in dict_of_diagnostics.items():
        columns += [k]
        unique_diags = list(set(unique_diags + [x["tene_pair"] for x in v]))

    # sort rows & columns
    unique_diags = sorted(unique_diags)
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

        # write first line
        print("{:45s}".format(" "), file=fout, end="")
        for col in columns:
            print("{:26s} {:26s} ".format("Temperature", "Density"), file=fout, end="")
        print("", file=fout)

        # write second line, i.e. column keys
        print("{:45s}".format("Slit Nr."), file=fout, end="")
        for col in columns:
            print("{:->6d}{:47s} ".format(col + offset, "-" * 47), file=fout, end="")
        print("", file=fout)

        # iterate for every diagnostic in unique_diags
        for diag in unique_diags:
            print("{:45s}".format(pair2str(diag)), file=fout, end="")
            for col in columns:
                entry = inDictOf(dict_of_diagnostics[col], diag)
                if entry is not None:
                    print(
                        "{:12.6e} {:12.6e} / {:12.6e} {:12.6e} ".format(
                            entry["sT"],
                            # err2scalar(entry["eT"], logger),
                            entry["eT"],
                            entry["sN"],
                            # err2scalar(entry["eN"], logger),
                            entry["eN"],
                        ),
                        file=fout,
                        end="",
                    )
                else:
                    print("{:53s} ".format(" "), file=fout, end="")
            print("", file=fout)
