import pyneb as pn
import numpy as np
import copy
from satellite import satdebug as db
from satellite import roman
import re
from dataclasses import dataclass
from typing import List, Tuple, Dict
import math
from typing import Optional

MIN_VALID_PERCENTAGE = 70.0  # e.g. require at least 70% non-NaN values


def inspectMcErr(err_array, min_percentage=MIN_VALID_PERCENTAGE, logger=None):
    """
    Convert an array of MC error realisations into a single scalar error.

    - If the fraction of non-NaN values is below `min_percentage`,
      we treat the diagnostic as failed and return (None, False).
    - Otherwise we compute the std dev over the non-NaN values and
      return (sigma, True).
    """
    arr = np.asarray(err_array, dtype=float)

    if arr.size == 0:
        logger.warning("inspectMcErr: empty error array; failing diagnostic")
        return None, False

    nan_mask = np.isnan(arr)
    n_total = arr.size
    n_nan = nan_mask.sum()
    n_valid = n_total - n_nan

    valid_percent = 100.0 * n_valid / n_total
    nan_percent = 100.0 - valid_percent

    if n_nan > 0:
        logger.debug(
            f"inspectMcErr: {n_nan}/{n_total} NaN values "
            f"({nan_percent:.1f}% NaN, {valid_percent:.1f}% valid)"
        )

    # Too many NaNs → fail
    if valid_percent < min_percentage:
        logger.warning(
            "inspectMcErr: too many NaN MC realisations "
            f"({nan_percent:.1f}% > allowed {100.0 - min_percentage:.1f}%). "
            "Marking diagnostic as failed."
        )
        return None, False

    # Enough valid values → pass
    # valid_vals = arr[~nan_mask]
    # sigma = float(np.std(valid_vals, ddof=1))  # ddof=1 → sample std dev

    return arr[~nan_mask], True


_DIAG_RE = re.compile(
    r"""^
        \[
          (?P<elem>[A-Z][a-z]*)          # element symbol: N, Ni, Fe, ...
          (?P<spec>[IVXLCDM]+)           # roman spectrum: I, II, III, ...
        \]
        \s*
        (?P<num>\d+(?:\.\d+)?)           # numerator wavelength
        \s*/\s*
        (?P<den>\d+(?:\.\d+)?)           # denominator wavelength
        (?P<plus>\+)?                    # optional '+'
        $
    """,
    re.VERBOSE,
)


@dataclass(frozen=True)
class DiagSpec:
    label: str
    element: str  # e.g. "Ni"
    spectrum_roman: str  # e.g. "II"
    num_wave: float
    den_wave: float
    den_plus: bool


def parse_diag_label(label: str) -> DiagSpec:
    m = _DIAG_RE.fullmatch(label.strip())
    if not m:
        raise ValueError(f"Unsupported diagnostic label format: {label!r}")
    return DiagSpec(
        label=label.strip(),
        element=m.group("elem"),
        spectrum_roman=m.group("spec"),
        num_wave=float(m.group("num")),
        den_wave=float(m.group("den")),
        den_plus=bool(m.group("plus")),
    )


def ratio_and_err(
    num: float, den: float, num_err: float, den_err: float
) -> Tuple[float, float]:
    """
    R = num/den
    σ_R ≈ R * sqrt( (σ_num/num)^2 + (σ_den/den)^2 )
    """
    if den <= 0 or num <= 0:
        return float("nan"), float("nan")
    R = num / den
    rel2 = 0.0
    if num_err is not None and num > 0:
        rel2 += (num_err / num) ** 2
    if den_err is not None and den > 0:
        rel2 += (den_err / den) ** 2
    return R, abs(R) * math.sqrt(rel2)


def findIntensities(spec: DiagSpec, fitsd: list, intensity_list: list, logger=None):
    # cstr eg "N2_6583", should be matched in fitsd and then the "fractional"
    # name returned, e.g. N2_6583.7A, to get intensity and error
    ion = spec.element
    spectrum = roman.roman2int(spec.spectrum_roman)
    line = spec.num_wave
    for entry in fitsd:
        if (
            ion == entry["element"]
            and spectrum == roman.roman2int(entry["spectrum"])
            and line == float(entry["atomic"])
        ):
            for ilist in intensity_list:
                if ilist["element_pn"] == entry["pnstr"]:
                    return (ilist["intensity"], ilist["intensity_err"])
    raise RuntimeError(f"[ERROR] Failed finding intensity for {spec}")


def get_line_from_my_dict(
    spec: DiagSpec, which: str, fitsd, intensity_list, logger=None
) -> Tuple[float, float]:
    # which is "num" or "den"
    if which == "num":
        wave = spec.num_wave
        plus = False
    else:
        wave = spec.den_wave
        plus = spec.den_plus

    i, ierr = findIntensities(spec, fitsd, intensity_list, logger)

    # TODO: replace with your own matcher.
    # Return (intensity, intensity_err) for that wave, or sum if plus=True.
    raise NotImplementedError


def computeTeNePairs(
    density_diagnostics: list,
    temperature_diagnostics: list,
    intensities: list,
    min_percentage,
    logger,
):
    """
    For every (T-diagnostic, n-diagnostic) pair:
      - compute observed ratios + ratio errors
      - call diags.getCrossTemDen using value_tem/value_den
    Returns a list of results dicts.
    """
    diags = pn.Diagnostics()

    # Register all diagnostics with PyNeb
    all_labels = list(
        dict.fromkeys(temperature_diagnostics + density_diagnostics)
    )  # preserve order, unique
    for lab in all_labels:
        diags.addDiag(lab)  # predefined labels in pn.diags_dict

    # Pre-parse
    T_specs = [parse_diag_label(l) for l in temperature_diagnostics]
    n_specs = [parse_diag_label(l) for l in density_diagnostics]

    results = []

    # Compute ratios once per diagnostic
    ratio_cache: Dict[str, Tuple[float, float]] = {}  # label -> (R, Rerr)

    def get_ratio(spec: DiagSpec) -> Tuple[float, float]:
        if spec.label in ratio_cache:
            return ratio_cache[spec.label]
        In, en = get_line(spec, "num")
        Id, ed = get_line(spec, "den")
        R, Rerr = ratio_and_err(In, Id, en, ed)
        ratio_cache[spec.label] = (R, Rerr)
        return R, Rerr

    # Cross-combine temperature and density diagnostics
    for t in T_specs:
        Rt, Rt_err = get_ratio(t)
        for d in n_specs:
            Rn, Rn_err = get_ratio(d)

            Te, Ne = diags.getCrossTemDen(
                t.label,
                d.label,
                value_tem=Rt,
                value_den=Rn,
            )

            results.append(
                {
                    "tem_diag": t.label,
                    "den_diag": d.label,
                    "R_tem": Rt,
                    "R_tem_err": Rt_err,
                    "R_den": Rn,
                    "R_den_err": Rn_err,
                    "Te": Te,
                    "Ne": Ne,
                }
            )

    return results


def computeTeNePairs_obsolete(
    density_diagnostics: list,
    tempterature_diagnostics: list,
    pnObs,
    pnErrObs,
    min_percentage,
    logger,
):
    def filterPairs(densityd, temperatured, pobs):
        user = [(td, dd) for td in temperatured for dd in densityd]
        diags = pn.Diagnostics()
        # construct all possible diagnostics from the given observation set
        diags.addDiagsFromObs(pobs)
        validLines = diags.getDiagLabels()
        # filter user list based on observation set
        return [(d[0], d[1]) for d in user if d[0] in validLines and d[1] in validLines]

    if len(filterPairs(density_diagnostics, tempterature_diagnostics, pnObs)) == 0:
        print("----------------------------------->>>>")
    db.debug_obs_labels(pnObs)

    tene_slit_dict = []
    diags = pn.Diagnostics()
    for pair in filterPairs(density_diagnostics, tempterature_diagnostics, pnObs):
        st, sn = diags.getCrossTemDen(pair[0], pair[1], obs=pnObs)
        et, en = diags.getCrossTemDen(pair[0], pair[1], obs=pnErrObs)
        et, okT = inspectMcErr(et, min_percentage, logger)
        en, okN = inspectMcErr(en, min_percentage, logger)
        tene_slit_dict.append(
            {
                "tene_pair": (pair[0], pair[1]),
                "sT": st,
                "sN": sn,
                "eT": et,
                "eN": en,
            }
        )

    return tene_slit_dict, not (okT and okN)


def err2scalar(err_array, logger):
    """
    Define the fucntion to produce a single (i.e. scalar) error value,
    when we have an array of error values (i.e. from Monte-Carlo
    simulations.
    """
    if np.isnan(err_array).any():
        logger.warning(
            "WARNING  array contains nan! {:} for computing error in temperature/density diagnostics".format(
                err_array
            )
        )
    return np.std(err_array[~np.isnan(err_array)])


def computeTeNePairs_obsolete(
    density_diagnostics: list, tempterature_diagnostics: list, pnObs, pnErrObs, logger
):
    def filterPairs(densityd, temperatured, pobs):
        user = [(td, dd) for td in temperatured for dd in densityd]
        diags = pn.Diagnostics()
        # construct all possible diagnostics from the given observation set
        diags.addDiagsFromObs(pobs)
        validLines = diags.getDiagLabels()
        # filter user list based on observation set
        return [(d[0], d[1]) for d in user if d[0] in validLines and d[1] in validLines]

    includes_nan = False
    tene_slit_dict = []
    diags = pn.Diagnostics()
    for pair in filterPairs(density_diagnostics, tempterature_diagnostics, pnObs):
        st, sn = diags.getCrossTemDen(pair[0], pair[1], obs=pnObs)
        et, en = diags.getCrossTemDen(pair[0], pair[1], obs=pnErrObs)
        if np.isnan(np.array([et, en])).any():
            logger.warning(
                f"Nan value(s) encountered while computing Te/Ne diagnostics!"
            )
            includes_nan = True
        tene_slit_dict.append(
            {
                "tene_pair": (pair[0], pair[1]),
                "sT": st,
                "sN": sn,
                "eT": et,
                "eN": en,
            }
        )

    return tene_slit_dict, includes_nan


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
                            err2scalar(entry["eT"], logger),
                            entry["sN"],
                            err2scalar(entry["eN"], logger),
                        ),
                        file=fout,
                        end="",
                    )
                else:
                    print("{:53s} ".format(" "), file=fout, end="")
            print("", file=fout)
