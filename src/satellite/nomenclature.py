from __future__ import annotations
import satellite.roman as sr
import pyneb as pn
import numpy as np
import re
from typing import Optional, Tuple, Union


def closestPyNebElement(element: str, atomic_number: int, logger=None) -> str:
    min_diff = 1000000
    pn_wl = None
    atomic_number = float(atomic_number)
    for wl in pn.LINE_LABEL_LIST[element]:
        if re.fullmatch("[0-9]*A", wl):
            diff = abs(float(wl[0:-1]) - atomic_number)
            if diff < min_diff:
                pn_wl = wl
                min_diff = diff
    if logger:
        logger.warning(
            "PyNeb is missing wavelength {:} for element {:}; using {:} instead".format(
                atomic_number, element, pn_wl
            )
        )

    return pn_wl, round(pn_wl[:-1]) if pn_wl.endswith("A") else round(pn_wl)


def unlisted(atom, spectrum_int, wavelength, logger=None):
    entries = pn.atomicData.getDataFile()[f"{atom}{spectrum_int}"]
    if pn.atomicData.getDataFile(f"{atom}{spectrum_int}", "atom") is None:
        if "atom" in entries:
            pn.atomicData.setDataFile(entries["atom"])
            if logger:
                logger.info(f"Loaded {entries['atom']}")
    if pn.atomicData.getDataFile(f"{atom}{spectrum_int}", "coll") is None:
        if "coll" in entries:
            pn.atomicData.setDataFile(entries["coll"])
            if logger:
                logger.info(f"Loaded {entries['coll']}")
    if pn.atomicData.getDataFile(f"{atom}{spectrum_int}", "rec") is None:
        if "rec" in entries:
            pn.atomicData.setDataFile(entries["rec"])
            if logger:
                logger.info(f"Loaded {entries['rec']}")
    elmnt = pn.Atom(atom, spectrum_int)
    waves = np.array(elmnt.lineList)
    target = float(wavelength)
    i = np.argmin(np.abs(waves - target))
    closest_wave = waves[i]
    if logger:
        logger.info(
            f"Closest line to {atom}{spectrum_int}_{wavelength} is line {closest_wave}"
        )
    return f"{atom}{spectrum_int}_{round(closest_wave)}A", round(closest_wave)


def objectIntensityPyNebCode(atom: str, spectrum: str, atomic_number: int, logger=None):
    """ """
    ## Spectrum as integer value
    try:
        spectrum = int(spectrum)
    except:
        spectrum = sr.roman2int(spectrum)

    ## The atom
    if atom in ["H", "He"]:
        pnatom = "{:}{:}r".format(atom, spectrum)
    else:
        pnatom = "{:}{:}".format(atom, spectrum)

    ## Transition Lines (if any)
    try:
        wls = pn.LINE_LABEL_LIST[pnatom]
        pnatomic = f"{atomic_number}A"
        if pnatomic not in wls:
            pnatomic, pnline = closestPyNebElement(pnatom, atomic_number)
        else:
            pnline = round(atomic_number)
        pnstr = "_".join([pnatom, pnatomic])
    except:
        wls = None
        if logger:
            logger.warning(
                f"Failed matching object {atom}/{spectrum} (aka {pnatom}) to PyNeb (see LINE_LABEL_LIST)"
            )
            logger.warning(f"Will try to match a specific data file ...")

    if wls is not None:
        return pnstr, pnline
    else:
        return unlisted(atom, spectrum, atomic_number, logger)


# Generic Wavelength for easy manipulation (parsing, etc)
WaveLike = Union[int, float, str]


def _wave_to_angstrom(w: WaveLike) -> float:
    """Parse wavelength to Angstrom.
    Accepts: 5007, 5007.0, '5007', '5007A', '10.5m', '500.7nm', 'Å' (treated as A).
    Returns the Wavelength as a floating point number
    """
    if isinstance(w, (int, float)):
        return float(w)

    s = str(w).strip().replace("Å", "A")
    m = re.fullmatch(r"(\d+(?:\.\d+)?)(?:\s*)(A|m|nm)", s)
    if m:
        val = float(m.group(1))
        unit = m.group(2)
        if unit == "A":
            return val
        if unit == "nm":
            return val * 10.0
        if unit == "m":  # microns in PyNeb labels
            return val * 1e4

    # Fallback: extract first float-looking token, assume Angstrom
    m2 = re.search(r"(\d+(?:\.\d+)?)", s)
    if not m2:
        raise ValueError(f"Could not parse wavelength from {w!r}")
    return float(m2.group(1))


def _parse_line_label_to_angstrom(wl_label: str) -> Optional[float]:
    """Parse a PyNeb line fragment like '5007A', '88.3m', '7319A+' into Angstrom."""
    s = wl_label.strip().replace("Å", "A")

    # Find number + unit somewhere inside (tolerate suffixes like '+' or '.')
    m = re.search(r"(\d+(?:\.\d+)?)(A|m)", s)
    if m:
        val = float(m.group(1))
        unit = m.group(2)
        return val * (1e4 if unit == "m" else 1.0)

    # If no explicit unit but has a number, assume Angstrom
    m2 = re.search(r"(\d+(?:\.\d+)?)", s)
    if m2:
        return float(m2.group(1))

    return None


def _pick_closest_label(
    wl_labels: list[str],
    target_ang: float,
) -> Tuple[str, float]:
    """Return (best_wl_label_fragment, best_wavelength_angstrom)."""
    best_label = None
    best_ang = None
    best_diff = float("inf")

    for lab in wl_labels:
        ang = _parse_line_label_to_angstrom(lab)
        if ang is None:
            continue
        diff = abs(ang - target_ang)
        if diff < best_diff:
            best_diff = diff
            best_label = lab  # wavelenth string
            best_ang = ang  # wavelength float

    if best_label is None or best_ang is None:
        raise ValueError("No parseable wavelengths found in candidate label list")

    return best_label, best_ang


def _try_load_any_data_files(ion: str, logger=None) -> bool:
    """Try to make PyNeb aware of some data files for `ion` by scanning available files."""
    try:
        files = pn.atomicData.getAllAvailableFiles(ion)
    except Exception:
        return False

    if not files:
        return False

    # Pick the first matching file per type (atom/coll/rec)
    def pick(pattern: str) -> Optional[str]:
        for f in files:
            if re.search(pattern, f, flags=re.IGNORECASE):
                return f
        return None

    atom_f = pick(r"_atom_|\batom\b")  # atomic
    coll_f = pick(r"_coll_|\bcoll\b")  # collision
    rec_f = pick(r"_rec_|\brec\b")  # recombination

    loaded_any = False
    for f in (atom_f, coll_f, rec_f):
        if f:
            try:
                pn.atomicData.setDataFile(f)
                loaded_any = True
                if logger:
                    logger.info(f"Loaded atomic data file: {f}")
            except Exception:
                pass

    return loaded_any


def best_pyneb_line(
    element: str,
    spectrum: Union[int, str],
    wavelength: WaveLike,
    *,
    logger=None,
    use_recomb_suffix_for_hhe: bool = True,
) -> Tuple[str, float]:
    """
    Given (element, spectrum, wavelength), return:
      (full_pyneb_line_label, chosen_wavelength_in_angstrom)

    Examples:
      best_pyneb_line("O", "III", 5006.8) -> ("O3_5007A", 5007.0)
      best_pyneb_line("H", 1, 4861)       -> ("H1r_4861A", 4861.0)
      best_pyneb_line("S", 4, "10.5m")    -> ("S4_10.5m", 105000.0)
    """
    # Normalize spectrum, aka make it an integer
    spec_int = (
        int(spectrum)
        if isinstance(spectrum, (int, np.integer))
        else (
            int(spectrum)
            if str(spectrum).strip().isdigit()
            else sr.roman2int(str(spectrum))
        )
    )

    # Get the wavelength as a floating point number (the user-declared one)
    target_ang = _wave_to_angstrom(wavelength)

    # Ion label for LINE_LABEL_LIST (H/He recombination labels are usually *r)
    if use_recomb_suffix_for_hhe and element in {"H", "He"}:
        ion_label = f"{element}{spec_int}r"
    else:
        ion_label = f"{element}{spec_int}"

    # Ion key for atomicData file lookup typically does NOT use the 'r' suffix
    ion_data_key = f"{element}{spec_int}"

    # 1) Try LINE_LABEL_LIST first
    if ion_label in pn.LINE_LABEL_LIST:
        wl_labels = list(pn.LINE_LABEL_LIST[ion_label])
        best_frag, best_ang = _pick_closest_label(wl_labels, target_ang)
        full = f"{ion_label}_{best_frag}"

        # warn if not exact fragment match
        # (exactness is fuzzy; treat within 1e-3 A as exact)
        if abs(best_ang - target_ang) > 1e-3 and logger:
            logger.warning(
                f"PyNeb has no exact {ion_label}_{wavelength}; using closest {full}"
            )
        return full, best_ang

    if logger:
        logger.warning(
            f"{ion_label} not found in pn.LINE_LABEL_LIST; trying atomic data files..."
        )

    # 2) Try to load any data files for the ion
    _try_load_any_data_files(ion_data_key, logger=logger)

    # 3) Instantiate Atom and use its computed transitions (collisionally excited lines)
    try:
        atom_obj = pn.Atom(element, spec_int)
        waves = np.asarray(atom_obj.lineList, dtype=float)
        if waves.size == 0:
            raise ValueError("Atom.lineList is empty")
        i = int(np.argmin(np.abs(waves - target_ang)))
        best_ang = float(waves[i])

        # Build a PyNeb-looking label fragment (Angstrom by default)
        best_frag = f"{int(round(best_ang))}A"
        full = f"{ion_data_key}_{best_frag}"

        if logger:
            logger.info(
                f"Closest Atom() line to {ion_data_key}_{wavelength} is {full} (≈{best_ang:.4f} A)"
            )
        return full, best_ang

    except Exception as e:
        raise ValueError(
            f"Could not match {element}{spec_int} at {wavelength!r}: "
            f"not in LINE_LABEL_LIST and Atom() fallback failed."
        ) from e
