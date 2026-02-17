from __future__ import annotations

import re, os
from typing import Optional, Tuple, Union
from pathlib import Path

import numpy as np
import pyneb as pn

import satellite.roman as sr

WaveLike = Union[int, float, str]


def wave_to_angstrom(w: WaveLike) -> float:
    if isinstance(w, (int, float)):
        return float(w)
    s = str(w).strip().replace("Å", "A")
    m = re.fullmatch(r"(\d+(?:\.\d+)?)(?:\s*)(A|nm|m)", s)
    if m:
        val = float(m.group(1))
        unit = m.group(2)
        if unit == "A":
            return val
        if unit == "nm":
            return val * 10.0
        if unit == "m":  # microns
            return val * 1e4
    # fallback: first number, assume Å
    m2 = re.search(r"(\d+(?:\.\d+)?)", s)
    if not m2:
        raise ValueError(f"Could not parse wavelength from {w!r}")
    return float(m2.group(1))


def label_to_angstrom(lbl: str) -> Optional[float]:
    s = lbl.strip().replace("Å", "A")
    m = re.search(r"(\d+(?:\.\d+)?)(A|m)", s)
    if not m:
        return None
    val = float(m.group(1))
    return val * (1e4 if m.group(2) == "m" else 1.0)


def closest_from_labels(ion_label: str, target_ang: float) -> Tuple[str, float]:
    labels = pn.LINE_LABEL_LIST[ion_label]
    best_lbl, best_ang, best_diff = None, None, float("inf")
    for lab in labels:
        ang = label_to_angstrom(lab)
        if ang is None:
            continue
        diff = abs(ang - target_ang)
        if diff < best_diff:
            best_lbl, best_ang, best_diff = lab, ang, diff
    if best_lbl is None or best_ang is None:
        raise ValueError(f"No parseable labels in LINE_LABEL_LIST[{ion_label}]")
    return best_lbl, best_ang


def format_ang_label(wang: float) -> str:
    # mimic the “one decimal when needed” feel you see in printTransition()
    if abs(wang - round(wang)) < 1e-6:
        return f"{int(round(wang))}A"
    s = f"{wang:.1f}".rstrip("0").rstrip(".")
    return f"{s}A"


# ---------- atomic data loading helper (strips '* ') ----------

_STATUS_PREFIX_RE = re.compile(r"^[*XDD]\s*", flags=re.IGNORECASE)


def pick_default(candidates):
    # prefer starred/default then name
    return sorted(candidates, key=lambda t: (not t[0], t[1]))[0][1]


def pick_chianti(candidates):
    chi = [c for c in candidates if c[1].lower().endswith(".chianti")]
    if not chi:
        return None
    return sorted(chi, key=lambda t: (not t[0], t[1]))[0][1]


def try_load_any_data_files(
    ion: str, logger=None, *, allow_chianti: bool = False
) -> bool:
    """
    Load one file per bucket (atom/coll/rec).
    If allow_chianti=False: never choose *.chianti files.
    If allow_chianti=True: prefer *.chianti files when present.
    """
    try:
        raw = pn.atomicData.getAllAvailableFiles(ion) or []
    except Exception:
        return False

    def clean(entry: str) -> tuple[bool, str]:
        s = entry.strip()
        is_default = s.startswith("*")
        s = _STATUS_PREFIX_RE.sub("", s).strip()
        return is_default, s

    def kind(fname: str):
        f = fname.lower()
        if "_atom_" in f or f.endswith("_atom.chianti") or "_atom." in f:
            return "atom"
        if "_coll_" in f or f.endswith("_coll.chianti") or "_coll." in f:
            return "coll"
        if "_rec_" in f or f.endswith(".func"):
            return "rec"
        return None

    buckets = {"atom": [], "coll": [], "rec": []}
    for e in raw:
        is_def, fn = clean(e)
        k = kind(fn)
        if k and fn:
            buckets[k].append((is_def, fn))

    loaded_any = False
    for k in ("atom", "coll", "rec"):
        if not buckets[k]:
            continue

        fn = None
        if allow_chianti:
            fn = pick_chianti(buckets[k])

        if fn is None:
            # default selection, but avoid *.chianti when not allowed
            pool = (
                buckets[k]
                if allow_chianti
                else [c for c in buckets[k] if not c[1].lower().endswith(".chianti")]
            )
            if not pool:
                continue
            fn = pick_default(pool)

        try:
            pn.atomicData.setDataFile(fn)
            loaded_any = True
            if logger:
                logger.info(f"Loaded {k} data for {ion}: {fn}")
        except Exception as ex:
            if logger:
                logger.warning(f"Failed loading {k} for {ion} from {fn}: {ex}")

    return loaded_any


# ---------- chianti stuff ----------
def has_valid_xuvtop() -> bool:
    """
    True iff XUVTOP is set and looks like a CHIANTI database root.
    (CHIANTI root typically contains 'masterlist' and 'VERSION'.)
    """
    x = os.environ.get("XUVTOP")
    if not x:
        return False
    p = Path(x).expanduser()
    return (
        p.exists()
        and p.is_dir()
        and (p / "masterlist").exists()
        and (p / "VERSION").exists()
    )


def try_create_atom_with_chianti(element: str, spec: int, logger=None):
    """
    Try to construct a CHIANTI-backed Atom.
    Returns ChiantiAtom instance or None.
    """
    try:
        import pyneb as pn

        chianti_atom = pn.utils.pn_chianti.ChiantiAtom(element, spec)

        if logger:
            logger.info(f"Loaded {element}{spec} from CHIANTI database")

        return chianti_atom

    except Exception as ex:
        if logger:
            logger.warning(f"Failed to load {element}{spec} from CHIANTI: {ex}")
        return None


def chianti_closest_wavelength_A(element: str, spec: int, target_ang: float) -> float:
    from ChiantiPy.core.ion import ion as ChiantiIon

    """
    Return closest transition wavelength (Å) from CHIANTI wgfa table.
    Requires XUVTOP to be set (CHIANTI DB root).
    """
    ion_name = f"{element.lower()}_{spec}"

    # T, Ne just to satisfy constructor; irrelevant for wgfa wavelengths
    ch = ChiantiIon(ion_name, temperature=1e4, eDensity=1e2)

    # Ensure radiative transitions table is loaded
    if not hasattr(ch, "Wgfa") or ch.Wgfa is None:
        ch.wgfaRead()

    waves = np.asarray(ch.Wgfa.get("wvl", []), dtype=float)
    waves = waves[np.isfinite(waves)]
    waves = waves[waves > 0.0]
    if waves.size == 0:
        raise RuntimeError(f"CHIANTI: no wavelengths found for {ion}")

    i = int(np.argmin(np.abs(waves - target_ang)))
    return float(waves[i])


def pyneb_atom_has_lines(atom_obj) -> bool:
    """
    True if Atom looks usable (has a non-empty numeric lineList).
    This catches the 'data not available but Atom still created' situation.
    """
    try:
        ll = getattr(atom_obj, "lineList", None)
        if ll is None:
            return False
        arr = np.asarray(ll, dtype=float)
        return arr.ndim == 1 and arr.size > 0 and np.isfinite(arr).all()
    except Exception:
        return False


# ---------- main function ----------


def ensure_pyneb_atom(element: str, spec: int, ion_key: str, logger=None):
    # 1) Try plain PyNeb Atom
    try:
        atom = pn.Atom(element, spec)
        if pyneb_atom_has_lines(atom):
            return atom
    except Exception:
        pass

    # 2) Try loading NON-chianti files first
    try_load_any_data_files(ion_key, logger=logger, allow_chianti=False)
    try:
        atom = pn.Atom(element, spec)
        if pyneb_atom_has_lines(atom):
            return atom
    except Exception:
        pass

    # 3) Only now try CHIANTI (if configured)
    if has_valid_xuvtop():
        try_load_any_data_files(ion_key, logger=logger, allow_chianti=True)
        try:
            atom = pn.Atom(element, spec)
            if pyneb_atom_has_lines(atom):
                return atom
        except Exception:
            pass

    return None


def best_pyneb_line(
    element: str,
    spectrum: Union[int, str],
    wavelength: WaveLike,
    logger=None,
) -> Tuple[str, float]:
    """
    Returns (label, closest_wavelength_A).

    Key behavior change vs your current version:
    - For non-(H,He), we *prefer* the Atom transition table (like printTransition),
      and only fall back to LINE_LABEL_LIST if Atom can't be constructed.
    """
    # spectrum -> int
    if isinstance(spectrum, (int, np.integer)):
        spec = int(spectrum)
    else:
        s = str(spectrum).strip()
        spec = int(s) if s.isdigit() else sr.roman2int(s)

    target_ang = wave_to_angstrom(wavelength)

    ion_key = f"{element}{spec}"  # for pn.Atom / data files
    ion_label = (
        f"{element}{spec}r" if element in {"H", "He"} else ion_key
    )  # for LINE_LABEL_LIST

    # H/He: use LINE_LABEL_LIST (recombination labels)
    if element in {"H", "He"} and ion_label in pn.LINE_LABEL_LIST:
        frag, best_ang = closest_from_labels(ion_label, target_ang)
        return f"{ion_label}_{frag}", best_ang

    # Other ions: try Atom first (this matches printTransition behavior)
    atom_obj = None

    # If not using chianti, this is the correct branch
    #
    # try:
    #    atom_obj = pn.Atom(element, spec)
    # except Exception:
    #    try_load_any_data_files(ion_key, logger=logger)
    #    atom_obj = pn.Atom(element, spec)
    #
    atom_obj = ensure_pyneb_atom(element, spec, ion_key, logger)
    if atom_obj is None:
        raise RuntimeError(
            f"Could not construct Atom for {ion_key}. "
            f"Try setting pn.atomicData.setDataFileDict('PYNEB_23_01') "
            f"and ensure CHIANTI *.chianti files are discoverable."
        )

    # Prefer getTransition if available; otherwise use lineList directly
    closest_ang = None
    if hasattr(atom_obj, "getTransition"):
        try:
            tr = atom_obj.getTransition(target_ang)  # some versions accept positional
        except TypeError:
            tr = atom_obj.getTransition(wave=target_ang)

        # robust extraction: if getTransition returns a dict/tuple, find the float wavelength inside
        if isinstance(tr, dict):
            for k in ("wave", "wavelength", "lambda"):
                if k in tr:
                    closest_ang = float(tr[k])
                    break
        elif isinstance(tr, (tuple, list)):
            # find first float-like item
            for x in tr:
                if isinstance(x, (int, float, np.floating)):
                    closest_ang = float(x)
                    break

    if closest_ang is None:
        waves = np.asarray(atom_obj.lineList, dtype=float)
        i = int(np.argmin(np.abs(waves - target_ang)))
        closest_ang = float(waves[i])

    frag = format_ang_label(closest_ang)

    logger.info(
        f"Best fit for line {element}{spec}_{target_ang} in pyneb is {ion_key}_{frag}"
    )
    return f"{ion_key}_{frag}", closest_ang
