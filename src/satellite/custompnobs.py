"""Alternative utilities to imitate pyneb's pn.Observation() setup with custom functions

2. Load file to create pn.Observations()s, i.e.
        logger.info("Setting-up PyNeb ...")
        # PyNeb stuff; PyNeb will read the 'test.dat' file (for the slit).
        sobs = pn.Observation()
        sobs.readData(
            "test.dat", fileFormat="lines_in_rows_err_cols", errIsRelative=False
        )
        sobs.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)
        sobs.extinction.law = ext_law
        sobs.correctData(normWave=4861.0)

3. And a similar one that will be used to estimate respective errors, using Monte Carlo simulations:
      eobs = pn.Observation()
            eobs.readData(
                "test.dat", fileFormat="lines_in_rows_err_cols", errIsRelative=False
            )
            eobs.addMonteCarloObs(N=monte_carlo_fake_obs)
            eobs.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)
            eobs.extinction.law = ext_law
            eobs.correctData(normWave=4861.0)

            RC = pn.RedCorr(E_BV=sobs.extinction.E_BV[0], R_V=pn_rv, law=ext_law)

            # Use Monte-Carlo simulations for E(B-V) and c(Hb) undertainties
            # 1. factor to convert E(B–V) to c(Hβ)
            RC_test = pn.RedCorr(E_BV=1.0, R_V=pn_rv, law=ext_law)
            f = RC_test.cHbeta
            # 2. Get uncertainty on E(B–V) from Monte Carlo results
            ebv_err = eobs.extinction.E_BV.std()
            # 3. Convert to c(Hβ) uncertainty
            chbeta_err = f * ebv_err
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import numpy as np
import pyneb as pn


_NUM_UNIT = re.compile(r"(\d+(?:\.\d+)?)(A|m)")


def label_to_angstrom(label: str) -> float:
    """
    Parse 'O3_5006.8A' or 'O3_88.3m' -> wavelength in Angstrom.
    I have never encountered these .[0-9]{1}m values but still, could potentially
    be there
    """
    if "_" not in label:
        raise ValueError(f"Bad line label (no '_'): {label!r}")
    wfrag = label.split("_", 1)[1]
    m = _NUM_UNIT.search(wfrag)
    if not m:
        raise ValueError(f"Cannot parse wavelength from label {label!r}")
    val = float(m.group(1))
    unit = m.group(2)
    return val * (1e4 if unit == "m" else 1.0)


@dataclass
class LineRow:
    label: str
    I_obs: float
    I_err: float
    wave_A: float  ## Note, this is a float for us, not an int as in Pyneb (!!)


def read_simple_lines_file(path: str) -> List[LineRow]:
    """
    Parses a 3-column file (aka a 'test.dat' file):
        LINE  value  err
        O3_5006.8A  49.0  2.0
    """
    rows: List[LineRow] = []
    with open(path, "r", encoding="utf-8") as f:
        header = f.readline()  # discard
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split()
            if len(parts) < 3:
                continue
            label = parts[0]
            I = float(parts[1])
            e = float(parts[2])
            rows.append(
                LineRow(label=label, I_obs=I, I_err=e, wave_A=label_to_angstrom(label))
            )
    return rows


def _merge_duplicates(rows: List[LineRow]) -> List[LineRow]:
    """
    Optional: combine duplicate labels by summing intensities
    and adding errors in quadrature.

    This should not happen but config does not explicitelly dissalow it, and we
    have no checks to guard against it.
    """
    acc: Dict[str, Tuple[float, float, float]] = {}
    for r in rows:
        if r.label not in acc:
            acc[r.label] = (r.I_obs, r.I_err**2, r.wave_A)
        else:
            I, e2, w = acc[r.label]
            acc[r.label] = (I + r.I_obs, e2 + r.I_err**2, w)  # keep first wavelength
    out = []
    for lab, (I, e2, w) in acc.items():
        out.append(LineRow(lab, I, float(np.sqrt(e2)), w))
    return out


def deredden_and_normalize(
    rows: List[LineRow],
    *,
    ext_law: str,
    R_V: float,
    ha_label: str = "H1r_6563A",
    hb_label: str = "H1r_4861A",
    r_theo: float = 2.85,
    norm_wave_A: float = 4861.0,
) -> Tuple[Dict[str, float], pn.RedCorr]:
    """
    Mimics:
      obs.def_EBV(label1=ha, label2=hb, r_theo=2.85)
      obs.extinction.law = ext_law
      obs.correctData(normWave=4861.0)

    but works with ANY labels (fractional Angstrom, Fe2, etc.).
    """
    by_label = {r.label: r for r in rows}
    if ha_label not in by_label or hb_label not in by_label:
        raise KeyError(f"Missing lines: need {ha_label} and {hb_label}")

    FHa = by_label[ha_label].I_obs
    FHb = by_label[hb_label].I_obs

    rc = pn.RedCorr(R_V=R_V, law=ext_law)
    # Exactly like the doc example: obs_over_theo = (FHa/FHb) / r_theo :contentReference[oaicite:4]{index=4}
    rc.setCorr(
        (FHa / FHb) / r_theo, label_to_angstrom(ha_label), label_to_angstrom(hb_label)
    )

    # Relative correction to Hβ: getCorr(wave, rel_wave=4861) :contentReference[oaicite:5]{index=5}
    corr: Dict[str, float] = {}
    for r in rows:
        corr[r.label] = float(rc.getCorr(r.wave_A, norm_wave_A))
    return corr, rc


def monte_carlo_errors(
    rows: List[LineRow],
    *,
    ext_law: str,
    R_V: float,
    N: int,  # number of Monte Carlo draws
    seed: Optional[int] = None,
    ha_label: str = "H1r_6563A",
    hb_label: str = "H1r_4861A",
    r_theo: float = 2.85,
    norm_wave_A: float = 4861.0,
    clip_negative: bool = True,
) -> Tuple[Dict[str, float], float]:
    """
    Mimics:
      eobs.addMonteCarloObs(N)
      eobs.def_EBV(...)
      eobs.correctData(...)
      ebv_err = eobs.extinction.E_BV.std()

    Returns:
      (rel_sigma_by_label, ebv_std)
      where rel_sigma = std(I_corr)/mean(I_corr) for each line.
    """
    rng = np.random.default_rng(seed)

    labels = [r.label for r in rows]
    waves = np.array([r.wave_A for r in rows], dtype=float)
    I0 = np.array([r.I_obs for r in rows], dtype=float)
    s0 = np.array([r.I_err for r in rows], dtype=float)

    # 1) Draw fake observations (this matches addMonteCarloObs)
    # draws[k, i] is the simulated intensity for line i in Monte Carlo realization k
    # drawn from a normal distribution N(I0[i], s0[i]**2)
    draws = rng.normal(loc=I0[None, :], scale=s0[None, :], size=(N, len(rows)))
    if clip_negative:
        draws = np.clip(draws, 0.0, None)

    # figure out which columns in draws correspond to Hα and Hβ, so the Monte
    # Carlo loop can grab the right simulated values for the Balmer decrement.
    idx = {lab: i for i, lab in enumerate(labels)}
    iHa = idx[ha_label]
    iHb = idx[hb_label]

    rc = pn.RedCorr(R_V=R_V, law=ext_law)

    Icorr = np.empty_like(draws)
    ebv = np.empty(N, dtype=float)

    for k in range(N):
        FHa = draws[k, iHa]
        FHb = draws[k, iHb]
        # avoid division blow-ups
        if FHb <= 0:
            ebv[k] = np.nan
            Icorr[k, :] = np.nan
            continue

        rc.setCorr(
            (FHa / FHb) / r_theo,
            label_to_angstrom(ha_label),
            label_to_angstrom(hb_label),
        )
        ebv[k] = float(np.atleast_1d(rc.E_BV)[0])

        # vectorized per-draw correction relative to Hβ :contentReference[oaicite:6]{index=6}
        # imitates: correctData(normWave=4861)
        # kinda: I(corr)[λ(i)] = I(obs)[λ(i)] * Corr(λ(i), λ(norm))
        cf = np.array(rc.getCorr(waves, norm_wave_A), dtype=float)
        Icorr[k, :] = draws[k, :] * cf
        Icorr[k, :] *= 100.0 / Icorr[k, iHb]

    # stats (ignore NaNs if any)
    mean = np.nanmean(Icorr, axis=0)
    std = np.nanstd(Icorr, axis=0)

    rel_sigma = {
        lab: float(std[i] / mean[i]) if mean[i] > 0 else np.nan
        for i, lab in enumerate(labels)
    }
    ebv_std = float(np.nanstd(ebv))
    return rel_sigma, ebv_std
