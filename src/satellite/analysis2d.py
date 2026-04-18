import copy
import datetime
import logging
import multiprocessing as mp
import os
import pickle
import sys
import tempfile
import warnings

import numpy as np
import pyneb as pn

import satellite.abundance as sb
import satellite.cfgio as cfgio
import satellite.custompnobs as sc
import satellite.fitsutils as fs
import satellite.icf as sf
import satellite.intensity as si
import satellite.ionic_abundancies as sa
import satellite.ratios as so
import satellite.satlogger as sl
import satellite.specific_slit as ss
import satellite.tene as st


def _cfg_bool(value, default=False):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _ordered_dict(dct: dict):
    return {key: dct[key] for key in sorted(dct.keys())}


def _relabel_output_headers_to_pixels(*filenames: str) -> None:
    for fn in filenames:
        if not fn or not os.path.exists(fn):
            continue
        with open(fn, "r") as fin:
            content = fin.read()
        content = content.replace("Slit Nr.", "Pixel Nr.")
        content = content.replace("\nSlit ", "\nPixel ")
        if content.startswith("Slit "):
            content = "Pixel " + content[len("Slit "):]
        with open(fn, "w") as fout:
            fout.write(content)


def _validate_shapes(fitsd: list):
    expected_shape = None
    for entry in fitsd:
        for key in ("fns", "fne"):
            arr = fs.loadFitsImageData(entry[key])
            if arr.ndim != 2:
                raise RuntimeError(
                    f"FITS image {entry[key]} is not 2-D (shape={arr.shape})"
                )
            if expected_shape is None:
                expected_shape = arr.shape
            elif arr.shape != expected_shape:
                raise RuntimeError(
                    f"Shape mismatch: {entry[key]} has shape {arr.shape}, expected {expected_shape}"
                )
    return expected_shape


def _get_pixel_ranges(shape: tuple, analysis_cfg: dict):
    nrows, ncols = shape
    row_start = int(analysis_cfg.get("row_start", 0))
    col_start = int(analysis_cfg.get("col_start", 0))
    row_step = int(analysis_cfg.get("row_step", 1))
    col_step = int(analysis_cfg.get("col_step", 1))
    row_end = analysis_cfg.get("row_end", None)
    col_end = analysis_cfg.get("col_end", None)
    row_end = nrows if row_end is None else min(int(row_end), nrows)
    col_end = ncols if col_end is None else min(int(col_end), ncols)

    if row_start < 0 or col_start < 0:
        raise RuntimeError("row_start and col_start must be >= 0")
    if row_step <= 0 or col_step <= 0:
        raise RuntimeError("row_step and col_step must be > 0")
    if row_end <= row_start or col_end <= col_start:
        raise RuntimeError("Requested pixel window is empty")

    return (
        range(row_start, row_end, row_step),
        range(col_start, col_end, col_step),
    )


def _load_all_fits_arrays(fitsd: list, logger):
    loaded = []
    for entry in fitsd:
        signal = np.array(fs.loadFitsImageData(entry["fns"]), dtype=float, copy=True)
        error = np.array(fs.loadFitsImageData(entry["fne"]), dtype=float, copy=True)
        np.nan_to_num(signal, copy=False)
        np.nan_to_num(error, copy=False)
        obj = copy.deepcopy(entry)
        obj["signal_data"] = signal
        obj["error_data"] = error
        loaded.append(obj)
        logger.info(
            f"Loaded FITS pair for {entry['element']}{entry['spectrum']}_{entry['atomic']}"
        )
    return loaded


def _clone_pixel_payload(fitsd_loaded: list, row: int, col: int) -> list:
    cpd = []
    for entry in fitsd_loaded:
        cloned = {
            key: value
            for key, value in entry.items()
            if key not in ("signal_data", "error_data")
        }
        cloned["sslit_sum"] = float(entry["signal_data"][row, col])
        cloned["eslit_sum"] = float(entry["error_data"][row, col])
        cpd.append(cloned)
    return cpd


def _run_pixel_pipeline(
    pixel_id: int,
    row: int,
    col: int,
    fitsd_loaded: list,
    ratios: list,
    density_diagnostics: list,
    temperature_diagnostics: list,
    ext_law: str,
    pn_atomic_data: str,
    monte_carlo_fake_obs: int,
    max_nan_in_diagnostics_allowed_percentage: int,
    pn_rv: float,
    energy_parameter: float,
    temp_dir: str,
    logger,
):
    cpd = _clone_pixel_payload(fitsd_loaded, row, col)

    tmp_line_file = os.path.join(temp_dir, f"pixel_{pixel_id:08d}.dat")
    si.makeIntensitiesDataFile(
        cpd,
        ss.reference_element,
        ["sslit_sum", "eslit_sum"],
        tmp_line_file,
        100e0,
        logger,
    )

    try:
        rows = sc.read_simple_lines_file(tmp_line_file)
        rows = sc._merge_duplicates(rows)
        corr_factor, rc_central = sc.deredden_and_normalize(
            rows, ext_law=ext_law, R_V=pn_rv, r_theo=2.85, norm_wave_A=4861.0
        )
        I_corr = {r.label: r.I_obs * corr_factor[r.label] for r in rows}

        nan_diagnostics = True
        times_nan_encountered = 0
        max_nan_tries = (
            max_nan_in_diagnostics_allowed_percentage * monte_carlo_fake_obs // 10
        )
        max_nan_tries = max(1, max_nan_tries)

        while nan_diagnostics and times_nan_encountered < max_nan_tries:
            rel_sigma, ebv_err = sc.monte_carlo_errors(
                rows, ext_law=ext_law, R_V=pn_rv, N=monte_carlo_fake_obs, seed=0
            )

            rc_test = pn.RedCorr(E_BV=1.0, R_V=pn_rv, law=ext_law)
            fac = rc_test.cHbeta
            chbeta = float(np.atleast_1d(rc_central.cHbeta)[0])
            chbeta_err = float(fac * ebv_err)

            global_intensities = {
                "intensities": [
                    {
                        "element_pn": r.label,
                        "wavelength_A": r.wave_A,
                        "intensity": I_corr[r.label],
                        "intensity_err": rel_sigma[r.label],
                    }
                    for r in rows
                ],
                "E_BV": float(np.atleast_1d(rc_central.E_BV)[0]),
                "E_BVError": ebv_err,
                "cHbeta": chbeta,
                "cHbetaError": chbeta_err,
                "fac": fac,
                "FHb": ss.findEntry(
                    ss.reference_element["element"],
                    ss.reference_element["spectrum"],
                    ss.reference_element["atomic"],
                    cpd,
                )["sslit_sum"]
                * energy_parameter,
                "FHb_error": ss.findEntry(
                    ss.reference_element["element"],
                    ss.reference_element["spectrum"],
                    ss.reference_element["atomic"],
                    cpd,
                )["eslit_sum"]
                * energy_parameter,
            }

            global_ratios = {}
            for ratio in ratios:
                global_ratios[ratio] = so.computeRatio(
                    ratio, cpd, global_intensities["intensities"], logger
                )

            global_tene, diagnostics_ok = st.computeTeNePairs(
                density_diagnostics,
                temperature_diagnostics,
                global_intensities,
                ss.WAVELENGTH_TOLERANCE_FOR_TENE,
                logger,
                mc_N=monte_carlo_fake_obs,
                seed=int(datetime.datetime.now().strftime("%Y%m%d%H%M%S")),
            )
            nan_diagnostics = not diagnostics_ok
            times_nan_encountered += 1

        if times_nan_encountered >= max_nan_tries and nan_diagnostics:
            raise RuntimeError(
                f"Failed computing non-nan diagnostics after {times_nan_encountered} tries"
            )

        global_ionic_abundancies = sa.computeIonicAbundancies(
            cpd,
            global_tene,
            global_intensities,
            1 - max_nan_in_diagnostics_allowed_percentage,
            logger,
            tol_A=ss.WAVELENGTH_TOLERANCE_FOR_TENE,
            mc_N=monte_carlo_fake_obs,
            seed=int(datetime.datetime.now().strftime("%Y%m%d%H%M%S")),
        )

        elemspec_abundancies = sb.computeAbundancies(cpd, global_ionic_abundancies, logger)
        global_icfs = sf.computeIcfsWithErrors(elemspec_abundancies, logger)
        global_element_abundancies = sf.ionicAbundance2elementAbundance(
            elemspec_abundancies, logger
        )

        return (
            global_intensities,
            global_ratios,
            global_tene,
            global_ionic_abundancies,
            global_icfs,
            global_element_abundancies,
        )
    finally:
        try:
            os.remove(tmp_line_file)
        except OSError:
            pass


G_ANALYSIS_CFG = None
G_FITS_LOADED = None
G_RATIOS = None
G_DENSITY_DIAGNOSTICS = None
G_TEMPERATURE_DIAGNOSTICS = None
G_REF_IDX = None
G_EXT_LAW = None
G_PN_ATOMIC_DATA = None
G_MONTE_CARLO = None
G_MAX_NAN_PCT = None
G_PN_RV = None
G_ENERGY_PARAMETER = None
G_WORKER_LOG_LEVEL = None
G_TEMP_DIR = None


def _split_rows_into_chunks(row_range: range, jobs: int) -> list:
    rows = list(row_range)
    if not rows:
        return []
    jobs = max(1, min(jobs, len(rows)))
    base = len(rows) // jobs
    rem = len(rows) % jobs
    chunks = []
    start = 0
    for idx in range(jobs):
        size = base + (1 if idx < rem else 0)
        chunks.append(rows[start : start + size])
        start += size
    return [chunk for chunk in chunks if chunk]


def _build_chunk_tasks(row_range: range, col_range: range, jobs: int) -> list:
    row_chunks = _split_rows_into_chunks(row_range, jobs)
    row_to_chunk = {}
    for chunk_idx, rows in enumerate(row_chunks):
        for row in rows:
            row_to_chunk[row] = chunk_idx

    chunk_tasks = [[] for _ in row_chunks]
    pixel_id = 0
    for row in row_range:
        chunk_idx = row_to_chunk[row]
        for col in col_range:
            chunk_tasks[chunk_idx].append((pixel_id, row, col))
            pixel_id += 1
    return chunk_tasks


def _init_worker(
    analysis_cfg: dict,
    fits_loaded: list,
    ratios: list,
    density_diagnostics: list,
    temperature_diagnostics: list,
    ref_idx: int,
    ext_law: str,
    pn_atomic_data: str,
    monte_carlo_fake_obs: int,
    max_nan_in_diagnostics_allowed_percentage: int,
    pn_rv: float,
    energy_parameter: float,
    log_level: int,
    temp_dir: str,
):
    global G_ANALYSIS_CFG
    global G_FITS_LOADED
    global G_RATIOS
    global G_DENSITY_DIAGNOSTICS
    global G_TEMPERATURE_DIAGNOSTICS
    global G_REF_IDX
    global G_EXT_LAW
    global G_PN_ATOMIC_DATA
    global G_MONTE_CARLO
    global G_MAX_NAN_PCT
    global G_PN_RV
    global G_ENERGY_PARAMETER
    global G_WORKER_LOG_LEVEL
    global G_TEMP_DIR

    G_ANALYSIS_CFG = analysis_cfg
    G_FITS_LOADED = fits_loaded
    G_RATIOS = ratios
    G_DENSITY_DIAGNOSTICS = density_diagnostics
    G_TEMPERATURE_DIAGNOSTICS = temperature_diagnostics
    G_REF_IDX = ref_idx
    G_EXT_LAW = ext_law
    G_PN_ATOMIC_DATA = pn_atomic_data
    G_MONTE_CARLO = monte_carlo_fake_obs
    G_MAX_NAN_PCT = max_nan_in_diagnostics_allowed_percentage
    G_PN_RV = pn_rv
    G_ENERGY_PARAMETER = energy_parameter
    G_WORKER_LOG_LEVEL = log_level
    G_TEMP_DIR = temp_dir


def _worker_process_chunk(args: tuple) -> str:
    chunk_idx, tasks, outdir = args
    logger = sl.setup_logger(
        f"2d_analysis_worker_{chunk_idx}",
        G_WORKER_LOG_LEVEL,
        None,
    )

    global_intensities = {}
    global_ratios = {}
    global_tene = {}
    global_ionic_abundancies = {}
    global_element_abundancies = {}
    global_icfs = {}
    failed_pixels = []

    require_positive_reference = _cfg_bool(
        G_ANALYSIS_CFG.get("require_positive_reference", False), default=False
    )

    for task_idx, (pixel_id, row, col) in enumerate(tasks):
        ref_signal = float(G_FITS_LOADED[G_REF_IDX]["signal_data"][row, col])
        ref_error = float(G_FITS_LOADED[G_REF_IDX]["error_data"][row, col])

        if not np.isfinite(ref_signal) or not np.isfinite(ref_error):
            failed_pixels.append((pixel_id, row, col, "non-finite reference signal/error"))
            continue

        if require_positive_reference and ref_signal <= 0.0:
            failed_pixels.append(
                (pixel_id, row, col, f"non-positive reference signal ({ref_signal})")
            )
            continue

        try:
            (
                pixel_intensities,
                pixel_ratios,
                pixel_tene,
                pixel_ionic_abundancies,
                pixel_icfs,
                pixel_element_abundancies,
            ) = _run_pixel_pipeline(
                pixel_id,
                row,
                col,
                G_FITS_LOADED,
                G_RATIOS,
                G_DENSITY_DIAGNOSTICS,
                G_TEMPERATURE_DIAGNOSTICS,
                G_EXT_LAW,
                G_PN_ATOMIC_DATA,
                G_MONTE_CARLO,
                G_MAX_NAN_PCT,
                G_PN_RV,
                G_ENERGY_PARAMETER,
                G_TEMP_DIR,
                logger,
            )
            global_intensities[pixel_id] = pixel_intensities
            global_ratios[pixel_id] = pixel_ratios
            global_tene[pixel_id] = pixel_tene
            global_ionic_abundancies[pixel_id] = pixel_ionic_abundancies
            global_icfs[pixel_id] = pixel_icfs
            global_element_abundancies[pixel_id] = pixel_element_abundancies
        except Exception as exc:
            failed_pixels.append((pixel_id, row, col, str(exc).replace("\n", " | ")))

        if task_idx > 0 and task_idx % 100 == 0:
            logger.info(
                f"Chunk {chunk_idx}: processed {task_idx + 1}/{len(tasks)}; successful={len(global_intensities)} failed={len(failed_pixels)}"
            )

    payload = {
        "intensities": global_intensities,
        "ratios": global_ratios,
        "tene": global_tene,
        "ionic_abundancies": global_ionic_abundancies,
        "icfs": global_icfs,
        "element_abundancies": global_element_abundancies,
        "failed_pixels": failed_pixels,
    }
    outfile = os.path.join(outdir, f"chunk_{chunk_idx:03d}.pkl")
    with open(outfile, "wb") as fout:
        pickle.dump(payload, fout, protocol=pickle.HIGHEST_PROTOCOL)
    return outfile


def _merge_chunk_pickles(chunk_files: list, logger):
    global_intensities = {}
    global_ratios = {}
    global_tene = {}
    global_ionic_abundancies = {}
    global_element_abundancies = {}
    global_icfs = {}
    failed_pixels = []

    for fn in sorted(chunk_files):
        logger.info(f"Merging {fn}")
        with open(fn, "rb") as fin:
            payload = pickle.load(fin)
        global_intensities.update(payload["intensities"])
        global_ratios.update(payload["ratios"])
        global_tene.update(payload["tene"])
        global_ionic_abundancies.update(payload["ionic_abundancies"])
        global_icfs.update(payload["icfs"])
        global_element_abundancies.update(payload["element_abundancies"])
        failed_pixels.extend(payload["failed_pixels"])

    return (
        _ordered_dict(global_intensities),
        _ordered_dict(global_ratios),
        _ordered_dict(global_tene),
        _ordered_dict(global_ionic_abundancies),
        _ordered_dict(global_icfs),
        _ordered_dict(global_element_abundancies),
        sorted(failed_pixels, key=lambda x: x[0]),
    )


def analysis2d(
    fitsd: list,
    analysis_cfg: dict,
    ratios: list,
    density_diagnostics: list,
    temperature_diagnostics: list,
    ext_law: str,
    pn_atomic_data: str,
    monte_carlo_fake_obs: int,
    max_nan_in_diagnostics_allowed_percentage: int,
    pn_rv: float,
    energy_parameter: float,
    intensities_out: str,
    ratios_out: str,
    diagnostics_out: str,
    abundancies_out: str,
    total_abundancies_out: str,
    pixel_map_out: str,
    failed_pixels_out: str,
    logger,
):
    if ext_law not in pn.RedCorr().getLaws():
        err_msg = f"ERROR. Invalid extinction given: {ext_law} "
        print(err_msg, file=sys.stderr)
        raise RuntimeError(err_msg)

    if pn_atomic_data not in pn.atomicData.getPredefinedDataFileDict().keys():
        err_msg = f"ERROR. Invalid PyNeb atomic data set: {pn_atomic_data} "
        print(err_msg, file=sys.stderr)
        raise RuntimeError(err_msg)

    pn.atomicData.includeDeprecatedPath()
    try:
        pn.atomicData.setDataFileDict(pn_atomic_data)
    except Exception:
        pn.atomicData.includeFitsPath()
        pn.atomicData.setDataFileDict(pn_atomic_data)

    if logger:
        pyneb_logger = logging.getLogger("pyneb")
        for handler in pyneb_logger.handlers[:]:
            pyneb_logger.removeHandler(handler)
        pyneb_logger.setLevel(logger.level)
        for handler in logger.handlers:
            pyneb_logger.addHandler(handler)
        pyneb_logger.propagate = False
        logging.captureWarnings(True)
        warnings.simplefilter("always")
        sys.stderr = sl.PyNebLogRedirector(logger, logging.WARNING)

    pn.log_.open_file("pyneblog.log")

    fitsd = ss.getIonTransmittionLines(fitsd, logger)
    image_shape = _validate_shapes(fitsd)
    fits_loaded = _load_all_fits_arrays(fitsd, logger)
    ref_idx = cfgio.indexOf(
        ss.reference_element["element"],
        ss.reference_element["spectrum"],
        ss.reference_element["atomic"],
        fits_loaded,
    )
    if ref_idx < 0:
        raise RuntimeError("Reference line H i 4861 was not found in the FITS list")

    row_range, col_range = _get_pixel_ranges(image_shape, analysis_cfg)
    total_pixels = len(row_range) * len(col_range)
    jobs = max(1, int(analysis_cfg.get("jobs", 1)))
    keep_chunks = _cfg_bool(analysis_cfg.get("keep_chunks", False), default=False)
    chunk_dir = analysis_cfg.get("chunk_dir", None) or None

    logger.info(
        f"Starting 2-D analysis over {total_pixels} pixels from shape {image_shape} with jobs={jobs}"
    )

    with open(pixel_map_out, "w") as fmap:
        print("# pixel_id row0 col0 y1 x1", file=fmap)
        pixel_id = 0
        for row in row_range:
            for col in col_range:
                print(f"{pixel_id} {row} {col} {row + 1} {col + 1}", file=fmap)
                pixel_id += 1

    chunk_tasks = _build_chunk_tasks(row_range, col_range, jobs)
    logger.info("Chunk sizes (pixels): %s", [len(chunk) for chunk in chunk_tasks])

    created_temp_dir = False
    if chunk_dir is None:
        chunk_dir = tempfile.mkdtemp(prefix="2d_analysis_")
        created_temp_dir = True
    else:
        os.makedirs(chunk_dir, exist_ok=True)

    _init_worker(
        analysis_cfg,
        fits_loaded,
        ratios,
        density_diagnostics,
        temperature_diagnostics,
        ref_idx,
        ext_law,
        pn_atomic_data,
        monte_carlo_fake_obs,
        max_nan_in_diagnostics_allowed_percentage,
        pn_rv,
        energy_parameter,
        logger.level,
        chunk_dir,
    )

    try:
        if jobs == 1 or len(chunk_tasks) <= 1:
            chunk_files = [_worker_process_chunk((0, chunk_tasks[0], chunk_dir))]
        else:
            available_methods = mp.get_all_start_methods()
            if "fork" in available_methods:
                ctx = mp.get_context("fork")
                with ctx.Pool(processes=len(chunk_tasks)) as pool:
                    chunk_files = pool.map(
                        _worker_process_chunk,
                        [(idx, tasks, chunk_dir) for idx, tasks in enumerate(chunk_tasks)],
                    )
            else:
                logger.warning(
                    "Multiprocessing start method 'fork' is not available; running 2-D analysis serially."
                )
                chunk_files = []
                for idx, tasks in enumerate(chunk_tasks):
                    chunk_files.append(_worker_process_chunk((idx, tasks, chunk_dir)))

        (
            global_intensities,
            global_ratios,
            global_tene,
            global_ionic_abundancies,
            global_icfs,
            global_element_abundancies,
            failed_pixels,
        ) = _merge_chunk_pickles(chunk_files, logger)

        with open(failed_pixels_out, "w") as ffail:
            print("# pixel_id row0 col0 reason", file=ffail)
            for pixel_id, row, col, reason in failed_pixels:
                print(f"{pixel_id} {row} {col} {reason}", file=ffail)

        if len(global_intensities) == 0:
            raise RuntimeError(
                "No pixels completed successfully. Check failed_pixels.dat and consider enabling require_positive_reference or reducing the region."
            )

        logger.info("Writing 2-D analysis output tables ...")
        si.printIntensities(global_intensities, intensities_out, logger)
        so.printRatios(global_ratios, ratios_out, logger)
        st.printDiagnostics(global_tene, diagnostics_out, logger)
        sa.printIonicAbundancies(global_ionic_abundancies, abundancies_out, logger)
        sf.printIcfs(
            global_icfs,
            global_element_abundancies,
            total_abundancies_out,
            logger,
        )
        _relabel_output_headers_to_pixels(
            intensities_out,
            ratios_out,
            diagnostics_out,
            abundancies_out,
            total_abundancies_out,
        )

        logger.info(
            f"2-D analysis done. Successful pixels: {len(global_intensities)} / {total_pixels}; failed pixels: {len(failed_pixels)}"
        )
        print("2-D analysis done!")
    finally:
        if 'chunk_files' in locals() and not keep_chunks:
            for fn in chunk_files:
                try:
                    os.remove(fn)
                except OSError:
                    pass
        if created_temp_dir and not keep_chunks:
            try:
                os.rmdir(chunk_dir)
            except OSError:
                pass
        pn.log_.close_file()
