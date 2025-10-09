import numpy as np
import pyneb as pn
import scipy.ndimage as nd
import copy
import sys, os
import re
import logging
import warnings

import satellite.fitsutils as fs
import satellite.intensity as si
import satellite.cfgio as sc
import satellite.roman as sr
import satellite.satlogger as sl
import satellite.tene as st
import satellite.ionic_abundancies as sa
import satellite.abundance as sb
import satellite.icf as sf
import satellite.ratios as so


def getFitsSlit(fits_fn: str, slit: dict, logger=None):
    mat = fs.rotate2d(
        fs.loadFitsImageData(fits_fn), -1.0 * slit["PA"], slit["y"] - 1, slit["x"] - 1
    )
    return fs.getVerticalSlit(
        mat, slit["y"] - 1, slit["x"] - 1, slit["w"], slit["h"], logger
    )


def extract_ion(label):
    """Example:
    t = "[OI] 5577/6300+"
    ion = extract_ion(t)
    print(ion)  # Output: "[OI]"
    """
    match = re.match(r"\[(.*?)\]", label)  # Find text inside square brackets
    return match.group(0) if match else None  # Return full [Ion] if found


reference_element = {"element": "H", "spectrum": "i", "atomic": 4861}


def findEntry(element, spectrum, atomic, _list):
    for entry in _list:
        if (
            entry["element"] == element
            and entry["spectrum"] == spectrum
            and entry["atomic"] == atomic
        ):
            return entry
    raise RuntimeError(
        f"ERROR Failed finding entry {element}{spectrum}_{atomic} in list!"
    )


def setName(slits, slit_idx, default_name):
    if len(slits) == 1:
        return default_name
    dirpath, basename = os.path.split(default_name)
    basename = f"{slit_idx:03d}_" + basename
    return os.path.join(dirpath, basename)


def radial_slit_analysis(
    fitsd: list,
    slits: list,
    ratios: list,
    density_diagnostics: list,
    tempterature_diagnostics: list,
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
    logger,
):

    # check that the extinction law is valid
    if ext_law not in pn.RedCorr().getLaws():
        err_msg = "ERROR. Invalid extinction given: {:} ".format(ext_law)
        print(err_msg, file=sys.stderr)
        raise RuntimeError(err_msg)

    # check that the atomic data set is valid
    if pn_atomic_data not in pn.atomicData.getPredefinedDataFileDict().keys():
        err_msg = "ERROR. Invalid PyNeb atomic data set: {:} ".format(pn_atomic_data)
        print(err_msg, file=sys.stderr)
        raise RuntimeError(err_msg)
    # (at least!) some of the default atomic data sets cannot be loaded unless we
    # include a 'deprecated' path.
    pn.atomicData.includeDeprecatedPath()
    try:
        pn.atomicData.setDataFileDict(pn_atomic_data)
    except:
        # if the error persists, also try this include ...
        pn.atomicData.includeFitsPath()
        pn.atomicData.setDataFileDict(pn_atomic_data)

    # attach PyNeb's logger to our, if we have one!
    if logger:
        pyneb_logger = logging.getLogger("pyneb")
        # Remove existing PyNeb handlers (Fixes duplicate log outputs)
        for handler in pyneb_logger.handlers[:]:
            pyneb_logger.removeHandler(handler)
        pyneb_logger.setLevel(logger.level)
        for handler in logger.handlers:
            pyneb_logger.addHandler(handler)
        pyneb_logger.propagate = False
        # Redirect PyNeb warnings (stderr) to our logger
        # Redirects warnings.warn() to logging.WARNING
        logging.captureWarnings(True)
        warnings.simplefilter("always")  # Ensure all warnings are captured
        sys.stderr = sl.PyNebLogRedirector(logger, logging.WARNING)

    pn.log_.open_file("pyneblog.log")

    global_intensities = {}
    global_ratios = {}
    global_tene = {}
    global_ionic_abundancies = {}
    global_element_abundancies = {}
    global_icfs = {}

    if len(slits) < 1:
        logger.error(f'No slits given; stopping now!')
        raise RuntimeError('No slits found to process!')
    if len(fitsd) < 1:
        logger.error(f'No fits given/found; stopping now!')
        raise RuntimeError('No FITS files found to process!')

    # for every slit
    for slit_idx, slit in enumerate(slits):
        print("Processing slit {:d}/{:d}".format(slit_idx + 1, len(slits)))

        # List of submatrices (whole slit), for every FITS file
        slit_subm = []
        for idx, fits in enumerate(fitsd):
            ar = getFitsSlit(fits["fns"], slit)
            np.nan_to_num(ar, False)
            ae = getFitsSlit(fits["fne"], slit)
            np.nan_to_num(ae, False)
            slit_subm.append((ar, ae))

        # how many rows per slit/submatrix ?
        rows = slit_subm[0][0].shape[0]
        if rows < 1:
            logger.error(f'Height of slit is {rows}! Cannot operate on zero rows, skipping slit')
            continue

        for row_idx in range(rows):
            print(f"Processing row {row_idx}/{rows}")
            # copy of dictionary
            cpd = copy.deepcopy(fitsd)

            # process every FITS, each row
            for idx, _ in enumerate(fitsd):
                # sum all elements of row (FITS)
                pics = slit_subm[idx]
                cpd[idx]["sslit_sum"] = slit_subm[idx][0][row_idx].sum()
                cpd[idx]["eslit_sum"] = slit_subm[idx][1][row_idx].sum()

            # compile the intensities data file (for PyNeb) and write the test.dat file.
            # TODO we do not need to pass the cpd list here. We can pass a more simple/small
            # list.
            si.makeIntensitiesDataFile(
                cpd, reference_element, ["sslit_sum", "eslit_sum"], "test.dat"
            )

            # PyNeb stuff; PyNeb will read the 'test.dat' file (for the slit).
            sobs = pn.Observation()
            sobs.readData(
                "test.dat", fileFormat="lines_in_rows_err_cols", errIsRelative=False
            )
            sobs.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)
            sobs.extinction.law = ext_law
            sobs.correctData(normWave=4861.0)

            # create Monte Carlo simulations untill TeNe diagnostics contains no nan
            nan_diagnostics = True
            times_nan_encountered = 0
            MAX_NAN_IN_DIAGNOSTICS_ALLOWED = (
                max_nan_in_diagnostics_allowed_percentage * monte_carlo_fake_obs // 100
            )
            # MAX_NAN_IN_DIAGNOSTICS_ALLOWED cannot be zeros, or else we won;t get into the loop.
            MAX_NAN_IN_DIAGNOSTICS_ALLOWED = max(1, MAX_NAN_IN_DIAGNOSTICS_ALLOWED)
            while (
                nan_diagnostics
                and times_nan_encountered < MAX_NAN_IN_DIAGNOSTICS_ALLOWED
            ):
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

                # Compute intensity for each FITS/atom; add to global dictionary for printing
                # later on.
                """Example: global_intensities[1] = {'intensities': [...], 'E_BV': rc.E_BV, 'cHbeta': rc.cHbeta} """
                global_intensities[row_idx] = {
                    "intensities": si.computeIntensities(
                        fitsd, sobs, eobs, RC, reference_element, logger
                    ),
                    "E_BV": RC.E_BV,
                    "E_BVError": ebv_err,
                    "cHbeta": RC.cHbeta,
                    "cHbetaError": chbeta_err,
                    "fac": f,
                    "FHb": findEntry(
                        reference_element["element"],
                        reference_element["spectrum"],
                        reference_element["atomic"],
                        cpd,
                    )["sslit_sum"]
                    * energy_parameter,
                    "FHb_error": findEntry(
                        reference_element["element"],
                        reference_element["spectrum"],
                        reference_element["atomic"],
                        cpd,
                    )["eslit_sum"]
                    * energy_parameter,
                }

                # Compute intensity ratios
                global_ratios[row_idx] = {}
                for ratio in ratios:
                    try:
                        global_ratios[row_idx][ratio] = so.computeRatio(
                            ratio, global_intensities[row_idx]["intensities"]
                        )
                    except:
                        logger.info("Skipping ratio {:}".format(ratio))

                # Compute diagnostics (te/ne pairs)
                global_tene[row_idx], nan_diagnostics = st.computeTeNePairs(
                    density_diagnostics, tempterature_diagnostics, sobs, eobs, logger
                )
                if nan_diagnostics:
                    logger.warning(
                        f"Encountered nan value in diagnostics; restarting computations for slit! ({times_nan_encountered}/{MAX_NAN_IN_DIAGNOSTICS_ALLOWED})"
                    )
                times_nan_encountered += 1

            if nan_diagnostics or (
                times_nan_encountered >= MAX_NAN_IN_DIAGNOSTICS_ALLOWED
            ):
                msg = f"ERROR. Failed computing non-nan diagnostics after {times_nan_encountered} tries. Giving up!"
                logger.error(msg)
                # raise RuntimeError(msg)
                continue

            # Compute Ionic Abundancies
            logger.debug("> Calling computeIonicAbundancies ...")
            global_ionic_abundancies[row_idx] = sa.computeIonicAbundancies(
                cpd, global_tene[row_idx], sobs, eobs, logger
            )

            # Compute abundancies per element, e.g.
            # [...
            # 'He2': (np.float64(0.025174496324910436), np.float64(0.00020749028502411882)),
            # 'Ar3': (np.float64(5.035899412128091e-07), np.float64(2.0404083173261792e-08)),
            # 'Cl3': (np.float64(2.166739769269426e-08), np.float64(9.154696103324478e-10)),
            #  'N1': (np.float64(2.4613015709276348e-08), np.float64(9.722314034253516e-09)),
            # ...]
            logger.debug("> Calling computeAbundancies ...")
            elemspec_abundancies = sb.computeAbundancies(
                cpd, global_ionic_abundancies[row_idx], logger
            )

            logger.debug("> Calling computeIcfsWithErrors ...")
            global_icfs[row_idx] = sf.computeIcfsWithErrors(
                elemspec_abundancies, logger
            )
            logger.debug("> Calling  ionicAbundance2elementAbundance ...")
            global_element_abundancies[row_idx] = sf.ionicAbundance2elementAbundance(
                elemspec_abundancies, logger
            )
            # ----> New part <----
            # global_icfs[slit_idx] = sf.computeManualIcfs(elemspec_abundancies, logger)

        ## <-- End Looping Slits --> ##

        print(f"All rows done for slit #{slit_idx} ... writing files ...")
        if len(slits) > 1:
            dirpath, basename = os.path.split(intensities_out)
        si.printIntensities(
            global_intensities, setName(slits, slit_idx, intensities_out), logger
        )
        so.printRatios(global_ratios, setName(slits, slit_idx, ratios_out), logger)
        st.printDiagnostics(
            global_tene, setName(slits, slit_idx, diagnostics_out), logger
        )
        sa.printIonicAbundancies(
            global_ionic_abundancies, setName(slits, slit_idx, abundancies_out), logger
        )
        sf.printIcfs(
            global_icfs,
            global_element_abundancies,
            setName(slits, slit_idx, total_abundancies_out),
            logger,
        )

    print("All done!")
    pn.log_.close_file()
