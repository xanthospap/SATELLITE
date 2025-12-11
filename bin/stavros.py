#!/usr/bin/env python

from satellite import version
from satellite import cfgio
from satellite import roman
from satellite import specific_slit
from satellite import angular_slit
from satellite import radial_slit
from satellite import satlogger
from satellite import plotters

import argparse
import os
import sys
import logging

satellite_version = "2.r1"


class myFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter
):
    pass


parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description="Spectroscopic Analysis Tool for intEgraL FieLd unIt daTacubEs",
    epilog=(
        """National Observatory of Athens,
    Institute for Astronomy, Astrophysics, Space Applications and Remote Sensing\n
    Send bug reports to:
    Stavros Akras, xanthos@mail.ntua.gr
    May, 2024"""
    ),
)

parser.add_argument(
    "-i",
    "--config-file",
    metavar="CONFIG",
    dest="config",
    default=None,
    required=True,
    help="Configuration file for analysis.",
)

parser.add_argument(
    "-a",
    "--atomic-data-set",
    metavar="ATOMIC_DATA_SET",
    dest="pn_atomic_data",
    default="PYNEB_21_01",
    required=False,
    help="Predefined atomic data set provided by PyNeb. For a complete list, see https://github.com/Morisset/PyNeb_devel/blob/master/docs/Notebooks/PyNeb_manual_3.ipynb",
)

parser.add_argument(
    "-e",
    "--extinction-law",
    metavar="PN_extinction",
    dest="pn_extinction",
    default="S79 H83 CCM89",
    required=False,
    help="Predefined extinction (reddening) correction. See https://notebook.community/Morisset/PyNeb_devel/docs/Notebooks/PyNeb_manual_5",
)

parser.add_argument(
    "-m",
    "--monte-carlo-simulations",
    metavar="MONTE_CARLO_FAKE_OBS",
    dest="monte_carlo_fake_obs",
    type=int,
    default=4,
    required=False,
    help="Number of Mone-Carlo fake observations for error estimation.",
)

parser.add_argument(
    "--max-allowed-mc-tries",
    metavar="MAX_ALLOWED_MC_TRIES",
    dest="max_allowed_mc_tries",
    type=int,
    default=15,
    required=False,
    help="Sometimes the Monte Carlo simulated obs result in values that are invalid, producing Nan values for various computations. This value, given as a percentage, allows the reproduction of Monte Carlo simulated values up to some limit w.r.t the MONTE_CARLO_FAKE_OBS. E.g. if MONTE_CARLO_FAKE_OBS=500 and MAX_ALLOWED_MC_TRIES=15, then we can try at maximum 15%% of 500 i.e. 75 times.",
)

parser.add_argument(
    "-RV",
    metavar="R_V",
    dest="pn_rv",
    type=float,
    default=3.1e0,
    required=False,
    help="Reddening correction RV parameter for Pyneb; note that R_V = AV/E_BV. Default value is 3.1",
)

parser.add_argument(
    "-EP",
    metavar="ENERGY_PARAMETER",
    dest="energy_parameter",
    type=float,
    default=1e0,
    required=False,
    help="Energy parameter value for scaling F(Hb).",
)

parser.add_argument(
    "--intensities-out",
    metavar="INTENSITIES_OUTPUT_FILE",
    dest="intensities_out",
    default="line_intensities.dat",
    required=False,
    help="Output file to write computed inensity results.",
)

parser.add_argument(
    "--ratios-out",
    metavar="RATIOS_OUTPUT_FILE",
    dest="ratios_out",
    default="line_ratios.dat",
    required=False,
    help="Output file to write computed (line) ratio results.",
)

parser.add_argument(
    "--diagnostics-out",
    metavar="DIAGNOSTICS_OUTPUT_FILE",
    dest="diagnostics_out",
    default="tn_diagnostics.dat",
    required=False,
    help="Output file to write computed temperature/density diagnostics results.",
)

parser.add_argument(
    "--abundances-out",
    metavar="IONIC_ABUNDANCIES_OUTPUT_FILE",
    dest="abundancies_out",
    default="abundancies.dat",
    required=False,
    help="Output file to write computed ionic abundancies results.",
)

parser.add_argument(
    "--total-abundances-out",
    metavar="TOTAL_ABUNDANCIES_OUTPUT_FILE",
    dest="total_abundancies_out",
    default="total_abundancies.dat",
    required=False,
    help="Output file to write computed total abundancies, including ICFs and DIMSs.",
)

parser.add_argument(
    "--verbose", action="store_true", dest="verbose", help="Verbose mode on"
)

parser.add_argument(
    "--missing-fits-is-error",
    action="store_true",
    dest="missing_fits_is_error",
    help="Fail with error if any of the input FITS file specified in the given config file is missing.",
)

parser.add_argument(
    "--bar-plots",
    action="store_true",
    dest="barplot",
    help="Plot using barplots (excluding total abundancies) instead of scatter plots. Only relevant if '--no-plots' is not set.",
)

parser.add_argument(
    "--log",
    metavar="LOG_FILE",
    dest="log_file",
    default=None,
    required=False,
    help="Write log to an output file.",
)

parser.add_argument(
    "--no-plots", action="store_true", dest="no_plots", help="Do no produce plotts."
)

if __name__ == "__main__":

    # parse cmd
    args = parser.parse_args()

    # setup a logger
    logger = satlogger.setup_logger("specific_slit", logging.DEBUG, args.log_file)

    # parse the config file
    config = cfgio.parseConfigInout(args.config)

    # check input FITS files
    fits_info, missing_files = cfgio.checkInputFits(
        cfgio.configFitsFileList(config), logger
    )
    if missing_files != []:
        err = "\n".join(missing_files)
        logger.warning("Missing FITS files: {:}".format(err))
        if args.missing_fits_is_error:
            sys.exit(1)

        # try:
    if cfgio.doSpecificSlitAnalysis(config):
        # Specific Slit Analysis
        specific_slit.specific_slit_analysis(
            fits_info,
            cfgio.configSpecificSlitAnalysis(config),
            cfgio.configElementRatiosList(config),
            cfgio.configDensityDiagnostics(config),
            cfgio.configTemperatureDiagnostics(config),
            args.pn_extinction,
            args.pn_atomic_data,
            args.monte_carlo_fake_obs,
            args.max_allowed_mc_tries,
            args.pn_rv,
            args.energy_parameter,
            args.intensities_out,
            args.ratios_out,
            args.diagnostics_out,
            args.abundancies_out,
            args.total_abundancies_out,
            "specific_slit_corners.dat",
            logger,
        )

    if cfgio.doAngularSlitAnalysis(config):
        # Angular Slit Analysis
        angular_slit.angular_slit_analysis(
            fits_info,
            cfgio.configAngularSlitAnalysis(config),
            cfgio.configElementRatiosList(config),
            cfgio.configDensityDiagnostics(config),
            cfgio.configTemperatureDiagnostics(config),
            args.pn_extinction,
            args.pn_atomic_data,
            args.monte_carlo_fake_obs,
            args.max_allowed_mc_tries,
            args.pn_rv,
            args.energy_parameter,
            args.intensities_out,
            args.ratios_out,
            args.diagnostics_out,
            args.abundancies_out,
            args.total_abundancies_out,
            "angular_slit_corners.dat",
            logger,
        )

    if cfgio.doRadialSlitAnalysis(config):
        # Radial Slit Analysis
        radial_slit.radial_slit_analysis(
            fits_info,
            cfgio.configRadialSlitAnalysis(config),
            cfgio.configElementRatiosList(config),
            cfgio.configDensityDiagnostics(config),
            cfgio.configTemperatureDiagnostics(config),
            args.pn_extinction,
            args.pn_atomic_data,
            args.monte_carlo_fake_obs,
            args.max_allowed_mc_tries,
            args.pn_rv,
            args.energy_parameter,
            args.intensities_out,
            args.ratios_out,
            args.diagnostics_out,
            args.abundancies_out,
            args.total_abundancies_out,
            "radial_slit_corners.dat",
            logger,
        )
    # except Exception as e:
    #    logger.error(f"Analysis failed! Error was: {e}")
    #    print(f"Analysis failed! Error was: {e}")
    #    sys.exit(9)

    # if needed draw the plots ...
    if not args.no_plots:
        print("Compiling plots ...")
        plotters.plotLineAbundancies(
            args.abundancies_out, "line_abundancies.pdf", args.barplot
        )
        plotters.plotLineIntensities(
            args.intensities_out, "line_intensities.pdf", args.barplot
        )
        plotters.plotLineRatios(args.ratios_out, "line_ratios.pdf", args.barplot)
        plotters.PlotTeNeDiagnostics(
            args.diagnostics_out, "temperature.pdf", "density.pdf", args.barplot
        )
        plotters.plotTotalAbundancies(
            args.total_abundancies_out, "total_abundancies.pdf"
        )
