#!/usr/bin/env python

from satellite import version
from satellite import cfgio
from satellite import roman
from satellite import specific_slit

import argparse
import os
import sys

satellite_version = '2.r1'

class myFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawTextHelpFormatter):
    pass

parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description=
    'Spectroscopic Analysis Tool for intEgraL FieLd unIt daTacubEs',
    epilog=('''National Observatory of Athens,
    Institute for Astronomy, Astrophysics, Space Applications and Remote Sensing\n
    Send bug reports to:
    Stavros Akras, xanthos@mail.ntua.gr
    May, 2024'''))

parser.add_argument(
    '-i',
    '--config-file',
    metavar='CONFIG',
    dest='config',
    default=None,
    required=True,
    help='Configuration file for analysis.')

parser.add_argument(
    '-a',
    '--atomic-data-set',
    metavar='ATOMIC_DATA_SET',
    dest='pn_atomics',
    default='PYNEB_21_01',
    required=False,
    help='Predefined atomic data set provided by PyNeb. For a complete list, see https://github.com/Morisset/PyNeb_devel/blob/master/docs/Notebooks/PyNeb_manual_3.ipynb')

parser.add_argument(
    '-e',
    '--extinction',
    metavar='PN_extinction',
    dest='pn_extinction',
    default='S79 H83 CCM89',
    required=False,
    help='Predefined extinction (reddening) correction. See https://notebook.community/Morisset/PyNeb_devel/docs/Notebooks/PyNeb_manual_5')

parser.add_argument(
    '--intensities-out',
    metavar='INTENSITIES_OUTPUT_FILE',
    dest='intensities_out',
    default='output_line_intensities',
    required=False,
    help='Output file to write computed inensity results.')

parser.add_argument(
    '--ratios-out',
    metavar='RATIOS_OUTPUT_FILE',
    dest='ratios_out',
    default='output_line_ratios',
    required=False,
    help='Output file to write computed ratio results.')

parser.add_argument(
    '--verbose',
    action='store_true',
    dest='verbose',
    help='Verbose mode on')


if __name__ == "__main__":

# parse cmd
    args = parser.parse_args()

# verbose print
    verboseprint = print if args.verbose else lambda *a, **k: None
    
    verboseprint('Satellite (module) version: {:}'.format(version.satellite_version()))
    verboseprint('Satellite (binary) version: {:}'.format(satellite_version))

# parse the config file
    config = cfgio.parseConfigInout(args.config)

# check input FITS files
    fits_info, missing_files = cfgio.checkInputFits(cfgio.configFitsFileList(config))
    if missing_files != []:
        err = '\n'.join(missing_files)
        print("[ERROR] Missing FITS files: {:}".format(err), file=sys.stderr)
        sys.exit(1)

# Specific Slit Analysis
    specific_slit.specific_slit_analysis(fits_info, cfgio.configSpecificSlitAnalysis(config), cfgio.configElementRatiosList(config), args.intensities_out, args.ratios_out)
