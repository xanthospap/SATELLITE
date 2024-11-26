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
    '--verbose',
    action='store_true',
    dest='verbose',
    help='Verbose mode on')


def checkInputFits(fits_info: dict):
    missing_files = []
    for idx, fits in enumerate(fits_info):
        file_is_missing = False
        fitsfn = fits['fn']
        if not os.path.isfile(fitsfn):
            print("[WRNNG] Missing Fits file {:}".format(fitsfn), file=sys.stderr)
            bn = os.path.basename(fitsfn)
            file_is_missing = True
# find spectrum in filename
            sstr = fits['spectrum']
# replace roman spectrum with int and see if file exists
            if bn.find(sstr) >= 0:
                gfn = bn.replace(sstr, str(roman.roman2int(sstr)), 1)
                guess = os.path.join(os.path.dirname(fitsfn), gfn)
                if os.path.isfile(guess):
# change the name in the return list
                    print("[DEBUG] Fits filename {:} is missing; using {:}".format(os.path.basename(fitsfn), gfn), file=sys.stderr)
                    fits_info[idx]['fn'] = guess
                    file_is_missing = False
            if file_is_missing and bn.find(sstr) >= 0:
                gfn = bn.replace(sstr, sstr.upper(), 1)
                guess = os.path.join(os.path.dirname(fitsfn), gfn)
                if os.path.isfile(guess):
                    print("[DEBUG] Fits filename {:} is missing; using {:}".format(os.path.basename(fitsfn), gfn), file=sys.stderr)
                    fits_info[idx]['fn'] = guess
                    file_is_missing = False
        if file_is_missing: missing_files.append(fitsfn)
    return fits_info, missing_files

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
    fits_info, missing_files = checkInputFits(cfgio.configFitsFileList(config))
    if missing_files != []:
        err = '\n'.join(missing_files)
        print("[ERROR] Missing FITS files: {:}".format(err), file=sys.stderr)
        sys.exit(1)

# Specific Slit Analysis
    slits = cfgio.configSpecificSlitAnalysis(config)
    specific_slit.specific_slit_analysis(fits_info, slits)

# get the FITS file list, with corresponding elements
#    fits_info = cfgio.configFitsFileList(config)
#    for idx, fits in enumerate(fits_info): 
#        print(idx, fits)
#        # processFits( fits['fn'] ) #TODO
#        # fits[idx]['values'] = # TODO
