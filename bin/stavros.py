#!/usr/bin/env python

#from satellite import cfgio as io
from satellite import version
from satellite import cfgio
import argparse

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

if __name__ == "__main__":

# parse cmd
    args = parser.parse_args()

# verbose print
    verboseprint = print if args.verbose else lambda *a, **k: None
    
    verboseprint('Satellite (module) version: {:}'.format(version.satellite_version()))
    verboseprint('Satellite (binary) version: {:}'.format(satellite_version))

# parse the config file
    config = cfgio.parseConfigInout(args.config)

# get the FITS file list, with corresponding elements
    fits_info = cfgio.configFitsFileList(config)
    for idx, fits in enumerate(fits_info): 
        processFits( fits['fn'] ) #TODO
        fits[idx]['values'] = # TODO

# ion abundances
