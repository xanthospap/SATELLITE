#! /usr/bin/python

import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

parser = argparse.ArgumentParser(
    formatter_class=RawTextHelpFormatter,
    description=
    'Python FITS viewer.',
    epilog=(r'''National Observatory of Athens,
    Institute for Astronomy, Astrophysics, Space Applications and Remote Sensing.
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    Jan, 2024'''))

parser.add_argument(
    '-f',
    '--fits-file',
    metavar='FITS_FILE',
    dest='fits',
    default=None,
    required=True,
    #type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'),
    help='FITS file to view.')

parser.add_argument(
    '-p',
    '--print-info',
    dest='show_info',
    action = 'store_true',
    help='Print basic FITS info to STDOUT.')

parser.add_argument(
    '--histogram',
    dest='hist',
    action = 'store_true',
    help='Append histogram to image.')

def printInfo(fitsfn, image_data):
    print('FITS filename : {:}'.format(fitsfn))
    print('Num pixels    : {:}x{:}'.format(image_data.shape[0], image_data.shape[1]))
    print('Min           : {:}'.format(np.min(image_data)))
    print('Max           : {:}'.format(np.max(image_data)))
    print('Mean          : {:}'.format(np.mean(image_data)))
    print('Stdev         : {:}'.format(np.std(image_data)))

def fitsPlot(fitsfn, append_histogram=False, show_info=False):
    ftsin = fits.open(fitsfn)
    image_data = ftsin[0].data
    ftsin.close()
    if show_info: printInfo(fitsfn, image_data)
    if not append_histogram:
        plt.imshow(image_data, cmap='gray')
        plt.colorbar()
        #plt.show()
    else:
        fig, axis = plt.subplots(1,2)
        axis[0].set_title(fitsfn)
        axis[0].imshow(image_data, cmap='gray')
        #axis[0].colorbar()
        axis[1].set_title('Histogram')
        axis[1].hist(image_data.flatten(), bins='auto')
    plt.show()

if __name__ == "__main__":
    args = parser.parse_args()
    fitsPlot(args.fits, args.hist, args.show_info)
