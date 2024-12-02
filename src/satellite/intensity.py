import pyneb as pn
import sys
from satellite import roman

def objectIntensityPyNebCode(atom: str, spectrum: str, atomic_number: int):
#
# A full list can be obtained as:
# import pyneb as pn
# pn.LINE_LABEL_LIST
#
    if atom == 'H' and spectrum == 'i':
        pnstr = '{:}_{:}A'.format('H1r', atomic_number)
    elif atom == 'He':
        pnstr = '{:}{:}r_{:}A'.format(atom, roman.roman2int(spectrum), atomic_number)
    else:
        pnstr = '{:}{:}_{:}A'.format(atom, roman.roman2int(spectrum), atomic_number)
# Validate
    try:
        pn.LINE_LABEL_LIST[pnstr.split('_')[0]].index(pnstr.split('_')[1])
        return pnstr
    except:
        print('[ERROR] Failed matching object {:}/{:}/{:} to PyNeb (see LINE_LABEL_LIST)'.format(atom, spectrum, atomic_number), file=sys.stderr)
# Didn't return, that's an error ...
    raise RuntimeError('[ERROR] Failed matching object {:}/{:}/{:} to PyNeb (see LINE_LABEL_LIST)'.format(atom, spectrum, atomic_number))

def makeIntensitiesDataFile(fitsd: list, denominator: int, value_key: str, fn: str, factor=100e0):
    # denominator reference value
    denom = fitsd[denominator][value_key]
    with open(fn, 'w') as fout:
        for obj in fitsd:
            pnlabel = objectIntensityPyNebCode(obj['element'], obj['spectrum'], obj['atomic'])
            print('{:} {:07.3f} {:07.3f}'.format(pnlabel, obj[value_key]/denom, obj[value_key]/denom/10.), file=fout)
    return fn
