import sys
import pyneb as pn

from satellite import cfgio
from satellite import astroflux

def makeIntensitiesDataFile(fitsd: list, reference_element: dict, value_keys: list, fn: str, factor=100e0):
    ref_index = cfgio.indexOf(reference_element['element'], reference_element['spectrum'], reference_element['atomic'], fitsd)
    ref_sval = fitsd[ref_index][value_keys[0]]
    ref_eval = fitsd[ref_index][value_keys[1]]
    with open(fn, 'w') as fout:
        print('LINE test err', file=fout)
        for obj in fitsd:
            pnlabel = cfgio.objectIntensityPyNebCode(obj['element'], obj['spectrum'], obj['atomic'])
            print('{:} {:+9e} {:+9e}'.format(pnlabel, factor*obj[value_keys[0]]/ref_sval, astroflux.fluxError(obj[value_keys[0]], obj[value_keys[1]], ref_sval, ref_eval)), file=fout)
                
    return fn
