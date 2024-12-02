import sys
import re
import pyneb as pn

from satellite import roman
from satellite import cfgio
from satellite import astroflux

def closestPyNebElement(key: str, atomic_number: int):
    min_diff = 1000000
    pn_wl = None
    for wl in pn.LINE_LABEL_LIST[key]:
        if re.fullmatch('[0-9]*A', wl):
            diff = abs(int(wl[0:-1]) - atomic_number)
            if diff < min_diff: 
                pn_wl = wl
                min_diff = diff
    print("[WRNNG] PyNeb is missing wavelength {:} for element {:}; using {:} instead".format(atomic_number, key, pn_wl))
    return pn_wl

def objectIntensityPyNebCode(atom: str, spectrum: str, atomic_number: int):
#
# A full list can be obtained as:
# import pyneb as pn
# pn.LINE_LABEL_LIST
#
    if atom in ['H', 'He']:
        pnatom = "{:}{:}r".format(atom, roman.roman2int(spectrum))
    else:
        pnatom = "{:}{:}".format(atom, roman.roman2int(spectrum))
    
    try:
        wls = pn.LINE_LABEL_LIST[pnatom]
    except:
        wls = None
        print('[ERROR] Failed matching object {:}/{:} (aka {:}) to PyNeb (see LINE_LABEL_LIST)'.format(atom, spectrum, pnatom), file=sys.stderr)
    if wls is None: 
        raise RuntimeError('[ERROR] Failed matching object {:}/{:}/{:} to PyNeb (see LINE_LABEL_LIST)'.format(atom, spectrum, atomic_number))

    pnatomic = '{:}A'.format(atomic_number)
    if pnatomic not in wls:
         pnatomic = closestPyNebElement(pnatom, atomic_number)

    pnstr = '_'.join([pnatom, pnatomic])

# Validate
    try:
        pn.LINE_LABEL_LIST[pnstr.split('_')[0]].index(pnstr.split('_')[1])
        return pnstr
    except:
        print('[ERROR] Failed matching object {:}/{:}/{:} to PyNeb (see LINE_LABEL_LIST)'.format(atom, spectrum, atomic_number), file=sys.stderr)
    raise RuntimeError('[ERROR] Failed matching object {:}/{:}/{:} to PyNeb (see LINE_LABEL_LIST)'.format(atom, spectrum, atomic_number))

#def getMatchingErrorFits(fits_obj: dict, fits_list: list):
#    for idx, obj in enumerate(fits_list):
#        print("\tfn={:}, type={:}".format(obj['fn'],obj['type']))
#        if obj['type'] == 'error':
#            matched = True
#            print("\tcomparing [{:}]vs[{:}] [{:}]vs[{:}] [{:}]vs[{:}]".format(obj['element'], fits_obj['element'], obj['spectrum'], fits_obj['spectrum'],obj['atomic'], fits_obj['atomic']))
#            for key in ['element', 'spectrum', 'atomic']:
#                if obj[key] != fits_obj[key]: matched = False
#            if matched:
#                return idx
#    atom, spectrum, atomic_number = fits_obj['element'], fits_obj['spectrum'], fits_obj['atomic']
#    raise RuntimeError('[ERROR] Failed cannot find a matching error FITS for element {:}/{:}/{:} '.format(atom, spectrum, atomic_number))

def makeIntensitiesDataFile(fitsd: list, reference_element: dict, value_keys: list, fn: str, factor=100e0):
    ref_index = cfgio.indexOf(reference_element['element'], reference_element['spectrum'], reference_element['atomic'], fitsd)
    ref_sval = fitsd[ref_index][value_keys[0]]
    ref_eval = fitsd[ref_index][value_keys[1]]
    with open(fn, 'w') as fout:
        print('LINE test err', file=fout)
        for obj in fitsd:
            pnlabel = objectIntensityPyNebCode(obj['element'], obj['spectrum'], obj['atomic'])
            print('{:} {:+9e} {:+9e}'.format(pnlabel, factor*obj[value_keys[0]]/ref_sval, astroflux.fluxError(obj[value_keys[0]], obj[value_keys[1]], ref_sval, ref_eval)), file=fout)
                
    return fn
