import satellite.roman as sr
import pyneb as pn
import re

def closestPyNebElement(key: str, atomic_number: int, logger=None):
    min_diff = 1000000
    pn_wl = None
    for wl in pn.LINE_LABEL_LIST[key]:
        if re.fullmatch('[0-9]*A', wl):
            diff = abs(int(wl[0:-1]) - atomic_number)
            if diff < min_diff:
                pn_wl = wl
                min_diff = diff
    if logger:
        logger.warning("PyNeb is missing wavelength {:} for element {:}; using {:} instead".format(
            atomic_number, key, pn_wl))
    return pn_wl

def objectIntensityPyNebCode(atom: str, spectrum: str, atomic_number: int, logger=None):
    #
    # A full list can be obtained as:
    # import pyneb as pn
    # pn.LINE_LABEL_LIST
    #
    try:
        spectrum = int(spectrum)
    except:
        spectrum = sr.roman2int(spectrum)
    if atom in ['H', 'He']:
        pnatom = "{:}{:}r".format(atom, spectrum)
    else:
        pnatom = "{:}{:}".format(atom, spectrum)

    try:
        wls = pn.LINE_LABEL_LIST[pnatom]
    except:
        wls = None
        if logger:
            logger.error(
                "Failed matching object {:}/{:} (aka {:}) to PyNeb (see LINE_LABEL_LIST)".format(atom, spectrum, pnatom))
    if wls is None:
        raise RuntimeError(
            '[ERROR] Failed matching object {:}/{:}/{:} to PyNeb (see LINE_LABEL_LIST)'.format(atom, spectrum, atomic_number))

    pnatomic = '{:}A'.format(atomic_number)
    if pnatomic not in wls:
        pnatomic = closestPyNebElement(pnatom, atomic_number)

    pnstr = '_'.join([pnatom, pnatomic])

# Validate
    try:
        pn.LINE_LABEL_LIST[pnstr.split('_')[0]].index(pnstr.split('_')[1])
        return pnstr
    except:
        if logger:
            logger.error(
                '[ERROR] Failed matching object {:}/{:}/{:} to PyNeb (see LINE_LABEL_LIST)'.format(atom, spectrum, atomic_number))
    raise RuntimeError(
        '[ERROR] Failed matching object {:}/{:}/{:} to PyNeb (see LINE_LABEL_LIST)'.format(atom, spectrum, atomic_number))

def satellite_str2pyneb_str(sstr: str, logger=None):
    parts = sstr.split('_')
    assert (len(parts) == 2)
    atomic_number = int(parts[1])
    g = re.fullmatch("([A-Za-z]*)([0-9])", parts[0].strip())
    if g:
        return objectIntensityPyNebCode(g[1], g[2], atomic_number, logger)
    g = re.fullmatch("([A-Za-z]*)([iIvV]*)(r+)", parts[0].strip())
    if g:
        return objectIntensityPyNebCode(g[1], sr.roman2int(g[2]), atomic_number, logger)
    if logger:
        logger.error("Failed resolving satellite string {:}".format(sstr))
    raise RuntimeError(
        "[ERROR] Failed resolving satellite string {:}".format(sstr))

