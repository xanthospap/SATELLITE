import satellite.roman as sr
import pyneb as pn
import numpy as np
import re

def closestPyNebElement(element: str, atomic_number: int, logger=None) -> str:
    """ Given an an elemnt and its atomic number, find the closest 
        corresponding emission line available in PyNeb

        Note
        ----
        The available emission lines are placed in pyneb.LINE_LABEL_LIST.
        E.g.
        >>> pn.LINE_LABEL_LIST['O']
        will give all available emission lines for 'O'.

        Parameters
        ----------
        element: str
            The element to consider (e.g. 'H', 'He', 'O')
        atomic_number: int
            The atomic number as an integer value (e.g. 5199, 6583)
        logger: logging
            This parameter can be set to an instance of logging to
            write warnings on a specified stream. Default value is
            None.
        
        Returns
        -------
        str:
            The closest emission line available in PyNeb, i.e. a line 
            for the same element, with a value as close as possible to 
            atomic_number.

        Example
        -------
        >>> closestPyNebElement('N1', 5199)   -> '5200A'
        >>> closestPyNebElement('N2', 6583)   -> '6584A'
        >>> closestPyNebElement('Cl3', 5517)  -> '5518A'
        >>> closestPyNebElement('He2r', 5412) -> '5411A'
    """
    min_diff = 1000000
    pn_wl = None
    for wl in pn.LINE_LABEL_LIST[element]:
        if re.fullmatch('[0-9]*A', wl):
            diff = abs(int(wl[0:-1]) - atomic_number)
            if diff < min_diff:
                pn_wl = wl
                min_diff = diff
    if logger:
        logger.warning("PyNeb is missing wavelength {:} for element {:}; using {:} instead".format(
            atomic_number, element, pn_wl))
    return pn_wl

def unlisted(atom, spectrum_int, wavelength, logger=None):
    
    # if logger: logger.info(f'Looking for a data file to load {atom}{spectrum_int}_{wavelength} ...')
    entries = pn.atomicData.getDataFile()[f"{atom}{spectrum_int}"]
    if 'atom' in entries:
        #if entries['atom'] not in pn.atomicData.getDataFile('{atom}{spectrum}'):
        pn.atomicData.setDataFile(entries['atom'])
        if logger: logger.info(f"Loaded {entries['atom']}")
    if 'coll' in entries:
        #if entries['coll'] not in pn.atomicData.getDataFile('{atom}{spectrum}'):
        pn.atomicData.setDataFile(entries['coll'])
        if logger: logger.info(f"Loaded {entries['coll']}")
    elmnt = pn.Atom(atom, spectrum_int)
    waves = np.array(elmnt.lineList)
    target = float(wavelength)
    i = np.argmin(np.abs(waves - target))
    closest_wave = waves[i]
    if logger: logger.info(f"Closest line to {atom}{spectrum_int}_{wavelength} is line {closest_wave}")
    return f'{atom}{spectrum_int}_{int(closest_wave)}'

def objectIntensityPyNebCode(atom: str, spectrum: str, atomic_number: int, logger=None):
    """
        A full list can be obtained as:
        import pyneb as pn
        pn.LINE_LABEL_LIST

        Example
        -------
        >>> objectIntensityPyNebCode('Ar', '3', 7136) -> 'Ar3_7136A'
        >>> objectIntensityPyNebCode('Cl', '3', 5517) -> 'Cl3_5518A'
        >>> objectIntensityPyNebCode('Cl', '3', 5538) -> 'Cl3_5538A'
        >>> objectIntensityPyNebCode('H', '1', 4861) -> 'H1r_4861A'
        >>> objectIntensityPyNebCode('H', '1', 6563) -> 'H1r_6563A'
        >>> objectIntensityPyNebCode('H', '1', 4861) -> 'H1r_4861A'
        >>> objectIntensityPyNebCode('He', '1', 5876) -> 'He1r_5876A'
        >>> objectIntensityPyNebCode('He', '1', 6678) -> 'He1r_6678A'
        >>> objectIntensityPyNebCode('He', '2', 5412) -> 'He2r_5411A'
        >>> objectIntensityPyNebCode('N', '1', 5199) -> 'N1_5200A'
        >>> objectIntensityPyNebCode('N', '2', 5755) -> 'N2_5755A'
        >>> objectIntensityPyNebCode('N', '2', 6548) -> 'N2_6548A'
        >>> objectIntensityPyNebCode('N', '2', 6583) -> 'N2_6584A'
        >>> objectIntensityPyNebCode('O', '1', 6300) -> 'O1_6300A'
        >>> objectIntensityPyNebCode('O', '2', 7320) -> 'O2_7320A'
        >>> objectIntensityPyNebCode('O', '2', 7330) -> 'O2_7330A'
        >>> objectIntensityPyNebCode('O', '3', 4959) -> 'O3_4959A'
        >>> objectIntensityPyNebCode('O', '3', 5007) -> 'O3_5007A'
        >>> objectIntensityPyNebCode('S', '2', 6716) -> 'S2_6716A'
    """
    ## Spectrum as integer value
    try:
        spectrum = int(spectrum)
    except:
        spectrum = sr.roman2int(spectrum)

    ## The atom 
    if atom in ['H', 'He']:
        pnatom = "{:}{:}r".format(atom, spectrum)
    else:
        pnatom = "{:}{:}".format(atom, spectrum)
    
    ## Transition Line
    try:
        wls = pn.LINE_LABEL_LIST[pnatom]
    except:
        wls = None
        if logger:
            logger.warning(
                f"Failed matching object {atom}/{spectrum} (aka {pnatom}) to PyNeb (see LINE_LABEL_LIST)")
            logger.warning(
                f"Will try to match a specific data file ...")

    if wls is not None:

        pnatomic = f'{atomic_number}A'
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
                    f'[ERROR] Failed matching object {atom}/{spectrum}/{atomic_number} to PyNeb (see LINE_LABEL_LIST)')
        raise RuntimeError(
            f'[ERROR] Failed matching object {atom}/{spectrum}/{atomic_number} to PyNeb (see LINE_LABEL_LIST)')
    else:
        return unlisted(atom, spectrum, atomic_number, logger)
