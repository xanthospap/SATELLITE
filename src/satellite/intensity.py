import sys
import pyneb as pn
import numpy as np

# from satellite import cfgio
from satellite import astroflux
import satellite.roman as sr
import satellite.cfgio as sc


def makeIntensitiesDataFile(fitsd: list, reference_element: dict, value_keys: list, fn: str, factor=100e0):
    """ Write a relative intensities(?) data file to be used by PyNeb.

        For each of the entries in fitsd, compute relative intensity as:
            factor * value_key[0] / ref_value_key 
        where ref_value_key is the value_key[0] of the reference element.
        Corresponding error values are also computed.

        Obviously, the reference_element must be included in fitsd.

        Parameters
        -----------
        fitsd: A list of dictionaries. Each dictionary must at least have:
            {'element': ..., 'spectrum': ..., 'atomic': ..., value_keys[0]: ..., value_keys[1]: ...}
        reference_element: The reference element for computations (including
            element name, spectrum and atomic number. Example:
            {'element': 'H', 'spectrum': 'i', 'atomic': 4861}
        value_keys: A list containing the two keys (within fitsd) for the 
            computations. First key is for the intensity(?) value, second is
            for the error. Example: ['sslit_sum', 'eslit_sum']
            So, the values of 'sslit_sum' and 'eslit_sum' will be extracted 
            from each entry of fitsd.
        fn: The filename of the file to be created/written. Example: 'test.dat'

        Returns:
        ---------
        str: The filename of the file created/written. The contents of the 
             file are like:
                LINE test err
                H1r_6563A +3.237649e+02 +1.733553e+00
                H1r_4861A +1.000000e+02 +5.440056e-01
                He1r_5876A +1.488744e+01 +8.607707e-02
                He1r_6678A +4.157329e+00 +3.087141e-02
                He2r_5411A +2.028148e+00 +2.442375e-02
                N1_5200A +5.773088e-03 +3.834214e-05
                N2_5755A +2.216069e-01 +1.060179e-02
                N2_6548A +2.255210e+00 +1.520736e-02
                N2_6584A +6.765631e+00 +4.562210e-02
                O1_6300A +1.515809e-02 +1.255417e-04
                O2_7320A +1.322103e+00 +1.535686e-02
    """
    ref_index = sc.indexOf(
        reference_element['element'], reference_element['spectrum'], reference_element['atomic'], fitsd)
    ref_sval = fitsd[ref_index][value_keys[0]]
    ref_eval = fitsd[ref_index][value_keys[1]]
    with open(fn, 'w') as fout:
        print('LINE test err', file=fout)
        for obj in fitsd:
            pnlabel = sc.objectIntensityPyNebCode(
                obj['element'], obj['spectrum'], obj['atomic'])
            print('{:} {:+9e} {:+9e}'.format(pnlabel, factor*obj[value_keys[0]]/ref_sval, astroflux.fluxError(
                obj[value_keys[0]], obj[value_keys[1]], ref_sval, ref_eval)), file=fout)
    return fn


def computeIntensities(fitsd: dict, pnObs, pnErrObs, pnRC, reference_element: dict, logger):
    """
    """
    ref_pnstr = sc.objectIntensityPyNebCode(
        reference_element['element'], reference_element['spectrum'], reference_element['atomic'], logger)
    iref = float(pnObs.getIntens()[ref_pnstr])
    eref = float(pnObs.getError()[ref_pnstr])
    intensities_list = []
    for idx, fits in enumerate(fitsd):
        scor = pnRC.getCorr(fits['atomic'], reference_element['atomic'])
        ecor = pnRC.getErrCorr(fits['atomic'], np.std(
            pnErrObs.extinction.E_BV), reference_element['atomic'])
        pnstr = sc.objectIntensityPyNebCode(
            fits['element'], fits['spectrum'], fits['atomic'], logger)
        iele = float(pnObs.getIntens()[pnstr])
        eele = float(pnObs.getError()[pnstr])
        err = np.sqrt(eele**2 + eref**2 + float(ecor/scor)**2)
        intensities_list.append({'element_pn': pnstr, 'element': '{:}{:}_{:}'.format(
            fits['element'], sr.roman2int(fits['spectrum']), fits['atomic']), 'intensity': iele, 'intensity_err': err})
    return intensities_list


def printIntensities(dict_of_intensities_list, fn, logger):
    """
    """
    element_format = 'element_pn'

    def searchElem(subdict, elem):
        for entry in subdict:
            if entry[element_format] == elem:
                return entry
        return None
    # extract all, unique elements; this is the sequence they will be
    # written at
    all_elements = list(set(d[element_format]
                        for sublist in dict_of_intensities_list.values() for d in sublist))
    num_slits = len(dict_of_intensities_list)
    with open(fn, 'w') as fout:
        print("{:15s}{:}".format('#', ''.join(["Slit {:5d}{:25s}".format(
            d, ' ') for d in range(num_slits)])), file=fout)
        print('{:}{:}'.format('#', '-'*(15-1+35*num_slits)), file=fout)
        for elem in all_elements:
            print('{:<15s}'.format(elem), end='', file=fout)
            for k, v in dict_of_intensities_list.items():
                elem_dict_slit = searchElem(v, elem)
                if elem_dict_slit is not None:
                    print('{:+.9e} \u00B1 {:.9e} '.format(
                        elem_dict_slit['intensity'], elem_dict_slit['intensity_err']), end='', file=fout)
                else:
                    print('{>14} \u00B1 {:>14} '.format(
                        'nan', 'nan'), end='', file=fout)
            print('', file=fout)
    return fn
