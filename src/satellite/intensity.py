import sys
import pyneb as pn
import numpy as np

from satellite import astroflux
import satellite.roman as sr
import satellite.cfgio as sc
import satellite.nomenclature as sn


def makeIntensitiesDataFile(
    fitsd: list, reference_element: dict, value_keys: list, fn: str, factor=100e0
) -> str:
    """Write a relative line intensities data file to be used by PyNeb.

    For each of the entries in fitsd, compute relative intensity as:
        factor * value_key[0] / ref_value_key
    where ref_value_key is the value_key[0] of the reference element.
    Corresponding error values are also computed.

    Obviously, the reference_element must be included in fitsd.

    Parameters
    ----------
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

    Returns
    -------
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

    Example
    -------
    >>> # note that the following dictionary holds more key/value pairs than needed,
    >>> # but this is not a problem!
    >>> fitsd = [{'element': 'H', 'spectrum': 'i', 'atomic': 6563, 'ref_tene': '[SIII] 6312/9069 [ClIII] 5538/5518', 'fns': '/foo/data/HI_6563s.fit', 'fne': '/foo/data/HI_6563e.fit', 'sslit_sum': 908805000.0, 'eslit_sum': 3384860.0},
        {'element': 'H', 'spectrum': 'i', 'atomic': 4861, 'ref_tene': '[SIII] 6312/9069 [ClIII] 5538/5518', 'fns': '/foo/data/HI_4861s.fit', 'fne': '/foo/data/HI_4861e.fit', 'sslit_sum': 280699000.0, 'eslit_sum': 1079765.0},
        {'element': 'He', 'spectrum': 'i', 'atomic': 5876, 'ref_tene': '[SIII] 6312/9069 [ClIII] 5538/5518', 'fns': '/foo/data/HeI_5876s.fit', 'fne': '/foo/data/HeI_5876e.fit', 'sslit_sum': 41788900.0, 'eslit_sum': 180384.7},
        {'element': 'He', 'spectrum': 'i', 'atomic': 6678, 'ref_tene': '[SIII] 6312/9069 [ClIII] 5538/5518', 'fns': '/foo/data/HeI_6678s.fit', 'fne': '/foo/data/HeI_6678e.fit', 'sslit_sum': 11669580.0, 'eslit_sum': 74122.61},
        {'element': 'He', 'spectrum': 'ii', 'atomic': 5412, 'ref_tene': '[SIII] 6312/9069 [ClIII] 5538/5518', 'fns': '/foo/data/HeII_5412s.fit', 'fne': '/foo/data/HeII_5412e.fit', 'sslit_sum': 5692990.0, 'eslit_sum': 64965.5},
        {'element': 'N', 'spectrum': 'i', 'atomic': 5199, 'ref_tene': '[SIII] 6312/9069 [ClIII] 5538/5518', 'fns': '/foo/data/N1_5199s.fit', 'fne': '/foo/data/N1_5199e.fit', 'sslit_sum': 16205.0, 'eslit_sum': 87.736},
        {'element': 'N', 'spectrum': 'ii', 'atomic': 5755, 'ref_tene': '[SIII] 6312/9069 [ClIII] 5538/5518', 'fns': '/foo/data/N2_5755s.fit', 'fne': '/foo/data/N2_5755e.fit', 'sslit_sum': 622048.3, 'eslit_sum': 29662.768},
        {'element': 'N', 'spectrum': 'ii', 'atomic': 6548, 'ref_tene': '[SIII] 6312/9069 [ClIII] 5538/5518', 'fns': '/foo/data/N2_6548s.fit', 'fne': '/foo/data/N2_6548e.fit', 'sslit_sum': 6330353.0, 'eslit_sum': 35060.0},
        {'element': 'O', 'spectrum': 'i', 'atomic': 6300, 'ref_tene': '[SIII] 6312/9069 [ClIII] 5538/5518', 'fns': '/foo/data/O1_6300s.fit', 'fne': '/foo/data/O1_6300e.fit', 'sslit_sum': 42548.598, 'eslit_sum': 312.079},
        {'element': 'O', 'spectrum': 'ii', 'atomic': 7320, 'ref_tene': '[SIII] 6312/9069 [ClIII] 5538/5518', 'fns': '/foo/data/O2_7320s.fit', 'fne': '/foo/data/O2_7320e.fit', 'sslit_sum': 3711130.0, 'eslit_sum': 40674.1},
        {'element': 'O', 'spectrum': 'iii', 'atomic': 5007, 'ref_tene': '[SIII] 6312/9069 [ClIII] 5538/5518', 'fns': '/foo/data/O3_5007s.fit', 'fne': '/foo/data/O3_5007e.fit', 'sslit_sum': 3239288800.0, 'eslit_sum': 13766130.0},
        {'element': 'S', 'spectrum': 'ii', 'atomic': 6716, 'ref_tene': '[SIII] 6312/9069 [ClIII] 5538/5518', 'fns': '/foo/data/S2_6716s.fit', 'fne': '/foo/data/S2_6716e.fit', 'sslit_sum': 1472219.0, 'eslit_sum': 29649.402},
        {'element': 'S', 'spectrum': 'iii', 'atomic': 9069, 'ref_tene': '[SIII] 6312/9069 [ClIII] 5538/5518', 'fns': '/foo/data/S3_9069s.fit', 'fne': '/foo/data/S3_9069e.fit', 'sslit_sum': 68055600.0, 'eslit_sum': 301202.0},
        {'element': 'Ar', 'spectrum': 'iii', 'atomic': 7136, 'ref_tene': '[SIII] 6312/9069 [ClIII] 5538/5518', 'fns': '/foo/data/Ar3_7136s.fit', 'fne': '/foo/data/Ar3_7136e.fit', 'sslit_sum': 47880200.0, 'eslit_sum': 177312.1},
        {'element': 'Cl', 'spectrum': 'iii', 'atomic': 5538, 'ref_tene': '[SIII] 6312/9069 [ClIII] 5538/5518', 'fns': '/foo/data/Cl3_5538s.fit', 'fne': '/foo/data/Cl3_5538e.fit', 'sslit_sum': 1412538.0, 'eslit_sum': 38494.0}
    ]
    >>> reference_element = {'element': 'H', 'spectrum': 'i', 'atomic': 4861}
    >>> value_keys = ['sslit_sum', 'eslit_sum']
    >>> # write line intensities to a file named 'foobar.dat'
    >>> makeIntensitiesDataFile(fitsd, reference_elements, value_keys, 'foobar.dat')
    """
    ref_index = sc.indexOf(
        reference_element["element"],
        reference_element["spectrum"],
        reference_element["atomic"],
        fitsd,
    )
    ref_sval = fitsd[ref_index][value_keys[0]]
    ref_eval = fitsd[ref_index][value_keys[1]]
    with open(fn, "w") as fout:
        print("LINE test err", file=fout)
        for obj in fitsd:
            pnlabel = sn.objectIntensityPyNebCode(
                obj["element"], obj["spectrum"], obj["atomic"]
            )
            print(
                "{:} {:+9e} {:+9e}".format(
                    pnlabel,
                    factor * obj[value_keys[0]] / ref_sval,
                    astroflux.fluxError(
                        obj[value_keys[0]], obj[value_keys[1]], ref_sval, ref_eval
                    ),
                ),
                file=fout,
            )
    return fn


def computeIntensities(
    fitsd: dict, pnObs, pnErrObs, pnRC, reference_element: dict, logger
) -> list:
    """Compute line intensities based on PyNeb Observation.

    for line L and reference line R:
    intensity(L) = pnObs.getIntens(L)
    intensity_error(L) = {pnObs.getError(L)**2 +
            pnObs.getError(R)**2 +
            (pnRc.getCorr(LR)/pnRC.getErrCorr())**2}^{1/2}

    Parameters
    ----------
    fitsd: A list of dictionaries. Each dictionary must at least have:
        {'element': ..., 'spectrum': ..., 'atomic': ...,  ...}
    pnObs: pyneb.Observation
        A pyneb Observation instance, for which we will call:
        pnObs.getIntens and pnObs.getError to compute intensities and
        intensity error, for every line (i.e. combination of element, spectrum
        and atomic for every entry of fitsd).
    pnErrObs: pyneb.Observation
        This pyneb Observation is only used for calling:
        pnErrObs.extinction.E_BV
    pnRC: pyneb.RedCor
        A pyneb.RedCorr instance, created using something like:
        RC = pn.RedCorr(E_BV=pnObs.extinction.E_BV[0], law=some_ext_law)
    reference_element: The reference element for computations (including
        element name, spectrum and atomic number. Example:
        {'element': 'H', 'spectrum': 'i', 'atomic': 4861}
    logger: logging
        This parameter can be set to an instance of logging to
        write warnings on a specified stream. Default value is
        None.

    Return
    ------
    list:
        A list where each entry is a dictionary, corresponding to one entry in
        the fitsd input dictionary.
        For every fitsd entry, we have a corresponding dictionary with keys:
        * 'element_pn'    -> name of line according to pyneb
        * 'element'       -> name of line according to convention
        * 'intensity'     -> intensity value
        * 'intensity_err' -> intensity error
        Example output:
    [
        {'element_pn': 'H1r_6563A', 'element': 'H1_6563', 'intensity': 287.7865254633664, 'intensity_err': 0.007634893902546808},
        {'element_pn': 'H1r_4861A', 'element': 'H1_4861', 'intensity': 100.0, 'intensity_err': 0.00769340097526913},
        {'element_pn': 'He1r_5876A', 'element': 'He1_5876', 'intensity': 13.756839258392082, 'intensity_err': 0.007939511805793675},
        {'element_pn': 'He1r_6678A', 'element': 'He1_6678', 'intensity': 3.6734657420080756, 'intensity_err': 0.009206951765251258},
        {'element_pn': 'He2r_5411A', 'element': 'He2_5412', 'intensity': 1.9310008831029002, 'intensity_err': 0.013214297830023747},
        {'element_pn': 'N1_5200A', 'element': 'N1_5199', 'intensity': 0.005599837374762486, 'intensity_err': 0.008585203219191666},
        {'element_pn': 'N2_5755A', 'element': 'N2_5755', 'intensity': 0.2062965614133243, 'intensity_err': 0.04814892901497873},
        {'element_pn': 'O3_4959A', 'element': 'O3_4959', 'intensity': 381.27855164329196, 'intensity_err': 0.007902639517996805},
        {'element_pn': 'O3_5007A', 'element': 'O3_5007', 'intensity': 1138.896872086618, 'intensity_err': 0.00790265027212728},
        {'element_pn': 'S3_9069A', 'element': 'S3_9069', 'intensity': 19.552643938487247, 'intensity_err': 0.008005853216039908},
        {'element_pn': 'Ar3_7136A', 'element': 'Ar3_7136', 'intensity': 14.74327561603046, 'intensity_err': 0.007625687045022606},
        {'element_pn': 'Cl3_5538A', 'element': 'Cl3_5538', 'intensity': 0.474866526139717, 'intensity_err': 0.028054418400823428}
    ]
    """
    ref_pnstr = sn.objectIntensityPyNebCode(
        reference_element["element"],
        reference_element["spectrum"],
        reference_element["atomic"],
        logger,
    )
    iref = float(pnObs.getIntens()[ref_pnstr])
    eref = float(pnObs.getError()[ref_pnstr])
    intensities_list = []
    for idx, fits in enumerate(fitsd):
        scor = pnRC.getCorr(fits["atomic"], reference_element["atomic"])
        ecor = pnRC.getErrCorr(
            fits["atomic"],
            np.std(pnErrObs.extinction.E_BV),
            reference_element["atomic"],
        )
        pnstr = sn.objectIntensityPyNebCode(
            fits["element"], fits["spectrum"], fits["atomic"], logger
        )
        iele = float(pnObs.getIntens()[pnstr])
        eele = float(pnObs.getError()[pnstr])
        err = np.sqrt(eele**2 + eref**2 + float(ecor / scor) ** 2)
        intensities_list.append(
            {
                "element_pn": pnstr,
                "element": "{:}{:}_{:}".format(
                    fits["element"], sr.roman2int(fits["spectrum"]), fits["atomic"]
                ),
                "intensity": iele,
                "intensity_err": err,
            }
        )
    return intensities_list


def printIntensities(dict_of_intensities, fn, logger):
    """Example: global_intensities[1] = {'intensities': [...], 'E_BV': rc.E_BV, 'cHbeta': rc.cHbeta}"""
    element_format = "element_pn"

    # given an (inner) dictionary of 'intensities' return the specific
    # dictionary of a line, if it exists, else None
    def inDictOf(intensitiesDict, line):
        for entry in intensitiesDict:
            if entry[element_format] == line:
                return entry
        return None

    # extract unique lines and store columns
    unique_lines = []
    columns = []
    for k, v in dict_of_intensities.items():
        aplist = []
        columns += [k]
        for entry in v["intensities"]:
            aplist += [entry["element_pn"]]
        unique_lines = list(set(unique_lines + aplist))

    # sort rows & columns
    unique_lines = sorted(unique_lines)
    columns = sorted(columns)

    # if columns are numeric values and start at 0, then add an 1 offset so they start from 1
    offset = ""
    try:
        [int(c) for c in columns]
        if columns[0] == 0:
            offset = 1
    except:
        pass

    with open(fn, "w") as fout:
        # write first line, i.e. column keys
        print("{:15s}".format("Slit Nr."), file=fout, end="")
        for col in columns:
            print("{:->6d}{:25s} ".format(col + offset, "-" * 25), file=fout, end="")
        print("", file=fout)

        # iterate for every line in unique_lines
        for line in unique_lines:
            print("{:15s}".format(line), file=fout, end="")
            for col in columns:
                entry = inDictOf(dict_of_intensities[col]["intensities"], line)
                if entry is not None:
                    print(
                        "{:15.9e} {:15.9e} ".format(
                            entry["intensity"], entry["intensity_err"]
                        ),
                        file=fout,
                        end="",
                    )
                else:
                    print("{:31s} ".format(" "), file=fout, end="")
            print("", file=fout)

        # cHbeta line
        print("{:15s}".format("cHbeta"), file=fout, end="")
        for col in columns:
            print(
                "{:15.9e} {:15.9e} ".format(
                    dict_of_intensities[col]["cHbeta"],
                    dict_of_intensities[col]["cHbetaError"],
                ),
                file=fout,
                end="",
            )
        # Total F(Hb) line
        print("\n{:15s}".format("F(Hb)"), file=fout, end="")
        for col in columns:
            print(
                "{:15.9e} {:15.9e} ".format(
                    dict_of_intensities[col]["FHb"],
                    dict_of_intensities[col]["FHb_error"],
                ),
                file=fout,
                end="",
            )
    return fn
