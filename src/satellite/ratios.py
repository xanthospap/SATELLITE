import numpy as np


def partialResolve(pstr: str) -> list:
    """

    Example
    -------
        >>> partialResolve("He2_5412") -> [1.0, 'He2_5412']
        >>> partialResolve("N2_6583") -> [1.0, 'N2_6583']
        >>> partialResolve("H1_6563") -> [1.0, 'H1_6563']
        >>> partialResolve("N2_6548+N2_6583") -> [1.0, 'N2_6548', 1.0, 'N2_6583']
        >>> partialResolve("H1_6563") -> [1.0, 'H1_6563']
        >>> partialResolve("N2_6548+N2_6583") -> [1.0, 'N2_6548', 1.0, 'N2_6583']
        >>> partialResolve("N2_5755") -> [1.0, 'N2_5755']
        >>> partialResolve("N2_6548+N2_6583") -> [1.0, 'N2_6548', 1.0, 'N2_6583']
        >>> partialResolve("O3_4959+O3_5007") -> [1.0, 'O3_4959', 1.0, 'O3_5007']
        >>> partialResolve("N1_5199") -> [1.0, 'N1_5199']
        >>> partialResolve("H1_4861") -> [1.0, 'H1_4861']
        >>> partialResolve("S2_6716+S2_6731") -> [1.0, 'S2_6716', 1.0, 'S2_6731']
        >>> partialResolve("H1_6563") -> [1.0, 'H1_6563']
        >>> partialResolve("S2_6716") -> [1.0, 'S2_6716']
        >>> partialResolve("S2_6731") -> [1.0, 'S2_6731']
        >>> partialResolve("S2_6716+S2_6731") -> [1.0, 'S2_6716', 1.0, 'S2_6731']
    """
    lstpl = []
    lst = []
    for i in pstr.split("+"):
        lstpl += [1e0, i.strip()]
    for idx, s in enumerate(lstpl):
        try:
            1 + s
            lst += [s]
        except:
            if "-" in s:
                lst += [s.split("-")[0]]
                for ss in s.split("-")[1:]:
                    lst += [-1, ss]
            else:
                lst += [s]
    return lst


def resolveRatioStr(rstr: str):
    """

    Example
    -------
    >>> resolveRatioStr("He1_5876/H1_6563") -> [1.0, 'He1_5876'], [1.0, 'H1_6563']
    >>> resolveRatioStr("He1_6678/H1_6563") -> [1.0, 'He1_6678'], [1.0, 'H1_6563']
    >>> resolveRatioStr("He2_4686/H1_4861") -> [1.0, 'He2_4686'], [1.0, 'H1_4861']
    >>> resolveRatioStr("He2_5412/H1_4861") -> [1.0, 'He2_5412'], [1.0, 'H1_4861']
    >>> resolveRatioStr("He1_5876/He2_4686") -> [1.0, 'He1_5876'], [1.0, 'He2_4686']
    >>> resolveRatioStr("He1_5876/He2_5412") -> [1.0, 'He1_5876'], [1.0, 'He2_5412']
    >>> resolveRatioStr("N2_6583/H1_6563") -> [1.0, 'N2_6583'], [1.0, 'H1_6563']
    >>> resolveRatioStr("N2_6548+N2_6583)/H1_6563") -> [1.0, 'N2_6548', 1.0, 'N2_6583'], [1.0, 'H1_6563']
    >>> resolveRatioStr("(N2_6548+N2_6583)/N2_5755") -> [1.0, 'N2_6548', 1.0, 'N2_6583'], [1.0, 'N2_5755']
    >>> resolveRatioStr("(N2_6548+N2_6583)/(O3_4959+O3_5007)") -> [1.0, 'N2_6548', 1.0, 'N2_6583'], [1.0, 'O3_4959', 1.0, 'O3_5007']
    >>> resolveRatioStr("N1_5199/H1_4861") -> [1.0, 'N1_5199'], [1.0, 'H1_4861']
    >>> resolveRatioStr("(S2_6716+S2_6731)/H1_6563") -> [1.0, 'S2_6716', 1.0, 'S2_6731'], [1.0, 'H1_6563']
    >>> resolveRatioStr("S2_6716/S2_6731") -> [1.0, 'S2_6716'], [1.0, 'S2_6731']
    >>> resolveRatioStr("(S2_6716+S2_6731)/(S3_6312+S3_9069)") -> [1.0, 'S2_6716', 1.0, 'S2_6731'], [1.0, 'S3_6312', 1.0, 'S3_9069']
    >>> resolveRatioStr("(O1_6300+O1_6363)/H1_6563") -> [1.0, 'O1_6300', 1.0, 'O1_6363'], [1.0, 'H1_6563']
    """
    [nom, denom] = rstr.split("/")
    nom = nom.lstrip("(").rstrip(")")
    denom = denom.lstrip("(").rstrip(")")
    return partialResolve(nom), partialResolve(denom)


def computeRatio(ratio: str, intensity_list: list, logger=None):
    def getIntensity(element):
        for obj in intensity_list:
            if obj["element"] == element:
                return obj["intensity"], obj["intensity_err"]
        if logger:
            logger.error(
                "Cannot find element {:} for ratio {:} in intensities list!".format(
                    element, ratio
                )
            )
        raise RuntimeError(
            "[ERROR] Cannot find element {:} for ratio {:} in intensities list!".format(
                element, ratio
            )
        )
        return None, None

    ar, par = resolveRatioStr(ratio)
    var = 0e0
    vpar = 0e0
    par1 = 0e0
    par2 = 0e0
    for idx in range(1, len(ar), 2):
        # val, err = getIntensity(sc.satellite_str2pyneb_str(ar[idx]))
        val, err = getIntensity(ar[idx])
        var += ar[idx - 1] * val
        par1 += ar[idx - 1] * (err / val * np.log(10))
    for idx in range(1, len(par), 2):
        # val, err = getIntensity(sc.satellite_str2pyneb_str(par[idx]))
        val, err = getIntensity(par[idx])
        vpar += par[idx - 1] * val
        par2 += par[idx - 1] * (err / val * np.log(10))
    return var / vpar, np.sqrt(par1**2 + par2**2)


def printRatios(dict_of_ratios, fn, logger) -> str:
    # extract unique lines and store columns
    unique_ratios = []
    columns = []
    for k, v in dict_of_ratios.items():
        columns += [k]
        unique_ratios = list(set(unique_ratios + [ky for ky in v]))

    # sort rows & columns
    unique_ratios = sorted(unique_ratios)
    columns = sorted(columns)

    # if columns are numeric values and start at 0, then add an 1 offset so they start from 1
    offset = ""
    try:
        [int(c) for c in columns]
        if columns[0] == 0:
            offset = 1
    except:
        pass

    def inDictOf(idxDct, ratio):
        return idxDct[ratio] if ratio in idxDct else None

    with open(fn, "w") as fout:
        # write first line, i.e. column keys
        print("{:40s}".format("Slit Nr."), file=fout, end="")
        for col in columns:
            print("{:->6d}{:25s} ".format(col + offset, "-" * 25), file=fout, end="")
        print("", file=fout)

        # iterate for every line in unique_lines
        for ratio in unique_ratios:
            print("{:40s}".format(ratio), file=fout, end="")
            for col in columns:
                entry = inDictOf(dict_of_ratios[col], ratio)
                if entry is not None:
                    print(
                        "{:15.7e} {:15.7e} ".format(
                            np.log10(entry[0]), np.log10(entry[1])
                        ),
                        file=fout,
                        end="",
                    )
                else:
                    print("{:31s} ".format(" "), file=fout, end="")
            print("", file=fout)

    return fn
