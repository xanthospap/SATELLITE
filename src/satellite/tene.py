import pyneb as pn
import numpy as np

MIN_VALID_PERCENTAGE = 70.0  # e.g. require at least 70% non-NaN values


def inspectMcErr(err_array, min_percentage=MIN_VALID_PERCENTAGE, logger=None):
    """
    Convert an array of MC error realisations into a single scalar error.

    - If the fraction of non-NaN values is below `min_percentage`,
      we treat the diagnostic as failed and return (None, False).
    - Otherwise we compute the std dev over the non-NaN values and
      return (sigma, True).
    """
    arr = np.asarray(err_array, dtype=float)

    if arr.size == 0:
        logger.warning("inspectMcErr: empty error array; failing diagnostic")
        return None, False

    nan_mask = np.isnan(arr)
    n_total = arr.size
    n_nan = nan_mask.sum()
    n_valid = n_total - n_nan

    valid_percent = 100.0 * n_valid / n_total
    nan_percent = 100.0 - valid_percent

    if n_nan > 0:
        logger.debug(
            f"inspectMcErr: {n_nan}/{n_total} NaN values "
            f"({nan_percent:.1f}% NaN, {valid_percent:.1f}% valid)"
        )

    # Too many NaNs → fail
    if valid_percent < min_percentage:
        logger.warning(
            "inspectMcErr: too many NaN MC realisations "
            f"({nan_percent:.1f}% > allowed {100.0 - min_percentage:.1f}%). "
            "Marking diagnostic as failed."
        )
        return None, False

    # Enough valid values → pass
    # valid_vals = arr[~nan_mask]
    # sigma = float(np.std(valid_vals, ddof=1))  # ddof=1 → sample std dev

    return arr[~nan_mask], True


def computeTeNePairs(
    density_diagnostics: list,
    tempterature_diagnostics: list,
    pnObs,
    pnErrObs,
    min_percentage,
    logger,
):
    def filterPairs(densityd, temperatured, pobs):
        user = [(td, dd) for td in temperatured for dd in densityd]
        diags = pn.Diagnostics()
        # construct all possible diagnostics from the given observation set
        diags.addDiagsFromObs(pobs)
        validLines = diags.getDiagLabels()
        # filter user list based on observation set
        return [(d[0], d[1]) for d in user if d[0] in validLines and d[1] in validLines]

    tene_slit_dict = []
    diags = pn.Diagnostics()
    for pair in filterPairs(density_diagnostics, tempterature_diagnostics, pnObs):
        st, sn = diags.getCrossTemDen(pair[0], pair[1], obs=pnObs)
        et, en = diags.getCrossTemDen(pair[0], pair[1], obs=pnErrObs)
        et, okT = inspectMcErr(et, min_percentage, logger)
        en, okN = inspectMcErr(en, min_percentage, logger)
        tene_slit_dict.append(
            {
                "tene_pair": (pair[0], pair[1]),
                "sT": st,
                "sN": sn,
                "eT": et,
                "eN": en,
            }
        )

    return tene_slit_dict, not (okT and okN)


def err2scalar(err_array, logger):
    """
    Define the fucntion to produce a single (i.e. scalar) error value,
    when we have an array of error values (i.e. from Monte-Carlo
    simulations.
    """
    if np.isnan(err_array).any():
        logger.warning(
            "WARNING  array contains nan! {:} for computing error in temperature/density diagnostics".format(
                err_array
            )
        )
    return np.std(err_array[~np.isnan(err_array)])


def computeTeNePairs_obsolete(
    density_diagnostics: list, tempterature_diagnostics: list, pnObs, pnErrObs, logger
):
    def filterPairs(densityd, temperatured, pobs):
        user = [(td, dd) for td in temperatured for dd in densityd]
        diags = pn.Diagnostics()
        # construct all possible diagnostics from the given observation set
        diags.addDiagsFromObs(pobs)
        validLines = diags.getDiagLabels()
        # filter user list based on observation set
        return [(d[0], d[1]) for d in user if d[0] in validLines and d[1] in validLines]

    includes_nan = False
    tene_slit_dict = []
    diags = pn.Diagnostics()
    for pair in filterPairs(density_diagnostics, tempterature_diagnostics, pnObs):
        st, sn = diags.getCrossTemDen(pair[0], pair[1], obs=pnObs)
        et, en = diags.getCrossTemDen(pair[0], pair[1], obs=pnErrObs)
        if np.isnan(np.array([et, en])).any():
            logger.warning(
                f"Nan value(s) encountered while computing Te/Ne diagnostics!"
            )
            includes_nan = True
        tene_slit_dict.append(
            {
                "tene_pair": (pair[0], pair[1]),
                "sT": st,
                "sN": sn,
                "eT": et,
                "eN": en,
            }
        )

    return tene_slit_dict, includes_nan


def printDiagnostics(dict_of_diagnostics, fn, logger):
    def pair2str(diag):
        return "{:}/{:}".format(diag[0], diag[1])

    def inDictOf(diags_list, diag):
        for d in diags_list:
            if d["tene_pair"] == diag:
                return d
        return None

    # extract unique diagnostics and store columns
    unique_diags = []
    columns = []
    for k, v in dict_of_diagnostics.items():
        columns += [k]
        unique_diags = list(set(unique_diags + [x["tene_pair"] for x in v]))

    # sort rows & columns
    unique_diags = sorted(unique_diags)
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

        # write first line
        print("{:45s}".format(" "), file=fout, end="")
        for col in columns:
            print("{:26s} {:26s} ".format("Temperature", "Density"), file=fout, end="")
        print("", file=fout)

        # write second line, i.e. column keys
        print("{:45s}".format("Slit Nr."), file=fout, end="")
        for col in columns:
            print("{:->6d}{:47s} ".format(col + offset, "-" * 47), file=fout, end="")
        print("", file=fout)

        # iterate for every diagnostic in unique_diags
        for diag in unique_diags:
            print("{:45s}".format(pair2str(diag)), file=fout, end="")
            for col in columns:
                entry = inDictOf(dict_of_diagnostics[col], diag)
                if entry is not None:
                    print(
                        "{:12.6e} {:12.6e} / {:12.6e} {:12.6e} ".format(
                            entry["sT"],
                            err2scalar(entry["eT"], logger),
                            entry["sN"],
                            err2scalar(entry["eN"], logger),
                        ),
                        file=fout,
                        end="",
                    )
                else:
                    print("{:53s} ".format(" "), file=fout, end="")
            print("", file=fout)
