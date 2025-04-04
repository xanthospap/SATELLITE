import pyneb as pn
import numpy as np


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


def computeTeNePairs(
    density_diagnostics: list, tempterature_diagnostics: list, pnObs, pnErrObs, logger
):
    def filterPairs(densityd, temperatured, pobs):
        user = [(td, dd) for td in temperatured for dd in densityd]
        diags = pn.Diagnostics()
        # construct all possible diagnostics from the give observation set
        diags.addDiagsFromObs(pobs)
        validLines = diags.getDiagLabels()
        # filter user list based on observation set
        return [(d[0], d[1]) for d in user if d[0] in validLines and d[1] in validLines]

    tene_slit_dict = []
    diags = pn.Diagnostics()
    for pair in filterPairs(density_diagnostics, tempterature_diagnostics, pnObs):
        st, sn = diags.getCrossTemDen(pair[0], pair[1], obs=pnObs)
        et, en = diags.getCrossTemDen(pair[0], pair[1], obs=pnErrObs)
        tene_slit_dict.append(
            {"tene_pair": (pair[0], pair[1]), "sT": st, "sN": sn, "eT": et, "eN": en}
        )
    return tene_slit_dict


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
