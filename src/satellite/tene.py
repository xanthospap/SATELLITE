import pyneb as pn
import numpy as np

def err2scalar(err_array):
    """ Define the fucntion to produce a single (i.e. scalar) error value, 
        when we have an array of error values (i.e. from Monte-Carlo 
        simulations.
    """
    if np.isnan(err_array).any():
        print('WARNING  array contains nan! {:}'.format(err_array))
    return np.std(err_array)


def computeTeNePairs_obsolete(density_diagnostics: list, tempterature_diagnostics: list, pnObs, pnErrObs, logger):
    print("--------------------------------------------------------computeTeNePairs-start");
    diags = pn.Diagnostics()
    # Register all diagnostics with PyNeb
    diags.addDiag(density_diagnostics+tempterature_diagnostics)
    # Retrieve the list of valid diagnostics recognized by PyNeb
    all_diags = diags.getDiagLabels()
    # Filter out only valid diagnostics
    valid_temp = [t for t in tempterature_diagnostics if t in all_diags]
    valid_dens = [d for d in density_diagnostics if d in all_diags]
    # Automatically generate valid (temperature, density) pairs
    valid_pairs = [(t, d) for t in valid_temp for d in valid_dens]
    print('Temperature Diagnostics: {:}'.format(valid_temp))
    print('Density     Diagnostics: {:}'.format(valid_dens))
    print('Valid       Diagnostics: {:}'.format(valid_pairs))
    # Iterate valid pairs and add results to dictionary
    tene_slit_dict = []
    for t, d in valid_pairs:
        try:
            st, sn = diags.getCrossTemDen(t, d, obs=pnObs)
            et, en = diags.getCrossTemDen(t, d, obs=pnErrObs)
            tene_slit_dict.append(
                {'tene_pair': (t, d), 'sT': st, 'sN': sn, 'eT': et, 'eN': en})
            # global_tene[slit_idx] = tene_slit_dict
        except:
            logger.info(
                f"Skipping Tem/Den pair {t} (Temp) ↔ {d} (Density)")
    print("--------------------------------------------------------computeTeNePairs-stop");
    return tene_slit_dict

def computeTeNePairs(density_diagnostics: list, tempterature_diagnostics: list, pnObs, pnErrObs, logger):
    print("--------------------------------------------------------computeTeNePairs-start");
    diags = pn.Diagnostics()
    tene_slit_dict = []
    for td in tempterature_diagnostics:
        for dd in density_diagnostics:
            try:
                st, sn = diags.getCrossTemDen(td, dd, obs=pnObs)
                et, en = diags.getCrossTemDen(td, dd, obs=pnErrObs)
                tene_slit_dict.append(
                    {'tene_pair': (td, dd), 'sT': st, 'sN': sn, 'eT': et, 'eN': en})
            except:
                logger.info(
                    f"Skipping Tem/Den pair {td} (Temp) ↔ {dd} (Density)")
    print("--------------------------------------------------------computeTeNePairs-stop");
    return tene_slit_dict

def printTempDens(dict_of_tene, fn, logger):
    num_slits = len(dict_of_tene)
    with open(fn, 'w') as fout:
        all_rows = []
        for idx, elist in dict_of_tene.items():
            for entry in elist:
                # each elist is in the following form:
                # {'tene_pair': (t,d), 'sT': st, 'sN': sn, 'eT': et, 'eN': en}
                if entry['tene_pair'] not in all_rows:
                    all_rows.append(entry['tene_pair'])
        print("{:45s}{:}".format('#', ''.join(["Slit {:5d}{:16s}{:26s}".format(
            d, ' ', ' ') for d in range(num_slits)])), file=fout)
        print("{:45s}{:11s}   {:11s} {:11s}   {:11s}".format(
            "Te/Ne Pair", "Tempature", "", "Density", ""), file=fout)
        print('-'*(45+26*2*num_slits), file=fout)
        for pair in all_rows:
            print("{:45s}".format('/'.join(pair)), end='', file=fout)
# value of TeNe pair (pair) for all slits ...
            for sj in range(num_slits):
                elist = dict_of_tene[sj]
                edict = next(
                    (item for item in elist if item["tene_pair"] == pair), None)
# TODO handle case where edict is None!!
                print("{:.5e} \u00B1 {:.5e} {:.5e} \u00B1 {:.5e} ".format(edict['sT'], err2scalar(
                    edict['eT']), edict['sN'], err2scalar(edict['eN'])), end='', file=fout)
# pair done !
            print('', file=fout)
    return fn
