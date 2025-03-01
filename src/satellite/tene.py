import pyneb as pn
import numpy as np

def computeTeNePairs(density_diagnostics: list, tempterature_diagnostics: list, pnObs, pnErrObs, logger):
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
    return tene_slit_dict
