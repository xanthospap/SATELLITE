# define_Te_Ne_for_null_spaxels_script.py:
# Defines the Te and Ne classes of arrays with zero values.
# (C) Stavros Akras

import numpy as np
from satellite import fluxes as flx

def nullTeNe(param_mod_name, param_model_values):

    index_step_ang = param_mod_name.index("position_angle_step")
    index_star_ang = param_mod_name.index("minimum_position_angle")
    index_end_ang = param_mod_name.index("maximum_position_angle")
    start_angle = param_model_values[index_star_ang]
    end_angle = param_model_values[index_end_ang]
    step_angle = param_model_values[index_step_ang]

    numbersize = int((end_angle - start_angle) / step_angle) + 1
    if numbersize > 0:
        Te = flx.TeNe(numbersize)
        Ne = flx.TeNe(numbersize)
        return Te, Ne


def nullTeNe_specific_slits():

    number_specific_slits = 10
    Te = flx.TeNe(number_specific_slits)
    Ne = flx.TeNe(number_specific_slits)
    return Te, Ne
