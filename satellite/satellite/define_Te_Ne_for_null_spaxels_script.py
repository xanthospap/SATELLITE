# define_Te_Ne_for_null_spaxels_script.py:
# Defines the Te and Ne classes of arrays with zero values.
# (C) Stavros Akras

import numpy as np


class TeNe:

    def __init__(self, numbersize):
        self.NIISII = np.zeros(numbersize)
        self.OISII = np.zeros(numbersize)
        self.OIISII = np.zeros(numbersize)
        self.OIIISII = np.zeros(numbersize)
        self.SIIISII = np.zeros(numbersize)
        self.OIIOII = np.zeros(numbersize)
        self.NIIOII = np.zeros(numbersize)
        self.OIOII = np.zeros(numbersize)
        self.OIIIClIII = np.zeros(numbersize)
        self.SIIIClIII = np.zeros(numbersize)
        self.OIIIArVI = np.zeros(numbersize)
        self.SIIIArVI = np.zeros(numbersize)
        self.NIIClIII = np.zeros(numbersize)


def nullTeNe(param_mod_name, param_model_values):

    index_step_ang = param_mod_name.index("position_angle_step")
    index_star_ang = param_mod_name.index("minimum_position_angle")
    index_end_ang = param_mod_name.index("maximum_position_angle")
    start_angle = param_model_values[index_star_ang]
    end_angle = param_model_values[index_end_ang]
    step_angle = param_model_values[index_step_ang]

    numbersize = int((end_angle - start_angle) / step_angle) + 1
    if numbersize > 0:
        Te = TeNe(numbersize)
        Ne = TeNe(numbersize)
        return Te, Ne


def nullTeNe_specific_slits():

    number_specific_slits = 10
    Te = TeNe(number_specific_slits)
    Ne = TeNe(number_specific_slits)
    return Te, Ne
