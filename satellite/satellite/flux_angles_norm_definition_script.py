# flux_angles_norm_definition_script.py:
# Defines the flux_angles_norm and flux_angles_norm_error arrays
# in case the rotation analysis module is desactivated and returns two classes of arrays with zero values.
# (C) Stavros Akras


class flux_angles_norm:

    def __init__(self):
        self.Ha_6563 = []
        self.Hb_4861 = []
        self.Hg_4340 = []
        self.Hd_4101 = []
        self.HeIa_5876 = []
        self.HeIb_6678 = []
        self.HeIIa_4686 = []
        self.HeIIb_5412 = []
        self.NIIa_5755 = []
        self.NIIb_6548 = []
        self.NIIc_6584 = []
        self.NI_5199 = []
        self.OIIIa_4363 = []
        self.OIIIb_4959 = []
        self.OIIIc_5007 = []
        self.OIIa_3727 = []
        self.OIIb_3729 = []
        self.OIIc_7320 = []
        self.OIId_7330 = []
        self.OIa_5577 = []
        self.OIb_6300 = []
        self.OIc_6363 = []
        self.SIIa_6716 = []
        self.SIIb_6731 = []
        self.SIIIa_6312 = []
        self.SIIIb_9069 = []
        self.ClIIIa_5517 = []
        self.ClIIIb_5538 = []
        self.ArIII_7136 = []
        self.ArIVa_4712 = []
        self.ArIVb_4740 = []
        self.CI_8727 = []
        self.CII_6461 = []
        self.NeIIIa_3868 = []
        self.NeIIIb_3967 = []


def flux_ang_norm_definition(param_mod_name, param_model_values):

    flux_an = flux_angles_norm()
    flux_an_error = flux_angles_norm()

    index_step_ang = param_mod_name.index("position_angle_step")
    index_star_ang = param_mod_name.index("minimum_position_angle")
    index_end_ang = param_mod_name.index("maximum_position_angle")
    start_angle = param_model_values[index_star_ang]
    end_angle = param_model_values[index_end_ang]
    step_angle = param_model_values[index_step_ang]

    ang = []
    kkangle = 0
    for i in range(start_angle, end_angle, step_angle):
        ang.append(0)
        ang[kkangle] = i

        flux_an.Ha_6563.append(-99)
        flux_an.Hb_4861.append(-99)
        flux_an.Hg_4340.append(-99)
        flux_an.Hd_4101.append(-99)
        flux_an.HeIa_5876.append(-99)
        flux_an.HeIb_6678.append(-99)
        flux_an.HeIIa_4686.append(-99)
        flux_an.HeIIb_5412.append(-99)
        flux_an.NIIa_5755.append(-99)
        flux_an.NIIb_6548.append(-99)
        flux_an.NIIc_6584.append(-99)
        flux_an.NI_5199.append(-99)
        flux_an.OIIIa_4363.append(-99)
        flux_an.OIIIb_4959.append(-99)
        flux_an.OIIIc_5007.append(-99)
        flux_an.OIIa_3727.append(-99)
        flux_an.OIIb_3729.append(-99)
        flux_an.OIIc_7320.append(-99)
        flux_an.OIId_7330.append(-99)
        flux_an.OIa_5577.append(-99)
        flux_an.OIb_6300.append(-99)
        flux_an.OIc_6363.append(-99)
        flux_an.SIIa_6716.append(-99)
        flux_an.SIIb_6731.append(-99)
        flux_an.SIIIa_6312.append(-99)
        flux_an.SIIIb_9069.append(-99)
        flux_an.ClIIIa_5517.append(-99)
        flux_an.ClIIIb_5538.append(-99)
        flux_an.ArIII_7136.append(-99)
        flux_an.ArIVa_4712.append(-99)
        flux_an.ArIVb_4740.append(-99)
        flux_an.CI_8727.append(-99)
        flux_an.CII_6461.append(-99)
        flux_an.NeIIIa_3868.append(-99)
        flux_an.NeIIIb_3967.append(-99)

        flux_an_error.Ha_6563.append(-99)
        flux_an_error.Hb_4861.append(-99)
        flux_an_error.Hg_4340.append(-99)
        flux_an_error.Hd_4101.append(-99)
        flux_an_error.HeIa_5876.append(-99)
        flux_an_error.HeIb_6678.append(-99)
        flux_an_error.HeIIa_4686.append(-99)
        flux_an_error.HeIIb_5412.append(-99)
        flux_an_error.NIIa_5755.append(-99)
        flux_an_error.NIIb_6548.append(-99)
        flux_an_error.NIIc_6584.append(-99)
        flux_an_error.NI_5199.append(-99)
        flux_an_error.OIIIa_4363.append(-99)
        flux_an_error.OIIIb_4959.append(-99)
        flux_an_error.OIIIc_5007.append(-99)
        flux_an_error.OIIa_3727.append(-99)
        flux_an_error.OIIb_3729.append(-99)
        flux_an_error.OIIc_7320.append(-99)
        flux_an_error.OIId_7330.append(-99)
        flux_an_error.OIa_5577.append(-99)
        flux_an_error.OIb_6300.append(-99)
        flux_an_error.OIc_6363.append(-99)
        flux_an_error.SIIa_6716.append(-99)
        flux_an_error.SIIb_6731.append(-99)
        flux_an_error.SIIIa_6312.append(-99)
        flux_an_error.SIIIb_9069.append(-99)
        flux_an_error.ClIIIa_5517.append(-99)
        flux_an_error.ClIIIb_5538.append(-99)
        flux_an_error.ArIII_7136.append(-99)
        flux_an_error.ArIVa_4712.append(-99)
        flux_an_error.ArIVb_4740.append(-99)
        flux_an_error.CI_8727.append(-99)
        flux_an_error.CII_6461.append(-99)
        flux_an_error.NeIIIa_3868.append(-99)
        flux_an_error.NeIIIb_3967.append(-99)

    return flux_an, flux_an_error, ang


def flux_ang_norm_definition_specificslits(param_mod_name, param_model_values):

    flux_an = flux_angles_norm()
    flux_an_error = flux_angles_norm()

    for i in range(0, 10):

        flux_an.Ha_6563.append(-99)
        flux_an.Hb_4861.append(-99)
        flux_an.Hg_4340.append(-99)
        flux_an.Hd_4101.append(-99)
        flux_an.HeIa_5876.append(-99)
        flux_an.HeIb_6678.append(-99)
        flux_an.HeIIa_4686.append(-99)
        flux_an.HeIIb_5412.append(-99)
        flux_an.NIIa_5755.append(-99)
        flux_an.NIIb_6548.append(-99)
        flux_an.NIIc_6584.append(-99)
        flux_an.NI_5199.append(-99)
        flux_an.OIIIa_4363.append(-99)
        flux_an.OIIIb_4959.append(-99)
        flux_an.OIIIc_5007.append(-99)
        flux_an.OIIa_3727.append(-99)
        flux_an.OIIb_3729.append(-99)
        flux_an.OIIc_7320.append(-99)
        flux_an.OIId_7330.append(-99)
        flux_an.OIa_5577.append(-99)
        flux_an.OIb_6300.append(-99)
        flux_an.OIc_6363.append(-99)
        flux_an.SIIa_6716.append(-99)
        flux_an.SIIb_6731.append(-99)
        flux_an.SIIIa_6312.append(-99)
        flux_an.SIIIb_9069.append(-99)
        flux_an.ClIIIa_5517.append(-99)
        flux_an.ClIIIb_5538.append(-99)
        flux_an.ArIII_7136.append(-99)
        flux_an.ArIVa_4712.append(-99)
        flux_an.ArIVb_4740.append(-99)
        flux_an.CI_8727.append(-99)
        flux_an.CII_6461.append(-99)
        flux_an.NeIIIa_3868.append(-99)
        flux_an.NeIIIb_3967.append(-99)

        flux_an_error.Ha_6563.append(-99)
        flux_an_error.Hb_4861.append(-99)
        flux_an_error.Hg_4340.append(-99)
        flux_an_error.Hd_4101.append(-99)
        flux_an_error.HeIa_5876.append(-99)
        flux_an_error.HeIb_6678.append(-99)
        flux_an_error.HeIIa_4686.append(-99)
        flux_an_error.HeIIb_5412.append(-99)
        flux_an_error.NIIa_5755.append(-99)
        flux_an_error.NIIb_6548.append(-99)
        flux_an_error.NIIc_6584.append(-99)
        flux_an_error.NI_5199.append(-99)
        flux_an_error.OIIIa_4363.append(-99)
        flux_an_error.OIIIb_4959.append(-99)
        flux_an_error.OIIIc_5007.append(-99)
        flux_an_error.OIIa_3727.append(-99)
        flux_an_error.OIIb_3729.append(-99)
        flux_an_error.OIIc_7320.append(-99)
        flux_an_error.OIId_7330.append(-99)
        flux_an_error.OIa_5577.append(-99)
        flux_an_error.OIb_6300.append(-99)
        flux_an_error.OIc_6363.append(-99)
        flux_an_error.SIIa_6716.append(-99)
        flux_an_error.SIIb_6731.append(-99)
        flux_an_error.SIIIa_6312.append(-99)
        flux_an_error.SIIIb_9069.append(-99)
        flux_an_error.ClIIIa_5517.append(-99)
        flux_an_error.ClIIIb_5538.append(-99)
        flux_an_error.ArIII_7136.append(-99)
        flux_an_error.ArIVa_4712.append(-99)
        flux_an_error.ArIVb_4740.append(-99)
        flux_an_error.CI_8727.append(-99)
        flux_an_error.CII_6461.append(-99)
        flux_an_error.NeIIIa_3868.append(-99)
        flux_an_error.NeIIIb_3967.append(-99)

    return flux_an, flux_an_error
