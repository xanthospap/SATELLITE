# TeNe_angles_script.py:
# Computes c(Hb), Te and Ne parameters for the 10 pseudo-slits in the specific slits analysis module using the PyNeb package (Luridiana et al. 2015)
# as well as the ionic/elemental abundances and ICFs.
# The results are also provided in an ascii file (PyNeb_output_total_abund_ICFs_specific_slits.txt
# and PyNeb_output_ionic_abund_specific_slits.txt and PyNeb_output_Te_and_Ne_specific_slits.txt)
# (C) Stavros Akras

from __future__ import print_function

import numpy as np
import pyneb as pn
import os, fnmatch
from numpy import inf

from satellite import ionicabundances_specific_slits_script as iasss
from satellite import element_abundances_ICFs_specific_slits_script as eaIsss
from satellite import norm_flux_error_script as nfes
from satellite import fluxes as flx

dct = {"H1r_6563A": 'Ha_6563',
       "H1r_4861A": 'Hb_4861',
       "H1r_4341A":
       'Hg_4340',
       "H1r_4102A":
       'Hd_4101',
       "He2r_4686A":
       'HeIIa_4686',
       "He2r_5411A":
       'HeIIb_5412',
       "He1r_5876A":
       'HeIa_5876',
       "He1r_6678A":
       'HeIb_6678',
       "S2_6716A":
       'SIIa_6716',
       "S2_6731A":
       'SIIb_6731',
       "S3_6312A":
       'SIIIa_6312',
       "S3_9069A":
       'SIIIb_9069',
       "N1_5200A":
       'NI_5199',
       "N2_5755A":
       'NIIa_5755',
       "N2_6548A":
       'NIIb_6548',
       "N2_6584A":
       'NIIc_6584',
       "O1_5577A":
       'OIa_5577',
       "O1_6300A":
       'OIb_6300',
       "O1_6364A":
       'OIc_6363',
       "O2_3726A":
       'OIIa_3727',
       "O2_3729A":
       'OIIb_3729',
       "O2_7319A+":
       'OIIc_7320',
       "O2_7330A+":
       'OIId_7330',
       "O3_4363A":
       'OIIIa_4363',
       "O3_4959A":
       'OIIIb_4959',
       "O3_5007A":
       'OIIIc_5007',
       "Cl3_5518A":
       'ClIIIa_5517',
       "Cl3_5538A":
       'ClIIIb_5538',
       "Ne3_3869A":
       'NeIIIa_3868',
       "Ne3_3968A":
       'NeIIIb_3967',
       "Ar3_7136A":
       'ArIII_7136',
       "Ar4_4711A":
       'ArIVa_4712',
       "Ar4_4740A":
       'ArIVb_4740',
       "C1_8728A":
       'CI_8727'
       }


def TeNe(flux_angles, flux_angles_error, line_names, lines_available,
         param_estimated, param_requered, param_mod_name, param_model_values):

    Te = flx.TeNe()
    Te_error = flx.TeNe()
    Te_fake = flx.TeNe()

    Ne = flx.TeNe()
    Ne_error = flx.TeNe()
    Ne_fake = flx.TeNe()

    ion_abun = flx.Ion_Abun()
    ion_abun_error = flx.Ion_Abun()

    elem_abun_KB = flx.Elem_Abun_KB()
    elem_abun_KB_error = flx.Elem_Abun_KB()
    elem_abun_KB_ratio = flx.Elem_Abun_KB_Ratio()
    elem_abun_KB_ratio_error = flx.Elem_Abun_KB_Ratio()

    elem_ICF_KB = flx.Elem_Abun_KB()
    elem_ICF_KB_error = flx.Elem_Abun_KB()

    elem_abun_DIMS = flx.Elem_Abun_KB()
    elem_abun_DIMS_error = flx.Elem_Abun_KB()
    elem_abun_DIMS_ratio = flx.Elem_Abun_KB_Ratio()
    elem_abun_DIMS_ratio_error = flx.Elem_Abun_KB_Ratio()

    elem_ICF_DIMS = flx.Elem_Abun_KB()
    elem_ICF_DIMS_error = flx.Elem_Abun_KB()

    ####################
    ### for all the angle in the rotation angle task###
    ####################
    pn.log_.open_file('my_log.txt')
    file10 = open('general_output_file.txt', 'a')

    index_atomicData = 12
    name_atomicdata = str(param_mod_name[index_atomicData])

    if (name_atomicdata != "Chianti"):
        if (name_atomicdata != "IRAF_09"):
            pn.atomicData.includeFitsPath()
            pn.atomicData.includeDeprecatedPath()
            pn.atomicData.setDataFileDict(name_atomicdata)
        else:
            pn.atomicData.includeFitsPath()
            pn.atomicData.setDataFileDict(name_atomicdata)
    else:
        pn.atomicData.addAllChianti()

    file6 = open('PyNeb_output_total_abund_ICFs_specific_slits.txt', 'w')
    print('He      ',
          'N_KB    ',
          'ICF(N)_KB',
          'N_DI    ',
          'ICF(N)_DI',
          'O_KB    ',
          'ICF(O)_KB',
          'O_DI    ',
          'ICF(O)_DI',
          'S_KB    ',
          'ICF(S)_KB',
          'S_DI    ',
          'ICF(S)_DI',
          'Ne_KB  ',
          'ICF(Ne)_KB',
          'Ne_DI  ',
          'ICF(Ne)_DI',
          'Cl_KB  ',
          'ICF(Cl)_KB',
          'Cl_DI  ',
          'ICF(Cl)_DI',
          'Ar_KB ',
          'ICF(Ar)_KB',
          'Ar_DI   ',
          'ICF(Ar)_DI',
          'NO_ratio_KB',
          'NO_ratio_DIMS',
          'SO_ratio_KB',
          'SO_ratio_DIMS',
          'NeO_ratio_KB',
          'NeO_ratio_DIMS',
          'ClO_ratio_KB',
          'ClO_ratio_DIMS',
          'ArO_ratio_KB',
          'ArO_ratio_DIMS',
          file=file6)

    file5 = open('PyNeb_output_ionic_abund_specific_slits.txt', 'w')
    print('He+(5876)',
          'He+(6678)',
          'He++(4686)',
          'He++(5412)',
          'N0(5199)',
          'N+(5755)',
          'N+(6548)',
          'N+(6584)',
          'O0(5577)',
          'O0(6300)',
          'O0(6363)',
          'O+(3727)',
          'O+(3729)',
          'O+(7320)',
          'O+(7330)',
          'O+++(4363)',
          'O+++(4959)',
          'O+++(5007)',
          'S+(6716)      ',
          'S+(6731)',
          'S++(6312)',
          'S++(9069)',
          'Ne++(3868)',
          'Ne++(3968)',
          'Cl++(5518)',
          'Cl++(5538)',
          'Ar++(7136)',
          'Ar+++(4712)',
          'Ar+++(4740)',
          file=file5)

    file4 = open('PyNeb_output_Te_and_Ne_specific_slits.txt', 'w')
    print("Te_NII",
          "Ne_SII ",
          "Te_OI",
          "Ne_SII",
          "Te_OII",
          "Ne_SII",
          "Te_OIII",
          "Ne_SII",
          "Te_SIII",
          "Ne_SII",
          "Te_OII",
          "Ne_OII",
          "Te_NII",
          "Ne_OII",
          "Te_OI ",
          "Ne_OII",
          "Te_OIII",
          "Ne_ClIII",
          "Te_SIII",
          "Ne_ClIII",
          "Te_OIII",
          "Ne_ArVI",
          "Te_SIII",
          "Ne_ArVI",
          "Te_NII",
          "Ne_ClIII",
          file=file4)

    for i in range(10):
        Te.append(0, 1)
        Ne.append(0, 1)
        Te_error.append(0, 1)
        Ne_error.append(0, 1)


        print("######################################",
              " specific slit number ",
              i,
              " ##############################",
              file=file4)
        print("######################################",
              " specific slit number ",
              i,
              " ##############################",
              file=file5)
        print("######################################",
              " specific slit number ",
              i,
              " ##############################",
              file=file6)

        file3 = open('test.dat', 'w')
        print("LINE test err", file=file3)

        def foo(element):
            print(">> calling foo with key={:}, value={:}".format(element, dct[element]))
            error = nfes.flux_Error(getattr(flux_angles, dct[element])[i],
                getattr(flux_angles_error, dct[element])[i],
                flux_angles.Hb_4861[i], flux_angles_error.Hb_4861[i])

            print('{:} {:03.3f} {:03.3f}'.format(
                element,
                getattr(flux_angles, dct[element])[i] * 100 / flux_angles.Hb_4861[i],
                error),
                file=file3)
        
        for elmnt in dct: foo(elmnt)

        file3.close()

        index_extlaw = 11
        name_ext = str(param_mod_name[index_extlaw])
        if name_ext == 'CCM89':
            name_ext = str('CCM89')
        elif name_ext == 'K76':
            name_ext = str('K76')
        elif name_ext == 'F99':
            name_ext = str('F99')
        elif name_ext == 'S79_H83_CCM89':
            name_ext = str('S79') + ' ' + str('H83') + ' ' + str('CCM89')
        elif name_ext == 'CCM89_Bal07':
            name_ext = str('CCM89') + ' ' + str('Bal07')
        elif name_ext == 'CCM89_oD94':
            name_ext = str('CCM89') + ' ' + str('oD94')
        elif name_ext == 'SM79_Gal':
            name_ext = str('SM79') + ' ' + str('Gal')
        elif name_ext == 'G03_LMC':
            name_ext = str('G03') + ' ' + str('LMC')
        elif name_ext == 'MCC99_FM90_LMC':
            name_ext = str('MCC99') + ' ' + str('FM90') + ' ' + str('LMC')
        #if name_ext=='F99_like': name_ext=str('F99') + ' ' + str('-') + ' ' + str('like')
        elif name_ext == 'F88_F99_LMC':
            name_ext = str('F88') + ' ' + str('F99') + ' ' + str('LMC')
        elif name_ext == 'No_correction':
            name_ext = str('No') + ' ' + str('correction')
        RV_ext = float(param_model_values[index_extlaw]) * 0.1

        pn.log_.level = -5  # set this to 3 to have more details
        # Adopt an extinction law
        extinction_law = name_ext

        # Define the data file
        obs_data = 'test.dat'

        ### Read and deredden observational data
        # define an Observation object and assign it to name 'obs'
        obs = pn.Observation()
        obs_fake = pn.Observation()

        # fill obs with data read from file obs_data, with lines varying across rows and a default percent error on line intensities
        obs_fake.readData(obs_data,
                          fileFormat='lines_in_rows_err_cols',
                          errIsRelative=False)  #, err_default=0.05)
        obs.readData(obs_data,
                     fileFormat='lines_in_rows_err_cols',
                     errIsRelative=False)  #, err_default=0.05)

        index_MC = param_mod_name.index("number_fake_obs_MC")
        obs_fake.addMonteCarloObs(N=param_model_values[index_MC])

        obs.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)
        obs_fake.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)

        print(
            "######################################################################################",
            file=file10)
        print(
            "######################################################################################",
            file=file10)
        print("PyNeb: cH(b) = {:.2f} +/- {:.2f}".format(
            obs.extinction.cHbeta[0], np.std(obs_fake.extinction.cHbeta)),
              " slit number=",
              i,
              file=file10)
        print("PyNeb: E(B-V) = {:.2f} +/- {:.2f}".format(
            obs.extinction.E_BV[0], np.std(obs_fake.extinction.E_BV)),
              " slit number=",
              i,
              file=file10)

        obs.extinction.law = extinction_law
        obs.extinction.R_V = RV_ext
        obs.correctData(normWave=4861.)

        obs_fake.extinction.law = extinction_law
        obs_fake.extinction.R_V = RV_ext
        obs_fake.correctData(normWave=4861.)
        # deredden data with fitspatric'law
        # instantiate the Diagnostics class
        diags = pn.Diagnostics()
        # include in diags the relevant line ratios

        ########################
        ### Density diagnostics
        ########################
        index_line1 = line_names.index("S2_6716s")
        index_line2 = line_names.index("S2_6731s")
        index_Ne_SII = "no"
        if lines_available[index_line1] == "yes" and lines_available[
                index_line2] == "yes":
            diags.addDiag([
                '[SII] 6731/6716',
            ])
            diags.addClabel('[SII] 6731/6716', '[SII]a')
            index_Ne_SII = "yes"

        index_line1 = line_names.index("O2_3727s")
        index_line2 = line_names.index("O2_3729s")
        index_Ne_OII = "no"
        if lines_available[index_line1] == "yes" and lines_available[
                index_line2] == "yes":
            diags.addDiag([
                '[OII] 3726/3729',
            ])
            index_Ne_OII = "yes"

        index_line1 = line_names.index("Cl3_5517s")
        index_line2 = line_names.index("Cl3_5538s")
        index_Ne_ClIII = "no"
        if lines_available[index_line1] == "yes" and lines_available[
                index_line2] == "yes":
            diags.addDiag([
                '[ClIII] 5538/5518',
            ])
            index_Ne_ClIII = "yes"

        index_line1 = line_names.index("Ar4_4712s")
        index_line2 = line_names.index("Ar4_4740s")
        index_Ne_ArVI = "no"
        if lines_available[index_line1] == "yes" and lines_available[
                index_line2] == "yes":
            diags.addDiag([
                '[ArIV] 4740/4711',
            ])
            index_Ne_ArVI = "yes"

        ###########################
        ### Temperature diagnostics
        ###########################

        index_line1 = line_names.index("N2_5755s")
        index_line2 = line_names.index("N2_6548s")
        index_line3 = line_names.index("N2_6583s")
        index_Te_NII = "no"
        if lines_available[index_line1] == "yes" and lines_available[
                index_line2] == "yes" and lines_available[index_line3] == "yes":
            diags.addDiag([
                '[NII] 5755/6584+',
            ])
            index_Te_NII = "yes"

        index_line1 = line_names.index("O3_4363s")
        index_line2 = line_names.index("O3_4959s")
        index_line3 = line_names.index("O3_5007s")
        index_Te_OIII = "no"
        if lines_available[index_line1] == "yes" and lines_available[
                index_line2] == "yes" and lines_available[index_line3] == "yes":
            diags.addDiag([
                '[OIII] 4363/5007+',
            ])
            index_Te_OIII = "yes"

        index_line1 = line_names.index("O1_6300s")
        index_line2 = line_names.index("O1_6363s")
        index_line3 = line_names.index("O1_5577s")
        index_Te_OI = "no"
        index_Te_OI2 = "no"
        if lines_available[index_line1] == "yes" and lines_available[
                index_line2] == "yes" and lines_available[index_line3] == "yes":
            diags.addDiag([
                '[OI] 5577/6300+',
            ])
            index_Te_OI = "yes"
        elif lines_available[index_line1] == "yes" and lines_available[
                index_line2] == "no" and lines_available[index_line3] == "yes":
            diags.addDiag([
                '[OI] 5577/6300',
            ])
            index_Te_OI2 = "yes"

        index_line1 = line_names.index("O2_3727s")
        index_line2 = line_names.index("O2_3729s")
        index_line3 = line_names.index("O2_7320s")
        index_line4 = line_names.index("O2_7330s")
        index_Te_OII = "no"
        if lines_available[index_line1] == "yes" and lines_available[
                index_line2] == "yes" and lines_available[
                    index_line3] == "yes" and lines_available[
                        index_line4] == "yes":
            diags.addDiag([
                '[OII] 3727+/7325+',
            ])
            index_Te_OII = "yes"

        index_line1 = line_names.index("S3_6312s")
        index_line2 = line_names.index("S3_9069s")
        #index_line3=line_names.index("S3_9592s")
        index_Te_SIII = "no"
        if lines_available[index_line1] == "yes" and lines_available[
                index_line2] == "yes":
            diags.addDiag([
                '[SIII] 6312/9069',
            ])
            index_Te_SIII = "yes"


#       elif lines_available[index_line1]=="yes" and lines_available[index_line2]=="yes" and lines_available[index_line3]=="yes":
#           diags.addDiag([
#                '[SIII] 6312/9200+',
#           ])
#           index_Te_SIII="yes"

#        # Create the emission maps to be compared to the observation data (some overkill here)
#        emisgrids = pn.getEmisGridDict(atom_list=diags.getUniqueAtoms(), den_max=1e6)

### Plot
# Create the contour plot as the intersection of tem-den emission maps with dereddened line ratios
#        diags.plot(emisgrids, obs)

# Place the title
#       plt.title(title)

# Display the plot
#        plt.show()

        print(
            "######################################################################################"
        )
        #The observed ratio can be automatically extracted from an Observation object named obs:
        if (index_Te_NII == "yes" and index_Ne_SII == "yes"):
            a, b = diags.getCrossTemDen('[NII] 5755/6584+',
                                        '[SII] 6731/6716',
                                        obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[NII] 5755/6584+',
                                          '[SII] 6731/6716',
                                          obs=obs)
            Te_fake.NIISII = a
            Ne_fake.NIISII = b

            Te.NIISII[i] = a1
            Ne.NIISII[i] = b1

            a[a == inf] = 0
            a = np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1 * 3 and i2 > a1 * 0.3]

            b[b == inf] = 0
            b = np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1 * 3 and i2 > b1 * 0.3]

            Te_error.NIISII[i] = np.std(np.nan_to_num(a))
            Ne_error.NIISII[i] = np.std(np.nan_to_num(b))
            print(
                'PyNeb: low zones (NII-SII): Te_NII  = {0:5.0f}, Ne_SII = {1:7.1f}  cm-1'
                .format(Te.NIISII[i], Ne.NIISII[i]),
                " slit number=",
                i,
                file=file10)

        if (index_Te_OI == "yes" and index_Te_OI2 == "no" and
                index_Ne_SII == "yes"):
            a, b = diags.getCrossTemDen('[OI] 5577/6300+',
                                        '[SII] 6731/6716',
                                        obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OI] 5577/6300+',
                                          '[SII] 6731/6716',
                                          obs=obs)
            Te_fake.OISII = a
            Ne_fake.OISII = b

            Te.OISII[i] = a1
            Ne.OISII[i] = b1

            a[a == inf] = 0
            a = np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1 * 3 and i2 > a1 * 0.3]

            b[b == inf] = 0
            b = np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1 * 3 and i2 > b1 * 0.3]

            Te_error.OISII[i] = np.std(np.nan_to_num(a))
            Ne_error.OISII[i] = np.std(np.nan_to_num(b))
            print(
                'PyNeb: low zones (OI-SII): Te_OI  = {0:5.0f},Ne_SII = {1:7.1f}  cm-1'
                .format(Te.OISII[i], Ne.OISII[i]),
                " slit number=",
                i,
                file=file10)

        if (index_Te_OI == "no" and index_Te_OI2 == "yes" and
                index_Ne_SII == "yes"):
            a, b = diags.getCrossTemDen('[OI] 5577/6300',
                                        '[SII] 6731/6716',
                                        obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OI] 5577/6300',
                                          '[SII] 6731/6716',
                                          obs=obs)
            Te_fake.OISII = a
            Ne_fake.OISII = b

            Te.OISII[i] = a1
            Ne.OISII[i] = b1

            a[a == inf] = 0
            a = np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1 * 3 and i2 > a1 * 0.3]

            b[b == inf] = 0
            b = np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1 * 3 and i2 > b1 * 0.3]

            Te_error.OISII[i] = np.std(np.nan_to_num(a))
            Ne_error.OISII[i] = np.std(np.nan_to_num(b))
            print(
                'PyNeb: low zones (OI-SII): Te_OI  = {0:5.0f},Ne_SII = {1:7.1f}  cm-1'
                .format(Te.OISII[i], Ne.OISII[i]),
                " slit number=",
                i,
                file=file10)

        if (index_Te_OII == "yes" and index_Ne_SII == "yes"):
            a, b = diags.getCrossTemDen('[OII] 3727+/7325+',
                                        '[SII] 6731/6716',
                                        obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OII] 3727+/7325+',
                                          '[SII] 6731/6716',
                                          obs=obs)
            Te_fake.OIISII = a
            Ne_fake.OIISII = b

            Te.OIISII[i] = a1
            Ne.OIISII[i] = b1

            a[a == inf] = 0
            a = np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1 * 3 and i2 > a1 * 0.3]

            b[b == inf] = 0
            b = np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1 * 3 and i2 > b1 * 0.3]

            Te_error.OIISII[i] = np.std(np.nan_to_num(a))
            Ne_error.OIISII[i] = np.std(np.nan_to_num(b))
            print(
                'PyNeb: low zones (OII-SII): Te_OII  = {0:5.0f},Ne_SII = {1:7.1f}  cm-1'
                .format(Te.OIISII[i], Ne.OIISII[i]),
                " slit number=",
                i,
                file=file10)

        if (index_Te_OIII == "yes" and index_Ne_SII == "yes"):
            a, b = diags.getCrossTemDen('[OIII] 4363/5007+',
                                        '[SII] 6731/6716',
                                        obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OIII] 4363/5007+',
                                          '[SII] 6731/6716',
                                          obs=obs)
            Te_fake.OIIISII = a
            Ne_fake.OIIISII = b

            Te.OIIISII[i] = a1
            Ne.OIIISII[i] = b1

            a[a == inf] = 0
            a = np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1 * 3 and i2 > a1 * 0.3]

            b[b == inf] = 0
            b = np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1 * 3 and i2 > b1 * 0.3]

            Te_error.OIIISII[i] = np.std(np.nan_to_num(a))
            Ne_error.OIIISII[i] = np.std(np.nan_to_num(b))
            print(
                'PyNeb: medium zones (OIII-SII): Te_OIII  = {0:5.0f},Ne_SII = {1:7.1f}  cm-1'
                .format(Te.OIIISII[i], Ne.OIIISII[i]),
                " slit number=",
                i,
                file=file10)

        if (index_Te_SIII == "yes" and index_Ne_SII == "yes"):
            a, b = diags.getCrossTemDen('[SIII] 6312/9069',
                                        '[SII] 6731/6716',
                                        obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[SIII] 6312/9069',
                                          '[SII] 6731/6716',
                                          obs=obs)
            Te_fake.SIIISII = a
            Ne_fake.SIIISII = b

            Te.SIIISII[i] = a1
            Ne.SIIISII[i] = b1

            a[a == inf] = 0
            a = np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1 * 3 and i2 > a1 * 0.3]

            b[b == inf] = 0
            b = np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1 * 3 and i2 > b1 * 0.3]

            Te_error.SIIISII[i] = np.std(np.nan_to_num(a))
            Ne_error.SIIISII[i] = np.std(np.nan_to_num(b))
            print(
                'PyNeb: medium zones (SIII-SII): Te_SIII  = {0:5.0f},Ne_SII = {1:7.1f}  cm-1'
                .format(Te.SIIISII[i], Ne.SIIISII[i]),
                " slit number=",
                i,
                file=file10)

        if (index_Te_OII == "yes" and index_Ne_OII == "yes"):
            a, b = diags.getCrossTemDen('[OII] 3727+/7325+',
                                        '[OII] 3726/3729',
                                        obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OII] 3727+/7325+',
                                          '[OII] 3726/3729',
                                          obs=obs)
            Te_fake.OIIOII = a
            Ne_fake.OIIOII = b

            Te.OIIOII[i] = a1
            Ne.OIIOII[i] = b1

            a[a == inf] = 0
            a = np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1 * 3 and i2 > a1 * 0.3]

            b[b == inf] = 0
            b = np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1 * 3 and i2 > b1 * 0.3]

            Te_error.OIIOII[i] = np.std(np.nan_to_num(a))
            Ne_error.OIIOII[i] = np.std(np.nan_to_num(b))
            print(
                'PyNeb: low zones (OII-OII): Te_OII  = {0:5.0f},Ne_OII = {1:7.1f}  cm-1'
                .format(Te.OIIOII[i], Ne.OIIOII[i]),
                " slit number=",
                i,
                file=file10)

        if (index_Te_NII == "yes" and index_Ne_OII == "yes"):
            a, b = diags.getCrossTemDen('[NII] 5755/6584+',
                                        '[OII] 3726/3729',
                                        obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[NII] 5755/6584+',
                                          '[OII] 3726/3729',
                                          obs=obs)
            Te_fake.NIIOII = a
            Ne_fake.NIIOII = b

            Te.NIIOII[i] = a1
            Ne.NIIOII[i] = b1

            a[a == inf] = 0
            a = np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1 * 3 and i2 > a1 * 0.3]

            b[b == inf] = 0
            b = np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1 * 3 and i2 > b1 * 0.3]

            Te_error.NIIOII[i] = np.std(np.nan_to_num(a))
            Ne_error.NIIOII[i] = np.std(np.nan_to_num(b))
            print(
                'PyNeb: low zones (NII-OII): Te_NII  = {0:5.0f},Ne_OII = {1:7.1f}  cm-1'
                .format(Te.NIIOII[i], Ne.NIIOII[i]),
                " slit number=",
                i,
                file=file10)

        if (index_Te_OI == "yes" and index_Te_OI2 == "no" and
                index_Ne_OII == "yes"):
            a, b = diags.getCrossTemDen('[OI] 5577/6300+',
                                        '[OII] 3726/3729',
                                        obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OI] 5577/6300+',
                                          '[OII] 3726/3729',
                                          obs=obs)
            Te_fake.OIOII = a
            Ne_fake.OIOII = b

            Te.OIOII[i] = a1
            Ne.OIOII[i] = b1

            a[a == inf] = 0
            a = np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1 * 3 and i2 > a1 * 0.3]

            b[b == inf] = 0
            b = np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1 * 3 and i2 > b1 * 0.3]

            Te_error.OIOII[i] = np.std(np.nan_to_num(a))
            Ne_error.OIOII[i] = np.std(np.nan_to_num(b))
            print(
                'PyNeb: low zones (OI-OII): Te_OI  = {0:5.0f},Ne_OII = {1:7.1f}  cm-1'
                .format(Te.OIOII[i], Ne.OIOII[i]),
                " slit number=",
                i,
                file=file10)

        if (index_Te_OI == "no" and index_Te_OI2 == "yes" and
                index_Ne_OII == "yes"):
            a, b = diags.getCrossTemDen('[OI] 5577/6300',
                                        '[OII] 3726/3729',
                                        obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OI] 5577/6300',
                                          '[OII] 3726/3729',
                                          obs=obs)
            Te_fake.OIOII = a
            Ne_fake.OIOII = b

            Te.OIOII[i] = a1
            Ne.OIOII[i] = b1

            a[a == inf] = 0
            a = np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1 * 3 and i2 > a1 * 0.3]

            b[b == inf] = 0
            b = np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1 * 3 and i2 > b1 * 0.3]

            Te_error.OIOII[i] = np.std(np.nan_to_num(a))
            Ne_error.OIOII[i] = np.std(np.nan_to_num(b))
            print(
                'PyNeb: low zones (OI-OII): Te_OI  = {0:5.0f},Ne_OII = {1:7.1f}  cm-1'
                .format(Te.OIOII[i], Ne.OIOII[i]),
                " slit number=",
                i,
                file=file10)

        if (index_Te_OIII == "yes" and index_Ne_ClIII == "yes"):
            a, b = diags.getCrossTemDen('[OIII] 4363/5007+',
                                        '[ClIII] 5538/5518',
                                        obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OIII] 4363/5007+',
                                          '[ClIII] 5538/5518',
                                          obs=obs)
            Te_fake.OIIIClIII = a
            Ne_fake.OIIIClIII = b

            Te.OIIIClIII[i] = a1
            Ne.OIIIClIII[i] = b1

            a[a == inf] = 0
            a = np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1 * 3 and i2 > a1 * 0.3]

            b[b == inf] = 0
            b = np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1 * 3 and i2 > b1 * 0.3]

            Te_error.OIIIClIII[i] = np.std(np.nan_to_num(a))
            Ne_error.OIIIClIII[i] = np.std(np.nan_to_num(b))
            print(
                'PyNeb: medium zones (OIII-ClIII): Te_OIII  = {0:5.0f},Ne_ClIII = {1:7.1f}  cm-1'
                .format(Te.OIIIClIII[i], Ne.OIIIClIII[i]),
                " slit number=",
                i,
                file=file10)

        if (index_Te_SIII == "yes" and index_Ne_ClIII == "yes"):
            a, b = diags.getCrossTemDen('[SIII] 6312/9069',
                                        '[ClIII] 5538/5518',
                                        obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[SIII] 6312/9069',
                                          '[ClIII] 5538/5518',
                                          obs=obs)
            Te_fake.SIIIClIII = a
            Ne_fake.SIIIClIII = b

            Te.SIIIClIII[i] = a1
            Ne.SIIIClIII[i] = b1

            a[a == inf] = 0
            a = np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1 * 3 and i2 > a1 * 0.3]

            b[b == inf] = 0
            b = np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1 * 3 and i2 > b1 * 0.3]

            Te_error.SIIIClIII[i] = np.std(np.nan_to_num(a))
            Ne_error.SIIIClIII[i] = np.std(np.nan_to_num(b))
            print(
                'PyNeb: medium zones (SIII-ClIII): Te_SIII  = {0:5.0f},Ne_ClIII = {1:7.1f}  cm-1'
                .format(Te.SIIIClIII[i], Ne.SIIIClIII[i]),
                " slit number=",
                i,
                file=file10)

        if (index_Te_OIII == "yes" and index_Ne_ArVI == "yes"):
            a, b = diags.getCrossTemDen('[OIII] 4363/5007+',
                                        '[ArIV] 4740/4711',
                                        obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OIII] 4363/5007+',
                                          '[ArIV] 4740/4711',
                                          obs=obs)
            Te_fake.OIIIArVI = a
            Ne_fake.OIIIArVI = b

            Te.OIIIArVI[i] = a1
            Ne.OIIIArVI[i] = b1

            a[a == inf] = 0
            a = np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1 * 3 and i2 > a1 * 0.3]

            b[b == inf] = 0
            b = np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1 * 3 and i2 > b1 * 0.3]

            Te_error.OIIIArVI[i] = np.std(np.nan_to_num(a))
            Ne_error.OIIIArVI[i] = np.std(np.nan_to_num(b))
            print(
                'PyNeb: medium zones (OIII-ArVI): Te_OIII  = {0:5.0f},Ne_ArVI = {1:7.1f}  cm-1'
                .format(Te.OIIIArVI[i], Ne.OIIIArVI[i]),
                " slit number=",
                i,
                file=file10)

        if (index_Te_SIII == "yes" and index_Ne_ArVI == "yes"):
            a, b = diags.getCrossTemDen('[SIII] 6312/9069',
                                        '[ArIV] 4740/4711',
                                        obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[SIII] 6312/9069',
                                          '[ArIV] 4740/4711',
                                          obs=obs)
            Te_fake.NIISII = a
            Ne_fake.NIISII = b

            Te.SIIIArVI[i] = a1
            Ne.SIIIArVI[i] = b1

            a[a == inf] = 0
            a = np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1 * 3 and i2 > a1 * 0.3]

            b[b == inf] = 0
            b = np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1 * 3 and i2 > b1 * 0.3]

            Te_error.SIIIArVI[i] = np.std(np.nan_to_num(a))
            Ne_error.SIIIArVI[i] = np.std(np.nan_to_num(b))
            print(
                'PyNeb: medium zones (SIII-ArVI): Te_SIII  = {0:5.0f},Ne_ArVI = {1:7.1f}  cm-1'
                .format(Te.SIIIArVI[i], Ne.SIIIArVI[i]),
                " slit number=",
                i,
                file=file10)

        if (index_Te_NII == "yes" and index_Ne_ClIII == "yes"):
            a, b = diags.getCrossTemDen('[NII] 5755/6584+',
                                        '[ClIII] 5538/5518',
                                        obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[NII] 5755/6584+',
                                          '[ClIII] 5538/5518',
                                          obs=obs)
            Te_fake.NIIClIII = a
            Ne_fake.NIIClIII = b

            Te.NIIClIII[i] = a1
            Ne.NIIClIII[i] = b1

            a[a == inf] = 0
            a = np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1 * 3 and i2 > a1 * 0.3]

            b[b == inf] = 0
            b = np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1 * 3 and i2 > b1 * 0.3]

            Te_error.NIIClIII[i] = np.std(np.nan_to_num(a))
            Ne_error.NIIClIII[i] = np.std(np.nan_to_num(b))
            print(
                'PyNeb: medium zones (NII-ClIII): Te_NII  = {0:5.0f},Ne_ClIII = {1:7.1f}  cm-1'
                .format(Te.NIIClIII[i], Ne.NIIClIII[i]),
                " slit number=",
                i,
                file=file10)

        print('{:06.0f}'.format(Te.NIISII[i]),
              '{:06.0f}'.format(Ne.NIISII[i]),
              '{:06.0f}'.format(Te.OISII[i]),
              '{:06.0f}'.format(Ne.OISII[i]),
              '{:06.0f}'.format(Te.OIISII[i]),
              '{:06.0f}'.format(Ne.OIISII[i]),
              '{:06.0f}'.format(Te.OIIISII[i]),
              '{:07.0f}'.format(Ne.OIIISII[i]),
              '{:07.0f}'.format(Te.SIIISII[i]),
              '{:06.0f}'.format(Ne.SIIISII[i]),
              '{:06.0f}'.format(Te.OIIOII[i]),
              '{:06.0f}'.format(Ne.OIIOII[i]),
              '{:06.0f}'.format(Te.NIIOII[i]),
              '{:06.0f}'.format(Ne.NIIOII[i]),
              '{:06.0f}'.format(Te.OIOII[i]),
              '{:06.0f}'.format(Ne.OIOII[i]),
              '{:07.0f}'.format(Te.OIIIClIII[i]),
              '{:08.0f}'.format(Ne.OIIIClIII[i]),
              '{:07.0f}'.format(Te.SIIIClIII[i]),
              '{:08.0f}'.format(Ne.SIIIClIII[i]),
              '{:07.0f}'.format(Te.OIIIArVI[i]),
              '{:07.0f}'.format(Ne.OIIIArVI[i]),
              '{:07.0f}'.format(Te.SIIIArVI[i]),
              '{:07.0f}'.format(Ne.SIIIArVI[i]),
              '{:07.0f}'.format(Te.NIIClIII[i]),
              '{:07.0f}'.format(Ne.NIIClIII[i]),
              file=file4)

        print('{:06.0f}'.format(Te_error.NIISII[i]),
              '{:06.0f}'.format(Ne_error.NIISII[i]),
              '{:06.0f}'.format(Te_error.OISII[i]),
              '{:06.0f}'.format(Ne_error.OISII[i]),
              '{:06.0f}'.format(Te_error.OIISII[i]),
              '{:06.0f}'.format(Ne_error.OIISII[i]),
              '{:06.0f}'.format(Te_error.OIIISII[i]),
              '{:07.0f}'.format(Ne_error.OIIISII[i]),
              '{:07.0f}'.format(Te_error.SIIISII[i]),
              '{:06.0f}'.format(Ne_error.SIIISII[i]),
              '{:06.0f}'.format(Te_error.OIIOII[i]),
              '{:06.0f}'.format(Ne_error.OIIOII[i]),
              '{:06.0f}'.format(Te_error.NIIOII[i]),
              '{:06.0f}'.format(Ne_error.NIIOII[i]),
              '{:06.0f}'.format(Te_error.OIOII[i]),
              '{:06.0f}'.format(Ne_error.OIOII[i]),
              '{:07.0f}'.format(Te_error.OIIIClIII[i]),
              '{:08.0f}'.format(Ne_error.OIIIClIII[i]),
              '{:07.0f}'.format(Te_error.SIIIClIII[i]),
              '{:08.0f}'.format(Ne_error.SIIIClIII[i]),
              '{:07.0f}'.format(Te_error.OIIIArVI[i]),
              '{:07.0f}'.format(Ne_error.OIIIArVI[i]),
              '{:07.0f}'.format(Te_error.SIIIArVI[i]),
              '{:07.0f}'.format(Ne_error.SIIIArVI[i]),
              '{:07.0f}'.format(Te_error.NIIClIII[i]),
              '{:07.0f}'.format(Ne_error.NIIClIII[i]),
              file=file4)

        ###########################
        ### calculate, print and save ionic abundances
        ##########################
        ion_abun, ion_abun_error, atom_abun_fake = iasss.ionicabundances(
            obs, obs_fake, Te, Te_error, Te_fake, Ne, Ne_error, Ne_fake, i,
            line_names, lines_available, ion_abun, ion_abun_error,
            param_estimated, param_requered, param_mod_name)

        print('{:5.3e}'.format(ion_abun.HeIa[i]),
              '{:5.3e}'.format(ion_abun.HeIb[i]),
              '{:5.4e}'.format(ion_abun.HeIIa[i]),
              '{:5.4e}'.format(ion_abun.HeIIb[i]),
              '{:5.2e}'.format(ion_abun.NI[i]),
              '{:5.2e}'.format(ion_abun.NIIa[i]),
              '{:5.2e}'.format(ion_abun.NIIb[i]),
              '{:5.2e}'.format(ion_abun.NIIc[i]),
              '{:5.2e}'.format(ion_abun.OIa[i]),
              '{:5.2e}'.format(ion_abun.OIb[i]),
              '{:5.2e}'.format(ion_abun.OIc[i]),
              '{:5.2e}'.format(ion_abun.OIIa[i]),
              '{:5.2e}'.format(ion_abun.OIIb[i]),
              '{:5.2e}'.format(ion_abun.OIIc[i]),
              '{:5.2e}'.format(ion_abun.OIId[i]),
              '{:5.4e}'.format(ion_abun.OIIIa[i]),
              '{:5.4e}'.format(ion_abun.OIIIb[i]),
              '{:5.4e}'.format(ion_abun.OIIIc[i]),
              '{:5.2e}'.format(ion_abun.SIIa[i]),
              '{:5.2e}'.format(ion_abun.SIIb[i]),
              '{:5.3e}'.format(ion_abun.SIIIa[i]),
              '{:5.3e}'.format(ion_abun.SIIIb[i]),
              '{:5.4e}'.format(ion_abun.NeIIIa[i]),
              '{:5.4e}'.format(ion_abun.NeIIIb[i]),
              '{:5.4e}'.format(ion_abun.ClIIIa[i]),
              '{:5.4e}'.format(ion_abun.ClIIIb[i]),
              '{:5.4e}'.format(ion_abun.ArIII[i]),
              '{:5.5e}'.format(ion_abun.ArIVa[i]),
              '{:5.5e}'.format(ion_abun.ArIVb[i]),
              file=file5)

        print('{:5.3e}'.format(ion_abun_error.HeIa[i]),
              '{:5.3e}'.format(ion_abun_error.HeIb[i]),
              '{:5.4e}'.format(ion_abun_error.HeIIa[i]),
              '{:5.4e}'.format(ion_abun_error.HeIIb[i]),
              '{:5.2e}'.format(ion_abun_error.NI[i]),
              '{:5.2e}'.format(ion_abun_error.NIIa[i]),
              '{:5.2e}'.format(ion_abun_error.NIIb[i]),
              '{:5.2e}'.format(ion_abun_error.NIIc[i]),
              '{:5.2e}'.format(ion_abun_error.OIa[i]),
              '{:5.2e}'.format(ion_abun_error.OIb[i]),
              '{:5.2e}'.format(ion_abun_error.OIc[i]),
              '{:5.2e}'.format(ion_abun_error.OIIa[i]),
              '{:5.2e}'.format(ion_abun_error.OIIb[i]),
              '{:5.2e}'.format(ion_abun_error.OIIc[i]),
              '{:5.2e}'.format(ion_abun_error.OIId[i]),
              '{:5.4e}'.format(ion_abun_error.OIIIa[i]),
              '{:5.4e}'.format(ion_abun_error.OIIIb[i]),
              '{:5.4e}'.format(ion_abun_error.OIIIc[i]),
              '{:5.2e}'.format(ion_abun_error.SIIa[i]),
              '{:5.2e}'.format(ion_abun_error.SIIb[i]),
              '{:5.3e}'.format(ion_abun_error.SIIIa[i]),
              '{:5.3e}'.format(ion_abun_error.SIIIb[i]),
              '{:5.4e}'.format(ion_abun_error.NeIIIa[i]),
              '{:5.4e}'.format(ion_abun_error.NeIIIb[i]),
              '{:5.4e}'.format(ion_abun_error.ClIIIa[i]),
              '{:5.4e}'.format(ion_abun_error.ClIIIb[i]),
              '{:5.4e}'.format(ion_abun_error.ArIII[i]),
              '{:5.5e}'.format(ion_abun_error.ArIVa[i]),
              '{:5.5e}'.format(ion_abun_error.ArIVb[i]),
              file=file5)

        ###########################
        ### calculate, print and save total elemental abundances and ICFs
        ##########################
        elem_abun_KB, elem_ICF_KB, elem_abun_DIMS, elem_ICF_DIMS, elem_abun_KB_error, elem_ICF_KB_error, elem_abun_DIMS_error, elem_ICF_DIMS_error, elem_abun_KB_ratio, elem_abun_KB_ratio_error, elem_abun_DIMS_ratio, elem_abun_DIMS_ratio_error = eaIsss.total_abundances_ICFs(
            ion_abun, ion_abun_error, atom_abun_fake, i, line_names,
            lines_available, elem_abun_KB, elem_ICF_KB, elem_abun_DIMS,
            elem_ICF_DIMS, elem_abun_KB_error, elem_ICF_KB_error,
            elem_abun_DIMS_error, elem_ICF_DIMS_error, elem_abun_KB_ratio,
            elem_abun_KB_ratio_error, elem_abun_DIMS_ratio,
            elem_abun_DIMS_ratio_error, param_mod_name)

        print('{:5.2e}'.format(elem_abun_KB.He[i]),
              '{:5.2e}'.format(elem_abun_KB.N[i]),
              '{:09.5f}'.format(elem_ICF_KB.N[i]),
              '{:5.2e}'.format(elem_abun_DIMS.N[i]),
              '{:09.5f}'.format(elem_ICF_DIMS.N[i]),
              '{:5.2e}'.format(elem_abun_KB.O[i]),
              '{:09.5f}'.format(elem_ICF_KB.O[i]),
              '{:5.2e}'.format(elem_abun_DIMS.O[i]),
              '{:09.5f}'.format(elem_ICF_DIMS.O[i]),
              '{:5.2e}'.format(elem_abun_KB.S[i]),
              '{:09.5f}'.format(elem_ICF_KB.S[i]),
              '{:5.2e}'.format(elem_abun_DIMS.S[i]),
              '{:09.5f}'.format(elem_ICF_DIMS.S[i]),
              '{:5.2e}'.format(elem_abun_KB.Ne[i]),
              '{:09.5f}'.format(elem_ICF_KB.Ne[i]),
              '{:5.2e}'.format(elem_abun_DIMS.Ne[i]),
              '{:09.5f}'.format(elem_ICF_DIMS.Ne[i]),
              '{:5.2e}'.format(elem_abun_KB.Cl[i]),
              '{:09.5f}'.format(elem_ICF_KB.Cl[i]),
              '{:5.2e}'.format(elem_abun_DIMS.Cl[i]),
              '{:09.5f}'.format(elem_ICF_DIMS.Cl[i]),
              '{:5.2e}'.format(elem_abun_KB.Ar[i]),
              '{:09.5f}'.format(elem_ICF_KB.Ar[i]),
              '{:5.2e}'.format(elem_abun_DIMS.Ar[i]),
              '{:09.5f}'.format(elem_ICF_DIMS.Ar[i]),
              '{:9.2e}'.format(elem_abun_KB_ratio.NO[i]),
              '{:10.2e}'.format(elem_abun_DIMS_ratio.NO[i]),
              '{:9.2e}'.format(elem_abun_KB_ratio.SO[i]),
              '{:10.2e}'.format(elem_abun_DIMS_ratio.SO[i]),
              '{:9.2e}'.format(elem_abun_KB_ratio.NeO[i]),
              '{:10.2e}'.format(elem_abun_DIMS_ratio.NeO[i]),
              '{:9.2e}'.format(elem_abun_KB_ratio.ClO[i]),
              '{:10.2e}'.format(elem_abun_DIMS_ratio.ClO[i]),
              '{:9.2e}'.format(elem_abun_KB_ratio.ClO[i]),
              '{:10.2e}'.format(elem_abun_DIMS_ratio.ClO[i]),
              file=file6)

        print('{:5.2e}'.format(elem_abun_KB_error.He[i]),
              '{:5.2e}'.format(elem_abun_KB_error.N[i]),
              '{:09.5f}'.format(elem_ICF_KB_error.N[i]),
              '{:5.2e}'.format(elem_abun_DIMS_error.N[i]),
              '{:09.5f}'.format(elem_ICF_DIMS_error.N[i]),
              '{:5.2e}'.format(elem_abun_KB_error.O[i]),
              '{:09.5f}'.format(elem_ICF_KB_error.O[i]),
              '{:5.2e}'.format(elem_abun_DIMS_error.O[i]),
              '{:09.5f}'.format(elem_ICF_DIMS_error.O[i]),
              '{:5.2e}'.format(elem_abun_KB_error.S[i]),
              '{:09.5f}'.format(elem_ICF_KB_error.S[i]),
              '{:5.2e}'.format(elem_abun_DIMS_error.S[i]),
              '{:09.5f}'.format(elem_ICF_DIMS_error.S[i]),
              '{:5.2e}'.format(elem_abun_KB_error.Ne[i]),
              '{:09.5f}'.format(elem_ICF_KB_error.Ne[i]),
              '{:5.2e}'.format(elem_abun_DIMS_error.Ne[i]),
              '{:09.5f}'.format(elem_ICF_DIMS_error.Ne[i]),
              '{:5.2e}'.format(elem_abun_KB_error.Cl[i]),
              '{:09.5f}'.format(elem_ICF_KB_error.Cl[i]),
              '{:5.2e}'.format(elem_abun_DIMS_error.Cl[i]),
              '{:09.5f}'.format(elem_ICF_DIMS_error.Cl[i]),
              '{:5.2e}'.format(elem_abun_KB_error.Ar[i]),
              '{:09.5f}'.format(elem_ICF_KB_error.Ar[i]),
              '{:5.2e}'.format(elem_abun_DIMS_error.Ar[i]),
              '{:09.5f}'.format(elem_ICF_DIMS_error.Ar[i]),
              '{:9.2e}'.format(elem_abun_KB_ratio_error.NO[i]),
              '{:10.2e}'.format(elem_abun_DIMS_ratio_error.NO[i]),
              '{:9.2e}'.format(elem_abun_KB_ratio_error.SO[i]),
              '{:10.2e}'.format(elem_abun_DIMS_ratio_error.SO[i]),
              '{:9.2e}'.format(elem_abun_KB_ratio_error.NeO[i]),
              '{:10.2e}'.format(elem_abun_DIMS_ratio_error.NeO[i]),
              '{:9.2e}'.format(elem_abun_KB_ratio_error.ClO[i]),
              '{:10.2e}'.format(elem_abun_DIMS_ratio_error.ClO[i]),
              '{:9.2e}'.format(elem_abun_KB_ratio_error.ClO[i]),
              '{:10.2e}'.format(elem_abun_DIMS_ratio_error.ClO[i]),
              file=file6)

        os.remove("test.dat")
    file4.close()
    file5.close()
    file10.close()

    pn.log_.close_file()

    return Te, Ne
