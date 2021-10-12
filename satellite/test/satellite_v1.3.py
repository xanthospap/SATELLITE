#!/usr/bin/env python

# SATELLITE, the Spectroscopic Analysis Tool for intEgraL fieLd unIt daTacubEs
# (C) 2021- Stavros Akras

# SATELLITE uses several software packages in Python: Matplotlib
# (Hunter 2007), NumPy (van der Walt et al. 2011), SciPy (Virtanen
# et al. 2020), seaborn (Waskom & the seaborn development team
# 2020) and AstroPy (Astropy Collaboration et al. 2013, 2018),
# as well as the PyNeb package developed by Luridiana V., Morisset C.,
# Shaw R. A., 2015, A&A, 573, A42

# SATELLITE is a free code: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation.

# It is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

from __future__ import print_function
import warnings
import datetime
import argparse
import sys

from satellite import read_input_script as read
from satellite import read_input_lines_parameters_script as read_ilps
from satellite import radial_analysis_script as ras
from satellite import rotate_line_fluxes_script as rotlfs
from satellite import specificPA_line_fluxes_script as sPAlfs
from satellite import analysis2D_script as a2Ds
from satellite import my_analysis2D_script as mya2Ds
from satellite import TeNe_angles_script as TeNeangles
from satellite import TeNe_specific_slits_script as TeNeslits
from satellite import TeNe_2D_script as TeNe2Ds
from satellite import TeNe_radial_script as TeNers
from satellite import flux_angles_norm_definition_script as fands
from satellite import define_Te_Ne_for_null_spaxels_script as dTeNenullspaxelss
from satellite import exclude_probematic_errors_script as epes
from satellite import physical_parameters_script as pps
from satellite import fluxes as flx

##  set the cmd parser
parser = argparse.ArgumentParser(
    description='Do some processing with planetary images ...'
    'by Stavros Akras',
    epilog='NOA - 2021')
parser.add_argument('-o',
                    '--output',
                    action='store',
                    required=False,
                    help='The output file',
                    metavar='OUT_FILE',
                    dest='out_file',
                    default='general_output_file.txt')
parser.add_argument(
    '--trace-ram-usage',
    action='store_true',
    help=
    'Trace RAM usage for debugging purposes; only available with python > 3.4',
    dest='trace_ram',
    default=False)

# parse cmd
args = parser.parse_args()

print('=========================================')
prog_start = datetime.datetime.now()
print("Start Program")
print(prog_start.strftime("%Y-%m-%d %H:%M:%S"))
print('=========================================')

if args.trace_ram:
    import tracemalloc

warnings.filterwarnings("ignore")
file10 = open(args.out_file, 'w')

#######################################################################################################
# read all input parameters given by the user.
#######################################################################################################
line_names, line_ext_error, lines_available, lines_radial, param_estimated, param_requered, param_mod_name, param_model_values, DD_name, DD_avail, DDxmin, DDxmax, DDymin, DDymax, par_plotname, par_plotymin, par_plotymax = read_ilps.read_input_lines_parameters(
)

### if new lines must be included in the code, they have to be added in the input file, read them in the read.read_input_script for the flux and error
### and finally define their tables in the flux classes

#######################################################################################################
# read line flux images
#######################################################################################################
flux, hdr = read.read_input_images(line_names, lines_available,
                                   param_model_values, "fluxes")

#######################################################################################################
# read line flux error images
#######################################################################################################
flux_error, hdrer = read.read_input_images(line_names, lines_available,
                                           param_model_values, "errors")

flux, flux_error = epes.exclude_problematic_error_values(
    flux, flux_error, param_mod_name, param_model_values)

index_size = param_mod_name.index("total_num_pixels_horiz")
index_pixel_scale = param_mod_name.index("pixel_scale")
sizex = param_model_values[index_size]
sizey = param_model_values[index_size]
pixscale = float(param_model_values[index_pixel_scale]) * 0.01
minx = miny = 0
print('x=', sizex, 'y=', sizey, 'pixel scale=', pixscale, file=file10)

#######################################################################################################
# calculated the exctinction coefficient based of the Ha,Hb,Hg and Hd line images
# and returns the table of the extinction.
#######################################################################################################
#datachb,datachg,datachd,datachber,datachger,datachder=ext.extinction(flux2D.Ha_6563,flux2D.Hb_4861,flux2D.Hg_4340,flux2D.Hd_4101,flux2D_error.Ha_6563,flux2D_error.Hb_4861,flux2D_error.Hg_4340,flux2D_error.Hd_4101,sizex,sizey,line_names,lines_available,param_estimated,param_requered)

#######################################################################################################
#convert the tables to fits images and save the corresponding images.
#######################################################################################################
#sfs.save_image(datachb,hdr,'c_Hb',param_estimated,param_requered)
#sfs.save_image(datachg,hdr,'c_Hg',param_estimated,param_requered)
#sfs.save_image(datachd,hdr,'c_Hd',param_estimated,param_requered)
#sfs.save_image(datachber,hdr,'c_Hb_error',param_estimated,param_requered)
#sfs.save_image(datachber,hdr,'c_Hg_error',param_estimated,param_requered)
#sfs.save_image(datachber,hdr,'c_Hd_error',param_estimated,param_requered)

#hs.histogram_plots(datachb,datachber,'c_Hb',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)

#dataposChb,dataposChg,dataposChd=[],[],[]

#dataposChb=pvs.positivevalues(datachb,flux2D.Ha_6563,flux2D.Hb_4861,sizex,sizey,-1)
#meanChb=mvs.meanvalue(dataposChb,0)
#stdCHb=svs.stdvalue(dataposChb,meanChb,0)

#dataposChg=pvs.positivevalues(datachg,flux2D.Ha_6563,flux2D.Hb_4861,sizex,sizey,-1)
#meanChg=mvs.meanvalue(dataposChg,0)
#stdCHg=svs.stdvalue(dataposChg,meanChg,0)

#dataposChd=pvs.positivevalues(datachd,flux2D.Ha_6563,flux2D.Hb_4861,sizex,sizey,-1)
#meanChd=mvs.meanvalue(dataposChd,0)
#stdCHd=svs.stdvalue(dataposChd,meanChd,0)

#Xb=3.6 #(Fitzpatrick law 1999)
#b_v=meanChb/(0.4*Xb)
#b_v=0.61*meanChb+0.024*meanChb*meanChb # from Pyneb
#print >>file10, "#######################################################################################"
#print >>file10, "average: c(Hb)","%.3f" % meanChb,"std_c(Hb)","%.3f" % stdCHb, "E(B-V)", "%.3f" % b_v, "with outliers"
#print >>file10, "average: c(Hg)","%.3f" % meanChg,"std_c(Hg)","%.3f" % stdCHg, "E(B-V)","with outliers"
#print >>file10, "average: c(Hd)","%.3f" % meanChd,"std_c(Hd)","%.3f" % stdCHd, "E(B-V)","with outliers"
#print >>file10, "#######################################################################################"
# keeps only the value of the table "> 0"
#dataposChb= [x for x in dataposChb if x > 0]

#meanchb,stdchb= ceos.estimates_without_outliers(dataposChb,datachb,sizex,sizey)
#b_v=meanChb/(0.4*Xb)
#print >>file10, "average: c(Hb)","%.3f" % meanchb,"std_c(Hb)","%.3f" % stdchb, "E(B-V)","%.3f" % b_v, "without outliers"
#print >>file10, "#######################################################################################"

#######################################################################################################
# calculate line fluxes, normalized lines relative to Hbeta, and a number of line ratios for radial slits with PA from X1 to X2 given by the user.
#######################################################################################################
index = param_estimated.index('angular_analysis_task')
if param_requered[index] == "yes":
    flux_angles, flux_angles_error, flux_angles_norm, flux_angles_norm_error, ratio_angles, ratio_angles_error, angles = rotlfs.rotate_line_fluxes(
        flux, flux_error, line_names, line_ext_error, lines_available,
        param_estimated, param_requered, param_mod_name, param_model_values,
        par_plotname, par_plotymin, par_plotymax)

    #######################################################################################################
    # use the derived fluxes and ratios as input to Pyneb and calculate paramteres such as
    # Te, Ne, ionic and chemical abudnances,  and ICFs
    #######################################################################################################
    Te_PA, Ne_PA = TeNeangles.TeNe(flux_angles, flux_angles_error, angles,
                                   line_names, lines_available, param_estimated,
                                   param_requered, param_mod_name,
                                   param_model_values, par_plotname,
                                   par_plotymin, par_plotymax)

#######################################################################################################
# radial analysis of various emission lines for a specific PA of a slit
#######################################################################################################
index = param_estimated.index('radial_analysis_task')
if param_requered[index] == "yes":
    flux_radial, flux_radial_error, distance_CS = ras.radial_analysis(
        flux, flux_error, line_names, line_ext_error, lines_available,
        lines_radial, param_mod_name, param_model_values, par_plotname,
        par_plotymin, par_plotymax)

    #######################################################################################################
    # use the derived fluxes and ratios as input to Pyneb and calculate paramteres such as
    # Te, Ne, ionic and chemical abudnances,  and ICFs
    #######################################################################################################
    Te_rad, Ne_rad = TeNers.TeNe(flux_radial, flux_radial_error, distance_CS,
                                 line_names, lines_available, param_estimated,
                                 param_requered, param_mod_name,
                                 param_model_values, par_plotname, par_plotymin,
                                 par_plotymax)

#######################################################################################################
# calculate line fluxes, normalized lines relative to Hbeta, and a number of line ratios for specific slits with PAs, widths and lengths given by the user.
#######################################################################################################
index = param_estimated.index('specific_slits_analysis_task')
if param_requered[index] == "yes":
    flux_spec_slit, flux_spec_slit_error, flux_spec_slit_norm, flux_spec_slit_norm_error, ratio_spec_slit, ratio_spec_slit_error = sPAlfs.specficPA_line_fluxes(
        flux, flux_error, line_names, line_ext_error, lines_available,
        lines_radial, param_estimated, param_requered, param_mod_name,
        param_model_values)

    #######################################################################################################
    # use the derived fluxes and ratios as input to Pyneb and calculate paramteres such as
    # Te, Ne, ionic and chemical abudnances,  and ICFs
    #######################################################################################################
    Te_slits, Ne_slits = TeNeslits.TeNe(flux_spec_slit, flux_spec_slit_error,
                                        line_names, lines_available,
                                        param_estimated, param_requered,
                                        param_mod_name, param_model_values)

#######################################################################################################
# 2D analysis and diagnostic diagrams
#######################################################################################################
index = param_estimated.index('2D_anaysis_task')
if param_requered[index] == "yes":

    index = param_estimated.index('angular_analysis_task')
    if param_requered[index] == "no":
        flux_angles_norm, flux_angles_norm_error, angles = fands.flux_ang_norm_definition(
            param_mod_name, param_model_values)
        Te_PA, Ne_PA = dTeNenullspaxelss.nullTeNe(param_mod_name,
                                                  param_model_values)

    index = param_estimated.index('specific_slits_analysis_task')
    if param_requered[index] == "no":
        flux_spec_slit_norm, flux_spec_slit_norm_error = fands.flux_ang_norm_definition_specificslits(
            param_mod_name, param_model_values)
        Te_slits, Ne_slits = dTeNenullspaxelss.nullTeNe_specific_slits(
            param_mod_name, param_model_values)

    a2Ds.analysis2D(flux, flux_error, flux_angles_norm, angles, line_names,
                    line_ext_error, lines_available, lines_radial,
                    param_estimated, param_requered, param_mod_name,
                    param_model_values, DD_name, DD_avail, DDxmin, DDxmax,
                    DDymin, DDymax, hdr, flux_spec_slit_norm,
                    flux_spec_slit_norm_error)

    #######################################################################################################
    # use the derived fluxes and ratios as input to Pyneb and calculate paramteres such as
    # Te, Ne, ionic and chemical abudnances,  and ICFs
    #######################################################################################################
    Te_2D, Ne_2D = TeNe2Ds.TeNe(flux, flux_error, line_names, line_ext_error,
                                lines_available, param_mod_name,
                                param_model_values, param_estimated,
                                param_requered, hdr, Te_PA, Ne_PA)

#######################################################################################################
# user can write his/her proper script for a 2D analysis and diagnostic diagrams
#######################################################################################################
index = param_estimated.index('My_2D_anaysis_task')
if param_requered[index] == "yes":
    mya2Ds.analysis2D(flux, flux_angles_norm, angles, param_mod_name,
                      param_model_values)

#######################################################################################################
# use the derived fluxes and ratios as input to Pyneb and calculate paramteres such as
# Te, Ne, ionic and chemical abudnances,  and ICFs
#######################################################################################################
#pps.analysiswithPyNeb(flux2D,flux2D_error,flux_angles,flux_angles_error,flux_spec_slit,flux_spec_slit_error,angles,line_names,lines_available,lines_radial,param_estimated,param_requered,param_mod_name,param_model_values,hdr)

#######################################################################################################
# user can write his/her proper PyNeb script
#######################################################################################################
index = param_estimated.index('My_physical_parameters_task')
if param_requered[index] == "yes":
    pps.analysiswithPyNeb(flux, flux_angles, flux_spec_slit, angles, line_names,
                          lines_available, lines_radial, param_estimated,
                          param_requered, param_mod_name, param_model_values,
                          hdr)

file10.close()

prog_end = datetime.datetime.now()
print('=========================================')
print("End Program")
print(prog_end.strftime("%Y-%m-%d %H:%M:%S"))
print("Elapsed time")
print((prog_end - prog_start))
print('=========================================')
