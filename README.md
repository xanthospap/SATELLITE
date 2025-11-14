# Introduction

The Spectroscopic Analysis Tool for intEgraL fieLd unIt daTacubEs (Satellite) is a 
newly developed code for the spectroscopic characterization of extended photo-ionised 
nebulae such as planetary nebulae, H II regions or galaxies in the optical regime 
observed with any integral field unit (IFU). SATELLITE has been written in python 
and was developed to be an automatic and user-friendly code. The user does not need 
any programming or coding knowledge.

The capabilities and performance of the satellite code (v1.0) were presented using 
the IFU data of the Abell~14 PN obtained with VIMOS@ESO. Satellite v1.3 has been 
applied to three more PNe: Hen 2-108 (VIMOS) (Miranda Marques et al. 2021, submitted), 
NGC 7009 (MUSE), and NGC 6778 (MUSE) (Akras et al. 2021, submitted).

Satellite carries out a spectroscopic analysis of extended ionized nebulae through 1D 
and 2D approach on a list of 35 emission line (the brightest and more frequently detected 
in ionized nebulae) via a number of pseudo-slits that simulate slit spectrometry and 2D 
emission line imaging. The analysis is performed in four different modules:

* (I) rotation analysis,
* (II) radial analysis,
* (III) specific slits analysis and
* (IV) 2D analysis

For each module, SATELLITE computes all the typically used nebular parameters and 
their uncertainties using as input information the available emission lines maps. 
For all four modules, the uncertainties of the line intensities as well as those 
of the nebular parameters are computed following a Monte Carlo approach considering 

* extinction coefficient (c(Hbeta)),
* electron temperatures and densities for different diagnostic lines,
* ionic and elemental abundances, abundances ratio relative to oxygen and the ionization correction factors,
* emission line ratios from a pre-defined list

The input parameters necessary to run the code are provided by the user in a single `yaml` file.
and they are related to:

* emission line and error maps (or additional error as percentage of its pixel flux. )
* the pixel scale of the IFU
* the interstellar extinction law
* number of replicate spectra
* the width, length, position angle (PA) and coordinates of the pseudo-slits for the 1D analysis
* the coordinates of the central star or the centre of the nebula
* atomic data
* Te and Ne from different diagnostic lines for the proper calculation of the ionic abundances
* the line ratios that the code will compute
* the emission line diagnostic diagrams that the code will construct
* the module or modules that the code will execute

The philosophy behind the development of the satellite code is, besides the unique 2D 
imaging spectroscopy that IFU technology provides, to carry out a detailed 1D spectroscopic 
analysis through a number of pseudo-slits that simulate slit spectrometry and emission 
line imaging in order to properly compare the results presented in previous studies. 

The description of each module as well as the input parameters and outcomes are described 
in more details in the following sections using as a representative example the analysis 
of the planetary nebulae NGC 7009 and its MUSE data from the Science Verification phase .

# Installation

To install the software type (in the `ROOT` folder): `pip install -e .`

# Running SATELLITE

## Main Script

The main executable can be found under the `bin` folder. The program come with 
an extensive help message that can be triggered via `-h` or `--help`. Here is 
the output:
```bash
usage: stavros.py [-h] -i CONFIG [-a ATOMIC_DATA_SET] [-e PN_extinction] [-m MONTE_CARLO_FAKE_OBS] [--max-allowed-mc-tries MAX_ALLOWED_MC_TRIES] [-RV R_V] [-EP ENERGY_PARAMETER]
                  [--intensities-out INTENSITIES_OUTPUT_FILE] [--ratios-out RATIOS_OUTPUT_FILE] [--diagnostics-out DIAGNOSTICS_OUTPUT_FILE]
                  [--abundances-out IONIC_ABUNDANCIES_OUTPUT_FILE] [--total-abundances-out TOTAL_ABUNDANCIES_OUTPUT_FILE] [--verbose] [--missing-fits-is-error] [--bar-plots]
                  [--log LOG_FILE] [--no-plots]

Spectroscopic Analysis Tool for intEgraL FieLd unIt daTacubEs

options:
  -h, --help            show this help message and exit
  -i, --config-file CONFIG
                        Configuration file for analysis. (default: None)
  -a, --atomic-data-set ATOMIC_DATA_SET
                        Predefined atomic data set provided by PyNeb. For a complete list, see https://github.com/Morisset/PyNeb_devel/blob/master/docs/Notebooks/PyNeb_manual_3.ipynb (default: PYNEB_21_01)
  -e, --extinction-law PN_extinction
                        Predefined extinction (reddening) correction. See https://notebook.community/Morisset/PyNeb_devel/docs/Notebooks/PyNeb_manual_5 (default: S79 H83 CCM89)
  -m, --monte-carlo-simulations MONTE_CARLO_FAKE_OBS
                        Number of Mone-Carlo fake observations for error estimation. (default: 4)
  --max-allowed-mc-tries MAX_ALLOWED_MC_TRIES
                        Sometimes the Monte Carlo simulated obs result in values that are invalid, producing Nan values for various computations. This value, given as a percentage, allows the reproduction of Monte Carlo simulated values up to some limit w.r.t the MONTE_CARLO_FAKE_OBS. E.g. if MONTE_CARLO_FAKE_OBS=500 and MAX_ALLOWED_MC_TRIES=15, then we can try at maximum 15% of 500 i.e. 75 times. (default: 15)
  -RV R_V               Reddening correction RV parameter for Pyneb; note that R_V = AV/E_BV. Default value is 3.1 (default: 3.1)
  -EP ENERGY_PARAMETER  Energy parameter value for scaling F(Hb). (default: 1.0)
  --intensities-out INTENSITIES_OUTPUT_FILE
                        Output file to write computed inensity results. (default: line_intensities.dat)
  --ratios-out RATIOS_OUTPUT_FILE
                        Output file to write computed (line) ratio results. (default: line_ratios.dat)
  --diagnostics-out DIAGNOSTICS_OUTPUT_FILE
                        Output file to write computed temperature/density diagnostics results. (default: tn_diagnostics.dat)
  --abundances-out IONIC_ABUNDANCIES_OUTPUT_FILE
                        Output file to write computed ionic abundancies results. (default: abundancies.dat)
  --total-abundances-out TOTAL_ABUNDANCIES_OUTPUT_FILE
                        Output file to write computed total abundancies, including ICFs and DIMSs. (default: total_abundancies.dat)
  --verbose             Verbose mode on (default: False)
  --missing-fits-is-error
                        Fail with error if any of the input FITS file specified in the given config file is missing. (default: False)
  --bar-plots           Plot using barplots (excluding total abundancies) instead of scatter plots. Only relevant if '--no-plots' is not set. (default: False)
  --log LOG_FILE        Write log to an output file. (default: None)
  --no-plots            Do no produce plotts. (default: False)

National Observatory of Athens,
    Institute for Astronomy, Astrophysics, Space Applications and Remote Sensing

    Send bug reports to:
    Stavros Akras, xanthos@mail.ntua.gr
    May, 2024
```

## Configuration File

To run SATELLITE you will need to provide a *configuration file* in `yaml` format. A 
sample config file can be found below:

```yaml
geometry:
data_list:
  prefix: /home/xanthos/Software/SATELLITE/data/
  suffix: fit
  # Pre-suffix for error images, example: Ar3_7136e.fit
  error_img_id: e
  # Pre-suffix for data images, example: Ar3_7136s.fit
  data_img_id: s
  element_list:
    - atom: H
      spectrums:
        - spectrum: i
          atomic_numbers:
            - atomic_number: 6563
              ref_TeNe: '[SIII] 6312/9069 [ClIII] 5538/5518'
              # if ref_type is not set, it defaults to 'CEL'
              ref_type: RL
            - atomic_number: 4861
              ref_TeNe: '[SIII] 6312/9069 [ClIII] 5538/5518'
              ref_type: RL
    - atom: He
      spectrums:
        - spectrum: i
          atomic_numbers:
            - atomic_number: 5876
              ref_TeNe: '[NII] 5755/6584+ [SII] 6731/6716'
              ref_type: RL
            - atomic_number: 6678
              ref_TeNe: '[NII] 5755/6584+ [SII] 6731/6716'
              ref_type: RL
        - spectrum: ii
          atomic_numbers:
            - atomic_number: 5412
              ref_TeNe: '[SIII] 6312/9069 [ClIII] 5538/5518'
              ref_type: RL
    - atom: Cl
      spectrums:
        - spectrum: iii
          atomic_numbers:
            - atomic_number: 5517
              ref_TeNe: '[SIII] 6312/9069 [ClIII] 5538/5518'
            - atomic_number: 5538
              ref_TeNe: '[SIII] 6312/9069 [ClIII] 5538/5518'
analysis:
  specific_slit_analysis:
    # to skip specific slit analysis set 'skip' to true
    skip: true
    slits:
      # For each slit enter: PA, width, height, xcrd, ycrd
      # (xcrd, ycrd) here is point of rotation, i.e. center of slit
      - 0, 24, 18, 145, 90
      - 0, 24, 20, 155, 129
  angular_slit_analysis:
    skip: true
    slits:
      # - start_angle, stop_angle, step_angle, width, height, xcrd, ycrd
      # (xcrd, ycrd) here is point of rotation, i.e. center-bottom point
      - 0, 360, 30, 24, 20, 155, 129
  radial_slit_analysis:
    skip: false
    slits: 
      # For each slit enter: PA, width, height, xcrd, ycrd
      # (xcrd, ycrd) here is point of rotation, i.e. center of slit
      - 0, 24, 2, 145, 90
      #- 0, 24, 20, 155, 129
  log_ratios:
    - He1_5876/H1_6563
    - He1_6678/H1_6563
    - He2_4686/H1_4861
    - He2_5412/H1_4861
    - (Cl3_5517+Cl3_5538)/H1_4861
    - C1_8727/H1_6563
    - C2_6461/H1_6563
  density_diagnostics:
    # Note: you'll have to quote lines containing square brackets
    - '[SII] 6731/6716'
    - '[OII] 3726/3729'
    - '[ClIII] 5538/5518'
    - '[ArIV] 4740/4711'
  temperature_diagnostics:
    - '[NII] 5755/6584+'
    - '[OIII] 4363/5007+'
    - '[OI] 5577/6300+'
    - '[OI] 5577/6300'
    - '[OII] 3727+/7325+'
    - '[SIII] 6312/9069'
```

## Supplementary Scripts

### ascii2html_tab.py

This script (located under `bin/` folder), can be used to transform total abundancies 
ascii output (from SATELLITE) to html table(s). It can be useful to turn obfuscated ascii 
tables to human-readable, user friendy tabulated data.

Running `ascii2html_tab.py` is pretty straight-forward; just provide a total abundancies 
data file (produced via SATELLITE) as command line argument.


```bash
usage: ascii2html_tab.py [-h] -b TOTAL_ABUNDANCIES

Spectroscopic Analysis Tool for intEgraL FieLd unIt daTacubEs

options:
  -h, --help            show this help message and exit
  -b, --total-abundancies TOTAL_ABUNDANCIES
                        Total abundancies data file. (default: None)

National Observatory of Athens,
    Institute for Astronomy, Astrophysics, Space Applications and Remote Sensing

    Send bug reports to:
    Stavros Akras, xanthos@mail.ntua.gr
    Nov, 2024

```

### satellite_plot.py

This script (located under `bin/` folder), can be used to plot SATELLITE results. 
Normally, each run of SATELLITE produces a list of (ascii) files as output. These, 
if not plotted at runtime, can be plotted later on using this program.

`satellite_plot.py` can plot any (or all) for SATELLITE result files. Users just need 
to specify the data files to be plotted as command line arguments.

A few plotting options are also available as command line parameters. Users needding 
to further customize plotting, can make minor source code changes.

```bash
usage: satellite_plot.py [-h] [-a LINE_ABUNDANCIES] [-i LINE_INTENCITIES] [-r LINE_RATIOS] [-t TENE_DIAGNOSTICS] [-b TOTAL_ABUNDANCIES] [--bar-plots]

Spectroscopic Analysis Tool for intEgraL FieLd unIt daTacubEs

options:
  -h, --help            show this help message and exit
  -a, --line-abundancies LINE_ABUNDANCIES
                        Line abundancies data file. (default: None)
  -i, --line-intensities LINE_INTENCITIES
                        Line intensities data file. (default: None)
  -r, --line-ratios LINE_RATIOS
                        Line ratios data file. (default: None)
  -t, --tene TENE_DIAGNOSTICS
                        Temperature/Density data file. (default: None)
  -b, --total-abundancies TOTAL_ABUNDANCIES
                        Total abundancies data file. (default: None)
  --bar-plots           Plot using barplots instead of scatter plots. (default: False)

National Observatory of Athens,
    Institute for Astronomy, Astrophysics, Space Applications and Remote Sensing

    Send bug reports to:
    Stavros Akras, xanthos@mail.ntua.gr
    Feb, 2024
```

### SATELLITE_corr_ArIVHeI4712_blended_lines.py

This script corrects the intensity of the `Ar IV λ4712` emission line for the contribution of the `HeI λ4711` emission line using the work by Benjamin et al. (1999)
(https://iopscience.iop.org/article/10.1086/306923/pdf). The script takes as input, single values. Therefore, it is suitable for the Specific Slit and Angular analysis modules.

In order to run the code, assign a value to the variables:
* `I_4712_obs` → observed intensity at λ4712 relative to Hb=100 (8th code line)
* `Te_obs` → electron temperature in K (9th code line)
* `ne_obs` → electron density in cm^-3 (10th code line)

*Developed at May 2024.*

### SATELLITE_corr_NIIOII_for_recombination.py

This script corrects the intensity of the `N II λ5755` and `O II λλ7320,7330` emission lines for the recombination contribution according the works by  Liu et al. (2000)
(https://ui.adsabs.harvard.edu/abs/2000MNRAS.312..585L) and Garcia-Rojas et al. (2021)
 (https://ui.adsabs.harvard.edu/abs/2022MNRAS.510.5444G). The script takes as input, single values. Therefore, it is suitable for the Specific Slit and Angular analysis modules.

In order to run the code, assign a value to the variables:
* `I_5755` → observed intensity at λ5755 relative to Hb=100 (8th code line)
* `I_5679` → observed intensity at λ5679 relative to Hb=100 (9th code line)
* `abund_N2_5679` → Ν+2 ionic abundance ORLs `N II` lines or from CELs `N III` lines (13th code line)
* `Te` → electron temperature in K (15th code line)
* `ne` → electron density in cm^-3 (16th code line)

* `I_7320` →  observed intensity at λ7320 relative to Hb=100 (51th code line)
* `I_7330` →  observed intensity at λ7330 relative to Hb=100 (52th code line)
* `I_4649` →  observed intensity at λ4649 relative to Hb=100 (55th code line)
* `I_4650` →  observed intensity at λ4650 relative to Hb=100 (56th code line)
* `abund_O2` →  O+2 ionic abundance (61th code line)
* `Te` → electron temperature in K (63th code line)
* `ne` → electron density in cm^-3 (64th code line)

*Developed Jul 2024.*

### SATELLITE_TeHeI_recombination.py

The script estimates electron temperature using the `He I` recombination lines according the works from Zhang et al. (2005A) (https://ui.adsabs.harvard.edu/abs/2005MNRAS.358..457Z), Zhang et al. (2005B) (https://ui.adsabs.harvard.edu/abs/2005RMxAC..23...15Z) and Mendez-Delgado et al. (2021) (https://ui.adsabs.harvard.edu/abs/2021MNRAS.502.1703M).  The script takes as input, single values. Therefore, it is suitable for the Specific Slit and Angular analysis modules.

In order to run the code, assign a value to the variables:
* `I_7281` → observed intensity at λ7281 relative to Hb=100 (17th code line)
* `I_6678` → observed intensity at λ6678 relative to Hb=100 (18th code line)
* `ne_input` → electron density in cm^-3 (63th code line)

The code also estimates electron temperature from Balmer jump at λ3645.

In order to run the code, assign a value to the variables:
* `Te` → electron temperature in K (75th code line)
* `ne` → electron density in cm^-3 (76th code line)
* `abund_He1` → He+ ionic abundance (77th code line)
* `abund_He2` → He+2 ionic abundance (78th code line)

The code also estimates electron temperature from Paschen jump at λ8200.

In order to run the code, assign a value to the variables:
* `ne` → electron density in cm^-3 (138th code line)
* `abund_He1` → He+ ionic abundance (139th code line)
* `abund_He2` → He+2 ionic abundance (140th code line)
* `jump` →  the Balmer jump as measured by the user in erg/s/cm^2
* `I_P11` → intensity in erg/s/cm^2 of Paschen 11 emission line

*Developed Jul 2024.*

### SATELLITE_PCA.py

The script uses emission line maps from IFU and applies a principal component analysis (PCA). The script takes as input, 2D maps (fits images). Therefore, it is suitable for the 2D module.

To run the code, assing the full path of a folder with all the emission lines in the following variables:
* image_names and data (8th and 14th lines)

*Developed Feb 2025.*

### SATELLITE_voronoi.py

The script takes a  2D map (fits_file) as an input and creates the corresponding voronoi map (output_fits) given a specific signal-to-noise ratio (based on the fits_file/error_fits ratio). The script takes as input, 2D maps (fits images). Therefore, it is suitable for the 2D module.

To run the code, assign a value to the following variables:  
* target_snr, the requested signal-to-noise ratio (48th line), and

the required input and outout files in the apply_voronoi_binning function: 
* input files: the fits_file (the map for the Voronoi application) and the error_fits (the error of the map) (111st line). 
* the output files: output_fits, output_values_fits, output_sum_fits (113th lines). 

** PS:The script does not always work properly. It depends on the given target_snr **

*Developed May 2025.*

### SATELLITE_TeHeI_2D.py

The script estimates electron temperature from the `He I` recombination lines based on the work by Mendez-Delgado et al. (2021) (https://ui.adsabs.harvard.edu/abs/2021MNRAS.502.1703M).  The script takes as input, 2D maps (fits images). Therefore, it is suitable for the 2D module.

To run the script provide 2D maps in the following variables: 
* `ratio.fits` (2D map of the ratio, 15th line)
* `ne_fits` (2D map of electron density fits, (16th  line)

*Developed Sep 2024.*

### SATELLITE_TePaschen_2D.py

The script estimates electron temperature from the Paschen jump at λ8200.

To run the script, assign a value to the variables:
* `den` →  electron density in cm^-3 (38th line)
* `Hep` →   He+ ionic abundance (39th line)
* `Hepp` →  He+2 ionic abundance (40th line)

and provide the 2D maps (fits files) in the following variables:
* `cont` (2D map cont_s_sized.fits, 41st line) 
* `HI` (Paschen 9 extinction corrected map,  42nd line)

*Developed Sep 2024.*

# Science behind SATELLITE

To read more on the scientific aspects underlying SATELLITE's analysis toolchain, you 
can visit the [original SATELLITE](https://github.com/StavrosAkras/SATELLITE) repository.
