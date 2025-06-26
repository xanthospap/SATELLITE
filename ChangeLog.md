# Version 1.3b1

- [x] PA angle (in config) now interpreted as clock-wise angle.
- [x] Precision/error estimates fixed.
- [x] Removed titles from all plots. Titles are now the y-axis labeles.
- [x] All y-axis ticks are now in scientific format.
- [x] Changed default matplotlib style-sheet to 'default'.
- [x] Optionally plot everything with barplots, except from total aundancies. In `bin/satellite_plot.py` bar plots are triggered with `--bar-plots` and in `bin/stavros.py` with `--bar-plots`.

# Version 1.2b1

- [x] Ratios are orinted as log10.
- [x] Slit indexing starts at 1 (not 1) for all output files.
- [x] Abundance not printed for H1r_4861A and H1r_6563A (abundancies.dat).
- [x] Total abundance not printed for element H (total_abundancies.dat).
- [x] Added RV as a command line parameter with a default value of 3.1. Help meesage now includes: `-RV R_V               Reddening correction RV parameter for Pyneb; note that R_V = AV/E_BV. Default value is 3.1 (default: 3.1)`.
- [x] Added uncertainty for cHbeta (line_intensities.dat).
- [x] Corrected uncertainty computation in ionic abundancies (abundancies.dat).
- [x] Added function computeIcfsWithErrors to compute ICFs uncertainties.
- [x] Total abundancies now contain uncertainties (total_abundancies.dat).
- [x] Added plotters for all output files (source at src/satellite/plotters.py). 

# Version 1.1b1

- [x] Added source code file nomenclature.py
- [x] Added source code ratios.py. All ratio handling moved there.
- [x] "DIMS14_32b" to "DIMS14_32".
- [x] Extinction law can be given as command line parameter: e.g. as -e "S79 H83 CCM89" or --extinction-law="S79 H83 CCM89".
- [x] cHbeta (and E_BV) now printed at intensities.out (last line).
- [x] All output files are now printed sorted alphabeticaly.
- [ ] Print F(Hb): **where?**
- [x] Allow changing of atomic data file; can be done via command line: e.g. as -a 'PYNEB_23_01' or --atomic-data-set='PYNEB_23_01'.
- [x] Allow setting number of monte carlo simulations via command line: e.g. as -m 12 or --monte-carlo-simulations=12.
- [x] New output for line intensities (line_intensities.dat).
- [x] New output for line ratios (line_ratios.dat).
- [x] Computation of Temperature/Density diagnostics changed. User supplied list is flitered based on observation set. This removes excess work and Pyneb warnings.
- [x] New output for ionic abundancies (abundancies.dat) 
- [ ] Optimize output for total_abundancies.dat
