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
