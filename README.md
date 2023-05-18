# research_codes

These are some scripts that I've written for seismology research. I'm making them available to those who may be able to use them. Please note that since they are research codes, they may contain bugs and/or require additioanl effort to adapt to your specific needs.

Empirical Green's Functions
- EGFcode.zip: This is a set of MatLab scripts for doing EGF analysis of earthquake waveforms. These scripts were written and used with a Matlab version before 2015 (2012, I think).

Magnitude of Completeness (Mc)
- Mc_analysis.m: This script can calculate Mc of a catalog with 3 different methods for a given year, span of years, or set of years.
- Mc_analysis_sequence.m: This script is similar to the above but can use specific dates for more control (e.g., to examine a specific earthquake sequence).
- Mc_analysis_grided.m: This script uses the maximum curvature method to calculate Mc for spatial bins for a given year, span of years, or set of years.
- import_scsn_files.m: This scripts converts SCSN yearly text files (https://service.scedc.caltech.edu/ftp/catalogs/SCEC_DC/) to the array format needed for the Mc scripts.

*Note: These scripts were written and tested with Matlab 2022. The grided version uses the mapping toolbox for plotting.

If you use these scripts for your research, please provide appropriate citation/acknowledgement.
