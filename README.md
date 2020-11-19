# Lauderdale_and_Cael_Exports_Draft_Manuscript
Model  code, output, and analysis routines (to be updated) looking at the impact of remineralization profile shape on the air-sea carbon balance, and characterizing structural uncertainty in the ocean's biological pump.

Generate all the coefficients for the six functional forms statistically fit to the reference power-law curve in three different ways using the MATLAB routine `profile_coefficients.m`. Coefficients are stored in `export_profile_coefficients.csv`, with the coefficients used in the manuscript supplied here.

Simulations were run with the "Checkpoint63m" version of MITgcm, substituting `pkg DIC` files `DIC_VARS.h`, `dic_readparms.F`, and `phos_flux.F` (which contains the different parameterizations), and adding a function for exponential integral, or upper incomplete gamma function `expint.F`. The functional forms can each be activated by setting `selectExportRemin` between 1--7 in `data.dic`, where:

• `selectExportRemin=1` for the simple exponential profile,

• `selectExportRemin=2` for the power-law curve ,

• `selectExportRemin=3` for the ballast profile,

• `selectExportRemin=4` for the rational profile,

• `selectExportRemin=5` for the double exponential profile,

• `selectExportRemin=6` for the stretched exponential profile,

• `selectExportRemin=7` for the gamma function profile.

The Jupyter Notebook `export_flux_analysis.ipynb` contains the analysis routines used to generate the figures.

Any questions or comments, please get in contact!

JML To Do:
- [X] file description
- [X] updated list of code for preformed tracers
- [ ] add MITgcm experiment input files
- [ ] add MITgcm ouput last model decadal average for ptracers, dic, and diagnostics
- [ ] update Jupyter notebook with relative paths to MWE model output uploaded here
