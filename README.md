# Lauderdale_and_Cael_Exports_Draft_Manuscript
<!---
[![DOI](https://zenodo.org/badge/207910435.svg)](https://zenodo.org/badge/latestdoi/207910435)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/seamanticscience/Lauderdale_etal_2020_PNAS?color=1b3370)
![GitHub last commit](https://img.shields.io/github/last-commit/seamanticscience/Lauderdale_etal_2020_PNAS?color=f44323)
![GitHub License](https://img.shields.io/github/license/seamanticscience/Lauderdale_etal_2020_PNAS?color=ffa500)
<a href="https://doi.org/10.1073/pnas.1917277117"><img src="http://img.shields.io/badge/paper%20link-doi:10.1073%2Fpnas.1917277117-lightgrey.svg" alt="Link to paper at https://doi.org/10.1073/pnas.1917277117"></a>
---!>

Model  code, output, and analysis routines (to be updated) looking at the impact of remineralization profile shape on the air-sea carbon balance, and characterizing structural uncertainty in the ocean's biological pump.

Generate all the coefficients for the six functional forms statistically fit to the reference power-law curve in three different ways using the MATLAB routine `profile_coefficients.m`. Coefficients are stored in `export_profile_coefficients.csv`, with the coefficients used in the manuscript supplied here.

Simulations were run with the "Checkpoint63m" version of MITgcm, substituting `pkg DIC` files `DIC_VARS.h`, `dic_readparms.F`, and `phos_flux.F` (which contains the different parameterizations), and adding a function for exponential integral, or upper incomplete gamma function `expint.F` [JML needs updated list of code for the preformed tracers]. The functional forms can each be activated by setting `selectExportRemin` between 1--7 in `data.dic`, where:

• `selectExportRemin=1` for the simple exponential profile,

• `selectExportRemin=2` for the power-law curve ,

• `selectExportRemin=3` for the ballast profile,

• `selectExportRemin=4` for the rational profile,

• `selectExportRemin=5` for the double exponential profile,

• `selectExportRemin=6` for the stretched exponential profile,

• `selectExportRemin=7` for the gamma function profile.

The Jupyter Notebook `export_flux_analysis.ipynb` contains the analysis routines used to generate the figures [JML needs updating with relative paths to MWE model output uploaded here.]

Any questions or comments, please get in contact!
