# Lauderdale and Cael "Impact of remineralization profile shape on the air-sea carbon balance" manuscript repository.
[![DOI](https://zenodo.org/badge/207910435.svg)](https://zenodo.org/badge/latestdoi/207910435)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/seamanticscience/Lauderdale_and_Cael_2021_GRL?color=1b3370)
![GitHub last commit](https://img.shields.io/github/last-commit/seamanticscience/Lauderdale_and_Cael_2021_GRL?color=f44323)
![GitHub License](https://img.shields.io/github/license/seamanticscience/Lauderdale_and_Cael_2021_GRL?color=ffa500)
<a href="https://doi.org/10.1029/2020GL091746"><img src="http://img.shields.io/badge/paper%20link-doi:10.1029%2F2020GL091746-lightgrey.svg" alt="Link to paper at https://doi.org/10.1029/2020GL091746"></a>

Model  code, output, and analysis routines looking at the impact of remineralization profile shape (i.e. functional form) on the air-sea carbon balance, and characterizing structural uncertainty in the ocean's biological pump. Read the pre-print on the [Earth and Space Science Open Archive](https://www.essoar.org/doi/abs/10.1002/essoar.10504824.1).

## Fitting the remineralization profiles:
Generate all the coefficients for the six functional forms statistically fit to the reference power-law curve in three different ways using the MATLAB routine `profile_coefficients.m`. Coefficients are stored in `export_profile_coefficients.csv`, with the coefficients used in the manuscript supplied here.

## Numerical model code and configuration:
Simulations were run with the "Checkpoint63m" version of MITgcm, substituting files in the `code_mods` directory. `pkg DIC` files `DIC_VARS.h`, `dic_readparms.F`, and `phos_flux.F` supply the different remineralization functions (with an added function for exponential integral, or upper incomplete gamma function `expint.F`). The functional forms can each be activated by setting `selectExportRemin` between 1--7 in `data.dic`, where:

• `selectExportRemin=1` for the simple exponential profile,

• `selectExportRemin=2` for the power-law curve,

• `selectExportRemin=3` for the ballast profile,

• `selectExportRemin=4` for the rational profile,

• `selectExportRemin=5` for the double exponential profile,

• `selectExportRemin=6` for the stretched exponential profile,

• `selectExportRemin=7` for the gamma function profile.

Preformed tracers, which aid in partitioning carbon and nutrients into "physically-preformed" and "biologically-regenerated" nutrients are included in the `code_mods` directory, by substituting the default `pkg GCHEM` files `GCHEM_OPTIONS.h` and `gchem_forcing_sep.F`, and adding `gchem_preformed_tracers.F`. Activate the tracers by increasing the number of tracers associated with `pkg PTRACERS`, and include definitions in `data.ptracers` file (see experiment input files for details). 

The same executable was used for each simulation. Coefficients for each remineralization function are passed using the four-value `KRemin` array in `data.dic` (filled with zeroes tot he right where needed).

## Numerical model simulations:
Each simulation is associated with its own folder in this repository, which includes input files and steady-state output as an average of the last 100 years of the siumulation. The naming convention is:

• `mar0.70`, `mar0.84`, or `mar0.98` for the reference simulations using the Martin Curve with _b_ values of 0.70, 0.84 (reference), and 0.98.

• `abs`, `rel`, or `efd` for fitting to the reference profile by minimizing absolute error, relative error, or e-folding depth of remineralization.

• `exp`, `bal`, `dbl`, `str`, `rat`, or `gam` for the remineralization profile used (simple exponential, ballast, double exponential, stretched exponential, rational, and gamma functions).

• `noflux` for the simulation where no particulate organic carbon is produced at the surface, and 100% of production is channeled instead to dissolved organic carbon that degrades with a timescale of 6 months.

## Analysis:
The Jupyter Notebook `export_flux_analysis.ipynb` contains the analysis routines used to generate the figures. Dependencies include `xarray`, `dask`, and `mitgcm_tools` ([get from github](https://github.com/seamanticscience/mitgcm_tools)).

Any questions or comments, please get in contact!
