# covid-uk

Stochastic age-structured model of SARS-nCoV-2 transmission for UK scenario projections.

## Quick start guide

### Installing dependencies for Mac OS

You will need to install gfortran binaries from here: https://github.com/fxcoudert/gfortran-for-macOS/releases

Once installed, run `gcc --version` in terminal to get your current version, e.g. `Target: x86_64-apple-darwin18.8.2.0`. Then run below in terminal to add library path for R:

`cd ~
mkdir .R
cd .R
echo FLIBS=-L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin18/8.2.0 -L/usr/local/gfortran/lib -lgfortran -lquadmath -lm >> Makevars
`

Finally, install nlopt: `brew install nlopt`

### Guide to files

Main parameter setting and model run script is in `UK.R` – there is option to set local path at top. Output collation and plotting functions are in `UK-view.R`. Underlying model code is in `covidm` folder.

### Reference

Davies NG et al. The effect of non-pharmaceutical interventions on COVID-19 cases, deaths and demand for hospital services in the UK: a modelling study. CMMID COVID-19 working group pre-print, 2020.