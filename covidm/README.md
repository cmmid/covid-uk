# covidm
Dynamic model of SARS-nCoV-2 transmission
Local (ngdavies) version from 30 Mar 2020, for UK scenario projections.


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

