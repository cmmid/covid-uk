# covidm.R
# main file to source for users of the covidm API.

# TODO
#  - recover better from user halting of execution of cm_fit
#  - recover better from crashes during cm_fit
# TODO
#  - make sure seq_len and seq_along are used instead of 1:length() for 0-length things

packageStartupMessage("Loading covidm.")

# Required libraries
suppressPackageStartupMessages({
    library(data.table)   # for data.table, an enhanced (and faster) data.frame
    library(ggplot2)      # for plotting
    library(Rcpp)         # for running the C++ model backend
    library(qs)           # for qsave and qread, faster equivalents of saveRDS and readRDS
    library(lubridate)    # for manipulating dates and times. NB requires stringr
    library(nloptr)       # for numerical optimization
    library(HDInterval)   # for summarizing results
})

# Settings for covidm
cm_option_ = function(name, default_value) {
    if (!exists(name)) {
        packageStartupMessage("Using default value for ", name, ": ", default_value)
    }
    return (get0(name, ifnotfound = default_value))
}

cm_path_           = cm_option_("cm_path", "~/Dropbox/nCoV/covidm/")
cm_force_rebuild_  = cm_option_("cm_force_rebuild", F)
cm_build_dir_      = cm_option_("cm_build_dir", paste0(cm_path_, "/build/"))
cm_build_verbose_  = cm_option_("cm_build_verbose", T)

# Attach R code
source(paste0(cm_path_, "/R/covidm_misc.R"))
source(paste0(cm_path_, "/R/covidm_run.R"))
source(paste0(cm_path_, "/R/covidm_params.R"))
source(paste0(cm_path_, "/R/covidm_interventions.R"))
source(paste0(cm_path_, "/R/covidm_fit.R"))
source(paste0(cm_path_, "/R/covidm_plot.R"))

# Build C++ code
packageStartupMessage("Attaching C++ code...")
sourceCpp(paste0(cm_path_, "/model_v1/corona.cpp"), 
          rebuild  = cm_force_rebuild_,
          cacheDir = cm_build_dir_,
          verbose  = cm_build_verbose_)
sourceCpp(paste0(cm_path_, "/fit_v1/fit.cpp"),
          rebuild  = cm_force_rebuild_, 
          cacheDir = cm_build_dir_,
          verbose  = cm_build_verbose_)

# Attach data
cm_matrices     = readRDS(paste0(cm_path_, "/data/all_matrices.rds"));
cm_populations  = readRDS(paste0(cm_path_, "/data/wpp2019_pop2020.rds"));
cm_structure_UK = readRDS(paste0(cm_path_, "/data/structure_UK.rds"));

