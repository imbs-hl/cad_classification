
source("../init.R", chdir = TRUE)

##============================================================================+
## load R packages ----
##============================================================================+
if (!require(pacman)) {
  install.packages("pacman")
}
pacman::p_load(batchtools)
pacman::p_load(R.utils)
pacman::p_load(checkmate)
pacman::p_load(mlr)
pacman::p_load(mlrMBO)
pacman::p_load(devtools)

##============================================================================+
## setup batchtools ----
##============================================================================+
slurmtmpl_file <- file.path(getwd(), "batchtools_config/slurm.batchtools.tmpl")
slurm_btconf_file <- file.path(getwd(), "batchtools_config/batchtools.slurm.R")
interactive_btconf_file <- file.path(getwd(), "batchtools_config/batchtools.interactive.R")
multicore_btconf_file <- file.path(getwd(), "batchtools_config/batchtools.multicore.R")

##============================================================================+
## custom functions ----
##============================================================================+
R.utils::sourceDirectory(file.path(scripts_dir, "functions"))

##============================================================================+
## learner parameter sets ----
##============================================================================+
R.utils::sourceDirectory(file.path(scripts_dir, "learner_parameter_sets"))

##============================================================================+
## custom packages ----
##============================================================================+
#devtools::install(pkg = "packages/gwasFilter")
#devtools::install(pkg = "packages/importanceFilter")
pacman::p_load(gwasFilter)
pacman::p_load(importanceFilter)
