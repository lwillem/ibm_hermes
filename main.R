############################################################################ #
# This file is part of the individual-based modelling framework called HERMES
#
# Goal: main script, i.e. workbench, to run individual-based simulations
#
# This script is distributed in the hope that it will be useful, but without
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #

#' Run a HERMES individual-based model scenario
#'
#' This script serves as a lightweight workbench that:
#' \itemize{
#'   \item loads core model functionality,
#'   \item configures scenario-specific parameters,
#'   \item executes the simulation,
#'   \item optionally runs regression tests.
#' }
#'


# Clear workspace
rm(list = ls())

# Load core model functions
source('lib/ibm_core.R')

# ------------------------------------------------------------------------ -
# PARAMETER CONFIGURATION ----
# ------------------------------------------------------------------------ -

# Central parameter object controlling all model behaviour
params <- get_default_parameters()

# Optional scenario-specific overrides
# print_model_parameters()
params$num_infected_seeds <- 10
params$bool_add_baseline  <- TRUE

# ------------------------------------------------------------------------ -
# RUN MODEL ----
# ------------------------------------------------------------------------ -

ibm_results <- run_ibm(params)


# ------------------------------------------------------------------------ -
# REGRESSION TESTING  ----
# ------------------------------------------------------------------------ -
# Optional regression testing for the baseline setting
run_ibm_regression_test()
