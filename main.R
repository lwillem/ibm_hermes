############################################################################ #
# This file is part of the individual-based modelling framework called HERMES
#
# Goal: main script, i.e. workbench, to run individual-based simulations
#
# This script is distributed in the hope that it will be useful,but without 
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #

# clear workspace
rm(list = ls())

# load script with functions
source('lib/ibm_core.R')

# load default parameters
param <- get_default_parameters()

# option to change model parameters, to list all options, run `print_model_parameters()`
# print_model_parameters()
# e.g. param$num_infected_seeds <- 10
# e.g. param$bool_add_baseline <- TRUE
param$num_infected_seeds <- 10
param$bool_add_baseline <- TRUE

# run the IBM
ibm_out <- run_ibm_location(param)

# regression testing
source('lib/ibm_test.R')
run_ibm_regression_test(overrule_param_check = TRUE)
