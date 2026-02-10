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
# params$num_infected_seeds <- 10
# params$bool_add_baseline  <- TRUE
# params$pop_size <- 1e4

# ------------------------------------------------------------------------ -
# RUN MODEL ----
# ------------------------------------------------------------------------ -

ibm_results <- run_ibm(params)

# output is stored in: 
ibm_results$params$output_dir

# explore population output
pop_data_file <- file.path(ibm_results$params$output_dir,'pop_data.rds')
pop_data <- readRDS(pop_data_file)
dim(pop_data)
head(pop_data)

# explore health states output
health_time_file <- file.path(ibm_results$params$output_dir,'health_time.rds')
health_time_data <- readRDS(health_time_file)
dim(health_time_data)
head(health_time_data$log_health)

# visualization of the evolution in health states over time
library(ggplot2)
ggplot(health_time_data$log_health, 
       aes(x = seq_len(nrow(health_time_data$log_health)), 
           y = health_time_data$log_health$I)) +
  geom_line() + 
  labs(
    x = "Time (in days)",
    y = "Proportion of infected individuals",
    title = "Evolution of the proportion of infected individuals"
  ) +
  theme_minimal()

# ------------------------------------------------------------------------ -
# REGRESSION TESTING  ----
# ------------------------------------------------------------------------ -
# Optional regression testing for the baseline setting
run_ibm_regression_test()

