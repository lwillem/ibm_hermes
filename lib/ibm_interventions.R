############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: design and define interventions implemented at different discrete 
#       time points during the unfolding of the epidemic
#
# This script is distributed in the hope that it will be useful, but without
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 sabrams, SIMID, UHASSELT, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #

# ------------------------------------------------------------------------ -
# DEPENDENCIES ----
# ------------------------------------------------------------------------ -

#' Run sequential HERMES individual-based model scenarios

# Clear workspace
rm(list = ls())

# Load core model functions
source('lib/ibm_core.R')

# ------------------------------------------------------------------------ -
# PHASE 1: INITIAL UNRESTRICTED DISEASE SPREAD ----
# ------------------------------------------------------------------------ -

# ------------------------------------------------------------------------ -
# PARAMETER CONFIGURATION (INITIAL UNRESTRICTED SPREAD) ----
# ------------------------------------------------------------------------ -

# Central parameter object controlling all model behaviour
params <- get_default_parameters()

# Optional scenario-specific overrides
print_model_parameters(params)

# New settings
params$output_dir <- "output/ibm_modified"
params$num_days <- 20
params$num_infected_seeds <- 1
params$bool_add_baseline  <- TRUE
params$general_mortality_rate <- NULL
params$pop_size <- 1e5

params$vaccine_coverage <- 0
params$vaccine_effectiveness <- 0

params$transmission_prob <- 0.08

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

# ------------------------------------------------------------------------ -
# PHASE 2: FULL LOCKDOWN IMPOSED OVER 7 DAY PERIOD ----
# ------------------------------------------------------------------------ -

# ------------------------------------------------------------------------ -
# PARAMETER CONFIGURATION 1 (LOCKDOWN IMPOSED OVER 7 DAY PERIOD) ----
# ------------------------------------------------------------------------ -

new_params1 <- params

new_params1$num_days <- 5
new_params1$num_contacts_community_day <- 5
new_params1$contact_prob_household <- 1.0
new_params1$contact_prob_school <- 0.1
new_params1$contact_prob_workplace <- 0.01

# ------------------------------------------------------------------------ -
# RUN IBM MODEL AGAIN ----
# ------------------------------------------------------------------------ -

new_ibm_results <- rerun_ibm(pop_data = pop_data, 
                             prev_log_health_matrix = health_time_data$log_health, 
                             params = new_params1, 
                             intervention = 1)

# output is stored in: 
new_ibm_results$params$output_dir

# explore population output
pop_data_file <- file.path(new_ibm_results$params$output_dir,'pop_data_intervention1.rds')
pop_data <- readRDS(pop_data_file)
dim(pop_data)
head(pop_data)

# explore health states output
health_time_file <- file.path(new_ibm_results$params$output_dir,'health_time_intervention1.rds')
health_time_data <- readRDS(health_time_file)
dim(health_time_data)
head(health_time_data$log_health)

# ------------------------------------------------------------------------ -
# PARAMETER CONFIGURATION 2 (LOCKDOWN IMPOSED OVER 7 DAY PERIOD) ----
# ------------------------------------------------------------------------ -

new_params2 <- new_params1

new_params1$num_days <- 15
new_params1$num_contacts_community_day <- 5
new_params1$contact_prob_household <- 1.0
new_params1$contact_prob_school <- 0.01
new_params1$contact_prob_workplace <- 0.01

# ------------------------------------------------------------------------ -
# RUN IBM MODEL AGAIN ----
# ------------------------------------------------------------------------ -

new_ibm_results <- rerun_ibm(pop_data = pop_data, 
                             prev_log_health_matrix = health_time_data$log_health, 
                             params = new_params2, 
                             intervention = 2)

# output is stored in: 
new_ibm_results$params$output_dir

# explore population output
pop_data_file <- file.path(new_ibm_results$params$output_dir,'pop_data_intervention2.rds')
pop_data <- readRDS(pop_data_file)
dim(pop_data)
head(pop_data)

# explore health states output
health_time_file <- file.path(new_ibm_results$params$output_dir,'health_time_intervention2.rds')
health_time_data <- readRDS(health_time_file)
dim(health_time_data)
head(health_time_data$log_health)

# ------------------------------------------------------------------------ -
# PHASE 3: EXIT STRATEGY IMPOSED OVER 7 DAY PERIOD ----
# ------------------------------------------------------------------------ -

params$num_days <- 10

# ------------------------------------------------------------------------ -
# RUN IBM MODEL AGAIN ----
# ------------------------------------------------------------------------ -

new_ibm_results <- rerun_ibm(pop_data = pop_data, 
                             prev_log_health_matrix = health_time_data$log_health, 
                             params = params, 
                             intervention = 3)

# output is stored in: 
new_ibm_results$params$output_dir

# explore population output
pop_data_file <- file.path(new_ibm_results$params$output_dir,'pop_data_intervention3.rds')
pop_data <- readRDS(pop_data_file)
dim(pop_data)
head(pop_data)

# explore health states output
health_time_file <- file.path(new_ibm_results$params$output_dir,'health_time_intervention3.rds')
health_time_data <- readRDS(health_time_file)
dim(health_time_data)
head(health_time_data$log_health)

# ------------------------------------------------------------------------ -
# PHASE 4: FULL VACCINATION IMPOSED OVER 14 DAY PERIOD ----
# ------------------------------------------------------------------------ -

# ------------------------------------------------------------------------ -
# PARAMETER CONFIGURATION 3 (VACCINATION IMPOSED OVER 7 DAY PERIOD) ----
# ------------------------------------------------------------------------ -

new_params3 <- params

new_params3$num_days <- 5
new_params3$vaccine_coverage <- 0.05
new_params3$vaccine_effectiveness <- 0.8

# ------------------------------------------------------------------------ -
# RUN IBM MODEL AGAIN ----
# ------------------------------------------------------------------------ -

new_ibm_results <- rerun_ibm(pop_data = pop_data, 
                             prev_log_health_matrix = health_time_data$log_health, 
                             params = new_params3, 
                             intervention = 4)

# output is stored in: 
new_ibm_results$params$output_dir

# explore population output
pop_data_file <- file.path(new_ibm_results$params$output_dir,'pop_data_intervention4.rds')
pop_data <- readRDS(pop_data_file)
dim(pop_data)
head(pop_data)

# explore health states output
health_time_file <- file.path(new_ibm_results$params$output_dir,'health_time_intervention4.rds')
health_time_data <- readRDS(health_time_file)
dim(health_time_data)
head(health_time_data$log_health)

# ------------------------------------------------------------------------ -
# PARAMETER CONFIGURATION 4 (VACCINATION IMPOSED OVER 7 DAY PERIOD) ----
# ------------------------------------------------------------------------ -

new_params4 <- params

new_params4$num_days <- 20
new_params4$vaccine_coverage <- 0.10
new_params4$vaccine_effectiveness <- 0.8

# ------------------------------------------------------------------------ -
# RUN IBM MODEL AGAIN ----
# ------------------------------------------------------------------------ -

new_ibm_results <- rerun_ibm(pop_data = pop_data, 
                             prev_log_health_matrix = health_time_data$log_health, 
                             params = new_params4, 
                             intervention = 5)

# output is stored in: 
new_ibm_results$params$output_dir

# explore population output
pop_data_file <- file.path(new_ibm_results$params$output_dir,'pop_data_intervention5.rds')
pop_data <- readRDS(pop_data_file)
dim(pop_data)
head(pop_data)

# explore health states output
health_time_file <- file.path(new_ibm_results$params$output_dir,'health_time_intervention5.rds')
health_time_data <- readRDS(health_time_file)
dim(health_time_data)
head(health_time_data$log_health)