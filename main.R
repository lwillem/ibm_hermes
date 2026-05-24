############################################################################ #
# This file is part of the individual-based modelling framework called HERMES
#
# Goal: main script, i.e. workbench, to run individual-based simulations
#
# This script is distributed in the hope that it will be useful, but without
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
#                    sabrams, SIMID, UHASSELT, UNIVERSITY OF ANTWERP, BELGIUM
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
print_model_parameters(params)
params$output_dir <- output/ibm_modified
params$num_days <- 20
params$num_infected_seeds <- 10
params$bool_add_baseline  <- TRUE
params$pop_size <- 1e5

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
# VISUALIZE MODEL OUTPUT ----
# ------------------------------------------------------------------------ -

# visualization of the evolution of susceptible individuals over time
library(ggplot2)
p1 <- ggplot(health_time_data$log_health, 
       aes(x = seq_len(nrow(health_time_data$log_health)), 
           y = health_time_data$log_health$S)) +
  geom_line() + 
  labs(
    x = "Time (in days)",
    y = "Proportion of susceptible individuals",
    title = "Evolution of the proportion of susceptible individuals"
  ) +
  theme_minimal()


fig_dir <- file.path(ibm_results$params$output_dir, 'figures')

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

# Save plot
ggsave(
  filename = file.path(fig_dir, "evolution_susceptible.png"),
  plot = p1,
  width = 6,
  height = 4
)

# visualization of the evolution of infected individuals over time
p2 = ggplot(health_time_data$log_health, 
       aes(x = seq_len(nrow(health_time_data$log_health)), 
           y = health_time_data$log_health$I)) +
  geom_line() + 
  labs(
    x = "Time (in days)",
    y = "Proportion of infected individuals",
    title = "Evolution of the proportion of infected individuals"
  ) +
  theme_minimal()

# Save plot
ggsave(
  filename = file.path(fig_dir, "evolution_infected.png"),
  plot = p2,
  width = 6,
  height = 4
)

# ------------------------------------------------------------------------ -
# DERIVE SEROLOGY ----
# ------------------------------------------------------------------------ -

# explore serological data output

# ------------------------------------------------------------------------ -
# REGRESSION TESTING  ----
# ------------------------------------------------------------------------ -
# Optional regression testing for the baseline setting
run_ibm_regression_test()

