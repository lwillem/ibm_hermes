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

# New settings
params$output_dir <- "output/ibm_modified"
params$num_days <- 20
params$num_infected_seeds <- 10
params$bool_add_baseline  <- TRUE
params$general_mortality_rate <- NULL
params$pop_size <- 1e4

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
# VISUALIZE INCUBATION PERIOD DISTRIBUTION ----
# ------------------------------------------------------------------------ -

# visualization of continuous antibody titer data 
p1 <- ggplot(subset(pop_data, !is.na(symptom_onset)), aes(x = symptom_onset)) + 
  geom_bar() + 
  labs(
    x = "Symptom onset (days since infection)",
    y = "Number of individuals",
    title = "Incubation period"
  ) +
  theme_minimal()

fig_dir <- file.path(ibm_results$params$output_dir, 'figures')

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

# Save plot
ggsave(
  filename = file.path(fig_dir, paste0("incubation_period", ".png")),
  plot = p1,
  width = 6,
  height = 4
)

# ------------------------------------------------------------------------ -
# DERIVE INCIDENCE DATA ----
# ------------------------------------------------------------------------ -

# Central parameter object for delayed incidence data generation
delay_params <- get_default_delay_parameters()

# Incidence data generation with delays
incidence_results <- obtain_incidence_data(pop_data, delay_params)

# explore contact tracing data output
incidence_data_file <- file.path(incidence_results$params$output_dir,'incidence_data.rds')
incidence_data <- readRDS(incidence_data_file)
dim(incidence_data)
head(incidence_data)

# ------------------------------------------------------------------------ -
# VISUALIZE INCIDENCE DATA ----
# ------------------------------------------------------------------------ -

# visualization of the incidence of cases in relation to time of reporting 
tab <- table(incidence_data$time_of_reporting)
incidence_plot_df <- data.frame(time = sort(unique(incidence_data$time_of_reporting)),
                                new_cases = as.vector(tab))

library(ggplot2)
p1 <- ggplot(incidence_plot_df, 
             aes(x = time, 
                 y = new_cases)) +
  geom_point() + 
  labs(
    x = "Time (in days since start of outbreak)",
    y = "Number of new cases",
    title = "Total reported incidence"
  ) +
  theme_minimal()

fig_dir <- file.path(incidence_results$params$output_dir, 'figures')

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

# Save plot
ggsave(
  filename = file.path(fig_dir, "total_reported_incidence_time.png"),
  plot = p1,
  width = 6,
  height = 4
)


# ------------------------------------------------------------------------ -
# DERIVE CONTACT TRACING DATA ----
# ------------------------------------------------------------------------ -

# Central parameter object for contact tracing
contact_tracing_params <- get_default_contact_tracing_parameters()

# Contact tracing data generation
contact_tracing_results <- sample_contact_tracing_data(pop_data, contact_tracing_params, verbose)

# explore contact tracing data output
contact_tracing_data_file <- file.path(contact_tracing_results$params$output_dir,'contact_tracing_data.rds')
contact_tracing_data <- readRDS(contact_tracing_data_file)
dim(contact_tracing_data)
head(contact_tracing_data)

# ------------------------------------------------------------------------ -
# DERIVE SEROLOGY ----
# ------------------------------------------------------------------------ -

# Central parameter object for serology generation
sero_params <- get_default_sero_parameters()

# Serological data generation
sero_results <- sample_serological_data(pop_data, sero_params, verbose)

# explore serological data output
sero_data_file <- file.path(sero_results$params$output_dir,'sero_data.rds')
sero_data <- readRDS(sero_data_file)
dim(sero_data)
head(sero_data)

# ------------------------------------------------------------------------ -
# VISUALIZE SEROLOGICAL SURVEY DATA ----
# ------------------------------------------------------------------------ -

# visualization of the seroprevalence by age 
sero_plot_df <- data.frame(age = sort(unique(round(sero_data$age))),
                           seroprev = tapply(sero_data$sero_status, 
                                              round(sero_data$age), 
                                              mean))

library(ggplot2)
p1 <- ggplot(sero_plot_df, 
             aes(x = age, 
                 y = seroprev)) +
  geom_point() + 
  labs(
    x = "Age (in years)",
    y = "Proportion of seropositive individuals",
    title = "Seroprevalence by age"
  ) +
  theme_minimal()

fig_dir <- file.path(sero_results$params$output_dir, 'figures')

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

# Save plot
ggsave(
  filename = file.path(fig_dir, paste0("seroprevalence_time",sero_results$params$sampling_time,".png")),
  plot = p1,
  width = 6,
  height = 4
)

# visualization of the seroprevalence by age group 
seroprev_by_age_group <- tapply(sero_data$sero_status, 
                                sero_data$age_group, 
                                mean)
npos_by_age_group <- tapply(sero_data$sero_status, 
                            sero_data$age_group, 
                            sum)
n_by_age_group <- tapply(sero_data$sero_status, 
                         sero_data$age_group, 
                         length)

z <- 1.96
denominator <- 1 + ((z**2)/n_by_age_group)
center_adj <- seroprev_by_age_group + (z**2/(2*n_by_age_group))
std_adj = sqrt((seroprev_by_age_group*(1 - seroprev_by_age_group)/n_by_age_group) + 
                 (z**2/(4*n_by_age_group**2)))

seroprev_by_age_group_ll = (center_adj - z*std_adj)/denominator
seroprev_by_age_group_ul = (center_adj + z*std_adj)/denominator

sd_seroprev_by_age_group = sqrt((seroprev_by_age_group*(1 - seroprev_by_age_group))/n_by_age_group)
seroprev_by_age_group_ll_asymp = seroprev_by_age_group - z*sd_seroprev_by_age_group
seroprev_by_age_group_ul_asymp = seroprev_by_age_group + z*sd_seroprev_by_age_group

sero_plot_df2 <- data.frame(age_group = names(seroprev_by_age_group),
                            seroprev = seroprev_by_age_group,
                            seroprev_ll = seroprev_by_age_group_ll,
                            seroprev_ul = seroprev_by_age_group_ul,
                            seroprev_ll_asymp = seroprev_by_age_group_ll_asymp,
                            seroprev_ul_asymp = seroprev_by_age_group_ul_asymp)

p2 <- ggplot(sero_plot_df2, 
             aes(x = factor(age_group), 
                 y = seroprev)) +
  geom_point() + 
  labs(
    x = "Age (in years)",
    y = "Proportion of seropositive individuals",
    title = "Seroprevalence by age group"
  ) +
  geom_errorbar(aes(ymin=seroprev_ll, ymax=seroprev_ul), width=.2,
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=seroprev_ll_asymp, ymax=seroprev_ul_asymp), width=.2, col = "orange",
                position=position_dodge(.9)) +
  theme_minimal() 

# Save plot
ggsave(
  filename = file.path(fig_dir, paste0("seroprevalence_age_group_time",sero_results$params$sampling_time,".png")),
  plot = p2,
  width = 6,
  height = 4
)

# visualization of continuous antibody titer data 
p3 <- ggplot(sero_data, aes(x = log_igg)) + 
  geom_histogram() + 
  labs(
    x = "Log IgG titer",
    y = "Number of seropositive individuals",
    title = "Antibody titer data"
  ) +
  theme_minimal()

# Save plot
ggsave(
  filename = file.path(fig_dir, paste0("antibody_titer_data_time",sero_results$params$sampling_time,".png")),
  plot = p3,
  width = 6,
  height = 4
)

# ------------------------------------------------------------------------ -
# REGRESSION TESTING  ----
# ------------------------------------------------------------------------ -
# Optional regression testing for the baseline setting
run_ibm_regression_test()

