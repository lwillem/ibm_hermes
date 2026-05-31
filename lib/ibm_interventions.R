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
params$num_days <- 15
params$num_infected_seeds <- 10
params$bool_add_baseline  <- TRUE
params$general_mortality_rate <- NULL
params$pop_size <- 1e5

params$vaccine_coverage <- 0
params$vaccine_effectiveness <- 0

params$transmission_prob <- 0.07

params$output_dir = "output/ibm_final"

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
# PHASE 2: FULL LOCKDOWN IMPOSED OVER 7 DAY PERIOD (DAY 15) ----
# ------------------------------------------------------------------------ -

# ------------------------------------------------------------------------ -
# PARAMETER CONFIGURATION 1 (LOCKDOWN IMPOSED OVER 7 DAY PERIOD) ----
# ------------------------------------------------------------------------ -

new_params1 <- params

new_params1$num_days <- 5
new_params1$num_infected_seeds <- 0
new_params1$num_contacts_community_day <- 3
new_params1$contact_prob_household <- 1.0
new_params1$contact_prob_school <- 0.05
new_params1$contact_prob_workplace <- 0.02

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

new_params2$num_days <- 10
new_params2$num_infected_seeds <- 0
new_params2$num_contacts_community_day <- 2
new_params2$contact_prob_household <- 1.0
new_params2$contact_prob_school <- 0.02
new_params2$contact_prob_workplace <- 0.02

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
# PHASE 3: EXIT STRATEGY IMPOSED OVER 7 DAY PERIOD (DAY 30) ----
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
# PHASE 4: FULL VACCINATION IMPOSED OVER 10 DAY PERIOD (DAY 40) ----
# ------------------------------------------------------------------------ -

# ------------------------------------------------------------------------ -
# PARAMETER CONFIGURATION 3 (VACCINATION IMPOSED OVER 10 DAY PERIOD) ----
# ------------------------------------------------------------------------ -

new_params3 <- params

new_params3$num_days <- 5
new_params3$num_infected_seeds <- 0

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
# PARAMETER CONFIGURATION 4 (VACCINATION IMPOSED OVER 10 DAY PERIOD) ----
# ------------------------------------------------------------------------ -

new_params4 <- params

new_params4$num_days <- 20
new_params4$num_infected_seeds <- 0
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

# ------------------------------------------------------------------------ -
# CREATE FINAL POPULATION DATA ----
# ------------------------------------------------------------------------ -

# ------------------------------------------------------------------------ -
# DATASET0: POPULATION DATA ----
# ------------------------------------------------------------------------ -

# store final population dataset
saveRDS(list(log_health = health_time_data$log_health, params = health_time_data$params),
        file = file.path(health_time_data$params$output_dir, "final_health_time.rds"))
saveRDS(pop_data,
        file = file.path(health_time_data$params$output_dir, "final_pop_data.rds"))

# explore population output
pop_data_file <- file.path(health_time_data$params$output_dir, "final_pop_data.rds")
final_pop_data <- readRDS(pop_data_file)

# explore health states output
health_time_file <- file.path(health_time_data$params$output_dir, "final_health_time.rds")
final_health_time_data <- readRDS(health_time_file)

# ------------------------------------------------------------------------ -
# POSTPROCESSING OF POPULATION DATA ----
# ------------------------------------------------------------------------ -

# ------------------------------------------------------------------------ -
# DATASET1: INCIDENCE DATA ----
# ------------------------------------------------------------------------ -

# Central parameter object for delayed incidence data generation
delay_params <- get_default_delay_parameters()
delay_params$output_dir <- params$output_dir

# Incidence data generation with delays
incidence_results <- obtain_incidence_data(final_pop_data, delay_params)

# explore contact tracing data output
incidence_data_file <- file.path(incidence_results$params$output_dir,'incidence_data.rds')
final_incidence_data <- readRDS(incidence_data_file)
head(final_incidence_data)

# store final contact tracing data output
saveRDS(final_incidence_data,
        file = file.path(delay_params$output_dir, "final_incidence_data.rds"))

# plot incidence data
tab <- table(final_incidence_data$time_of_symptom_onset)
final_incidence_plot_df <- data.frame(time = sort(unique(final_incidence_data$time_of_symptom_onset)),
                                new_cases = as.vector(tab),
                                cumulative_cases = cumsum(as.vector(tab)))

library(ggplot2)
p1 <- ggplot(final_incidence_plot_df, 
             aes(x = time, 
                 y = new_cases)) +
  geom_point() + 
  labs(
    x = "Time (in days since start of outbreak)",
    y = "Number of new cases",
    title = "Incidence"
  )

fig_dir <- file.path(incidence_results$params$output_dir, 'final_figures')

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

# Save plot
ggsave(
  filename = file.path(fig_dir, "final_incidence_time.png"),
  plot = p1,
  width = 6,
  height = 4
)

p2 <- ggplot(final_incidence_plot_df, 
             aes(x = time, 
                 y = cumulative_cases)) +
  geom_point() + 
  labs(
    x = "Time (in days since start of outbreak)",
    y = "Cumulative number of cases",
    title = "Cumulative incidence"
  ) 

fig_dir <- file.path(incidence_results$params$output_dir, 'final_figures')

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

# Save plot
ggsave(
  filename = file.path(fig_dir, "final_cumulative_incidence_time.png"),
  plot = p2,
  width = 6,
  height = 4
)

# ------------------------------------------------------------------------ -
# DATASET2: CONTACT TRACING DATA ----
# ------------------------------------------------------------------------ -

# Central parameter object for contact tracing
contact_tracing_params <- get_default_contact_tracing_parameters()
contact_tracing_params$output_dir <- params$output_dir

# Contact tracing data generation
contact_tracing_results <- sample_contact_tracing_data(final_pop_data, contact_tracing_params, verbose)

# explore contact tracing data output
contact_tracing_data_file <- file.path(contact_tracing_results$params$output_dir,'contact_tracing_data.rds')
final_contact_tracing_data <- readRDS(contact_tracing_data_file)
head(final_contact_tracing_data)

# store final contact tracing data output
saveRDS(final_contact_tracing_data,
        file = file.path(contact_tracing_params$output_dir, "final_contact_tracing_data.rds"))

# ------------------------------------------------------------------------ -
# DATASET3: SEROLOGY (PRIOR TO VACCINATION STARTED AT DAY 40) ----
# ------------------------------------------------------------------------ -

# Central parameter object for serology generation
sero_params1 <- get_default_sero_parameters()
sero_params1$output_dir <- params$output_dir

sero_params2 <- sero_params1
sero_params2$sampling_time <- 20

sero_params3 <- sero_params1
sero_params3$sampling_time <- 30

# Serological data generation
sero_results1 <- sample_serological_data(final_pop_data, sero_params1, verbose)
final_sero_data1 <- sero_results1$sero_data 
final_sero_data1$time_of_sampling <- sero_params1$sampling_time

sero_results2 <- sample_serological_data(final_pop_data, sero_params2, verbose)
final_sero_data2 <- sero_results2$sero_data 
final_sero_data2$time_of_sampling <- sero_params2$sampling_time

sero_results3 <- sample_serological_data(final_pop_data, sero_params3, verbose)
final_sero_data3 <- sero_results3$sero_data 
final_sero_data3$time_of_sampling <- sero_params3$sampling_time

# store final contact tracing data output
saveRDS(final_sero_data1,
        file = file.path(sero_params1$output_dir, "final_sero_data_time10.rds"))
saveRDS(final_sero_data2,
        file = file.path(sero_params2$output_dir, "final_sero_data_time20.rds"))
saveRDS(final_sero_data3,
        file = file.path(sero_params2$output_dir, "final_sero_data_time30.rds"))

head(final_sero_data1)

# visualization of the seroprevalence by age group 
final_sero_data = rbind(final_sero_data1, final_sero_data2, final_sero_data3)

seroprev_by_age_group <- tapply(final_sero_data$sero_status, 
                                list(final_sero_data$age_group, final_sero_data$time_of_sampling), mean)
npos_by_age_group <- tapply(final_sero_data$sero_status, 
                            list(final_sero_data$age_group, final_sero_data$time_of_sampling), sum)
n_by_age_group <- tapply(final_sero_data$sero_status, 
                         list(final_sero_data$age_group, final_sero_data$time_of_sampling), length)

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

sero_plot_df2 <- data.frame(age_group = rep(rownames(seroprev_by_age_group), 3),
                            seroprev = as.vector(seroprev_by_age_group),
                            seroprev_ll = as.vector(seroprev_by_age_group_ll),
                            seroprev_ul = as.vector(seroprev_by_age_group_ul),
                            seroprev_ll_asymp = as.vector(seroprev_by_age_group_ll_asymp),
                            seroprev_ul_asymp = as.vector(seroprev_by_age_group_ul_asymp),
                            time = rep(colnames(seroprev_by_age_group), 
                                       each = length(unique(final_sero_data$age_group))))

p1 <- ggplot(sero_plot_df2, 
             aes(x = factor(age_group), 
                 y = seroprev,
                 color = time)) +
  geom_point(position=position_dodge(.9)) + 
  labs(
    x = "Age group",
    y = "Proportion of seropositive individuals",
    title = "Seroprevalence by age group over time"
  ) +
  geom_errorbar(aes(ymin=seroprev_ll, ymax=seroprev_ul), width=.2,
                position=position_dodge(.9)) 

# Save plot
ggsave(
  filename = file.path(fig_dir, "final_seroprevalence_age_group_time.png"),
  plot = p1,
  width = 6,
  height = 4
)