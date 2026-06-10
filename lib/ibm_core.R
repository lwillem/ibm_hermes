############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: core HERMES individual-based model functions
#
# This file contains the main execution kernel of the HERMES IBM and
# a small number of internal helper functions. All epidemiological
# behaviour is controlled via the params list.
#
# This script is distributed in the hope that it will be useful, but without
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
#                    sabrams, SIMID, UHASSELT, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #

# ------------------------------------------------------------------------ -
# DEPENDENCIES ----
# ------------------------------------------------------------------------ -

library(progress)
library(extraDistr)

source('lib/ibm_population.R')
source('lib/ibm_parameters.R')
source('lib/ibm_plot.R')
source('lib/ibm_test.R')

source('lib/ibm_pop_death.R')

source('lib/ibm_serology.R')
source('lib/ibm_sero_parameters.R')
source('lib/ibm_contact_tracing.R')
source('lib/ibm_contact_tracing_parameters.R')
source('lib/ibm_reporting_delays.R')
source('lib/ibm_delay_parameters.R')

# ------------------------------------------------------------------------ -
# INTERNAL HELPERS ----
# ------------------------------------------------------------------------ -

#' Compute transmission probabilities
#'
#' Internal helper to compute per-contact transmission probabilities for
#' different social contexts. This function performs only deterministic
#' transformations of parameters and introduces no stochasticity.
#'
#' @param params Model parameter list.
#'
#' @return Named list with transmission probabilities.
compute_transmission_probs <- function(params) {

  prob_contact_community <- 1 - exp(
    -params$num_contacts_community_day / params$pop_size
  )

  list(
    prob_community = prob_contact_community * params$transmission_prob,
    prob_household = params$contact_prob_household * params$transmission_prob,
    prob_school    = params$contact_prob_school    * params$transmission_prob,
    prob_workplace = params$contact_prob_workplace * params$transmission_prob
  )
}

#' Check susceptibility status
#'
#' @param health_vector Character vector of health states.
#' @param states Vector with health state options
#' @return Logical vector indicating susceptibility, wich an be partial in case of vaccination
is_susceptible <- function(health_vector, states) {
  health_vector %in% c(states$S, states$V)
}

# ------------------------------------------------------------------------ -
# MAIN MODEL FUNCTION ----
# ------------------------------------------------------------------------ -

#' Run the individual-based model
#'
#' Executes a stochastic individual-based epidemic model with household,
#' school, workplace, and community transmission.
#'
#' @param params List of model parameters.
#' @param default_params Comparing results with default scenario.
#' @param verbose Logical; if TRUE, plots and summaries are produced.
#'
#' @return A list with elements \code{log_health} and \code{params}.
#' @export
run_ibm <- function(params, default_params = NULL, verbose = TRUE) {

  ## Defensive checks ----
  # -------------------------- -
  if (params$num_infected_seeds > params$pop_size) {
    warning("Population size smaller than number of infected seeds")
    return(NULL)
  }

  if (any(unlist(params[grepl('num|size|age|prob', names(params))]) < 0)) {
    warning("Negative parameter values are not allowed")
    return(NULL)
  }

  if (!is.logical(unlist(params[grepl('bool', names(params))]))) {
    warning("All 'bool_*' parameters must be logical")
    return(NULL)
  }

  if (!dir.exists(params$output_dir)) {
    dir.create(params$output_dir, recursive = TRUE)
  }

  ## Initialisation ----
  # ------------------------- -
  time_start <- Sys.time()
  set.seed(params$rng_seed)

  # initialize the population
  pop_data <- create_population_matrix(params)

  transmission_probs <- compute_transmission_probs(params)

  # define health states
  states <- data.frame(S = "S", I = "I", R = "R",  V = "V", D = "D")
  
  # vaccination
  id_vaccinated <- sample(params$pop_size,
                          params$pop_size * params$vaccine_coverage)
  pop_data$health[id_vaccinated] <- states$V
  pop_data$time_of_vaccination[id_vaccinated] <- 0
    
  # seed infections
  seed_ids <- sample(which(pop_data$health == states$S),
                     params$num_infected_seeds)
  pop_data$health[seed_ids] <- states$I
  pop_data$time_of_infection[seed_ids] <- 0

  # recovery and mortality probabilities
  prob_recovery <- 1 - exp(-1 / params$num_days_infected)
  prob_hospitalized <- params$prob_hospital*(1 - exp(-1 / params$num_days_before_hospitalized))
  prob_mortality_disease <- 1 - exp(-params$disease_mortality_rate)
  if (!is.null(params$general_mortality_rate)){
    prob_mortality_general <- 1 - exp(-params$general_mortality_rate)}
  
  # natural mortality
  if (is.null(params$general_mortality_rate)){
    prob_mortality_general <- rep(0, length(prob_mortality_disease))
    for (id in 1:params$pop_size){
      if (pop_data$sex[id] == "Male"){
        pop_data$age_natural_death[id] <- gen_pop(params$date, pop_data$age[id], mlt)
      } else {
        pop_data$age_natural_death[id] <- gen_pop(params$date, pop_data$age[id], flt)
      }
    }
    
    if (sum(is.na(pop_data$age_natural_death)) > 0) {
      warning("Population death times are not generated correctly")
      return(NULL)
    }
  }
  
  # log health matrix
  log_health_matrix <- matrix(
    NA, nrow = params$pop_size, ncol = params$num_days
  )

  ## Model loop ----
  # ------------------------ -
  pb <- progress_bar$new(
    format = paste0("Run ", basename(params$output_dir),
                    ": [:bar] :percent ETA: :eta"),
    total = params$num_days, clear = FALSE, width = 60
  )

  for (day in seq_len(params$num_days)) {

    pb$tick()

    is_infected <- pop_data$health == states$I
    is_death <- pop_data$health == states$D
    infected_ids <- which(is_infected)

    for (i in infected_ids) {

      prob_infection <- is_susceptible(pop_data$health, states) *
        transmission_probs$prob_community +
        (is_susceptible(pop_data$health, states) &
           pop_data$hh_id[i] == pop_data$hh_id) *
        transmission_probs$prob_household

      if (!is.na(pop_data$classroom_id[i])) {
        in_class <- is_susceptible(pop_data$health, states) &
          !is.na(pop_data$classroom_id) &
          pop_data$classroom_id[i] == pop_data$classroom_id
        prob_infection[in_class] <-
          prob_infection[in_class] + transmission_probs$prob_school
      }

      if (!is.na(pop_data$workplace_id[i])) {
        at_work <- is_susceptible(pop_data$health, states) &
          !is.na(pop_data$workplace_id) &
          pop_data$workplace_id[i] == pop_data$workplace_id
        prob_infection[at_work] <-
          prob_infection[at_work] + transmission_probs$prob_workplace
      }

      prob_infection[pop_data$health == states$V] <-
        prob_infection[pop_data$health == states$V] *
        (1 - params$vaccine_effectiveness)

      new_infections <- rbinom(
        params$pop_size, size = 1, prob = prob_infection
      ) == 1
      
      prob_symptoms <- params$prop_symptoms
      symptoms <- rbinom(params$pop_size, size = 1, prob = prob_symptoms) == 1
      symptom_onsets <- rdgamma(
        params$pop_size, scale = params$scale_symptom_onset, 
        shape = params$shape_symptom_onset
      )

      pop_data$health[new_infections] <- states$I
      pop_data$infector[new_infections] <- i
      pop_data$infector_id[new_infections] <- rownames(pop_data)[i]
      pop_data$infector_age[new_infections] <- pop_data$age[i]
      pop_data$time_of_infection[new_infections] <- day
      pop_data$secondary_cases[i] <-
        pop_data$secondary_cases[i] + sum(new_infections)
      pop_data$generation_interval[new_infections] <-
        day - pop_data$time_of_infection[i]
      
      pop_data$symptom_onset[new_infections & symptoms] <- symptom_onsets[new_infections & symptoms] 
      pop_data$time_of_symptom_onset[new_infections & symptoms] <- pop_data$symptom_onset[new_infections & symptoms] + day 
    }

    # Hospitalization
    is_symptoms <- (pop_data$health != states$R) & (pop_data$health != states$D) & (pop_data$time_of_symptom_onset < day)
    is_hospitalized <- !is.na(pop_data$time_of_hospitalization)
      
    hospitalized <- is_symptoms & (!is_hospitalized) &
      rbinom(params$pop_size, size = 1, prob_hospitalized) == 1
    pop_data$time_of_hospitalization[hospitalized] <- day
    
    # Recovery (both within and outside hospital)
    recovered <- is_infected &
      rbinom(params$pop_size, 1, prob_recovery) == 1
    pop_data$health[recovered] <- states$R

    # In case general mortality rates are provided, no distinction between
    # natural and disease-related mortality is made; in case of absence
    # of general mortality rates (but population death times) the death times
    # solely represent disease-related mortality. 
    prob_mortality <- prob_mortality_general[pop_data$age] +
      is_infected * prob_mortality_disease[pop_data$age]
    deaths <- rbinom(params$pop_size, 1, prob_mortality) == 1
    pop_data$health[!is_death & deaths] <- states$D
    pop_data$time_of_death[!is_death & deaths] <- day

    if (is.null(params$general_mortality_rate)){
      current_ages <- pop_data$age + (day/365)
      natural_deaths <- (!is_death) & ((current_ages > pop_data$age_natural_death) == 1)
      pop_data$health[natural_deaths] <- states$D
      pop_data$time_of_natural_death[natural_deaths] <- day
    }
    
    log_health_matrix[, day] <- pop_data$health
  }
  
  ## Correct times of symptom onset not to lie beyond death times ---
  #------------------------------------------------------------------ -
  pop_data$time_of_symptom_onset[!is.na(pop_data$time_of_symptom_onset) & 
                                   pop_data$time_of_symptom_onset > pop_data$time_of_death] <- NA
  pop_data$time_of_symptom_onset[!is.na(pop_data$time_of_symptom_onset) & 
                                   pop_data$time_of_symptom_onset > pop_data$time_of_natural_death] <- NA

  ## Output ----
  # -------------------------- -
  log_health <- as.data.frame(sapply(
    states,
    function(s) colSums(log_health_matrix == s) / params$pop_size
  ))

  if (verbose) {

    baseline <- if (params$bool_add_baseline) {
      run_ibm_default(default_params, verbose = FALSE)
    } else {
      NA
    }

    par(mfrow = params$plot_mfrow)
    plot_health_states(log_health, params, out_baseline = baseline)
    plot_secondary_cases(pop_data, params)
    plot_generation_interval(pop_data, params)
    plot_transmission_matrix(pop_data)
    par(mfrow = c(1, 1))

    print_model_parameters(params)
    print_model_results(log_health, time_start, baseline)
  }

  saveRDS(list(log_health = log_health, params = params),
          file = file.path(params$output_dir, "health_time.rds"))
  saveRDS(pop_data,
          file = file.path(params$output_dir, "pop_data.rds"))

  return(list(log_health = log_health, params = params))
}

# BASELINE RUNNER ----
# ---------------------------- -

#' Run baseline model with default parameters
#'
#' @param verbose Logical to disable all output 
run_ibm_default <- function(default_params = NULL, verbose = FALSE) {

  if (is.null(default_params)) default_params <- get_default_parameters()
  default_params$bool_add_baseline <- FALSE
  default_params$bool_show_demographics <- FALSE
  default_params$output_dir <-
    paste0(default_params$output_dir, "_default")

  run_ibm(default_params, verbose = verbose)
}

# function to print the model results
print_model_results <- function(log_health,time_start,out_baseline=NA){
  
  bool_add_baseline <- !any(is.na(out_baseline))
  if(bool_add_baseline){
    # default epidemic characteristics
    default_ti <- paste0('   [baseline: ',round((out_baseline$log_health$I[length(out_baseline$log_health$I)] +
                                                   out_baseline$log_health$R[length(out_baseline$log_health$R)])*100,digits=0),'%]')
    default_pp <- paste0('   [baseline: ',round(max(out_baseline$log_health$I)*100,digits=0),'%]')
    default_pd <- paste0('    [baseline: ',which(out_baseline$log_health$I == max(out_baseline$log_health$I))[1],']')
  }
  
  # print total incidence
  print('-------------')
  print('MODEL RESULTS')
  
  print(paste0('total incidence: ',round((log_health$I[length(log_health$I)] + log_health$R[length(log_health$R)])*100,digits=0),'%',
               ifelse(bool_add_baseline,default_ti,'')))
  
  # print peak details
  print(paste0('Peak prevalence: ',round(max(log_health$I)*100,digits=0),'%',
               ifelse(bool_add_baseline,default_pp,'')))
  print(paste0('Peak day:        ',which(log_health$I == max(log_health$I))[1], 
               ifelse(bool_add_baseline,default_pd,'')))
  
  # print total run time
  total_time <- as.double(Sys.time() - time_start,unit='secs')
  print(paste0('Total run time:  ',round(total_time,digits=0),'s'))
  
}

# ------------------------------------------------------------------------ -
# MODEL FUNCTION BASED ON PREVIOUS INITIALIZATION ----
# ------------------------------------------------------------------------ -

#' Run the individual-based model from previous initialization and state
#'
#' Executes a stochastic individual-based epidemic model with household,
#' school, workplace, and community transmission.
#'
#' @param pop_data Generated population data. 
#' @param prev_log_health_matrix Matrix with previous health status data.
#' @param params List of model parameters.
#' @param intervention ID for the intervention
#' @param verbose Logical; if TRUE, plots and summaries are produced.
#'
#' @return A list with elements \code{log_health} and \code{params}.
#' @export
rerun_ibm <- function(pop_data, prev_log_health_matrix, params, intervention, 
                      default_params = NULL, verbose = TRUE) {
  
  ## Initialisation ----
  # ------------------------- -
  time_start <- Sys.time()
  set.seed(params$rng_seed)
  
  # number of previously simulated days
  nprev_days <- nrow(prev_log_health_matrix)
  message(paste0("Number of previously simulated days is ", nprev_days))
  
  # consider (potentially) new transmission parameters
  transmission_probs <- compute_transmission_probs(params)
  
  # define health states
  states <- data.frame(S = "S", I = "I", R = "R", V = "V", D = "D")
  
  # additional seed infections
  is_still_susceptible <- which(pop_data$health == states$S)
  if (length(is_still_susceptible) > params$num_infected_seeds){
    seed_ids <- sample(is_still_susceptible,
                       params$num_infected_seeds)
    pop_data$health[seed_ids] <- states$I
    pop_data$time_of_infection[seed_ids] <- nprev_days 
  } else {
    message("No additional seeding performed")
  }

  # vaccination among those still alive and not yet vaccinated before
  # coverage expressed in terms of total population
  # people move to V only in case they are susceptible
  is_eligible <- which(pop_data$health != states$V & pop_data$health != states$D)
  if (length(is_eligible) > round(length(is_eligible) * params$vaccine_coverage, 0)){
    id_vaccinated <- sample(is_eligible,
                            round(length(is_eligible) * params$vaccine_coverage, 0))
    pop_data$health[id_vaccinated] <- states$V
    pop_data$time_of_vaccination[id_vaccinated] <- nprev_days
  } else {
    message("All available individuals vaccinated")
    id_vaccinated <- is_eligible
    pop_data$health[id_vaccinated] <- states$V
    pop_data$time_of_vaccination[id_vaccinated] <- nprev_days
  }
  
  # recovery and mortality probabilities
  prob_recovery <- 1 - exp(-1 / params$num_days_infected)
  prob_hospitalized <- params$prob_hospital*(1 - exp(-1 / params$num_days_before_hospitalized))
  prob_mortality_disease <- 1 - exp(-params$disease_mortality_rate)
  if (!is.null(params$general_mortality_rate)){
    prob_mortality_general <- 1 - exp(-params$general_mortality_rate)}
  
  # natural mortality (if is.null() = T - natural death times have been generated before)
  if (is.null(params$general_mortality_rate)){
    prob_mortality_general <- rep(0, length(prob_mortality_disease))
  }
  
  # log health matrix (health status for new days)
  log_health_matrix <- matrix(
    NA, nrow = params$pop_size, ncol = params$num_days
  )
  
  ## Model loop ----
  # ------------------------ -
  pb <- progress_bar$new(
    format = paste0("Run ", basename(params$output_dir),
                    ": [:bar] :percent ETA: :eta"),
    total = params$num_days, clear = FALSE, width = 60
  )
  
  for (day in (nprev_days + seq_len(params$num_days))) {
 
    pb$tick()
    
    is_infected <- pop_data$health == states$I
    is_death <- pop_data$health == states$D
    infected_ids <- which(is_infected)

    for (i in infected_ids) {
      
      prob_infection <- is_susceptible(pop_data$health, states) *
        transmission_probs$prob_community +
        (is_susceptible(pop_data$health, states) &
           pop_data$hh_id[i] == pop_data$hh_id) *
        transmission_probs$prob_household
      
      if (!is.na(pop_data$classroom_id[i])) {
        in_class <- is_susceptible(pop_data$health, states) &
          !is.na(pop_data$classroom_id) &
          pop_data$classroom_id[i] == pop_data$classroom_id
        prob_infection[in_class] <-
          prob_infection[in_class] + transmission_probs$prob_school
      }
      
      if (!is.na(pop_data$workplace_id[i])) {
        at_work <- is_susceptible(pop_data$health, states) &
          !is.na(pop_data$workplace_id) &
          pop_data$workplace_id[i] == pop_data$workplace_id
        prob_infection[at_work] <-
          prob_infection[at_work] + transmission_probs$prob_workplace
      }
      
      prob_infection[pop_data$health == states$V] <-
        prob_infection[pop_data$health == states$V] *
        (1 - params$vaccine_effectiveness) * is.na(pop_data$time_of_infection[pop_data$health == states$V])

      new_infections <- rbinom(
        params$pop_size, size = 1, prob = prob_infection
      ) == 1
      
      prob_symptoms <- params$prop_symptoms
      symptoms <- rbinom(params$pop_size, size = 1, prob = prob_symptoms) == 1
      symptom_onsets <- rdgamma(
        params$pop_size, scale = params$scale_symptom_onset, 
        shape = params$shape_symptom_onset
      )
  
      pop_data$health[new_infections] <- states$I
      pop_data$infector[new_infections] <- i
      pop_data$infector_id[new_infections] <- rownames(pop_data)[i]
      pop_data$infector_age[new_infections] <- pop_data$age[i]
      pop_data$time_of_infection[new_infections] <- day
      pop_data$secondary_cases[i] <-
        pop_data$secondary_cases[i] + sum(new_infections)
      pop_data$generation_interval[new_infections] <-
        day - pop_data$time_of_infection[i]
      
      pop_data$symptom_onset[new_infections & symptoms] <- symptom_onsets[new_infections & symptoms] 
      pop_data$time_of_symptom_onset[new_infections & symptoms] <- pop_data$symptom_onset[new_infections & symptoms] + day 
    }
    
    # Hospitalization
    is_symptoms <- (pop_data$health != states$R) & (pop_data$health != states$D) & (pop_data$time_of_symptom_onset < day)
    is_hospitalized <- !is.na(pop_data$time_of_hospitalization)
    
    hospitalized <- is_symptoms & (!is_hospitalized) &
      rbinom(params$pop_size, size = 1, prob_hospitalized) == 1
    pop_data$time_of_hospitalization[hospitalized] <- day
    
    # Recovery (both within and outside hospital)
    recovered <- is_infected &
      rbinom(params$pop_size, 1, prob_recovery) == 1
    pop_data$health[recovered] <- states$R
    
    # In case general mortality rates are provided, no distinction between
    # natural and disease-related mortality are made; in case of absence
    # of general mortality rates (but population death times) the death times
    # solely represent disease-related mortality. 
    prob_mortality <- prob_mortality_general[pop_data$age] +
      is_infected * prob_mortality_disease[pop_data$age]
    deaths <- rbinom(params$pop_size, 1, prob_mortality) == 1
    pop_data$health[!is_death & deaths] <- states$D
    pop_data$time_of_death[!is_death & deaths] <- day
    
    if (is.null(params$general_mortality_rate)){
      current_ages <- pop_data$age + (day/365)
      natural_deaths <- (!is_death) & ((current_ages > pop_data$age_natural_death) == 1)
      pop_data$health[natural_deaths] <- states$D
      pop_data$time_of_natural_death[natural_deaths] <- day
    }
    
    log_health_matrix[, day - nprev_days] <- pop_data$health
  }
  
  ## Correct times of symptom onset not to lie beyond death times ---
  #------------------------------------------------------------------ -
  pop_data$time_of_symptom_onset[!is.na(pop_data$time_of_symptom_onset) & 
                                   pop_data$time_of_symptom_onset > pop_data$time_of_death] <- NA
  pop_data$time_of_symptom_onset[!is.na(pop_data$time_of_symptom_onset) & 
                                   pop_data$time_of_symptom_onset > pop_data$time_of_natural_death] <- NA
  
  ## Output ----
  # -------------------------- -
  log_health <- as.data.frame(sapply(
    states,
    function(s) colSums(log_health_matrix == s) / params$pop_size
  ))
  
  # extended log_health matrix
  log_health = rbind(prev_log_health_matrix, log_health)
  
  if (verbose) {
    
    baseline <- if (params$bool_add_baseline) {
      run_ibm_default(default_params, verbose = FALSE)
    } else {
      NA
    }
    
    plot_params <- params
    plot_params$num_days <- nprev_days + params$num_days 
      
    par(mfrow = params$plot_mfrow)
    plot_health_states(log_health, params, out_baseline = baseline)
    plot_secondary_cases(pop_data, plot_params)
    plot_generation_interval(pop_data, plot_params)
    plot_transmission_matrix(pop_data)
    par(mfrow = c(1, 1))
    
    print_model_parameters(params)
    print_model_results(log_health, time_start, baseline)
  }
  
  saveRDS(list(log_health = log_health, params = params),
          file = file.path(params$output_dir, paste0("health_time_intervention", intervention, ".rds")))
  saveRDS(pop_data,
          file = file.path(params$output_dir, paste0("pop_data_intervention", intervention,".rds")))
  
  return(list(log_health = log_health, params = params))
}