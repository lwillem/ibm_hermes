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
############################################################################ #

# ------------------------------------------------------------------------ -
# DEPENDENCIES ----
# ------------------------------------------------------------------------ -

library(progress)

source('lib/ibm_population.R')
source('lib/ibm_parameters.R')
source('lib/ibm_plot.R')
source('lib/ibm_test.R')

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
#' @param verbose Logical; if TRUE, plots and summaries are produced.
#'
#' @return A list with elements \code{log_health} and \code{params}.
#' @export
run_ibm <- function(params, verbose = TRUE) {

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

  pop_data <- create_population_matrix(params)

  transmission_probs <- compute_transmission_probs(params)

  # define health states
  states <- data.frame(S ="S", I = "I", R = "R",  V = "V", D = "D")
  
  # Vaccination
  id_vaccinated <- sample(params$pop_size,
                          params$pop_size * params$vaccine_coverage)
  pop_data$health[id_vaccinated] <- states$V

  # Seed infections
  seed_ids <- sample(which(pop_data$health == states$S),
                     params$num_infected_seeds)
  pop_data$health[seed_ids] <- states$I
  pop_data$time_of_infection[seed_ids] <- 0

  # Recovery and mortality probabilities
  prob_recovery <- 1 - exp(-1 / params$num_days_infected)
  prob_mortality_general <- 1 - exp(-params$general_mortality_rate)
  prob_mortality_disease <- 1 - exp(-params$disease_mortality_rate)

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

      pop_data$health[new_infections] <- states$I
      pop_data$infector[new_infections] <- i
      pop_data$infector_age[new_infections] <- pop_data$age[i]
      pop_data$time_of_infection[new_infections] <- day
      pop_data$secondary_cases[i] <-
        pop_data$secondary_cases[i] + sum(new_infections)
      pop_data$generation_interval[new_infections] <-
        day - pop_data$time_of_infection[i]
    }

    recovered <- is_infected &
      rbinom(params$pop_size, 1, prob_recovery) == 1
    pop_data$health[recovered] <- states$R

    prob_mortality <- prob_mortality_general[pop_data$age] +
      is_infected * prob_mortality_disease[pop_data$age]
    deaths <- rbinom(params$pop_size, 1, prob_mortality) == 1
    pop_data$health[deaths] <- states$D

    log_health_matrix[, day] <- pop_data$health
  }

  
  ## Output ----
  # -------------------------- -
  log_health <- as.data.frame(sapply(
    states,
    function(s) colSums(log_health_matrix == s) / params$pop_size
  ))

  if (verbose) {

    baseline <- if (params$bool_add_baseline) {
      run_ibm_default(verbose = FALSE)
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
run_ibm_default <- function(verbose = FALSE) {

  default_params <- get_default_parameters()
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

