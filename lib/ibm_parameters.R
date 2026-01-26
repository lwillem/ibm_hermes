############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: define, print, and compare model parameters
#
# This script is distributed in the hope that it will be useful, but without
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #

#' Default HERMES model parameters
#'
#' Creates and returns a list containing all parameters controlling
#' demography, transmission, mortality, and output behaviour of the IBM.
#'
#' @return A named list of model parameters.
#' @export
get_default_parameters <- function() {

  params <- list(
    # population and time
    pop_size           = 2000,
    num_days           = 50,
    num_infected_seeds = 3,

    # vaccination
    vaccine_coverage        = 0.2,
    vaccine_effectiveness   = 0.8,

    # disease natural history
    num_days_infected = 7,
    transmission_prob = 0.1,

    # demography
    ages_adult              = 19:80,
    adult_age_tolerance     = 5,
    child_age_tolerance     = 4,
    household_age_gap       = 22:35,

    # schools
    target_school_size = 350,
    target_school_ages = 3:18,

    # contacts
    num_contacts_community_day = 4,
    contact_prob_household     = 0.9,
    contact_prob_school        = 0.5,
    contact_prob_workplace     = 0.1,

    # workplaces
    target_workplace_size = 10,
    target_workplace_ages = 19:65,

    # random number generation
    rng_seed = 2026,

    # mortality (by age)
    general_mortality_rate = c(rep(0, 60), seq(0, 0.01, length = 30)),
    disease_mortality_rate = c(rep(0.01, 60), seq(0.01, 0.7, length = 30)),

    # visualisation and output
    bool_show_demographics = TRUE,
    bool_add_baseline      = FALSE,
    bool_return_prevelance = FALSE,

    output_dir = "output/ibm_flu",
    plot_mfrow = c(2, 2)
  )

  return(params)
}

#' Print available parameter names
#'
#' Convenience function to inspect configurable parameters.
print_parameter_names <- function() {
  print(names(get_default_parameters()))
}

#' Print model parameters and values
#'
#' @param params Parameter list as used by the model.
print_model_parameters <- function(params) {
  print("MODEL PARAMETERS")
  for (param_name in names(params)) {
    values <- unlist(params[[param_name]])
    if (length(values) > 10) {
      values <- c(values[1:9], "...")
    }
    print(paste0(param_name, ": ", paste(values, collapse = ",")))
  }
}

#' Compare two parameter sets
#'
#' @param params1 First parameter list.
#' @param params2 Second parameter list.
#' @param param_exclude Character vector of parameter names to ignore.
#' @param verbose Logical; print changed parameters if TRUE.
#'
#' @return Logical indicating whether parameters are equal.
are_parameters_equal <- function(params1,
                                 params2,
                                 param_exclude = c("output_dir"),
                                 verbose = FALSE) {

  equal <- TRUE

  for (name in names(params1)) {
    if (!name %in% param_exclude) {
      if (!identical(unlist(params1[[name]]), unlist(params2[[name]]))) {
        if (verbose) message("Parameter changed: ", name)
        equal <- FALSE
      }
    }
  }

  return(equal)
}
