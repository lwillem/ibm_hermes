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

#' Default parameters related to reporting delay distribution
#'
#' Creates and returns a list containing all parameters controlling
#' the reporting delay process.
#'
#' @return A named list of model parameters.
#' @export
get_default_delay_parameters <- function() {

  params <- list(
    # delay parameters
    mean = 2,
    max_delay = 6,
    
    # missingness in symptomatic cases
    prob_missing = 0.1,
    
    # visualisation and output
    bool_show_demographics = TRUE,
    bool_add_baseline      = FALSE,
    bool_return_prevelance = FALSE,

    output_dir = "output/ibm_modified",
    plot_mfrow = c(2, 2)
  )

  return(params)
}

#' Print available parameter names
#'
#' Convenience function to inspect configurable parameters.
print_delay_parameter_names <- function() {
  print(names(get_default_delay_parameters()))
}
