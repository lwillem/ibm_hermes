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

#' Default sampling parameters related to the contact tracing data
#'
#' Creates and returns a list containing all parameters controlling
#' contact tracing.
#'
#' @return A named list of model parameters.
#' @export
get_default_contact_tracing_parameters <- function() {
  
  params <- list(
    # number of infector-infectee pairs
    n              = 200,
    
    # contact tracing window
    start_contact_tracing  = 0,
    end_contact_tracing  = 30,
    
    # probable time of contact
    cens_prob = 0.6,
    
    # sampling design
    design = "chains",
    
    # random number generation
    rng_seed = 072026,
    
    # visualisation and output
    bool_show_seroprevelance = FALSE,
    
    output_dir = "output/ibm_modified",
    plot_mfrow = c(1, 1)
  )
  
  return(params)
}

#' Print available parameter names
#'
#' Convenience function to inspect configurable parameters.
print_contact_tracing_parameter_names <- function() {
  print(names(get_default_contact_tracing_parameters()))
}
