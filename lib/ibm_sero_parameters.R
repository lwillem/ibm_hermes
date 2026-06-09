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

#' Default serological sampling parameters
#'
#' Creates and returns a list containing all parameters controlling
#' cross-sectional serological survey sampling.
#'
#' @return A named list of model parameters.
#' @export
get_default_sero_parameters <- function() {

  params <- list(
    # serological survey sample size
    n              = 1000,
    sampling_time  = 10,
    
    # log-antibody titers (in IU/mL)
    peak    = 5,
    decay   = 0.00,
    sigma1   = 0.5,
    sigma2   = 0.75,
    seroneg_offset = 0.25,
    LLOD    = log(5),
      
    # random number generation
    rng_seed = 062026,

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
print_sero_parameter_names <- function() {
  print(names(get_default_sero_parameters()))
}
