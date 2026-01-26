############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: to run regression testing for the modelling framework
# 
# This script is distributed in the hope that it will be useful,but without 
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #

#' Run an internal logic check for the IBM model
#'
#' This function runs the IBM with a fixed set of parameters and
#' compares the output against a stored reference run. It is meant
#' as a lightweight regression test to detect unintended changes
#' in model parameters or outputs.
#'
#' @param update_reference Logical. If TRUE, the current model output
#'   is stored as the new reference output.
#' @param overrule_param_check to overrule the behaviour that model output
#' is not checked if parameters differ
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Load default model parameters
#'   \item Apply fixed testing settings
#'   \item Run the IBM model
#'   \item Optionally update or compare against a reference output
#' }
#'
#' @return Invisibly returns TRUE if the model ran successfully.
#'
#' @examples
#' run_ibm_regression_test()
#' update_ibm_reference()
run_ibm_regression_test <- function(update_reference = FALSE,
                                    overrule_param_check = TRUE) {

  message("Start IBM regression test")
  
  # Run the IBM with default parameters
  ibm_output <- run_ibm_default(verbose = FALSE)
  param <- ibm_output$param
  
  message("Model execution successfully completed")

  reference_file <- paste0(param$output_dir,"_reference/ibm_out_baseline.rds")

  # Update reference output if requested
  if (update_reference) {

    if (!dir.exists(dirname(reference_file))) {
      dir.create(dirname(reference_file), recursive = TRUE)
    }

    saveRDS(ibm_output, reference_file)
    message("New reference model output stored: ", reference_file)
    return(invisible(TRUE))
  }

  # Compare against existing reference
  if (!file.exists(reference_file)) {
    message(
      "Model reference not available. ",
      "Run update_ibm_reference() to create one."
    )
    return(invisible(TRUE))
  }

  reference_output <- readRDS(reference_file)

  # Check parameter equality (exclude a specific subet)
  param_equal <- are_parameters_equal(ibm_output$param,
                                      reference_output$param,
                                      verbose = TRUE)
  
  # Only compare outputs if parameters are identical (this can be overruled)
  if (!overrule_param_check && !param_equal) {
    message("Model output not tested because parameters differ – consider to update the IBM reference output using `update_ibm_reference()`")
    return(invisible(TRUE))
  }

  if (!identical(
    unlist(ibm_output$log_health),
    unlist(reference_output$log_health)
  )) {
    message("Model output changed")
  } else {
    message("Model output unchanged – regression test complete")
  }

  invisible(TRUE)
}

#' Update the IBM reference output
#'
#' Convenience wrapper around \code{run_ibm_regression_test()} to
#' overwrite the stored reference model output.
#'
#' @return Invisibly returns TRUE.
update_ibm_reference <- function() {
  run_ibm_regression_test(update_reference = TRUE)
}
