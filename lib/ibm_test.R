############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: regression testing for the IBM
#
# This script is distributed in the hope that it will be useful, but without
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #

#' Run IBM regression test
#'
#' @param update_reference Logical; overwrite stored reference if TRUE.
#' @param overrule_param_check Logical; compare outputs even if parameters differ.
#' @return Invisibly returns TRUE.
run_ibm_regression_test <- function(update_reference = FALSE,
                                    overrule_param_check = TRUE) {

  message("Start IBM regression test")

  model_output <- run_ibm_default(verbose = FALSE)
  params <- model_output$params

  reference_file <- paste0(params$output_dir,
                           "_reference/ibm_out_baseline.rds")

  if (update_reference) {
    dir.create(dirname(reference_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(model_output, reference_file)
    message("Reference output updated")
    return(invisible(TRUE))
  }

  if (!file.exists(reference_file)) {
    message("Reference output missing – create one with update_ibm_reference()")
    return(invisible(TRUE))
  }

  reference_output <- readRDS(reference_file)

  params_equal <- are_parameters_equal(
    model_output$params,
    reference_output$params,
    verbose = TRUE
  )

  if (!overrule_param_check && !params_equal) {
    message("Parameters differ – output comparison skipped")
    return(invisible(TRUE))
  }

  if (!identical(unlist(model_output$log_health),
                 unlist(reference_output$log_health))) {
    message("Model output changed")
  } else {
    message("Model output unchanged – regression test passed")
  }

  invisible(TRUE)
}

#' Update regression test reference output
update_ibm_reference <- function() {
  run_ibm_regression_test(update_reference = TRUE)
}
