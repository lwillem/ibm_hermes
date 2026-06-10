############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: retrieving contact tracing data from the generated data
#
# This script is distributed in the hope that it will be useful, but without
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 sabrams, SIMID, UHASSELT, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #

# ------------------------------------------------------------------------ -
# MAIN DATA GENERATING FUNCTION ----
# ------------------------------------------------------------------------ -

#' Contact-tracing data
#'
#' @param pop_data Population data with infection histories and times of death.
#' @param contact_tracing_params List of parameters governing contact tracing.
#' @param verbose Logical; if TRUE, plots and summaries are produced.
#'
#' @return A list with elements \code{log_sero}.
#' @export
sample_contact_tracing_data <- function(pop_data, contact_tracing_params, verbose) {
  
  ## Reproducibility ----
  # ------------------------- -
  set.seed(contact_tracing_params$rng_seed)
  
  ## Defensive checks ----
  # -------------------------- -
  is_available = (!is.na(pop_data$time_of_infection)) & (!is.na(pop_data$time_of_symptom_onset)) & 
    (pop_data$time_of_symptom_onset > contact_tracing_params$start_contact_tracing) &
    (pop_data$time_of_symptom_onset < contact_tracing_params$end_contact_tracing)
  
  if (contact_tracing_params$n > sum(is_available)){
    warning("Contact tracing sample size is too large in view of available infectee-infector pairs")
    return(NULL)
  }
  
  ## Subset of available data from individuals that are infected and experienced symptoms
  ## within the specified contact tracing window 
  #-------------------------------------------------------------------------------------- -
  ss_pop_data = pop_data[is_available, ]

  ## Select infector-infectee pairs ----
  # ---------------------------------- -
  if (contact_tracing_params$design == "independent"){
    infectee_dat <- ss_pop_data[sample(nrow(ss_pop_data), size = contact_tracing_params$n), ]
    infector_dat <- pop_data[infectee_dat$infector, ]
  }
  
  if (contact_tracing_params$design == "chains"){
    index_cases = sample(rownames(ss_pop_data), size = contact_tracing_params$n)
    
    for (chain_id in 1:contact_tracing_params$n){
      if (chain_id == 1){
        infectee_dat <- ss_pop_data[index_cases[chain_id], ]
        secondary_cases <- ss_pop_data[ss_pop_data$infector_id == index_cases[chain_id], ]
        infectee_dat <- rbind(infectee_dat, secondary_cases)
      } else {
        infectee_dat <- rbind(infectee_dat, ss_pop_data[index_cases[chain_id], ])
        secondary_cases <- ss_pop_data[ss_pop_data$infector_id == index_cases[chain_id], ]
        infectee_dat <- rbind(infectee_dat, secondary_cases)
      }
    }

    infector_dat <- pop_data[infectee_dat$infector, ]
  }
  
  ## Output ----
  # -------------------------- -
  contact_tracing_data <- data.frame(time = contact_tracing_params$end_contact_tracing,
                                     infector_id = rownames(pop_data)[infectee_dat$infector], 
                                     infectee_id = rownames(infectee_dat), 
                                     infector_symptom_onset = infector_dat$time_of_symptom_onset, 
                                     infectee_symptom_onset = infectee_dat$time_of_symptom_onset,
                                     infector_age = infector_dat$age, 
                                     infectee_age = infectee_dat$age)
  
  saveRDS(contact_tracing_data,
          file = file.path(contact_tracing_params$output_dir, "contact_tracing_data.rds"))
  
  return(list(contact_tracing_data = contact_tracing_data, params = contact_tracing_params))
}