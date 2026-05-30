############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: Generate serological survey data at a specific cross-sectional time
#
# This script is distributed in the hope that it will be useful, but without
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 sabrams, SIMID, UHASSELT, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #

# ------------------------------------------------------------------------ -
# MAIN DATA GENERATING FUNCTION ----
# ------------------------------------------------------------------------ -

#' Cross-sectional serological survey data
#'
#' @param pop_data Population data with infection histories and times of death.
#' @param sero_params List of parameters governing serological survey sampling.
#' @param verbose Logical; if TRUE, plots and summaries are produced.
#'
#' @return A list with elements \code{log_sero}.
#' @export
sample_serological_data <- function(pop_data, sero_params, verbose) {
  
  ## Reproducibility ----
  # ------------------------- -
  set.seed(sero_params$rng_seed)
  
  ## Defensive checks ----
  # -------------------------- -
  is_alive = is.na(pop_data$time_of_death) | (pop_data$time_of_death > sero_params$sampling_time)
  if (sero_params$n > sum(is_alive)){
    warning("Serological survey sample size is too large")
    return(NULL)
  }
    
  ## Subset data from individuals that are alive
  #--------------------------------------------- -
  ss_pop_data = pop_data[is_alive, ]
  
  ## Select individuals if still alive ----
  # ------------------------------------- -
  sample_ids <- sample(nrow(ss_pop_data), size = sero_params$n)
  ind_ids <- rownames(ss_pop_data)[sample_ids]
  
  ## Age at sampling ----
  # ------------------------ -
  age_at_sampling <- ss_pop_data$age[sample_ids] + sero_params$sampling_time
  
  ## Serological data ---
  # ------------------------ - 
  age_group <- vector(length = sero_params$n)
  status <- vector(length = sero_params$n)
  ab_titer <- vector(length = sero_params$n)
  toi <- vector(length = sero_params$n)
  tso <- vector(length = sero_params$n)
  tsi <- vector(length = sero_params$n)
  tsso <- vector(length = sero_params$n)
  
  age_group_id <- findInterval(age_at_sampling, seq(0, 100, 10))
  age_group_vec <- c("[0, 10)", "[10, 20)", "[20, 30)", "[30, 40)",
                     "[40, 50)", "[50, 60)", "[60, 70)", "[70, 80)",
                     "[80, 90)", "[90, 100)", "100+")
    
  for (i in 1:length(sample_ids)){
    age_group[i] <- age_group_vec[age_group_id[i]]
    toi[i] = ss_pop_data[sample_ids[i], ]$time_of_infection
    tso[i] = ss_pop_data[sample_ids[i], ]$time_of_symptom_onset
    tsi[i] = ifelse(is.na(toi[i]), -10, max(age_at_sampling[i] - toi[i], 0))
    tsso[i] = ifelse(is.na(tso[i]), -10, max(age_at_sampling[i] - tso[i], 0))
    status[i] = as.numeric(tsi[i] >= 0)
    ab_titer[i] = status[i]*max(sero_params$LLOD, 
                                rnorm(1, mean = sero_params$peak*exp(-sero_params$decay*tsi[i]), 
                                              sd = sero_params$sigma))
  }

  ## Output ----
  # -------------------------- -
  sero_data <- data.frame(ind_id = sample_ids, age = age_at_sampling, age_group = age_group,
                          time_of_infection = toi, 
                          time_since_infection = ifelse(tsi < 0, NA, tsi/365),
                          time_since_symptom_onset = ifelse(tsso < 0, NA, tsso/365),
                          sero_status = status, log_igg = ab_titer)

  saveRDS(sero_data,
          file = file.path(sero_params$output_dir, "sero_data.rds"))
  
  return(list(sero_data = sero_data, params = sero_params))
}