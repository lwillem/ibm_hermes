############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: retrieving contact tracing data from the generated data
#
# This script is distributed in the hope that it will be useful, but without
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 sabrams, SIMID, UHASSELT, UNIVERSITY OF ANTWERP, BELGIUM
#                    hputter, LUMC, THE NETHERLANDS
############################################################################ #

# ------------------------------------------------------------------------ -
# MAIN DATA GENERATING FUNCTION ----
# ------------------------------------------------------------------------ -

#' Reporting delays
#'
#' @param pop_data Population data with infection histories and times of death.
#' @param delay_params List of parameters for the delay distribution.
#' @param current_date The date at which the incidence is reported
#' @param verbose Logical; if TRUE, plots and summaries are produced.
#'
#' @return A list with elements \code{log_sero}.
#' @export
obtain_incidence_data <- function(pop_data, delay_params, current_date, verbose = TRUE) {
  
  ## Reproducibility ----
  # ------------------------- -
  set.seed(delay_params$rng_seed)
  
  ## Defensive checks ----
  # -------------------------- -
  if (delay_params$max_delay < 0){
    warning("Maximum delay is negative")
    return(NULL)
  }
  
  ## Define incidence data ----
  #--------------------------- -
  incidence_data <- pop_data
  
  ## Reporting delays (truncation at truncation limit (in days) ----
  # -------------------------------------------------------------- -
  report_distribution <- dpois(0:delay_params$max_delay, lambda = delay_params$mean)
  trunc_report_distribution <- report_distribution / sum(report_distribution)
  
  whsymptoms <- which(!is.na(pop_data$time_of_symptom_onset))
  nsymptoms <- length(whsymptoms)
  
  # Reporting delays for symptomatic cases
  rep_delay <- sample(0:delay_params$max_delay, size = nsymptoms, replace = TRUE, 
                      prob = trunc_report_distribution)

  # Reported time of infection will be in a separate column
  incidence_data$time_of_reporting <- NA
  incidence_data$time_of_reporting[whsymptoms] <- incidence_data$time_of_symptom_onset[whsymptoms] + 
    rep_delay
  
  # Some of the symptomatic cases are missing altogether (i.e., reporting date set to NA)
  rep_miss <- sample(0:1, size = nsymptoms, replace = TRUE, 
                     prob = c(1-delay_params$prob_missing, delay_params$prob_missing))
  
  incidence_data$time_of_reporting[whsymptoms][which(rep_miss == 1)] <- NA
  
  if (verbose) {
    plot(incidence_data$time_of_infection, incidence_data$time_of_reporting,
         xlab = "Time of infection", ylab = "Time of reporting")
  }
  
  ## Calculate time stamp
  #-------------------------- -
  time_stamp <- current_date
  
  ## Output ----
  # -------------------------- -
  incidence_data <- subset(incidence_data, !is.na(incidence_data$time_of_reporting))
  incidence_data <- subset(incidence_data, incidence_data$time_of_symptom_onset < current_date)
  incidence_data <- subset(incidence_data, 
                           select = -c(health, infector, infector_id, infector_age,
                                       time_of_infection, time_of_natural_death,
                                       age_natural_death, symptom_onset, 
                                       generation_interval, secondary_cases))
  ind_id <- rownames(incidence_data)
  
  incidence_data <- data.frame(ind_id = ind_id,
                               subset(incidence_data,
                                      select = c("time_of_symptom_onset",
                                                 "time_of_reporting",
                                                 "time_of_hospitalization")),
                               time = time_stamp)
  rownames(incidence_data)
  saveRDS(incidence_data,
          file = file.path(delay_params$output_dir, "incidence_data.rds"))
  
  return(list(incidence_data = incidence_data, params = delay_params))
}
