############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: retrieving mortality data from the generated data
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
#' @param current_date Current date of evaluation
#'
#' @return A list with elements selected from population data \code{demo_data}.
#' @export
obtain_vaccination_data <- function(pop_data, current_date) {
  
  ## Calculate individual ids
  #-------------------------- -
  ind_id <- rownames(pop_data)
  
  ## Output ----
  # -------------------------- -
  vaccination_data <- data.frame(ind_id = ind_id, 
                               time_of_vaccination = pop_data$time_of_vaccination,
                               time = current_date)
  vaccination_data <- subset(vaccination_data, !is.na(time_of_vaccination))
  vaccination_data <- subset(vaccination_data, time_of_vaccination < current_date)
  
  return(vaccination_data)
}
