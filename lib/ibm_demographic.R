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
#'
#' @return A list with elements selected from population data \code{demo_data}.
#' @export
obtain_demographic_data <- function(pop_data) {
  
  ## Calculate individual ids
  #-------------------------- -
  ind_id <- rownames(pop_data)
  
  ## Output ----
  # -------------------------- -
  demographic_data <- data.frame(ind_id = ind_id, 
                                 subset(pop_data, 
                                        select = c("age", "sex", "hh_id", 
                                                      "classroom_id", 
                                                      "workplace_id")),
                                 time = rep(0, length(ind_id)))
  
  return(list(demographic_data = demographic_data))
}
