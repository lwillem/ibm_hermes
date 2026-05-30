############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: define natural mortality in the population.
#
# This script is distributed in the hope that it will be useful, but without
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 hputter, ldewreede, LUMC, THE NETHERLANDS
############################################################################ #

# ------------------------------------------------------------------------ -
# DEPENDENCIES ----
# ------------------------------------------------------------------------ -

source('lib/ibm_pop_death_data.R')

# ------------------------------------------------------------------------ -
# MAIN DATA GENERATING FUNCTION ----
# ------------------------------------------------------------------------ -

#' Population mortality
#'
#' @param date Calendar time  
#' @param age Current age of an individual
#' @param tbl Population life table (for males or females)
#'
#' @return A list with natural death time \code{death_time}.
#' @export

gen_pop <- function(date, age, tbl){
  
  ## Select relevant cohort ----
  # ------------------------------ -
  yr <- floor(date - age)
  
  if (yr > 2023){
    gpyr <- subset(tbl, Year == 2023)
  } else{
    gpyr <- subset(tbl, Year == yr)
  }

  fage <- floor(age)                  # age rounded below
  gpyrr <- subset(gpyr, Age >= fage)  # remaining ages
  
  ## Find t such that exp(-H(age + t)) = U, with U random Un(0, 1) ----
  # -------------------------------------------------------------------- -
  thr <- -log(runif(1))
  x <- gpyrr$Age
  x[1] <- age                   # replace first instance (fage) by actual age
  mx0 <- gpyrr$mx[-nrow(gpyrr)] # hazard rates until age 110
  mx1 <- gpyrr$mx[nrow(gpyrr)]  # hazard rate from age 110 onwards
  y <- c(0, cumsum(mx0 * diff(x)))
  idx <- max(which(y <= thr))
  x0 <- x[idx]
  if (x0 < 110) {
    y0 <- y[idx]
    y1 <- y[idx + 1]
    x1 <- x[idx + 1]
    m <- (y1 - y0) / (x1 - x0)
    res <- x0 + (thr - y0) / m
  } else {                      # going beyond 110, assuming constant rate mx1
    y0 <- y[idx]
    m <- mx1
    res <- x0 + (thr - y0) / m
  }
  
  ## Output ----
  # -------------------------- -
  return(res)
}
