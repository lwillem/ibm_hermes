############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: create and visualise a synthetic population
#
# This script is distributed in the hope that it will be useful, but without
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #

#' Create a synthetic population
#'
#' @param params Model parameter list.
#' @return Data frame describing individuals and their attributes.
create_population_matrix <- function(params) {

  num_households <- ceiling(params$pop_size / 2)

  adult1 <- data.frame(
    age = sample(params$ages_adult, num_households, replace = TRUE),
    hh_id = seq_len(num_households),
    hh_role = 1
  )

  adult2 <- data.frame(
    age = adult1$age +
      sample(-params$adult_age_tolerance:params$adult_age_tolerance,
             num_households, replace = TRUE),
    hh_id = seq_len(num_households),
    hh_role = 2
  )

  child1 <- data.frame(
    age = pmax(adult1$age, adult2$age) -
      sample(params$household_age_gap, num_households, replace = TRUE),
    hh_id = seq_len(num_households),
    hh_role = 3
  )

  child2 <- data.frame(
    age = child1$age -
      sample(seq_len(params$child_age_tolerance), num_households, replace = TRUE),
    hh_id = seq_len(num_households),
    hh_role = 4
  )

  children <- rbind(child1, child2)
  children <- children[children$age > 0 &
                         children$age < min(params$ages_adult), ]

  pop_data <- rbind(adult1, adult2, children)

  if (nrow(pop_data) > params$pop_size) {
    pop_data <- pop_data[sample(seq_len(nrow(pop_data)), params$pop_size), ]
  }

  num_schools <- ceiling(sum(pop_data$age %in% params$target_school_ages) /
                           params$target_school_size)
  pop_data <- set_schools(pop_data, num_schools, params$target_school_ages)

  num_workplaces <- ceiling(sum(pop_data$age %in% params$target_workplace_ages) /
                              params$target_workplace_size)
  pop_data <- set_workplaces(pop_data, num_workplaces, params$target_workplace_ages)

  pop_data <- data.frame(
    pop_data,
    health = "S",
    infector = NA,
    time_of_infection = NA,
    generation_interval = 0,
    secondary_cases = 0,
    stringsAsFactors = FALSE
  )
  
  if (params$bool_show_demographics) {
    plot_population_histograms(pop_data)
  }

  # remove househole role
  pop_data$hh_role <- NULL
  
  return(pop_data)
}


#' Explore demography distributions of the given synthetic population
#'
#' @param pop_data Model population matrix.
#' @param opt_mfrow mfrow parameter, default c(2,4)
plot_population_histograms <- function(pop_data, opt_mfrow = c(2,4)){
  
  # create a figure with 8 subplots
  par(mfrow = opt_mfrow)
  
  # get max age
  pop_age_max <- max(pop_data$age)
  hist(pop_data$age,-1:pop_age_max,main='total population',xlab='age')
  hist(pop_data$age[pop_data$hh_role==1],-1:pop_age_max,main='adult 1',xlab='age')
  hist(pop_data$age[pop_data$hh_role==2],-1:pop_age_max,main='adult 2',xlab='age')
  hist(pop_data$age[pop_data$hh_role==3],-1:pop_age_max,main='child 1',xlab='age')
  hist(pop_data$age[pop_data$hh_role==4],-1:pop_age_max,main='child 2',xlab='age')
  hist(table(pop_data$hh_id),main='household size',xlab='household size')
  
  # check class and workplace size
  if(any(!is.na(pop_data$classroom_id))) hist(table(pop_data$classroom_id),xlab='Size',main='School class size')
  if(any(!is.na(pop_data$workplace_id))) hist(table(pop_data$workplace_id),xlab='Size',main='Worplace size')
}


#' Define school ids for the given synthetic population
#'
#' @param pop_data model population matrix.
#' @param num_schools number of schools to set
#' @param target_school_ages ages to include in school setting
set_schools <- function(pop_data, num_schools, target_school_ages){
  if(num_schools > 0){
    # eg. 'class3_1' is the 1th classroom with 3-year olds children
    pop_data$classroom_id <- paste0('class', pop_data$age, '_', sample(num_schools, nrow(pop_data), replace = TRUE))
    
    # set 'classroom_id' for infants and adults outside the target ages to 'NA' (=none)
    boolean_school_pop    <- pop_data$age %in% target_school_ages
    pop_data$classroom_id[!boolean_school_pop] <- NA
  } else {
    pop_data$classroom_id <- NA
  }
  # return the result
  return(pop_data)
}

#' Define workplace id for the given synthetic population
#'
#' @param pop_data model population matrix.
#' @param num_workplaces number of workplaces to set
#' @param target_school_ages ages to include in school setting, hence excluded from workplace
set_workplaces <- function(pop_data, num_workplaces, target_workplace_ages){
  if(num_workplaces > 0){
    # sample a workplace for each individual
    pop_data$workplace_id <- sample(num_workplaces, nrow(pop_data), replace = TRUE)
    
    # set 'workplace_id' for individuals with an age outside target_workplace_ages to 'NA' (= no workplace)
    boolean_workplace_pop <- pop_data$age %in% target_workplace_ages
    pop_data$workplace_id[!boolean_workplace_pop] <- NA    
  } else {
    pop_data$workplace_id <- NA
  }
  # return the result
  return(pop_data)
}
