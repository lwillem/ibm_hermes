############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: functions to create and visualise a model population
#
# This script is distributed in the hope that it will be useful,but without 
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #

#' @title Create a synthetic population with explicit households
#'
#' @description This function creates a population with households
#'
#' @param pop_size  the final population size
#' @param num_schools the number of schools (which contains one class per age group)
#' @param num_workplaces the number of workplaces in the population
#'
#' @keywords external
#' @export
#pop_size <- 1e4
create_population_matrix <- function(param)
{
  # define a target number of households
  # start with households of size 4, with 20% extra (so we can delete some of them later)
  num_hh <- ceiling(param$pop_size / 2)
  
  ## create adult 1
  pop_data_adult1 <- data.frame(age = sample(param$ages_adult, num_hh, replace = TRUE),
                                hh_id = 1:num_hh,
                                member_id = 1)
  
  ## create adult 2: based on age gap with adult 1
  pop_data_adult2 <- data.frame(age = pop_data_adult1$age + sample(-param$adult_age_tolerance:param$adult_age_tolerance, num_hh, replace = TRUE),
                                hh_id = 1:num_hh,
                                member_id = 2)
  
  ## create child 1, based on age gap with youngest parent
  pop_data_child1 <- data.frame(age = pmax(pop_data_adult1$age,pop_data_adult2$age) - sample(param$household_age_gap, num_hh, replace = TRUE),
                                hh_id = 1:num_hh,
                                member_id = 3)
  
  ## create child 2: based on age gap with sibling
  pop_data_child2 <- data.frame(age = pop_data_child1$age - sample(1:param$child_age_tolerance, num_hh, replace = TRUE),
                                hh_id = 1:num_hh,
                                member_id = 4)
  
  # combine child data
  pop_data_child <- rbind(pop_data_child1,
                          pop_data_child2)
  
  # remove child ages < 0 or above 18
  pop_data_child <- pop_data_child[pop_data_child$age > 0  & 
                                     pop_data_child$age < min(param$ages_adult),]
  
  # create the population
  pop_data         <- rbind(pop_data_adult1,
                            pop_data_adult2,
                            pop_data_child)  # start from empty matrix
  dim(pop_data)
  
  # if needed, select all individuals within the given population size
  if(nrow(pop_data) > param$pop_size){
    pop_data <- pop_data[sample(1:nrow(pop_data),param$pop_size),]
  }
  
  # add health state: susceptible
  pop_data <- data.frame(pop_data,
                         health              = 'S',           # column to store the health state
                         infector            = NA,            # column to store the source of infection
                         time_of_infection   = NA,            # column to store the time of infection
                         generation_interval = 0,             # column to store the generation interval
                         secondary_cases     = 0,             # column to store the number of secondary cases
                         stringsAsFactors    = FALSE)
  
  # set school classes by age and number of schools
  num_schools <- ceiling(sum(pop_data$age %in% param$target_school_ages) / param$target_school_size)
  pop_data    <- set_schools(pop_data, num_schools, param$target_school_ages)
  
  # set workplace 
  num_workplaces <- ceiling(sum(pop_data$age %in% param$target_workplace_ages) / param$target_workplace_size)
  pop_data <- set_workplaces(pop_data, num_workplaces, param$target_school_ages)
  
  # option to plot demographics
  if(param$bool_show_demographics){
    plot_population_histograms(pop_data)
  }
  
  return(pop_data)
  
} # end function


plot_population_histograms <- function(pop_data, opt_mfrow = c(2,4)){
  
  # create a figure with 8 subplots
  par(mfrow = opt_mfrow)
  
  # get max age
  pop_age_max <- max(pop_data$age)
  hist(pop_data$age,-1:pop_age_max,main='total population',xlab='age')
  hist(pop_data$age[pop_data$member_id==1],-1:pop_age_max,main='adult 1',xlab='age')
  hist(pop_data$age[pop_data$member_id==2],-1:pop_age_max,main='adult 2',xlab='age')
  hist(pop_data$age[pop_data$member_id==3],-1:pop_age_max,main='child 1',xlab='age')
  hist(pop_data$age[pop_data$member_id==4],-1:pop_age_max,main='child 2',xlab='age')
  hist(table(pop_data$hh_id),main='household size',xlab='household size')
  
  # check class and workplace size
  if(any(!is.na(pop_data$classroom_id))) hist(table(pop_data$classroom_id),xlab='Size',main='School class size')
  if(any(!is.na(pop_data$workplace_id))) hist(table(pop_data$workplace_id),xlab='Size',main='Worplace size')
}

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

## workplaces
set_workplaces <- function(pop_data, num_workplaces, target_school_ages){
  if(num_workplaces > 0){
    # sample a workplace for each individual
    pop_data$workplace_id <- sample(num_workplaces, nrow(pop_data), replace = TRUE)
    
    # set 'workplace_id' for children to 'NA' (= no workplace)
    boolean_workplace_pop <- pop_data$age <= max(target_school_ages)
    pop_data$workplace_id[boolean_workplace_pop] <- NA    
  } else {
    pop_data$workplace_id <- NA
  }
  # return the result
  return(pop_data)
}

