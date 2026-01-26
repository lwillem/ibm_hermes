############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: to load the default parameters and interact with the user on parameter values.
# 
# This script is distributed in the hope that it will be useful,but without 
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #

# load the default parameters into the global environment (for development)
get_default_parameters <- function(){
  
param <- list(pop_size           = 2000,  # population size
              num_days           = 50,    # time horizon
              num_infected_seeds = 3,     # initial number of infections
              vaccine_coverage   = 0.2,   # initial proportion in vaccine state
              vaccine_effectiveness = 0.8, # vaccine-related protection against infection
              
              num_days_infected  = 7, # disease parameter
              transmission_prob  = 0.1, # transmission dynamics

              # demographic parameters
              ages_adult = 19:80,
              adult_age_tolerance     = 5,     # age tolerance between adults (both ways)
              child_age_tolerance     = 4,     # age tolerance between children (1 to 'child_age_tolerance)
              household_age_gap       = 22:35, # min age gap between adults and youngest child
              
              # school settings
              target_school_size = 350,     # one class per age group in each school
              target_school_ages = c(3:18), # define the age groups in school
              
              # social contact parameters
              num_contacts_community_day = 4,    # average number of "effective contacts" per day in the general community 
              contact_prob_household     = 0.9,  # probability for an "effective contact" at home (1 = fully connected)
              contact_prob_school        = 0.5,  # probability for an "effective contact" at school 
              contact_prob_workplace     = 0.1,  # probability for an "effective contact" at work 
              
              # workplace settings
              target_workplace_size = 10,
              target_workplace_ages = c(19:65),
              
              # random number generator settings
              rng_seed = 2026,
              
              # mortality parameters, by age
              general_mortality_rate = c(rep(0,60),seq(0,0.01,length = 30)),
              disease_mortality_rate = c(rep(0.01,60),seq(0.01,0.7,length = 30)),
              
              #  visualisation parameters
              bool_show_demographics  = TRUE,   # option to show the demography figures
              bool_add_baseline       = FALSE,  # option to add the prevalence with default param
              bool_return_prevelance  = FALSE,  # option to return the prevalence (and stop)

              # output_tag
              output_dir = 'output/ibm_flu',
              plot_mfrow = c(2,2)
              )
  
  return(param)
}

# function to print the parameter options
print_parameter_names <- function(){
  
  param <- get_default_parameters()
  
  print(names(param))
  
}

print_model_parameters <- function(param){
  print('MODEL PARAMETERS')
  # loop over the given parameters, add name & value
  for(i_param in names(param)){
    p_values <- unlist(param[i_param])
    if(length(p_values)>10){
      p_values <- c(p_values[1:9],'...')
    } 
    print(paste0(i_param,': ',paste(p_values,collapse = ',')))
  }
}

# Check parameter equality (exclude a specific subet)
are_parameters_equal <- function(param1, param2, param_exclude = c('output_dir'), verbose = FALSE){
  
  params_equal <- TRUE
  
  for (param_name in names(param1)) {
    if(!param_name %in% param_exclude)
      if (!identical(
        unlist(param1[[param_name]]),
        unlist(param2[[param_name]])
      )) {
        if(verbose) message("Parameter changed: ", param_name)
        params_equal <- FALSE
      }
  }
  
  return(params_equal)
}

