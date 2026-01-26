############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: script that load all required functions and scripts and contains the main modelling kernel.
#
# This script is distributed in the hope that it will be useful,but without 
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #

# load libraries 
library(progress) # progress bar

# load help function
source('lib/ibm_population.R')
source('lib/ibm_test.R')
source('lib/ibm_parameters.R')
source('lib/ibm_plot.R')

# main function to run the individual-based model based on social contact locations
run_ibm <- function(param, verbose = TRUE){
  
  ######################################################### #
  # DEFENSIVE CHECKS  ----
  ######################################################### #
  
  if(param$num_infected_seeds > param$pop_size){
    warning("ERROR: population size < number of infected seeds")
    return(NULL)
  }

  if(any(unlist(param[grepl('num',names(param)) | grepl('size',names(param)) |
                        grepl('age',names(param)) | grepl('prob',names(param))]) < 0)){
    warning("ERROR: negative values not allowed as function parameter")
    return(NULL)
  }
   
  if(!is.logical(unlist(param[grepl('bool',names(param))]))){
    warning("ERROR: 'bool_show_demographics', 'bool_add_baseline' and 'bool_return_prevelance' should be a boolean")
     return(NULL)
  }
  
  if(!dir.exists(param$output_dir)){
    dir.create(param$output_dir, recursive = TRUE)
  }
  
  ######################################################### #
  # INITIALIZE POPULATION & MODEL PARAMETERS  ----
  ######################################################### #
  
  # save start time
  time_start <- Sys.time()
  
  # initialize random number generator
  set.seed(param$rng_seed)
  
  # create a population matrix with:
  #   - age             the age of each individual
  #   - household_id    the household index of each individual
  #   - member_id       the household member index of each individual
  pop_data              <- create_population_matrix(param)
  
  # set contact and transmission parameters
  contact_prob_community         <- 1-exp(-param$num_contacts_community_day / param$pop_size)  # rate to probability
  transmission_prob_community    <- contact_prob_community * param$transmission_prob
  transmission_prob_household    <- param$contact_prob_household * param$transmission_prob
  transmission_prob_school       <- param$contact_prob_school    * param$transmission_prob
  transmission_prob_workplace    <- param$contact_prob_workplace * param$transmission_prob
    
  # set vaccine coverage (= fully protected)
  id_vaccinated                  <- sample(param$pop_size,param$pop_size*param$vaccine_coverage)
  pop_data$health[id_vaccinated] <- 'V'
  
  # introduce infected individuals in the population
  id_infected_seeds                             <- sample(which(pop_data$health=='S'),param$num_infected_seeds)
  pop_data$health[id_infected_seeds]            <- 'I'
  pop_data$time_of_infection[id_infected_seeds] <- 0
  
  # set recovery parameters
  recovery_rate        <- 1/param$num_days_infected
  recovery_probability <- 1-exp(-recovery_rate)      # convert rate to probability
  
  # set general mortality probability
  general_mortality_probability <- 1-exp(-param$general_mortality_rate)      # convert rate to probability
  
  # set disease-related mortality probability
  disease_mortality_probability <- 1-exp(-param$disease_mortality_rate)      # convert rate to probability
  
  # create matrix to log health states: one row per individual, one column per time step
  log_pop_data  <- matrix(NA,nrow=param$pop_size,ncol=param$num_days)
  
  ####################################### #
  # RUN THE MODEL        ----
  ####################################### #
  
  # init progress bar
  pb <- progress_bar$new(format = paste0("Run ",basename(param$output_dir),": [:bar] :percent ETA: :eta"),
                         total = param$num_days, clear = FALSE, width= 60)
  
  # LOOP OVER ALL DAYS
  for(day_i in 1:param$num_days)
  {
    
    # option to advance progress bar
    pb$tick()

    # step 2: identify infected individuals
    boolean_infected <- pop_data$health == 'I'   # = boolean TRUE/FALSE
    ind_infected     <- which(boolean_infected)  # = indices
    num_infected     <- length(ind_infected)     # = number
    
    # step 4: loop over all infected individuals
    p <- ind_infected[1]
    for(p in ind_infected)
    {
      
      # new infections are possible in the community and household
      transmission_prob_all <- is_susceptible(pop_data$health) * transmission_prob_community +
                               (is_susceptible(pop_data$health) & pop_data$hh_id[p]  == pop_data$hh_id)  * transmission_prob_household

      # add probability when in same class
      if(!is.na(pop_data$classroom_id[p])){
          flag_classroom <- is_susceptible(pop_data$health) & !is.na(pop_data$classroom_id) & pop_data$classroom_id[p] == pop_data$classroom_id
          transmission_prob_all[flag_classroom] <-  transmission_prob_all[flag_classroom] + transmission_prob_school
      }

      # add probability when at same workplace
      if(!is.na(pop_data$workplace_id[p])){
        flag_workplace <- is_susceptible(pop_data$health) & !is.na(pop_data$workplace_id) & pop_data$workplace_id[p] == pop_data$workplace_id
        transmission_prob_all[flag_workplace] <-  transmission_prob_all[flag_workplace] + transmission_prob_workplace
      }

      # account for vaccine-related protection
      transmission_prob_all[pop_data$health == 'V'] <- transmission_prob_all[pop_data$health == 'V'] * (1 - param$vaccine_effectiveness)
      
      # sample given the obtained probability
      flag_new_infection <- rbinom(param$pop_size, size = 1, prob = transmission_prob_all) == 1

      # mark new infected individuals
      pop_data$health[flag_new_infection] <- 'I'
      
      # log transmission details
      pop_data$infector[flag_new_infection]             <- p
      pop_data$infector_age[flag_new_infection]         <- pop_data$age[p]
      pop_data$time_of_infection[flag_new_infection]    <- day_i
      pop_data$secondary_cases[p]                       <- pop_data$secondary_cases[p] + sum(flag_new_infection)
      pop_data$generation_interval[flag_new_infection]  <- day_i - pop_data$time_of_infection[p]
    
    }
    
    # step 5: identify newly recovered individuals
    new_recovered <- boolean_infected & rbinom(param$pop_size, size = 1, prob = recovery_probability)
    pop_data$health[new_recovered] <- 'R'
    
    # step 6: identify newly deaths (can overrule recovery)
    pop_mortality_probability <- general_mortality_probability[pop_data$age] + 
                                  boolean_infected * disease_mortality_probability[pop_data$age]
    new_deaths <- rbinom(param$pop_size, size = 1, prob = pop_mortality_probability) == 1
    pop_data$health[new_deaths] <- 'D'
    
    # step 7: log population health states
    log_pop_data[,day_i] <- pop_data$health
    
  } # end for-loop for each day
  
  ####################################### #
  # PLOT RESULTS  ----
  ####################################### #
  
  # reformat the log matrix with one row per individual and one column per time step
  # 'colSums' = sum per column
  states <- c("S", "I", "R", "V", "D")
  log_health <- as.data.frame(sapply(
    states,
    \(x) colSums(log_pop_data == x) / param$pop_size
  ))
  
  
  # PRINT RESULTS AND PARAMETERS
  if(verbose) {

    # option to add the baseline prevalence, if requested and param are not default
    if(param$bool_add_baseline){
      out_baseline <- run_ibm_default()
    } else{
      out_baseline <- NA
    }
      
    # update figure configuration (use panels?)
    par(mfrow = param$plot_mfrow)
    
    # plot model output
    plot_health_states(log_health, param, out_baseline = out_baseline)
    plot_secondary_cases(pop_data,param)
    plot_generation_interval(pop_data,param)
    plot_transmission_matrix(pop_data = pop_data)
    
    # set back the default mfrow
    par(mfrow=c(1,1))
    
    # print model parameters
    print_model_parameters(param)
    
    # print model results
    print_model_results(log_health = log_health,
                        time_start = time_start,
                        out_baseline = out_baseline)
  }
  
  # save model results
  ibm_out <- list(log_health = log_health,
                  param = param)
  saveRDS(ibm_out, file = paste0(param$output_dir,'/health_time.rds'))
  saveRDS(pop_data, file = paste0(param$output_dir,'/pop_data.rds'))
  
  # return health states over time
  return(list(log_health = log_health, 
              param = param))
}

# help function to define whether a health state relates to susceptibility
is_susceptible <- function(vector_health){
  return(vector_health == 'S' | vector_health == 'V')
}

# main function to run the individual-based model based on social contact locations
run_ibm_default <- function(param, verbose = TRUE){
  
  
  # get default parameters 
  default_param <- get_default_parameters()
  
  # disable the 'add_baseline' option
  default_param$bool_add_baseline      <- FALSE
  default_param$bool_show_demographics <- FALSE
  default_param$output_dir             <- paste0(default_param$output_dir,'_default')
  
  # re-run the model with default parameters
  model_out <- run_ibm(param = default_param, verbose = FALSE)
  
  return(model_out)
}
  
  
  
# function to print the model results
print_model_results <- function(log_health,time_start,out_baseline=NA){
  
  bool_add_baseline <- !any(is.na(out_baseline))
  if(bool_add_baseline){
    # default epidemic characteristics
    default_ti <- paste0('   [baseline: ',round((out_baseline$log_health$I[length(out_baseline$log_health$I)] +
                                                   out_baseline$log_health$R[length(out_baseline$log_health$R)])*100,digits=0),'%]')
    default_pp <- paste0('   [baseline: ',round(max(out_baseline$log_health$I)*100,digits=0),'%]')
    default_pd <- paste0('    [baseline: ',which(out_baseline$log_health$I == max(out_baseline$log_health$I))[1],']')
  }
  
  # print total incidence
  print('-------------')
  print('MODEL RESULTS')
  
  print(paste0('total incidence: ',round((log_health$I[length(log_health$I)] + log_health$R[length(log_health$R)])*100,digits=0),'%',
               ifelse(bool_add_baseline,default_ti,'')))
  
  # print peak details
  print(paste0('Peak prevalence: ',round(max(log_health$I)*100,digits=0),'%',
               ifelse(bool_add_baseline,default_pp,'')))
  print(paste0('Peak day:        ',which(log_health$I == max(log_health$I))[1], 
               ifelse(bool_add_baseline,default_pd,'')))
  
  # print total run time
  total_time <- as.double(Sys.time() - time_start,unit='secs')
  print(paste0('Total run time:  ',round(total_time,digits=0),'s'))
  
}



