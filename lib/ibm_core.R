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

# main function to run the individual-based model based on social contact locations
run_ibm_location <- function(param, 
                             verbose = TRUE){
  
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
  pb <- progress_bar$new(format = "  run all days [:bar] :percent ETA: :eta",
                         total = param$num_days, clear = FALSE, width= 60)
  
  # LOOP OVER ALL DAYS
  for(day_i in 1:param$num_days)
  {
    
    # option to advance progress bar
    if(verbose) {
      pb$tick()
    }
    
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
  
  if(param$bool_return_prevelance){
    return(list(log_health = log_health, 
                      param = param))
  }
  
  # change figure configuration => 4 sub-plots
  if(param$bool_single_plot) {
    par(mfrow=c(2,2))
  }
  
  # plot health states over time
  plot(log_health$S,
       type='l',
       xlab='Time (days)',
       ylab='Population fraction',
      # main='location-specific IBM',
       ylim=c(0,1),
       lwd=2)
  lines(log_health$I,  col=2,lwd=2)
  lines(log_health$R,  col=3,lwd=2)
  lines(log_health$V,  col=4,lwd=2)
  lines(log_health$D,  col=5,lwd=2)
  
  legend('top',legend=c('S','I','R','V','D'),col=1:5,lwd=2,ncol=5,cex=0.7,bg='white',
         inset = c(0, -0.3), xpd = NA) # push legend above the plot)
  
  # option to add the baseline prevalence, if requested and param are not default
  if(param$bool_add_baseline){

    # get default parameters 
    default_param <- get_default_parameters()
    
    # if parameters did not change, do not add baseline
    param_equal <- are_parameters_equal(param, default_param)
    
    if(param_equal){
      out_baseline <- NA
    } else{
      
    # disable the 'add_baseline' option
    default_param$bool_add_baseline      <- FALSE
    default_param$bool_return_prevelance <- TRUE
    default_param$bool_show_demographics <- FALSE
    
    # re-run the model with default parameters
    out_baseline <- run_ibm_location(param = default_param)
    
    # add output to the line plot
    lines(out_baseline$log_health$I, col=2,lwd=2,lty=2)
    legend('topright',legend=c('I (baseline)'),col=2,lwd=2,lty=3,cex=0.7,bg='white')
    }
    
  } else{
    out_baseline <- NA
  }
  
  if(all(is.na(pop_data$secondary_cases))){
    pop_data$secondary_cases <- -1
  }

  # plot secondary cases
  boxplot(secondary_cases ~ time_of_infection, data=pop_data,
          xlab='time of infection (day)',
          ylab='secondary cases',
          main='secondary cases',
          ylim=c(0,max(10,pop_data$secondary_cases)),
          xlim=c(0,param$num_days),
          xaxt='n')
  axis(1,seq(0,param$num_days,5))
  
  # plot generation interval
  if(all(is.na(pop_data$generation_interval))){
    pop_data$generation_interval <- -1
  }
  boxplot(generation_interval ~ time_of_infection, data=pop_data,
          xlab='time of infection (day)',
          ylab='generation interval (days)',
          main='generation interval',
          ylim=c(0,max(10,pop_data$generation_interval)),
          xlim=c(0,param$num_days),
          xaxt='n')
  axis(1,seq(0,param$num_days,5))
  
  ####################################### #
  # TRANSMISSION MATRIX    ----
  ####################################### #
  # use 3-year age classes
  max_age                 <- max(pop_data$age)
  age_cat                 <- seq(1,max_age,3)
  transmission_age_matrix <- table(cut(pop_data$infector_age,age_cat,right=F),cut(pop_data$age,age_cat,right=F),dnn = list('age infector','age contact'))
  
  # create plot title with some driving parameters
  plot_title <- paste('num_cnt_community',param$num_contacts_community_day,' || ',
                      'p_cnt_household', param$contact_prob_household, '\n',
                      'p_cnt_school',param$contact_prob_school,' || ',
                      'p_cnt_workplace',param$contact_prob_workplace, '\n',
                      'p_transmission',param$transmission_prob)
  
  if(max(transmission_age_matrix)>25){
    plot_breaks <- c(0:4,seq(5,25,10),max(c(30,transmission_age_matrix)))
  } else{
    plot_breaks <- pretty(transmission_age_matrix)
  }
    
  # plot the matrix with color coding
  plot_transmission_matrix(mij = transmission_age_matrix,
                           breaks = plot_breaks,
                           main = plot_title,
                           num.digits = NA)
  
  # PRINT PARAMETERS AND RESULTS ----
  if(verbose) {
  print('MODEL PARAMETERS')
  # loop over the given parameters, add name & value
  for(i_param in names(param)){
    p_values <- unlist(param[i_param])
    if(length(p_values)>10){
      p_values <- c(p_values[1:9],'...')
    } 
      print(paste0(i_param,': ',paste(p_values,collapse = ',')))
  }
  
  # print total incidence
  print_model_results(log_health = log_health,
                      time_start = time_start,
                      out_baseline = out_baseline)
  }
  
  # set back the default par(mfrow)
  par(mfrow=c(1,1))
  
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

plot_transmission_matrix <- function(mij,
                                     breaks = NA,
                                     main = 'Transmission matrix',
                                     num.digits = NA){
  
  # if no breaks given, use uniform values
  if(all(is.na(breaks))){
    breaks <- pretty(mij)
  } 
  
  # add additional starting break, to include the minimum value as starting point (in the legend)
  breaks <- c(min(breaks)-1,breaks)
  
  # set midpoints
  midpoints <- matrix(
    breaks[-length(breaks)] + diff(breaks) / 2,
    nrow = 1, ncol = length(breaks) - 1
  )

  # set colors
  num.colors <- length(breaks)-1
  redc <- heat.colors(num.colors)
  
  # get plot region for matrix and legend based on current graphical parameters
  # note: based on layout from fields::imagePlot
  char.size     <- par()$cin[1] / par()$din[1] # get text character size
  offset        <- char.size * par()$mar[4] # space between legend and main plot
  legend.width  <- 1
  legend.mar    <- 8.1
  legend.shrink <- 0.9
  cex.lab       <- 1.2
  cex.axis      <- 0.8
  cex.text      <- 1
  
  # set legends' plot region
  legend_plot_region <- par()$plt
  legend_plot_region[2] <- 1 - (legend.mar * char.size)
  legend_plot_region[1] <- legend_plot_region[2] - (legend.width * char.size)

  # account for legend.shrink
  pr <- (legend_plot_region[4] - legend_plot_region[3]) * ((1 - legend.shrink) / 2)
  legend_plot_region[4] <- legend_plot_region[4] - pr
  legend_plot_region[3] <- legend_plot_region[3] + pr

  # set main matrix' plot region
  main_plot_region    <- par()$plt
  main_plot_region[2] <- min(main_plot_region[2], legend_plot_region[1] - offset)

  # defensive check for main and legends' plot region
  dp <- legend_plot_region[2] - legend_plot_region[1]
  legend_plot_region[1] <- min(main_plot_region[2] + offset, legend_plot_region[1])
  legend_plot_region[2] <- legend_plot_region[1] + dp

  # store old graphical parameters, and initiate the ones for the main plot
  old.par <- par(no.readonly = TRUE)
  par(plt = main_plot_region)

  # add image plot
  image(mij,
        xlab = 'Age infector',
        ylab = 'Age contact',
        main = main,
        cex.lab = 1.2,
        breaks = breaks,
        col = redc,
        xaxt = "n",
        yaxt = "n")
  
  # add axis labels
  plt_ticks <- seq(0, 1, length = nrow(mij))
  axis(2, at = plt_ticks, labels = c(colnames(mij)), cex.axis = cex.axis, tick = FALSE, las = 1, mgp = c(3, 0.3, 0))
  axis(1, at = plt_ticks, labels = c(colnames(mij)), cex.axis = cex.axis, tick = FALSE, las = 2, mgp = c(3, 0.3, 0))

  # add numeric values if num.digits != NA and cex.text > 0
  if (!is.na(num.digits) && !is.na(cex.text) && cex.text > 0) {
    # format results (rounding/scientific)
    if (any(max(mij, na.rm = TRUE) > 1)) {
      mij <- round(mij, digits = num.digits)
    } else {
      mij <- format(mij, digits = num.digits)
    }
    
    # get grid centers and add values
    e_grid <- expand.grid(plt_ticks, plt_ticks)
    text(e_grid, labels = mij, cex = cex.text)
  }
  
  # set graphical parameters for the legend
  par(new = TRUE, pty = "m", plt = legend_plot_region, err = -1)

  # include legend bar with axis
  image(x = 1:2, y = 1:length(breaks), z = midpoints,
        xaxt = "n", yaxt = "n", xlab = "",
        ylab = "", col = redc,
        breaks = breaks)
  #axis(side = 4, at = 1:length(breaks[-1]), labels= breaks[-1], mgp = c(3, 1, 0), las = 2)
  mtext("count", side = 2, cex = 0.8)  
  breaks_label <- breaks
  if(any(diff(breaks)>1)){
    sel_bin <- c(FALSE,diff(breaks) > 1)
    breaks_label_bin <- breaks_label
    breaks_label_bin[-1] <- paste(breaks_label[-length(breaks_label)], breaks_label[-1], sep=' - ')
    breaks_label[sel_bin] <- breaks_label_bin[sel_bin]
  }
  axis(side = 4, at = 1:length(breaks[-1])+0.5, labels = breaks_label[-1], cex.axis = cex.axis, mgp = c(3, 0.3, 0), las = 2, tick = FALSE)
  
  # restore original graphical parameters
  par(old.par)
  
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



