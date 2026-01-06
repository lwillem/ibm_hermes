############################################################################ #
# This file is part of the SIMID course material
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) 2026 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #


library(progress)

# load the default parameters into the global environment (for development)
get_default_parameters <- function(){
  
  attach(list(pop_size = 1000,  # population size
              num_days = 50,    # time horizon
              num_infected_seeds = 3, # initial infections
              vaccine_coverage = 0.1, # vaccine state
              
              num_days_infected  = 7, # disease parameter
              transmission_prob  = 0.1, # transmission dynamics

              # school settings
              # note: we model (abstract) school contacts in our simulation
              #num_schools            = 2,          # number of classes per age group
              target_school_size     = 350,
              target_school_ages     = c(3:18),
              
              # social contact parameters
              num_contacts_community_day = 4,    # average number of "effective contacts" per day in the general community 
              contact_prob_household     = 1,    # probability for an "effective contact" at home (1 = fully connected)
              contact_prob_school        = 0.5,  # probability for an "effective contact" at school 
              contact_prob_workplace     = 0.1,  # probability for an "effective contact" at work 
              
              #num_workplaces = 100,
              target_workplace_size = 10,
              target_workplace_ages = c(19:65),
              
              rng_seed = 2020,
              
              #  visualisation parameters
              bool_show_demographics  = TRUE,   # option to show the demography figures (location IBM only)
              add_baseline            = FALSE,  # option to add the prevalence with default param
              return_prevelance       = FALSE,  # option to return the prevalence (and stop)
              bool_single_plot        = TRUE    # option to specify the plot layout using sub-panels (or not)
              ))
}

# function to print the model results
print_model_results <- function(log_i,log_r,time_start,out_baseline=NA){
  
  add_baseline <- !any(is.na(out_baseline))
  if(add_baseline){
    # default epidemic characteristics
    default_ti <- paste0('   [baseline: ',round((out_baseline$log_i[length(out_baseline$log_i)] +
                                                   out_baseline$log_r[length(out_baseline$log_r)])*100,digits=0),'%]')
    default_pp <- paste0('   [baseline: ',round(max(out_baseline$log_i)*100,digits=0),'%]')
    default_pd <- paste0('    [baseline: ',which(out_baseline$log_i == max(out_baseline$log_i))[1],']')
  }
 
  # print total incidence
  print('-------------')
  print('MODEL RESULTS')
  
  print(paste0('total incidence: ',round((log_i[length(log_i)] + log_r[length(log_r)])*100,digits=0),'%',
               ifelse(add_baseline,default_ti,'')))
  
  # print peak details
  print(paste0('Peak prevalence: ',round(max(log_i)*100,digits=0),'%',
               ifelse(add_baseline,default_pp,'')))
  print(paste0('Peak day:        ',which(log_i == max(log_i))[1], 
               ifelse(add_baseline,default_pd,'')))
  
  # print total run time
  total_time <- as.double(Sys.time() - time_start,unit='secs')
  print(paste0('Total run time:  ',round(total_time,digits=0),'s'))
  
}

# function to run the individual-based model based on social contact locations
run_ibm_location <- function(pop_size              = 2000,     # population size                         
                             num_days              = 50,       # number of days to simulate (time step = one day) 
                             num_infected_seeds    = 3,        # initial number of infected individuals   
                             vaccine_coverage      = 0 ,       # vaccine coverage [0,1]                
                             rng_seed              = as.numeric(format(Sys.time(),'%S')), # random number seed = current time (seconds)
                             
                             # schools =  number of classes per age group
                             target_school_size     = 350,
                             target_school_ages     = c(3:18),
                             
                             # workplaces = total number of workplaces
                             target_workplace_size = 10,
                             target_workplace_ages = 18:65,
                             
                             # social contact parameters
                             num_contacts_community_day = 4,    # average number of "effective contacts" per day in the general community 
                             contact_prob_household     = 1,    # probability for an "effective contact" at home (1 = fully connected)
                             contact_prob_school        = 0.5,  # probability for an "effective contact" at school 
                             contact_prob_workplace     = 0.3,  # probability for an "effective contact" at work  
                             
                             # disease parameters
                             num_days_infected     = 7,        # average number of days individuals are infected/infectious   
                             transmission_prob     = 0.1,      # transmission probability per social contact                  
                             
                             # visualisation parameters
                             bool_show_demographics = TRUE,    # option to show the demography figures
                             add_baseline           = FALSE,   # option to add the prevalence with default param
                             return_prevelance      = FALSE,   # option to return the prevalence (and stop)
                             bool_single_plot       = TRUE     # option to specify the plot layout using sub-pabels (or not)
                             ){
  
  ######################################################### #
  # DEFENSIVE CHECKS  ----
  ######################################################### #
  
  if(num_infected_seeds > pop_size){
    warning("ERROR: population size < number of infected seeds")
    return(NULL)
  }

  if(any(c(pop_size,num_days,num_infected_seeds,vaccine_coverage,target_school_size,
           target_school_ages,target_workplace_size,target_workplace_ages,
           num_contacts_community_day,contact_prob_household,contact_prob_school,
           contact_prob_workplace,num_days_infected,transmission_prob)<0)){
    warning("ERROR: negative values not allowed as function parameter")
    return(NULL)
  }
  
  if(!is.logical(c(bool_show_demographics,add_baseline,return_prevelance))){
    warning("ERROR: 'bool_show_demographics', 'add_baseline' and 'return_prevelance' should be a boolean")
    return(NULL)
  }
  
  ######################################################### #
  # INITIALIZE POPULATION & MODEL PARAMETERS  ----
  ######################################################### #
  
  # save start time
  time_start <- Sys.time()
  
  # initialize random number generator
  set.seed(rng_seed)
  
  # create a population matrix with:
  #   - age             the age of each individual
  #   - household_id    the household index of each individual
  #   - member_id       the household member index of each individual
  pop_data              <- create_population_matrix(pop_size, 
                                                    target_school_size, 
                                                    target_school_ages, 
                                                    target_workplace_size, 
                                                    target_workplace_ages,
                                                    bool_show_demographics)
  # set contact and transmission parameters
  contact_prob_community         <- 1-exp(-num_contacts_community_day / pop_size)  # rate to probability
  transmission_prob_community    <- contact_prob_community * transmission_prob
  transmission_prob_household    <- contact_prob_household * transmission_prob
  transmission_prob_school       <- contact_prob_school    * transmission_prob
  transmission_prob_workplace    <- contact_prob_workplace * transmission_prob
    
  # set vaccine coverage
  id_vaccinated                  <- sample(pop_size,pop_size*vaccine_coverage)
  pop_data$health[id_vaccinated] <- 'V'
  
  # introduce infected individuals in the population
  id_infected_seeds                             <- sample(which(pop_data$health=='S'),num_infected_seeds)
  pop_data$health[id_infected_seeds]            <- 'I'
  pop_data$time_of_infection[id_infected_seeds] <- 0
  
  # set recovery parameters
  recovery_rate        <- 1/num_days_infected
  recovery_probability <- 1-exp(-recovery_rate)      # convert rate to probability
  
  # create matrix to log health states: one row per individual, one column per time step
  log_pop_data  <- matrix(NA,nrow=pop_size,ncol=num_days)
  
  
  ####################################### #
  # RUN THE MODEL        ----
  ####################################### #
  
  # init progress bar
  pb <- progress_bar$new(format = "  run all days [:bar] :percent eta: :eta",
                         total = num_days, clear = FALSE, width= 60)
  
  # LOOP OVER ALL DAYS
  for(day_i in 1:num_days)
  {
    
    # advance progress bar
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
      transmission_prob_all <- (pop_data$health == 'S') * transmission_prob_community +
                               ((pop_data$health == 'S') & pop_data$hh_id[p]  == pop_data$hh_id)  * transmission_prob_household

      # add probability when in same class
      if(!is.na(pop_data$classroom_id[p])){
          flag_classroom <- pop_data$health == 'S' & !is.na(pop_data$classroom_id) & pop_data$classroom_id[p] == pop_data$classroom_id
          transmission_prob_all[flag_classroom] <-  transmission_prob_all[flag_classroom] + transmission_prob_school
      }

      # add probability when at same workplace
      if(!is.na(pop_data$workplace_id[p])){
        flag_workplace <- pop_data$health == 'S' & !is.na(pop_data$workplace_id) & pop_data$workplace_id[p] == pop_data$workplace_id
        transmission_prob_all[flag_workplace] <-  transmission_prob_all[flag_workplace] + transmission_prob_workplace
      }

      # sample given the obtained probability
      flag_new_infection <- rbinom(pop_size, size = 1, prob = transmission_prob_all) == 1

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
    new_recovered <- boolean_infected & rbinom(pop_size, size = 1, prob = recovery_probability)
    pop_data$health[new_recovered] <- 'R'
    
    # step 6: log population health states
    log_pop_data[,day_i] <- pop_data$health
    
  } # end for-loop for each day
  
  ####################################### #
  # PLOT RESULTS  ----
  ####################################### #
  # reformat the log matrix with one row per individual and one column per time step
  # 'colSums' = sum per column
  log_s <- colSums(log_pop_data == 'S')  / pop_size
  log_i <- colSums(log_pop_data == 'I')  / pop_size
  log_r <- colSums(log_pop_data == 'R')  / pop_size
  log_v <- colSums(log_pop_data == 'V')  / pop_size
  
  if(return_prevelance){
    return(data.frame(log_i=log_i,log_r=log_r))
  }
  
  # change figure configuration => 4 sub-plots
  if(bool_single_plot) {
    par(mfrow=c(2,2))
  }
  
  # plot health states over time
  plot(log_s,
       type='l',
       xlab='Time (days)',
       ylab='Population fraction',
       main='location-specific IBM',
       ylim=c(0,1),
       lwd=2)
  lines(log_i,  col=2,lwd=2)
  lines(log_r,  col=3,lwd=2)
  lines(log_v,  col=4,lwd=2)
  
  legend('top',legend=c('S','I','R','V'),col=1:4,lwd=2,ncol=2,cex=0.7,bg='white')
  
  if(add_baseline){
    out_baseline <- run_ibm_location(rng_seed=rng_seed, return_prevelance = T, bool_show_demographics=F)
    lines(out_baseline$log_i,  col=2,lwd=2,lty=2)
    legend('topright',legend=c('I (baseline)'),col=2,lwd=2,lty=3,cex=0.7,bg='white')
  } else{
    out_baseline <- NA
  }
  
  if(all(is.na(pop_data$secondary_cases))){
    pop_data$secondary_cases <- -1
  }

  boxplot(secondary_cases ~ time_of_infection, data=pop_data,
          xlab='time of infection (day)',
          ylab='secondary cases',
          main='secondary cases',
          ylim=c(0,10),
          xlim=c(0,num_days),
          xaxt='n')
  axis(1,seq(0,num_days,5))
  
  if(all(is.na(pop_data$generation_interval))){
    pop_data$generation_interval <- -1
  }
  boxplot(generation_interval ~ time_of_infection, data=pop_data,
          xlab='time of infection (day)',
          ylab='generation interval (days)',
          main='generation interval',
          ylim=c(0,10),
          xlim=c(0,num_days),
          xaxt='n')
  axis(1,seq(0,num_days,5))
  
  ####################################### #
  # TRANSMISSION MATRIX    ----
  ####################################### #
  # use 3-year age classes
  max_age                 <- max(pop_data$age)
  age_cat                 <- seq(1,max_age,3)
  transmission_age_matrix <- table(cut(pop_data$infector_age,age_cat,right=F),cut(pop_data$age,age_cat,right=F),dnn = list('age infector','age contact'))
  
  # create plot title with some driving parameters
  plot_title <- paste('num_cnt_community',num_contacts_community_day,' || ',
                      'p_cnt_household', contact_prob_household, '\n',
                      'p_cnt_school',contact_prob_school,' || ',
                      'p_cnt_workplace',contact_prob_workplace, '\n',
                      'p_transmission',transmission_prob)
  
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
  # collect possible parameter names
  all_param <- c('pop_size','num_days' ,'num_infected_seeds','vaccine_coverage',
                 'rng_seed', 'num_days_infected','transmission_prob',
                 'num_contacts_community_day','contact_prob_household','contact_prob_school','contact_prob_workplace',
                 'num_schools','target_school_ages','num_workplaces','bool_show_demographics'
  )
  
  print('MODEL PARAMETERS')
  # loop over the given parameter names, if present, add name & value
  for(i_param in all_param){
    if(exists(i_param)){
      print(paste0(i_param,': ',paste(get(i_param),collapse = ',')))
    }
  }
  
  # print total incidence
  print_model_results(log_i = log_i,
                      log_r = log_r,
                      time_start=time_start,
                      out_baseline = out_baseline)

  # set back the defaul par(mfrow)
  par(mfrow=c(1,1))
}

#' @title Create a synthetic population with households
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
create_population_matrix <- function(pop_size, 
                                     target_school_size, 
                                     target_school_ages, 
                                     target_workplace_size, 
                                     target_workplace_ages,
                                     bool_show_demographics = TRUE)
{
  ## demographic parameters
  ages_adult <- 19:80
  ages_child <- 1:18
  adult_age_tolerance     <- -5:5    # age tolerance between adults
  child_age_tolerance     <- 1:4    # age tolerance between children
  household_age_gap       <- 22:35     # min age gap between adults and youngest child
  
  # start with households of size 4, with 20% extra (so we can delete some of them later)
  #num_hh <- ceiling(pop_size / 4 * 1.20)
  num_hh <- ceiling(pop_size / 2)
  
  ## create adult 1
  pop_data_adult1 <- data.frame(age = sample(ages_adult, num_hh, replace = TRUE),
                                hh_id = 1:num_hh,
                                member_id = 1)
  
  ## create adult 2: based on age gap with adult 1
  pop_data_adult2 <- data.frame(age = pop_data_adult1$age + sample(adult_age_tolerance, num_hh, replace = TRUE),
                                hh_id = 1:num_hh,
                                member_id = 2)
  
  ## create child 1, based on age gap with youngest parent
  pop_data_child1 <- data.frame(age = pmax(pop_data_adult1$age,pop_data_adult2$age) - sample(household_age_gap, num_hh, replace = TRUE),
                                hh_id = 1:num_hh,
                                member_id = 3)
  
  ## create child 2: based on age gap with sibling
  pop_data_child2 <- data.frame(age = pop_data_child1$age - sample(child_age_tolerance, num_hh, replace = TRUE),
                                hh_id = 1:num_hh,
                                member_id = 4)
  
  # combine child data
  pop_data_child <- rbind(pop_data_child1,
                          pop_data_child2)
  
  # remove child ages < 0 or above 18
  pop_data_child <- pop_data_child[pop_data_child$age > 0  & 
                                     pop_data_child$age < min(ages_adult),]
  
  # create the population
  pop_data         <- rbind(pop_data_adult1,
                            pop_data_adult2,
                            pop_data_child)  # start from empty matrix
  dim(pop_data)
  
  # if needed, select all individuals within the given population size
  if(nrow(pop_data) > pop_size){
    pop_data <- pop_data[sample(1:nrow(pop_data),pop_size),]
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
  num_schools <- ceiling(sum(pop_data$age %in% target_school_ages) / target_school_size)
  pop_data    <- set_schools(pop_data, num_schools, target_school_ages)
  
  # set workplace 
  num_workplaces <- ceiling(sum(pop_data$age %in% target_workplace_ages) / target_workplace_size)
  pop_data <- set_workplaces(pop_data, num_workplaces, target_school_ages)
  
  # option to plot demographics
  if(bool_show_demographics){
    # create a figure with 8 subplots
    par(mfrow=c(2,4))
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
  
  return(pop_data)
  
} # end function

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



