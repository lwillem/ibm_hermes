############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: to plot model output
#
# Copyright (C) 2026 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #

plot_health_states <- function(log_health, param, out_baseline = NA){
  
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
  if(any(!is.na(out_baseline))){
    
      # add baseline output to the line plot
      lines(out_baseline$log_health$I, col=2,lwd=2,lty=2)
      legend('topright',legend=c('I (baseline)'),col=2,lwd=2,lty=3,cex=0.7,bg='white')
  }
  
}

plot_secondary_cases <- function(pop_data,param){
  
  # defensive check
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
}


plot_generation_interval <- function(pop_data, param){
  
  # defensive check
  if(all(is.na(pop_data$generation_interval))){
    pop_data$generation_interval <- -1
  }
  
  # plot generation interval
  boxplot(generation_interval ~ time_of_infection, data=pop_data,
          xlab='time of infection (day)',
          ylab='generation interval (days)',
          main='generation interval',
          ylim=c(0,max(10,pop_data$generation_interval)),
          xlim=c(0,param$num_days),
          xaxt='n')
  axis(1,seq(0,param$num_days,5))
  
}

plot_transmission_matrix <- function(pop_data,
                                     breaks = NA,
                                     main = 'Transmission matrix',
                                     num.digits = NA){
  
  # use 3-year age classes
  max_age <- max(pop_data$age)
  age_cat <- seq(1,max_age,3)
  mij     <- table(cut(pop_data$infector_age,age_cat,right=F),
                   cut(pop_data$age,age_cat,right=F),
                   dnn = list('age infector','age contact'))
  
  if(max(mij)>25){
    plot_breaks <- c(0:4,seq(5,25,10),max(c(30,mij)))
  } else{
    plot_breaks <- pretty(mij)
  }
  
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

