############################################################################ #
# This file is part of the individual-based model framework called HERMES
#
# Goal: plotting utilities for IBM output
#
# This script is distributed in the hope that it will be useful, but without
# any warranty; See the LICENCE.txt for more details.
#
# Copyright (C) 2026 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
############################################################################ #

#' Plot population health states over time
#'
#' @param log_health Data frame with time series of health states.
#' @param params Model parameters.
#' @param out_baseline Optional baseline output for comparison.
plot_health_states <- function(log_health, params, out_baseline = NA) {

  plot(log_health$S, type = "l", xlab = "Time (days)",
       ylab = "Population fraction", ylim = c(0, 1), lwd = 2)
  lines(log_health$I, col = 2, lwd = 2)
  lines(log_health$R, col = 3, lwd = 2)
  lines(log_health$V, col = 4, lwd = 2)
  lines(log_health$D, col = 5, lwd = 2)

  legend("top", legend = c("S", "I", "R", "V", "D"),
         col = 1:5, lwd = 2, ncol = 5, cex = 0.7, bg = "white",
         inset = c(0, -0.3), xpd = NA)

  if (any(!is.na(out_baseline))) {
    lines(out_baseline$log_health$I, col = 2, lwd = 2, lty = 2)
    legend("topright", legend = "I (baseline)", col = 2,
           lwd = 2, lty = 2, cex = 0.7, bg = "white")
  }
}

#' Plot secondary cases over time
plot_secondary_cases <- function(pop_data, params) {

  if (all(is.na(pop_data$secondary_cases))) {
    pop_data$secondary_cases <- -1
  }

  boxplot(secondary_cases ~ time_of_infection, data = pop_data,
          xlab = "Time of infection (day)",
          ylab = "Secondary cases",
          main = "Secondary cases",
          ylim = c(0, max(10, pop_data$secondary_cases)),
          xlim = c(0, params$num_days),
          xaxt = "n")
  axis(1, seq(0, params$num_days, 5))
}

#' Plot generation interval distribution
plot_generation_interval <- function(pop_data, params) {

  if (all(is.na(pop_data$generation_interval))) {
    pop_data$generation_interval <- -1
  }

  boxplot(generation_interval ~ time_of_infection, data = pop_data,
          xlab = "Time of infection (day)",
          ylab = "Generation interval (days)",
          main = "Generation interval",
          ylim = c(0, max(10, pop_data$generation_interval)),
          xlim = c(0, params$num_days),
          xaxt = "n")
  axis(1, seq(0, params$num_days, 5))
}

#' Plot transmission matrix
#' 
plot_transmission_matrix2 <- function(pop_data,
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


#' Plot transmission matrix
#'
#' Creates a heatmap-style age-by-age transmission matrix based on infector
#' ages and infected individuals' ages observed in a simulation run.
#'
#' The matrix is constructed by binning ages into fixed-width age groups
#' (default: 3-year bins) and counting transmissions from each infector-age
#' bin to each contact-age bin.
#'
#' The plotting layout reserves space on the right for a custom color legend.
#' This implementation is inspired by the layout logic in `fields::imagePlot`,
#' but is implemented using base R only.
#'
#' @param pop_data Data frame produced by the IBM containing at least `age`: age 
#'   of each individual and `infector_age`: age of the infector for each infected 
#'   individual
#' @param breaks Numeric vector of breakpoints for color bins. If `NA` (default),
#'   pretty breaks are derived from the matrix counts.
#' @param main Character title of the plot.
#' @param num_digits Integer; number of digits to print inside cells.
#'   If `NA` (default), no numbers are printed.
#' @param age_bin_width Integer; width of age bins (years). Default is 3.
plot_transmission_matrix <- function(pop_data,
                                     breaks = NA,
                                     main = "Transmission matrix",
                                     num_digits = NA,
                                     age_bin_width = 3) {
  
  # Build transition matrix 
  max_age <- max(pop_data$age, na.rm = TRUE)
  
  # Age bin edges (left-closed, right-open via right = FALSE)
  age_breaks <- seq(1, max_age, by = age_bin_width)
  
  transmission_counts <- table(
    cut(pop_data$infector_age, age_breaks, right = FALSE),
    cut(pop_data$age,         age_breaks, right = FALSE),
    dnn = list("age infector", "age contact")
  )
  
  # If user did not supply breaks, derive breaks from data
  if (all(is.na(breaks))) {
    breaks <- pretty(transmission_counts)
  }
  
  # Add an extra "below minimum" break so the minimum is represented in legend
  breaks <- c(min(breaks) - 1, breaks)
  
  # Midpoints used for legend raster
  midpoints <- matrix(
    breaks[-length(breaks)] + diff(breaks) / 2,
    nrow = 1, ncol = length(breaks) - 1
  )
  
  # Colours (kept as base R heat palette to preserve previous output appearance)
  n_colours <- length(breaks) - 1
  colours <- heat.colors(n_colours)
  
  # Set layout for main matrix
  # NOTE: Layout logic based on fields::imagePlot style (implemented manually)
  char_size  <- par()$cin[1] / par()$din[1]     # text character size in plot coords
  offset     <- char_size * par()$mar[4]        # space between main plot and legend
  legend_w   <- 1
  legend_mar <- 8.1
  legend_shrink <- 0.9
  
  cex_lab  <- 1.2
  cex_axis <- 0.8
  cex_text <- 1
  
  legend_region <- par()$plt
  legend_region[2] <- 1 - (legend_mar * char_size)
  legend_region[1] <- legend_region[2] - (legend_w * char_size)
  
  # Shrink legend vertically
  pad <- (legend_region[4] - legend_region[3]) * ((1 - legend_shrink) / 2)
  legend_region[4] <- legend_region[4] - pad
  legend_region[3] <- legend_region[3] + pad
  
  main_region <- par()$plt
  main_region[2] <- min(main_region[2], legend_region[1] - offset)
  
  # Defensive adjustment so legend region remains valid
  region_width <- legend_region[2] - legend_region[1]
  legend_region[1] <- min(main_region[2] + offset, legend_region[1])
  legend_region[2] <- legend_region[1] + region_width
  
  # Preserve old graphical parameters and switch to main plot region
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  par(plt = main_region)
  
  # Start image
  image(
    transmission_counts,
    xlab = "Age infector",
    ylab = "Age contact",
    main = main,
    cex.lab = cex_lab,
    breaks = breaks,
    col = colours,
    xaxt = "n",
    yaxt = "n"
  )
  
  # Axis labels use bin labels from the table dimension names
  axis_ticks <- seq(0, 1, length.out = nrow(transmission_counts))
  axis(2, at = axis_ticks, labels = colnames(transmission_counts),
       cex.axis = cex_axis, tick = FALSE, las = 1, mgp = c(3, 0.3, 0))
  axis(1, at = axis_ticks, labels = colnames(transmission_counts),
       cex.axis = cex_axis, tick = FALSE, las = 2, mgp = c(3, 0.3, 0))
  
  # Optionally print numeric values in each cell
  if (!is.na(num_digits) && !is.na(cex_text) && cex_text > 0) {
    
    label_matrix <- transmission_counts
    if (any(max(transmission_counts, na.rm = TRUE) > 1)) {
      label_matrix <- round(label_matrix, digits = num_digits)
    } else {
      label_matrix <- format(label_matrix, digits = num_digits)
    }
    
    grid_centres <- expand.grid(axis_ticks, axis_ticks)
    text(grid_centres, labels = label_matrix, cex = cex_text)
  }
  
  # Add color-bar legend on right-hand side 
  par(new = TRUE, pty = "m", plt = legend_region, err = -1)
  
  image(
    x = 1:2,
    y = seq_len(length(breaks)),
    z = midpoints,
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    col = colours,
    breaks = breaks
  )
  
  mtext("count", side = 2, cex = 0.8)
  
  # Create nicer labels when breaks have gaps
  break_labels <- breaks
  if (any(diff(breaks) > 1)) {
    is_wide_bin <- c(FALSE, diff(breaks) > 1)
    ranged <- break_labels
    ranged[-1] <- paste(break_labels[-length(break_labels)], break_labels[-1], sep = " - ")
    break_labels[is_wide_bin] <- ranged[is_wide_bin]
  }
  
  axis(
    side = 4,
    at = seq_len(length(breaks) - 1) + 0.5,
    labels = break_labels[-1],
    cex.axis = cex_axis,
    mgp = c(3, 0.3, 0),
    las = 2,
    tick = FALSE
  )
}



