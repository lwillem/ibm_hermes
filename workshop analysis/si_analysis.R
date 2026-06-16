# Load the stored contact tracing data
#---------------------------------------
contact_tracing_dat = readRDS(file = "contact_tracing_data.rds")
head(contact_tracing_dat)

# Load the demographic data
#---------------------------
demographic_dat = readRDS(file = "demographic_data.rds")
head(demographic_dat)

# Combine the different data sources
#------------------------------------
combined_dat = merge(contact_tracing_dat, demographic_dat, by = "ind_id")
combined_dat$s = combined_dat$ind_symptom_onset - combined_dat$reported_contact_symptom_onset
combined_dat$sl = combined_dat$ind_symptom_onset - (combined_dat$reported_contact_symptom_onset + 1)
combined_dat$su = (combined_dat$ind_symptom_onset + 1) - combined_dat$reported_contact_symptom_onset

# Look at the directionality of contacts/exposures
#--------------------------------------------------
combined_dat$infectee <- combined_dat$ind_id
combined_dat$infector <- combined_dat$reported_contact

combined_dat$time_of_contact_upper < combined_dat$ind_symptom_onset
combined_dat$time_of_contact_upper < combined_dat$reported_contact_symptom_onset

# Installation of relevant packages for analyses
#------------------------------------------------
#library("TranStat")
library("EpiLPS")
library("mitey")
library("EpiDelays")

# Some exploration of naive estimation of the serial interval distribution
#--------------------------------------------------------------------------
si_data <- data.frame(xl = pmax(0,combined_dat$sl), xr = combined_dat$su)
SIfit <- nonparfit2(x = na.omit(si_data),  Bboot = 100)

par_fit <- parfitml(x = na.omit(si_data), family = "gamma", ci = "npboot")
summary(par_fit)

xgrid <- seq(0.01, 10, 0.01)
plot_df <- data.frame(x = xgrid,
                      y = dgamma(x = xgrid, 
                                 shape = par_fit$parfit$par1$point,
                                 rate = par_fit$parfit$par2$point))

# Nonparametric and parametric estimation of the serial interval distribution
#-----------------------------------------------------------------------------
library(ggplot2)
ggplot(plot_df, aes(x = x, y = y)) +
  geom_line()

## Helper function
##-----------------
nonparfit2 <- function (x, Bboot = 1000, pgbar = TRUE) 
{
  tic <- proc.time()
  n <- nrow(x)
  nc <- ncol(x)
  if (nc != 2) {
    stop("Data frame must have 2 columns")
  }
  dfck <- kerdata_check(x = x)
  if (dfck$result == "fail") {
    stop(dfck$message)
  }
  xl <- x[, 1]
  xr <- x[, 2]
  ninv <- 1/n
  pfeats <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
  pointestim <- function(tl, tr) {
    tmid <- 0.5 * (tl + tr)
    tw <- tr - tl
    Fhat <- function(t) ninv * sum((t - tl)/tw * (t >= tl & 
                                                    t <= tr) + (t > tr))
    tord <- sort(c(tl, tr))
    Fhattord <- sapply(tord, Fhat)
    qfun <- function(p) {
      pcub <- which(p <= Fhattord)[1]
      t1 <- tord[pcub - 1]
      t2 <- tord[pcub]
      Fhatt1 <- Fhat(t1)
      Fhatt2 <- Fhat(t2)
      dFinv <- 1/(Fhatt2 - Fhatt1)
      val <- (t1 * (Fhatt2 - p) + t2 * (p - Fhatt1)) * 
        dFinv
      return(val)
    }
    mu <- mean(tmid)
    sd <- sqrt(mean((tl^2 + tl * tr + tr^2)/3) - mu^2)
    qp <- sapply(pfeats, qfun)
    o <- list(mu = mu, sd = sd, qp = qp)
    return(o)
  }
  npp <- pointestim(xl, xr)
  npfeat <- c(npp$mu, npp$sd^2, npp$sd, npp$qp)
  if (isTRUE(pgbar)) {
    cat(paste0("Nonparametric fit \n", "Bootstrap progress (Bboot=", 
               Bboot, "): \n"))
    progbar <- utils::txtProgressBar(min = 1, max = Bboot, 
                                     initial = 1, style = 3, char = "*")
  }
  fboot <- matrix(0, nrow = Bboot, ncol = length(npfeat))
  for (b in 1:Bboot) {
    xb <- kerboot(x)
    bootest <- pointestim(tl = xb[, 1], tr = xb[, 2])
    fboot[b, ] <- c(bootest$mu, bootest$sd^2, bootest$sd, 
                    bootest$qp)
    if (isTRUE(pgbar)) 
      utils::setTxtProgressBar(progbar, b)
  }
  if (isTRUE(pgbar)) 
    close(progbar)
  feats <- c("mean", "var", "sd", paste0("q", c(1, 5, 25, 50, 
                                                75, 95, 99)))
  delayfit <- stats::setNames(vector("list", length(feats)), 
                              feats)
  delayfit <- kerstats(slist = delayfit, pestim = npfeat, boot = fboot,
                       se = NULL)
  toc <- proc.time() - tic
  o <- list(n = n, Bboot = Bboot, delayfit = delayfit, censtype = "single", 
            elapsed = toc[3])
  attr(o, "class") <- "nonparfit"
  return(o)
}

