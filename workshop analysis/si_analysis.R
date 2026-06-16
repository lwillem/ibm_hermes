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

# Installation of relevant packages for analyses
#------------------------------------------------
#library("TranStat")
library("EpiLPS")
library("mitey")
library("EpiDelays")

# Some exploration of naive estimation of the serial interval distribution
#--------------------------------------------------------------------------
si_data <- data.frame(xl = combined_dat$sl, xr = combined_dat$su)
SIfit <- nonparfit(x = na.omit(si_data),  Bboot = 100)

parfitml(x = na.omit(si_data), family = "gamma", ci = "npboot")