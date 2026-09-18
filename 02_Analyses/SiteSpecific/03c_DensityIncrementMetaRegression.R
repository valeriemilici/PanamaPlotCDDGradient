#Initialize Workspace ----------------------------------------------
rm(list =ls())
# Load libraries
library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(here)
library(metafor)
library(stringr)
library(patchwork)
library(tidyverse)

#Load Data
sites <- read.csv("2023data/SiteData_Gradient_Jeffs.csv")
sites[3,1] <- "ANAM"
#Load rAME output from Marginal Effects section
load(here("02_Analyses/(r)AMEOutput/mortality_increment_sitespecific.Rdata"))

# Separate rAME model output by density increment
rAME_1 <- rAME %>%
  filter(change == "one")%>%
  rename(BOSQUE = b.sp)
rAME_4 <- rAME %>%
  filter(change == "four")%>%
  rename(BOSQUE = b.sp)
rAME_8 <- rAME %>%
  filter(change == "eight")%>%
  rename(BOSQUE = b.sp)

# Join marginal effects and MAP data
rAME_1 <- left_join(rAME_1, sites, by = "BOSQUE")
rAME_4 <- left_join(rAME_4, sites, by = "BOSQUE")
rAME_8 <- left_join(rAME_8, sites, by = "BOSQUE")

# Prep data for metaregression --------
dat_meta_1 = metafor::escalc(measure = "GEN", # Set measure to generic, which passes the observed effect sizes or outcomes via the yi argument and the corresponding sampling variances via the vi argument (or the standard errors via the sei argument) to the function.
                             yi = estimate, # observed outcomes
                             sei = std.error, # standard errors
                             slab = BOSQUE, # label for site
                             data = rAME_1)
dat_meta_4 = metafor::escalc(measure = "GEN", # Set measure to generic, which passes the observed effect sizes or outcomes via the yi argument and the corresponding sampling variances via the vi argument (or the standard errors via the sei argument) to the function.
                             yi = estimate, # observed outcomes
                             sei = std.error, # standard errors
                             slab = BOSQUE, # label for site
                             data = rAME_4)
dat_meta_8 = metafor::escalc(measure = "GEN", # Set measure to generic, which passes the observed effect sizes or outcomes via the yi argument and the corresponding sampling variances via the vi argument (or the standard errors via the sei argument) to the function.
                             yi = estimate, # observed outcomes
                             sei = std.error, # standard errors
                             slab = BOSQUE, # label for site
                             data = rAME_8)

#save for plotting
save(dat_meta_1, file = "02_Analyses/SiteSpecific/output/dat_meta_1")
save(dat_meta_4, file = "02_Analyses/SiteSpecific/output/dat_meta_4")
save(dat_meta_8, file = "02_Analyses/SiteSpecific/output/dat_meta_8")
# Run metaregressions against MAP --------

metamod_1 = metafor::rma(yi = yi,
                         vi = vi,
                         mods = ~ scale(MAP),
                         method = "REML",
                         data = dat_meta_1)

metamod_4 = metafor::rma(yi = yi,
                         vi = vi,
                         mods = ~ scale(MAP),
                         method = "REML",
                         data = dat_meta_4)

metamod_8 = metafor::rma(yi = yi,
                         vi = vi,
                         mods = ~ scale(MAP),
                         method = "REML",
                         data = dat_meta_8)
## Examine model output
summary(metamod_1) #B = -0.076, se = 0.048, p = 0.11
summary(metamod_4) #B = 0.012, se = 0.055, p = 0.77
summary(metamod_8) #B =  0.056, se = 0.06, p = 0.35

#save output
save(metamod_1, file = "02_Analyses/SiteSpecific/output/metamod_1")
save(metamod_4, file = "02_Analyses/SiteSpecific/output/metamod_4")
save(metamod_8, file = "02_Analyses/SiteSpecific/output/metamod_8")
# Run metaregressions agaist Dry Deficit --------

metamod_dd_1 = metafor::rma(yi = yi,
                         vi = vi,
                         mods = ~ scale(DRY_DEFICIT),
                         method = "REML",
                         data = dat_meta_1)

metamod_dd_4 = metafor::rma(yi = yi,
                         vi = vi,
                         mods = ~ scale(DRY_DEFICIT),
                         method = "REML",
                         data = dat_meta_4)

metamod_dd_8 = metafor::rma(yi = yi,
                         vi = vi,
                         mods = ~ scale(DRY_DEFICIT),
                         method = "REML",
                         data = dat_meta_8)
## Examine model output
summary(metamod_dd_1) #B = -0.08, se = 0.04, p = 0.07
summary(metamod_dd_4) #B = 0.0, se = 0.05, p = 0.98
summary(metamod_dd_8) #B = 0.045, se = 0.06, p = 0.44

#save output
save(metamod_dd_1, file = "02_Analyses/SiteSpecific/output/metamod_dd_1")
save(metamod_dd_4, file = "02_Analyses/SiteSpecific/output/metamod_dd_4")
save(metamod_dd_8, file = "02_Analyses/SiteSpecific/output/metamod_dd_8")


