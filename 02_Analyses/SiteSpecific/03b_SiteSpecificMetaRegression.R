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

#Load Data
sites <- read.csv("2023data/SiteData_Gradient_Jeffs.csv")
sites[3,1] <- "ANAM"
#Load rAME output from Marginal Effects section
load(here("02_Analyses/(r)AMEOutput/mortality_sitespecific.Rdata"))

AME <- read.csv("02_Analyses/(r)AMEOutput/AME_sitespecific.csv")
rAME <- read.csv("02_Analyses/(r)AMEOutput/rAME_sitespecific.csv")
# Subset down to different measures of change
# rAME ---------
rAME_eq <- rAME %>%
  filter(change == "equilibrium") %>%
  rename(BOSQUE = b.sp)

rAME_in <- rAME %>%
  filter(change == "invasion")%>%
  rename(BOSQUE = b.sp)

# Join marginal effects and MAP data ------
rAME_eq <- left_join(rAME_eq, sites, by = "BOSQUE")
rAME_in <- left_join(rAME_in, sites, by = "BOSQUE")

# Reformat data for model fitting -------
# rAME ---------
dat_meta_eq = metafor::escalc(measure = "GEN", # Set measure to generic, which passes the observed effect sizes or outcomes via the yi argument and the corresponding sampling variances via the vi argument (or the standard errors via the sei argument) to the function.
                           yi = estimate, # observed outcomes
                           sei = std.error, # standard errors
                           slab = BOSQUE, # label for site
                           data = rAME_eq)

dat_meta_in = metafor::escalc(measure = "GEN", # Set measure to generic, which passes the observed effect sizes or outcomes via the yi argument and the corresponding sampling variances via the vi argument (or the standard errors via the sei argument) to the function.
                              yi = estimate, # observed outcomes
                              sei = std.error, # standard errors
                              slab = BOSQUE, # label for site
                              data = rAME_in)
save(dat_meta_eq, file = "02_Analyses/SiteSpecific/output/dat_meta_eq")
save(dat_meta_in, file = "02_Analyses/SiteSpecific/output/dat_meta_in")
# Fit model (MAP) ---------
# rAME ------------
#Equilibrium
metamod_eq = metafor::rma(yi = yi,
                       vi = vi,
                       mods = ~ scale(MAP),
                       method = "REML",
                       data = dat_meta_eq)

summary(metamod_eq)
# Teach yourself to interpret all of the metafor outputs so you can 
# judge whether these models are ok. At the momment this feels like magic.

save(metamod_eq, file = "02_Analyses/SiteSpecific/output/metamod_eq")

# Invasion

metamod_in = metafor::rma(yi = yi,
                          vi = vi,
                          mods = ~ scale(MAP),
                          method = "REML",
                          data = dat_meta_in)

summary(metamod_in)

save(metamod_in, file = "02_Analyses/SiteSpecific/output/metamod_in")

# Fit model (Dry Deficit) ---------
# rAME -----------
# Equilibrium
metamod2_eq = metafor::rma(yi = yi,
                       vi = vi,
                       mods = ~ scale(DRY_DEFICIT),
                       method = "REML",
                       data = dat_meta_eq)

summary(metamod2_eq)

save(metamod2_eq, file = "02_Analyses/SiteSpecific/output/metamod2_eq")

# Invasion
metamod2_in = metafor::rma(yi = yi,
                           vi = vi,
                           mods = ~ scale(DRY_DEFICIT),
                           method = "REML",
                           data = dat_meta_in)

summary(metamod2_in)

save(metamod2_in, file = "02_Analyses/SiteSpecific/output/metamod2_in")

##Move & Delete everything below this point to the appropriate place

# Bind predictions to dataframe ------
plotpreds <- cbind(pred, predict(object = metamod_7, newmods = pred$MAP))
plotpreds3 <- cbind(pred, predict(object = metamod_3, newmods = pred$MAP))
plotpreds0 <- cbind(pred, predict(object = metamod_0, newmods = pred$MAP))

plotpreds2 <- cbind(pred2, predict(object = metamod2_eq, newmods = pred2$DRY_DEFICIT))


plotpreds4 <- cbind(pred2, predict(object = metamod2_in, newmods = pred2$DRY_DEFICIT))

plotpreds5 <- cbind(pred, predict(object = metamod_eq_A, newmods = pred$MAP))
plotpreds6 <- cbind(pred2, predict(object = metamod2_eq_A, newmods = pred2$DRY_DEFICIT))

plotpreds7 <- cbind(pred, predict(object = metamod_in_A, newmods = pred$MAP))
plotpreds8 <- cbind(pred2, predict(object = metamod2_in_A, newmods = pred2$DRY_DEFICIT))
# Set abundance values for x axis
MAP_x_axis <- quantile(pred$MAP)
Def_x_axis <- quantile(pred2$DRY_DEFICIT)

# Plot prediction ----------
## Equilibrium rAME (MAP + DEF) ---------
MAP_Site_0 <- ggplot( ) + 
  geom_ribbon(data = plotpreds0, aes(x = MAP, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_0, aes( x = as.numeric(MAP), y = estimate,
                                          ymin = estimate - std.error,
                                          ymax = estimate + std.error)) +
  geom_line(data = plotpreds0, aes(x = MAP, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Mean Annual Precipitation", y = "rAME, ConDens = 0") + 
  scale_x_continuous(breaks = MAP_x_axis,
                     labels = MAP_x_axis) + 
  theme_bw(15)
MAP_Site_3 <- ggplot( ) + 
  geom_ribbon(data = plotpreds3, aes(x = MAP, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_3, aes( x = as.numeric(MAP), y = estimate,
                                          ymin = estimate - std.error,
                                          ymax = estimate + std.error)) +
  geom_line(data = plotpreds3, aes(x = MAP, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Mean Annual Precipitation", y = "rAME, ConDens = 3") + 
  scale_x_continuous(breaks = MAP_x_axis,
                     labels = MAP_x_axis) + 
  theme_bw(15) 


MAP_Site_7 <- ggplot( ) + 
  geom_ribbon(data = plotpreds, aes(x = MAP, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_7, aes( x = as.numeric(MAP), y = estimate,
                             ymin = estimate - std.error,
                             ymax = estimate + std.error)) +
  geom_line(data = plotpreds, aes(x = MAP, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Mean Annual Precipitation", y = "rAME, ConDens = 7") + 
  scale_x_continuous(breaks = MAP_x_axis,
                     labels = MAP_x_axis) + 
  theme_bw(15) 

MAP_Site_0 + MAP_Site_3 + MAP_Site_7

MAP_Site_3
MAP_Site_7

ggsave(plot = MAP_Site_R, filename = "figures/MAP_Site_eq_R.jpg")

Def_Site_R <- ggplot(plotpreds2, aes(x = DRY_DEFICIT, y = pred)) + 
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_eq, aes( x = as.numeric(DRY_DEFICIT), y = estimate,
                                           ymin = estimate - std.error,
                                           ymax = estimate + std.error)) +
  geom_line() + 
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Dry Season Deficit", y = "Equilibrium rAME") + 
  scale_x_continuous(breaks = Def_x_axis,
                     labels = Def_x_axis) + 
  theme_bw(15) 

Def_Site_R

ggsave(plot = Def_Site_R, filename = "figures/DEF_Site_eq_R.jpg")

## Invasion rAME (MAP + DEF) --------

MAP_Site_in_R <- ggplot( ) + 
  geom_ribbon(data = plotpreds3, aes(x = MAP, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_in, aes( x = as.numeric(MAP), y = estimate,
                                           ymin = estimate - std.error,
                                           ymax = estimate + std.error)) +
  geom_line(data = plotpreds3, aes(x = MAP, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Mean Annual Precipitation", y = "Invasion rAME") + 
  scale_x_continuous(breaks = MAP_x_axis,
                     labels = MAP_x_axis) + 
  theme_bw(15) 

MAP_Site_in_R

ggsave(plot = MAP_Site_in_R, filename = "figures/MAP_Site_in_R.jpg")

Def_Site_in_R <- ggplot(plotpreds4, aes(x = DRY_DEFICIT, y = pred)) + 
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_in, aes( x = as.numeric(DRY_DEFICIT), y = estimate,
                                           ymin = estimate - std.error,
                                           ymax = estimate + std.error)) +
  geom_line() + 
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Dry Season Deficit", y = "Invasion rAME") + 
  scale_x_continuous(breaks = Def_x_axis,
                     labels = Def_x_axis) + 
  theme_bw(15) 

Def_Site_in_R

ggsave(plot = Def_Site_in_R, filename = "figures/DEF_Site_in_R.jpg")

## Equilibrium AME (MAP + DEF) ------
MAP_Site_eq_A <- ggplot( ) + 
  geom_ribbon(data = plotpreds5, aes(x = MAP, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_eq_A, aes( x = as.numeric(MAP), y = estimate,
                                           ymin = estimate - std.error,
                                           ymax = estimate + std.error)) +
  geom_line(data = plotpreds5, aes(x = MAP, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Mean Annual Precipitation", y = "Equilibrium AME") + 
  scale_x_continuous(breaks = MAP_x_axis,
                     labels = MAP_x_axis) + 
  theme_bw(15) 

MAP_Site_eq_A

ggsave(plot = MAP_Site_eq_A, filename = "figures/MAP_Site_eq_A.jpg")

Def_Site_eq_A <- ggplot(plotpreds6, aes(x = DRY_DEFICIT, y = pred)) + 
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_eq_A, aes( x = as.numeric(DRY_DEFICIT), y = estimate,
                                           ymin = estimate - std.error,
                                           ymax = estimate + std.error)) +
  geom_line() + 
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Dry Season Deficit", y = "Equilibrium AME") + 
  scale_x_continuous(breaks = Def_x_axis,
                     labels = Def_x_axis) + 
  theme_bw(15) 

Def_Site_eq_A

ggsave(plot = Def_Site_eq_A, filename = "figures/DEF_Site_eq_A.jpg")
## Invasion AME (MAP + DEF) --------

MAP_Site_in_A <- ggplot( ) + 
  geom_ribbon(data = plotpreds7, aes(x = MAP, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_in_A, aes( x = as.numeric(MAP), y = estimate,
                                           ymin = estimate - std.error,
                                           ymax = estimate + std.error)) +
  geom_line(data = plotpreds7, aes(x = MAP, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Mean Annual Precipitation", y = "Invasion AME") + 
  scale_x_continuous(breaks = MAP_x_axis,
                     labels = MAP_x_axis) + 
  theme_bw(15) 

MAP_Site_in_A

ggsave(plot = MAP_Site_in_A, filename = "figures/MAP_Site_in_A.jpg")

Def_Site_in_A <- ggplot(plotpreds8, aes(x = DRY_DEFICIT, y = pred)) + 
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_in_A, aes( x = as.numeric(DRY_DEFICIT), y = estimate,
                                           ymin = estimate - std.error,
                                           ymax = estimate + std.error)) +
  geom_line() + 
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Dry Season Deficit", y = "Invasion AME") + 
  scale_x_continuous(breaks = Def_x_axis,
                     labels = Def_x_axis) + 
  theme_bw(15) 

Def_Site_in_A

ggsave(plot = Def_Site_in_A, filename = "figures/DEF_Site_in_A.jpg")

