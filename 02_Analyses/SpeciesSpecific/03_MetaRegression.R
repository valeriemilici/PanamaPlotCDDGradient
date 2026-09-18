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

#Load Data
sites <- read.csv("2023data/SiteData_Gradient_Jeffs.csv")
sites[3,1] <- "ANAM"
sites <- sites %>% mutate(MAP_s = (MAP - mean(MAP))/sd(MAP),
                          DD_s = (DRY_DEFICIT - mean(DRY_DEFICIT))/sd(DRY_DEFICIT))

#Load rAME output from Marginal Effects section
rAME_nl <- read.csv("02_Analyses/(r)AMEOutput/rAMEmortnl.csv") 

## Data Manipulation --------------------


# Subset down to just equilibrium change and invasion change
rAME.eq <- rAME_nl %>%
  filter(change == "equilibrium") %>%
 filter(estimate <= quantile(estimate, 0.95)) %>% filter(estimate >= quantile(estimate, 0.025)) %>%
filter(std.error <= quantile(std.error, 0.92)) %>%
 filter(std.error >= quantile(std.error, 0.025)) %>% 
  mutate(BOSQUE = str_split_i(b.sp, "_", 1)) %>%
  mutate(spp = str_split_i(b.sp, "_", -1))
# note: filtering the UL of the estimate and error to be 0.95 gets a model that 
# maximizes the data it uses while removing unusually large estimate (i.e., 22). 
# Changing the value from 0.975 to 0.95 does NOT affect the estimates of
# the relationships (or their significance) in the metaregression.

# I could go further. Santa Rita_RINOSQ has very large std.error

rAME.in <- rAME_nl %>%
  filter(change == "invasion") %>%
  filter(estimate <= quantile(estimate, 0.95)) %>%
  filter(estimate >= quantile(estimate, 0.025)) %>%
  filter(std.error <= quantile(std.error, 0.92)) %>%
  filter(std.error >= quantile(std.error, 0.025)) %>%
  mutate(BOSQUE = str_split_i(b.sp, "_", 1))  

# Join marginal effects and MAP data
rAME.eq <- left_join(rAME.eq, sites, by = "BOSQUE")
rAME.in <- left_join(rAME.in, sites, by = "BOSQUE")
#head(rAME.in, n = 10)

# Reformat data for model fitting ------
dat_meta_eq = metafor::escalc(measure = "GEN", # Set measure to generic, which passes the observed effect sizes or outcomes via the yi argument and the corresponding sampling variances via the vi argument (or the standard errors via the sei argument) to the function.
                           yi = estimate, # observed outcomes
                           sei = std.error, # standard errors
                           slab = b.sp, # label for species
                           data = rAME.eq)

dat_meta_in = metafor::escalc(measure = "GEN", # Set measure to generic, which passes the observed effect sizes or outcomes via the yi argument and the corresponding sampling variances via the vi argument (or the standard errors via the sei argument) to the function.
                              yi = estimate, # observed outcomes
                              sei = std.error, # standard errors
                              slab = b.sp, # label for species
                              data = rAME.in)
# Fit model (MAP) -------
metamod.eq = metafor::rma.mv(yi = yi,
                       V = vi,
                       mods = ~ MAP_s,
                       random = c(~1|BOSQUE, ~1|spp), #add random intercept for site
                       method = "REML",
                       data = dat_meta_eq)
# add random intercept for site

metamod.in = metafor::rma(yi = yi,
                          vi = vi,
                          mods = ~ MAP_s,
                          method = "REML",
                          data = dat_meta_in)
summary(metamod.eq)
summary(metamod.in)

forest(metamod.eq, 
       header = "Species x Site", 
       xlab = "Stabilizing CDD",
       order = MAP)

baujat(metamod.eq) 
#89 and 99 are very influential; Sherman_BROSUT & Sherman_TOVOLO
#Sherman_BROSUT has a reliable estimate of CNDD
#Sherman_TOVOLO has a reliable estimate of CPDD
#their influence opposes one another. Removing them from the MR slightly changes
# the model output, but does not alter the interpretation.

### Plot the models ---------
# Generate a prediction dataframe
pred <- expand_grid(MAP_s = seq(min(dat_meta_eq$MAP_s, na.rm = TRUE), 
                              max(dat_meta_eq$MAP_s, na.rm = TRUE),
                              length.out = 50) )
# Bind predictions to dataframe
pred_eq <- cbind(pred, predict(object = metamod.eq, newmods = pred$MAP_s))
#pred_in <- cbind(pred, predict(object = metamod.in, newmods = pred$MAP_s))

# Add unscaled MAP column for plotting
#scale(sites$MAP) #mean = 2358.75, sd = 543.74
pred_eq <- mutate(pred_eq, MAP = (MAP_s * 543.74) + 2358.75)
#pred_in <- mutate(pred_in, MAP = (MAP_s * 543.74) + 2358.75)

# Set abundance values for x axis
MAP_x_axis <- quantile(pred_eq$MAP)

# Plot prediction
## Equilibrium MAP Figure -------------------

#add label to different point types
rAME.eq <- rAME.eq %>% mutate(pointtype = if_else(spp == "seedling", "triangle", "point"))

MAP_eq <- ggplot() + 
  geom_ribbon(pred_eq, mapping = aes(x = MAP, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "goldenrod1") + 
  geom_line(pred_eq, mapping = aes(x = MAP, y = pred)) + 
  geom_jitter(rAME.eq, mapping = aes(x = MAP, y = estimate,
                                     size = 1/std.error,
                                     shape = point, color = pointtype)) +
  scale_alpha_manual(values = c(point = 0.5, triangle = 1)) +
  scale_color_manual(values = c(point = "black", triangle = "darkgrey")) +
  
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Mean Annual Precipitation (mm)", y = "Stabilizing CDD") + 
  scale_x_continuous(breaks = as.numeric(format(round(MAP_x_axis, 0))),
                     labels = as.numeric(format(round(MAP_x_axis, 0)))) + 
  ylim(-1,1) +
  theme_classic(16) +
  theme(legend.position = "none")

MAP_eq

ggsave(plot = MAP_eq, filename = "03_Figures/MS_Figs/MAP_eq.jpg")

## Invasion MAP Figure -----------------------------
MAP_in <- ggplot() + 
  geom_ribbon(pred_in, mapping = aes(x = MAP, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "steelblue2", alpha = 0.75) + 
  geom_line(pred_in, mapping = aes(x = MAP, y = pred)) + 
  geom_jitter(rAME.in, mapping = aes(x = MAP, y = estimate), alpha = 0.4) +
  #geom_pointrange(rAME.ss.in, mapping = aes( x = MAP, y = estimate,
   #                                          ymin = estimate - std.error,
    #                                         ymax = estimate + std.error),
     #             shape = 15, size = 0.8) +
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Mean Annual Precipitation (mm)", y = "Invasion rAME") + 
  scale_x_continuous(breaks = as.numeric(format(round(MAP_x_axis, 0))),
                     labels = as.numeric(format(round(MAP_x_axis, 0)))) + 
  theme_bw(15)  

MAP_in

ggsave(plot = MAP_in, filename = "figures/MAP_in.jpg")

## Equilibrium DD Figure  ------------------
DD_eq <- ggplot() + 
  geom_ribbon(pred_eq_dd, mapping = aes( x = DD, y = pred, ymin = ci.lb, ymax = ci.ub), fill = "steelblue2", alpha = 0.75) + 
  geom_line(pred_eq_dd, mapping = aes(x = DD, y = pred)) + 
  geom_jitter(rAME.eq, mapping = aes(x = DRY_DEFICIT, y = estimate), alpha = 0.4) +
 # geom_pointrange(rAME.ss.eq, mapping = aes( x = DRY_DEFICIT, y = estimate,
  #                                           ymin = estimate - std.error,
   #                                          ymax = estimate + std.error),
    #              shape = 15, size = 0.8) +
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Dry Season Deficit (mm)", y = " Equilibrium rAME") + 
  scale_x_continuous(breaks = as.numeric(format(round(DD_x_axis, 0))),
                     labels = as.numeric(format(round(DD_x_axis, 0)))) + 
  theme_bw(15) 

DD_eq

ggsave(plot = DD_eq, filename = "figures/DD_eq.jpg")

## Invasion DD Figure  ------------------
DD_in <- ggplot() + 
  geom_ribbon(pred_in_dd, mapping = aes( x = DD, y = pred, ymin = ci.lb, ymax = ci.ub), fill = "steelblue2", alpha = 0.75) + 
  geom_line(pred_in_dd, mapping = aes(x = DD, y = pred)) + 
  geom_jitter(rAME.in, mapping = aes(x = DRY_DEFICIT, y = estimate), alpha = 0.4) +
  geom_pointrange(rAME.ss.in, mapping = aes( x = DRY_DEFICIT, y = estimate,
                                             ymin = estimate - std.error,
                                             ymax = estimate + std.error),
                  shape = 15, size = 0.8) +
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Dry Season Deficit (mm)", y = " Equilibrium rAME") + 
  scale_x_continuous(breaks = as.numeric(format(round(DD_x_axis, 0))),
                     labels = as.numeric(format(round(DD_x_axis, 0)))) + 
  theme_bw(15) 

DD_in

ggsave(plot = DD_in, filename = "figures/DD_in.jpg")

# Average Marginal Effects --------------------------------------------------

## Data Manipulation --------------------
# Subset down to just equilibrium change and invasion change
AME.eq <- AME %>%
  filter(change == "equilibrium") %>%
  filter(std.error <= 5) %>%
  mutate(BOSQUE = str_split_i(b.sp, "_", -1)) 

AME.eq <- left_join(AME.eq, sites, by = "BOSQUE")

# Reformat data for model fitting ------
dat_meta_eq_AME = metafor::escalc(measure = "GEN", # Set measure to generic, which passes the observed effect sizes or outcomes via the yi argument and the corresponding sampling variances via the vi argument (or the standard errors via the sei argument) to the function.
                              yi = estimate, # observed outcomes
                              sei = std.error, # standard errors
                              slab = b.sp, # label for species
                              data = AME.eq)
# Fit model (MAP) -------
metamod.AME.eq = metafor::rma(yi = yi,
                          vi = vi,
                          mods = ~ MAP_s,
                          method = "REML",
                          data = dat_meta_eq_AME)

summary(metamod.AME.eq)
# seems to be just as mediocre as rAME

## Equilibrium MAP Figure -------------------
# Generate a prediction dataframe
pred_AME <- expand_grid(MAP_s = seq(min(dat_meta_eq_AME$MAP_s, na.rm = TRUE), 
                                max(dat_meta_eq_AME$MAP_s, na.rm = TRUE),
                                length.out = 50))

pred_eq_AME <- cbind(pred_AME, predict(object = metamod.AME.eq, newmods = pred_AME$MAP_s))
pred_eq_AME <- mutate(pred_eq_AME, MAP = (MAP_s * 543.74) + 2358.75)
MAP_x_axis <- quantile(pred_eq_AME$MAP)


MAP_eq_AME <- ggplot() + 
  geom_ribbon(pred_eq_AME, mapping = aes(x = MAP, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "steelblue2", alpha = 0.75) + 
  geom_line(pred_eq_AME, mapping = aes(x = MAP, y = pred)) + 
  geom_jitter(AME.eq, mapping = aes(x = MAP, y = estimate), alpha = 0.4) +
  #geom_pointrange(rAME.ss.eq, mapping = aes( x = MAP, y = estimate,
  #                                          ymin = estimate - std.error,
  #                                         ymax = estimate + std.error),
  #             shape = 15, size = 0.8) +
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Mean Annual Precipitation (mm)", y = "Equilibrium AME") + 
  scale_x_continuous(breaks = as.numeric(format(round(MAP_x_axis, 0))),
                     labels = as.numeric(format(round(MAP_x_axis, 0)))) + 
  theme_bw(15) 

MAP_eq_AME

