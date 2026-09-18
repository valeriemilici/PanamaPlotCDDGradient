#Initialize Workspace ----------------------------------------------
rm(list =ls())

library(tidyverse)
library(ggplot2)
library(patchwork)
library(metafor)
library(gridExtra)

# Load Models
load("02_Analyses/SiteSpecific/output/metamod_eq")
load("02_Analyses/SiteSpecific/output/metamod_in")
load("02_Analyses/SiteSpecific/output/metamod2_eq")
load("02_Analyses/SiteSpecific/output/metamod2_in")

#load data
load("02_Analyses/SiteSpecific/output/dat_meta_eq")
load("02_Analyses/SiteSpecific/output/dat_meta_in")

sites <- read.csv("2023data/SiteData_Gradient_Jeffs.csv")

# Generate a prediction dataframe ---------
MAPpred <- data_frame(MAP = seq(min(sites$MAP, na.rm = TRUE), 
                              max(sites$MAP, na.rm = TRUE),
                              length.out = 50),
                       s_MAP = seq(min(scale(sites$MAP)),
                                   max(scale(sites$MAP)),
                                   length.out = 50))

DEFpred <- data_frame(DRY_DEFICIT = seq(min(sites$DRY_DEFICIT, na.rm = T),
                                       max(sites$DRY_DEFICIT, na.rm = T),
                                       length.out = 50),
                      s_DEF = seq(min(scale(sites$DRY_DEFICIT)),
                                  max(scale(sites$DRY_DEFICIT)),
                                  length.out = 50))
# Set abundance values for x axis
MAP_x_axis <- quantile(MAPpred$MAP)
Def_x_axis <- quantile(DEFpred$DRY_DEFICIT)


##MAP ----
#Forest Plot
png(file = "03_Figures/MS_Figs/EQ_MAP_forest.png",
    width =1800, height = 2400, res = 300)
forest(metamod_eq, 
       header = "Site", 
       xlab = "rAME",
       order = MAP)
dev.off()

#Prediction Plot
pred1 <- cbind(MAPpred, predict(object = metamod_eq, newmods = MAPpred$s_MAP))

MAP_Site_EQ <- ggplot( ) + 
  geom_ribbon(data = pred1, aes(x = MAP, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "goldenrod1", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_eq, aes( x = as.numeric(MAP), y = estimate,
                                          ymin = estimate - std.error,
                                          ymax = estimate + std.error)) +
  geom_line(data = pred1, aes(x = MAP, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Mean Annual Precipitation", y = "rAME") + 
  scale_x_continuous(breaks = MAP_x_axis,
                     labels = MAP_x_axis) + 
  theme_classic(12)

ggsave(MAP_Site_EQ, filename = "03_Figures/MS_Figs/MAP_Site_Eq.png",
       width = 2600, height = 1800, units = "px")

#Forest Plot
png(file = "03_Figures/MS_Figs/IN_MAP_forest.png",
    width =720, height = 2400, res = 1080)

Forest_InMAP <- forest(metamod_in, 
       header = "Site", 
       xlab = "rAME",
       order = MAP)
dev.off()

#Prediction Plot
pred2 <- cbind(MAPpred, predict(object = metamod_in, newmods = MAPpred$s_MAP))

MAP_Site_IN <- ggplot( ) + 
  geom_ribbon(data = pred2, aes(x = MAP, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "goldenrod1", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_in, aes( x = as.numeric(MAP), y = estimate,
                                           ymin = estimate - std.error,
                                           ymax = estimate + std.error)) +
  geom_line(data = pred2, aes(x = MAP, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Mean Annual Precipitation", y = "rAME") + 
  scale_x_continuous(breaks = MAP_x_axis,
                     labels = MAP_x_axis) + 
  theme_classic(12)

ggsave(MAP_Site_IN, filename = "03_Figures/MS_Figs/MAP_Site_IN.png",
       width = 2600, height = 1800, units = "px")

## Dry Deficit ----
#Forest Plot
png(file = "03_Figures/MS_Figs/Eq_DEF_forest.png",
    width =1800, height = 2400, res = 300)
forest(metamod2_eq, 
       header = "Site", 
       xlab = "rAME",
       order = DRY_DEFICIT)
dev.off()

#Prediction Plot
pred3 <- cbind(DEFpred, predict(object = metamod2_eq, newmods = DEFpred$s_DEF))

DEF_Site_EQ <- ggplot( ) + 
  geom_ribbon(data = pred3, aes(x = DRY_DEFICIT, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_eq, aes( x = as.numeric(DRY_DEFICIT), y = estimate,
                                           ymin = estimate - std.error,
                                           ymax = estimate + std.error)) +
  geom_line(data = pred3, aes(x = DRY_DEFICIT, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Dry Season Deficit", y = "rAME") + 
  scale_x_continuous(breaks = Def_x_axis,
                     labels = Def_x_axis) + 
  theme_bw(12)

ggsave(DEF_Site_EQ, filename = "03_Figures/MS_Figs/DEF_Site_EQ.png",
       width = 2600, height = 1800, units = "px")

#Forest Plot
png(file = "03_Figures/MS_Figs/IN_DEF_forest.png",
    width =1800, height = 2400, res = 300)
Forest_InDef <- forest(metamod2_in, 
       header = "Site", 
       xlab = "rAME",
       order = DRY_DEFICIT) #the relationship is marginally negative

dev.off()

#Prediction Plot
pred4 <- cbind(DEFpred, predict(object = metamod2_in, newmods = DEFpred$s_DEF))

DEF_Site_IN <- ggplot( ) + 
  geom_ribbon(data = pred4, aes(x = DRY_DEFICIT, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_in, aes( x = as.numeric(DRY_DEFICIT), y = estimate,
                                           ymin = estimate - std.error,
                                           ymax = estimate + std.error)) +
  geom_line(data = pred4, aes(x = DRY_DEFICIT, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Dry Season Deficit", y = "rAME") + 
  scale_x_continuous(breaks = Def_x_axis,
                     labels = Def_x_axis) + 
  theme_bw(12)

ggsave(DEF_Site_IN, filename = "03_Figures/MS_Figs/DEF_Site_IN.png",
       width = 2600, height = 1800, units = "px")

