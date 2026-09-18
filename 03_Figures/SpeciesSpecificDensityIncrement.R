#Initialize Workspace ----------------------------------------------
rm(list =ls())

library(tidyverse)
library(ggplot2)
library(patchwork)
library(metafor)
library(gridExtra)

# Load Models

load(file = "02_Analyses/SpeciesSpecific/output/metamod_1")
load(file = "02_Analyses/SpeciesSpecific/output/metamod_4")
load(file = "02_Analyses/SpeciesSpecific/output/metamod_8")

#load data
load(file = "02_Analyses/SpeciesSpecific/output/dat_meta_1")
load(file = "02_Analyses/SpeciesSpecific/output/dat_meta_4")
load(file = "02_Analyses/SpeciesSpecific/output/dat_meta_8")

sites <- read.csv("2023data/SiteData_Gradient_Jeffs.csv")

# Generate a prediction dataframe ---------
MAPpred <- data_frame(MAP = seq(min(sites$MAP, na.rm = TRUE), 
                                max(sites$MAP, na.rm = TRUE),
                                length.out = 50),
                      s_MAP = seq(min(scale(sites$MAP)),
                                  max(scale(sites$MAP)),
                                  length.out = 50))

# Set abundance values for x axis
MAP_x_axis <- quantile(MAPpred$MAP)

### Combined plot of MAP relationship
#Prediction Plot
pred1 <- cbind(MAPpred, predict(object = metamod_1, newmods = MAPpred$s_MAP))
pred4 <- cbind(MAPpred, predict(object = metamod_4, newmods = MAPpred$s_MAP))
pred8 <- cbind(MAPpred, predict(object = metamod_8, newmods = MAPpred$s_MAP))

MAP1 <- ggplot( ) + 
  geom_ribbon(data = pred1, aes(x = MAP, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_1, aes( x = as.numeric(MAP), y = estimate,
                                          ymin = estimate - std.error,
                                          ymax = estimate + std.error),
                  position = position_dodge2(width = 40)) +
  geom_line(data = pred1, aes(x = MAP, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  ylim(c(-0.8,1.1)) +
  labs(x = "Mean Annual Precipitation",
       y = "Conspecific Seedling Density Effect \n on Mortality (rAME)") + 
  scale_x_continuous(breaks = MAP_x_axis,
                     labels = MAP_x_axis) + 
  ggtitle("Density: 0 - 1")+
  theme_classic(12) +
  theme(plot.title=element_text( vjust = -5, hjust = 0.05))



MAP4 <- ggplot( ) + 
  geom_ribbon(data = pred4, aes(x = MAP, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_4, aes( x = as.numeric(MAP), y = estimate,
                                          ymin = estimate - std.error,
                                          ymax = estimate + std.error),
                  position = position_dodge2(width = 40)) +
  geom_line(data = pred4, aes(x = MAP, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  ylim(c(-0.8, 1.1)) +
  labs(x = "Mean Annual Precipitation",
       y = "Conspecific Seedling Density Effect \n on Mortality (rAME)") + 
  scale_x_continuous(breaks = MAP_x_axis,
                     labels = MAP_x_axis) + 
  ggtitle("3 - 4")+
  theme_classic(12) +
  theme(plot.title=element_text( vjust = -5, hjust = 0.05))

MAP8 <- ggplot( ) + 
  geom_ribbon(data = pred8, aes(x = MAP, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_8, aes( x = as.numeric(MAP), y = estimate,
                                          ymin = estimate - std.error,
                                          ymax = estimate + std.error),
                  position = position_dodge2(width = 40)) +
  geom_line(data = pred8, aes(x = MAP, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
 ylim(c(-0.8, 1.1)) +
  labs(x = "Mean Annual Precipitation",
       y = "Conspecific Seedling Density Effect \n on Mortality (rAME)") + 
  scale_x_continuous(breaks = MAP_x_axis,
                     labels = MAP_x_axis) + 
  ggtitle("7 - 8")+
  theme_classic(12) +
  theme(plot.title=element_text( vjust = -5, hjust = 0.05))


MAP_DENS <- MAP1 + MAP4 + MAP8 +
  plot_annotation(tag_levels = "A") +
  plot_layout(axis_titles = "collect")


MAP_DENS

ggsave(MAP_DENS, filename = "03_Figures/MS_Figs/Dens_Increment_SS_MAP.png",
       width = 10, height = 5, units = "in")
