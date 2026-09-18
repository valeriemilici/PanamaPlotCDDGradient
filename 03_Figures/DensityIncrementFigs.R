#Initialize Workspace ----------------------------------------------
rm(list =ls())

library(tidyverse)
library(ggplot2)
library(patchwork)
library(metafor)
library(gridExtra)

# Load Models

load(file = "02_Analyses/SiteSpecific/output/metamod_1")
load(file = "02_Analyses/SiteSpecific/output/metamod_4")
load(file = "02_Analyses/SiteSpecific/output/metamod_8")

load(file = "02_Analyses/SiteSpecific/output/metamod_dd_1")
load(file = "02_Analyses/SiteSpecific/output/metamod_dd_4")
load(file = "02_Analyses/SiteSpecific/output/metamod_dd_8")

#load data
load(file = "02_Analyses/SiteSpecific/output/dat_meta_1")
load(file = "02_Analyses/SiteSpecific/output/dat_meta_4")
load(file = "02_Analyses/SiteSpecific/output/dat_meta_8")

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
                                           ymax = estimate + std.error)) +
  geom_line(data = pred1, aes(x = MAP, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  ylim(c(-0.5, 0.75)) +
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
                                          ymax = estimate + std.error)) +
  geom_line(data = pred4, aes(x = MAP, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  ylim(c(-0.5, 0.75)) +
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
                                          ymax = estimate + std.error)) +
  geom_line(data = pred8, aes(x = MAP, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  ylim(c(-0.5, 0.75)) +
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

ggsave(MAP_DENS, filename = "03_Figures/MS_Figs/Dens_Increment_MAP.png",
       width = 10, height = 4, units = "in")

### Combined plot of DEF relationship
#Prediction Plot
pred1dd <- cbind(DEFpred, predict(object = metamod_dd_1, newmods = DEFpred$s_DEF))
pred4dd <- cbind(DEFpred, predict(object = metamod_dd_4, newmods = DEFpred$s_DEF))
pred8dd <- cbind(DEFpred, predict(object = metamod_dd_8, newmods = DEFpred$s_DEF))

DEF1 <- ggplot( ) + 
  geom_ribbon(data = pred1dd, aes(x = DRY_DEFICIT, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_1, aes( x = as.numeric(DRY_DEFICIT), y = estimate,
                                          ymin = estimate - std.error,
                                          ymax = estimate + std.error)) +
  geom_line(data = pred1dd, aes(x = DRY_DEFICIT, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
 ylim(c(-0.4, 0.75)) +
  labs(x = "Dry Season Deficit", y = "rAME") + 
  scale_x_continuous(breaks = Def_x_axis,
                     labels = Def_x_axis) + 
  theme_bw(12)



DEF4 <- ggplot( ) + 
  geom_ribbon(data = pred4dd, aes(x = DRY_DEFICIT, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_4, aes( x = as.numeric(DRY_DEFICIT), y = estimate,
                                          ymin = estimate - std.error,
                                          ymax = estimate + std.error)) +
  geom_line(data = pred4dd, aes(x = DRY_DEFICIT, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  ylim(c(-0.4, 0.75)) +
  labs(x = "Dry Season Deficit", y = "rAME") + 
  scale_x_continuous(breaks = Def_x_axis,
                     labels = Def_x_axis) + 
  theme_bw(12)

DEF8 <- ggplot( ) + 
  geom_ribbon(data = pred8dd, aes(x = DRY_DEFICIT, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_8, aes( x = as.numeric(DRY_DEFICIT), y = estimate,
                                          ymin = estimate - std.error,
                                          ymax = estimate + std.error)) +
  geom_line(data = pred8dd, aes(x = DRY_DEFICIT, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
   ylim(c(-0.4, 0.75)) +
  labs(x = "Dry Season Deficit", y = "rAME") + 
  scale_x_continuous(breaks = Def_x_axis,
                     labels = Def_x_axis) + 
  theme_bw(12)

DEF_DENS <- DEF1 + DEF4 + DEF8 +
  plot_annotation(tag_levels = "A") +
  plot_layout(axis_titles = "collect")


DEF_DENS

ggsave(DEF_DENS, filename = "03_Figures/MS_Figs/Dens_Increment_DEF.png",
       width = 10, height = 4, units = "in")
## If I decide to create another plot of slope change across density see below:
# Create DF of slopes and error across models ------

slope <- function(mod){
  
  out <- data.frame(est = mod$b[2,1],
                    lwr = mod$ci.lb[2],
                    upr = mod$ci.ub[2])
  
  return(out)
}

#MAP models
modlist <- list(metamod_0, metamod_1, metamod_2, metamod_3, metamod_4,
                metamod_5, metamod_6, metamod_7)
predout <- data.frame(apply(sapply(modlist, slope), 1, unlist)) %>%
  mutate(condens = paste(c(0:7),rep("-", 8),c(1:8)),
         dens = c(1:8)) 
#DD models
modlist_dd <- list(metamod_dd_0, metamod_dd_1, metamod_dd_2, metamod_dd_3,
                   metamod_dd_4, metamod_dd_5, metamod_dd_6, metamod_dd_7)

predout_dd <- data.frame(apply(sapply(modlist_dd, slope), 1, unlist)) %>%
  mutate(condens = paste(c(0:7),rep("-", 8),c(1:8)),
         dens = c(1:8))


# Plot the relationships ----

ggplot(predout, mapping = aes(x = dens, y = est)) +
  geom_ribbon(mapping = aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  geom_line() +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
  theme_classic()

ggplot(predout_dd, mapping = aes(x = dens, y = est)) +
  geom_ribbon(mapping = aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  geom_line() +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
  theme_classic()