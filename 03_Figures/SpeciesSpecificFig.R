#Initialize Workspace ----------------------------------------------
rm(list =ls())

library(tidyverse)
library(ggplot2)
library(patchwork)
library(metafor)
library(gridExtra)

# Load Models

load("02_Analyses/SpeciesSpecific/metamod.spp.in")

#load data
load("02_Analyses/SpeciesSpecific/dat_meta_in")

sites <- read.csv("2023data/SiteData_Gradient_Jeffs.csv")
dat <- read.csv("2023data/2023censusdat.csv")

#data for forest plot species column
dat1 <- dat %>% dplyr::select(SPP, GENUS_SPECIES) %>%
  group_by(SPP) %>%
  slice(1) %>%
  ungroup()
dat_meta_in <- dat_meta_in %>% mutate(SPP = str_split_i(b.sp, "_", -1)) 
dat_meta_in <- left_join(dat_meta_in, dat1, by = "SPP")

# Generate a prediction dataframe ---------
MAPpred <- data_frame(MAP = seq(min(sites$MAP, na.rm = TRUE), 
                                max(sites$MAP, na.rm = TRUE),
                                length.out = 50),
                      s_MAP = seq(min(scale(sites$MAP)),
                                  max(scale(sites$MAP)),
                                  length.out = 50))
# Set abundance values for x axis
MAP_x_axis <- quantile(MAPpred$MAP)

### The Plots --------
#create the forest plot
png(file = "03_Figures/MS_Figs/IN_SPP_forest.png",
    width =2700, height = 2700, res = 300)
forest(metamod.in, 
       slab = NA,
       ilab = cbind(dat_meta_in$BOSQUE, dat_meta_in$GENUS_SPECIES),
       xlim = c(-3,2.5),
       ilab.xpos = c(-2.6, -1.6),
       header = F, 
       xlab = "rAME",
       order = MAP)
text(c(-2.6,-1.6, 1.9), metamod.in$k+2,
     c("Site", "Species", "Est. [95% CIs]"), cex=1.6, font=2)
dev.off()
#overall this looks good. Ideally, I get the species column to be italicized,
#but that's been too much of a struggle and I want to move on to the next plot.

#Prediction Plot
pred2 <- cbind(MAPpred, predict(object = metamod.in, newmods = MAPpred$s_MAP))

MAP_SPP_IN <- ggplot( ) + 
  geom_ribbon(data = pred2, aes(x = MAP, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "steelblue2", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_in, aes( x = as.numeric(MAP), y = estimate,
                                           ymin = estimate - std.error,
                                           ymax = estimate + std.error),
                  position = position_dodge2(width = 25)) +
  geom_line(data = pred2, aes(x = MAP, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Mean Annual Precipitation", y = "rAME") + 
  scale_x_continuous(breaks = MAP_x_axis,
                     labels = MAP_x_axis) + 
  theme_bw(12)
  

MAP_SPP_IN

ggsave(MAP_SPP_IN, filename = "03_Figures/MS_Figs/MAP_SPP_IN.png",
       width = 4, height = 4, units = "in")


