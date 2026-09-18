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
load(here("02_Analyses/(r)AMEOutput/mortality.Rdata")) #species specific rAME
rAME_ss <- read.csv("02_Analyses/(r)AMEOutput/rAME_sitespecific.csv")#site specific
     

rAME.spp.in <- rAME %>% filter(change == "invasion")  %>%
  mutate(BOSQUE = str_split_i(b.sp, "_", 1),
         SPP = str_split_i(b.sp, "_", -1)) %>%
  dplyr::select(-7)
rAME.site.in <- rAME_ss %>% filter(change == "invasion") %>%
  dplyr::rename("BOSQUE" = "b.sp") %>%
  mutate(SPP = "Comm. Ave.") %>%
  dplyr::select(-1)


rAME.spp.in <- left_join(rAME.spp.in, sites, by = "BOSQUE")
rAME.site.in <- left_join(rAME.site.in, sites, by = "BOSQUE")


#combine the estimates
rAME.combo.in <- rbind(rAME.spp.in, rAME.site.in)

# Reformat data for model fitting ------
dat_meta_in = metafor::escalc(measure = "GEN", # Set measure to generic, which passes the observed effect sizes or outcomes via the yi argument and the corresponding sampling variances via the vi argument (or the standard errors via the sei argument) to the function.
                              yi = estimate, # observed outcomes
                              sei = std.error, # standard errors
                              slab = BOSQUE, # label for species
                              data = rAME.combo.in)


metamod.in = metafor::rma(yi = yi,
                          vi = vi,
                          mods = ~ MAP_s,
                          method = "REML",
                          data = dat_meta_in)

summary(metamod.in)

#forest plot
png(file = "03_Figures/MS_Figs/IN_COMBO_forest.png",
    width =2700, height = 2700, res = 300)
forest(metamod.in, 
       slab = NA,
       ilab = cbind(dat_meta_in$BOSQUE, dat_meta_in$SPP),
       xlim = c(-3,2.5),
       ilab.xpos = c(-2.6, -1.6),
       header = F, 
       xlab = "rAME",
       order = MAP)
text(c(-2.6,-1.6, 1.9), metamod.in$k+2,
     c("Site", "Species", "Est. [95% CIs]"), cex=1.2, font=2)
dev.off()
# Generate a prediction dataframe ---------
MAPpred <- data_frame(MAP = seq(min(sites$MAP, na.rm = TRUE), 
                                max(sites$MAP, na.rm = TRUE),
                                length.out = 50),
                      s_MAP = seq(min(scale(sites$MAP)),
                                  max(scale(sites$MAP)),
                                  length.out = 50))
# Set abundance values for x axis
MAP_x_axis <- quantile(MAPpred$MAP)

#Prediction Plot
pred2 <- cbind(MAPpred, predict(object = metamod.in, newmods = MAPpred$s_MAP))

dat_meta_in2 <- dat_meta_in %>% mutate(type = ifelse(SPP == "Comm. Ave.", "Site", "Sp"))

MAP_COMBO_IN <- ggplot( ) + 
  geom_ribbon(data = pred2, aes(x = MAP, y = pred, ymin = ci.lb, ymax = ci.ub),
              fill = "darkgrey", alpha = 0.75) + 
  geom_pointrange(data = dat_meta_in2, aes( x = as.numeric(MAP), y = estimate,
                                           ymin = estimate - std.error,
                                           ymax = estimate + std.error,
                                           color = type),
                  position = position_dodge2(width = 25)) +
  scale_color_manual(values = c("red", "black")) +
  
  geom_line(data = pred2, aes(x = MAP, y = pred)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "Mean Annual Precipitation", y = "rAME") + 
  scale_x_continuous(breaks = MAP_x_axis,
                     labels = MAP_x_axis) + 
  theme_bw(12)+
theme(legend.position = "none") 
MAP_COMBO_IN

ggsave(plot = MAP_COMBO_IN, "03_Figures/MS_Figs/MAP_COMBO_IN.png",
       width = 4, height = 3, units = "in")
