#look for weird patterns in the data for the species x site combos that 
#don't model well.
rm(list =ls())
library(tidyverse)
library(ggplot2)
library(car)
library(lme4)

#load data
dat <- read.csv("2025data/2025censusdat.csv")

dat$census <- factor(dat$CENSO)
dat$spp <- factor(dat$SPP)
dat$plot <- factor(dat$PLOT)
dat$height <- as.numeric(dat$ht.std)
dat$interval <- as.numeric(dat$TIME_SINCE_LAST_CENSUS)
dat$status <- as.numeric(if_else(dat$STATUS == 1, 0, 1)) #dead == 1
dat$b.spp <- str_c(dat$BOSQUE, "_", dat$spp)
dat$con_dens.p <- as.numeric(dat$con_dens.p)
dat$tot_dens.p <- as.numeric(dat$tot_dens.p)


# If I model with GLMM does con dens have an effect?
mort <- glmer(status ~  tot_dens.p + con_dens.p *scale(MAP) +
                  (1|CENSO) +  (1 + con_dens.p|spp) + (1+ tot_dens.p|spp) + 
                (1|plot) + (1|BOSQUE),
                family = binomial(link = "cloglog"),
              offset = log(interval),
                glmerControl(optimizer = "bobyqa",
                             optCtrl = list(maxfun = 10000)),
                data = filter(dat, !is.na(height)))

anova(mort)
summary(mort)
