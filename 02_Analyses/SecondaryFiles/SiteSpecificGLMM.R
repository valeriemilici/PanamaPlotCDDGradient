###Initialize Workspace
rm(list =ls())
# Load libraries
library(tidyverse)
library(lme4)
library(lmerTest)

# Load Data
dat <- read.csv("2023data/2023censusdat.csv")

# Split the data into individual sites ----------------------
met <- dat %>% filter(BOSQUE == "Metropolitano")
anam <- dat %>% filter(BOSQUE == "ANAM")
charco <- dat %>% filter(BOSQUE == "Charco")
oleo <- dat %>% filter(BOSQUE == "Oleoducto")
pp <- dat %>% filter(BOSQUE == "Panama Pacifico")
bv <- dat %>% filter(BOSQUE == "Buena Vista")
sr <- dat %>% filter(BOSQUE == "Santa Rita")
sherman <- dat %>% filter(BOSQUE == "Sherman")

# filter each site to only include species for which there is variation in condens
spvar <- function(dat){
  dat1 <- dat %>% group_by(SPP) %>%
    summarize(spvar = var(con_dens)) %>%
    filter(spvar > 0)
  
  sp <- dat1$SPP
  
datvar <- dat %>% filter(SPP %in% sp) %>%
  group_by(SPP) %>%
  mutate(height = scale(ALT)) %>%
  ungroup()


}

#filter each site
met2 <- spvar(met)
anam2<-spvar(anam)
charco2<- spvar(charco)
oleo2<- spvar(oleo)
pp2<- spvar(pp)
bv2<- spvar(bv)
sr2<- spvar(sr)
sherman2<- spvar(sherman)
dat2<- spvar(dat)

### Model each site -------------------------------

met.surv <- glmer(SURV_NEXT ~ scale(con_dens) + scale(tot_dens) +
                    scale(height[,1]) +
                   (1|SPP)+ (1|PLOT) + (1|ID),
                 data = filter(met2, STATUS == 1),
                 family = binomial(link = "cloglog"),
                 glmerControl(optimizer = "bobyqa",
                              optCtrl = list(maxfun = 10000)),
                 offset = log(TIME_SINCE_LAST_CENSUS))
anam.surv <- glmer(SURV_NEXT ~ scale(con_dens) + scale(tot_dens) +
                    scale(height[,1]) +
                    (1|SPP)+ (1|PLOT) + (1|ID),
                  data = filter(anam2, STATUS == 1),
                  family = binomial(link = "cloglog"),
                  glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 10000)),
                  offset = log(TIME_SINCE_LAST_CENSUS))
charco.surv <- glmer(SURV_NEXT ~ scale(con_dens) + scale(tot_dens) +
                    scale(height[,1]) +
                    (1|SPP)+ (1|PLOT) + (1|ID),
                  data = filter(charco2, STATUS == 1),
                  family = binomial(link = "cloglog"),
                  glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 10000)),
                  offset = log(TIME_SINCE_LAST_CENSUS))
oleo.surv <- glmer(SURV_NEXT ~ scale(con_dens) + scale(tot_dens) +
                    scale(height[,1]) +
                    (1|SPP)+ (1|PLOT) + (1|ID),
                  data = filter(oleo2, STATUS == 1),
                  family = binomial(link = "cloglog"),
                  glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 10000)),
                  offset = log(TIME_SINCE_LAST_CENSUS))
bv.surv <- glmer(SURV_NEXT ~ scale(con_dens) + scale(tot_dens) +
                    scale(height[,1]) +
                    (1|SPP)+ (1|PLOT) + (1|ID),
                  data = filter(bv2, STATUS == 1),
                  family = binomial(link = "cloglog"),
                  glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 10000)),
                  offset = log(TIME_SINCE_LAST_CENSUS))
pp.surv <- glmer(SURV_NEXT ~ scale(con_dens) + scale(tot_dens) +
                    scale(height[,1]) +
                    (1|SPP)+ (1|PLOT) + (1|ID),
                  data = filter(pp2, STATUS == 1),
                  family = binomial(link = "cloglog"),
                  glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 10000)),
                  offset = log(TIME_SINCE_LAST_CENSUS))
sr.surv <- glmer(SURV_NEXT ~ scale(con_dens) + scale(tot_dens) +
                    scale(height[,1]) +
                    (1|SPP)+ (1|PLOT) + (1|ID),
                  data = filter(sr2, STATUS == 1),
                  family = binomial(link = "cloglog"),
                  glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 10000)),
                  offset = log(TIME_SINCE_LAST_CENSUS))
sherman.surv <- glmer(SURV_NEXT ~ scale(con_dens) + scale(tot_dens) +
                    scale(height[,1]) +
                    (1|SPP)+ (1|PLOT) + (1|ID),
                  data = filter(sherman2, STATUS == 1),
                  family = binomial(link = "cloglog"),
                  glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 10000)),
                  offset = log(TIME_SINCE_LAST_CENSUS))
dat.surv <- glmer(STATUS ~ scale(MAP)*(scale(con_dens) + scale(tot_dens)) +
                        scale(height[,1]) +
                        (1|SPP) + (1|PLOT/ID) + (1|BOSQUE),
                      data = filter(dat2, STATUS == 1),
                      family = binomial(link = "cloglog"),
                      glmerControl(optimizer = "bobyqa",
                                   optCtrl = list(maxfun = 10000)),
                      offset = log(TIME_SINCE_LAST_CENSUS))
#this offset term is wrong. It needs to be "Time to next census" because the 
#model as written is forward thinking.

summary(met.surv)
summary(anam.surv)
summary(charco.surv)
summary(pp.surv)
summary(sherman.surv)
summary(dat.surv)
