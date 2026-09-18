###Initialize Workspace --------------------------------------------------------
rm(list =ls())
# Load libraries
library(tidyr)
library(dplyr)
library(lubridate) # For calculating census intervals


#load data
dat <- read.csv("2025data/cleaned seedling census - excluding post dry season censuses - 2025-04-27.csv")
sites <- read.csv("2023data/SiteData_Gradient_Jeffs.csv")

#Remove observations with missing information
dat1 <- dat %>% filter(EXCLUDE_SPECIES != 1) %>%
  #remove Lianas
  filter(LIANA == 0) %>%
  #add columns for survival/mortality next census
  arrange(CENSO) %>% group_by(PLOT, TAG) %>% 
  mutate(SURV_NEXT = if_else(lead(STATUS == 1), 1,0)) %>%
  mutate(MORT_NEXT = if_else(lead(STATUS == 0), 1, 0))
  

#calculate conspecific and total density in each plot
Cdens <- dat1 %>% group_by(PLOT, CENSO, SPP) %>%
  tally() %>%
  mutate(con_dens = n - 1)

Cdens <- Cdens[,-4]

Tdens <- dat1 %>% group_by(PLOT, CENSO) %>%
  tally() %>%
  mutate(tot_dens = n - 1) 

Tdens <- Tdens[,-3]

#merge back together
dat2 <- left_join(dat1, Cdens, join_by(PLOT, CENSO, SPP))
dat3 <- left_join(dat2, Tdens, join_by(PLOT, CENSO))

#create het dens column
dat4 <- dat3 %>% mutate(het_dens = tot_dens - con_dens) %>%
  group_by(PLOT, ID) %>%
  arrange(CENSO, .by_group = T) %>%
  mutate(con_dens.p = lag(con_dens),
         tot_dens.p = lag(tot_dens)) %>%
  ungroup() %>%
  #remove each individual's first census (can't be modeled)
  filter(FIRST_CENSUS != 1) %>%
  #remove seedlings that fall outside of the height threshold
  filter(!is.na(ALT_LAST_CENSUS)) %>%
  filter(ALT_LAST_CENSUS >= 200) %>%
  filter(ALT_LAST_CENSUS <= 1200) %>%
  #standardize height within species
  ungroup() %>%
  group_by(SPP) %>%
  mutate(ht.std = scale(ALT_LAST_CENSUS))

# rename one site to match dat
sites[3,1] <- "ANAM"

dat4 <- left_join(dat4, sites, by = "BOSQUE")

write.csv(dat4, "2025data/2025censusdat.csv")
