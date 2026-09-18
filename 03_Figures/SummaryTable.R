### Script used to calculate values to put in a data summary table at the start
### of the results section 

### Initialize Workspace -------------------------------------------------------
rm(list =ls())

library(tidyverse)

## Data

dat <- read.csv("2023data/2023censusdat.csv")

#make sure the data are the appropriate type and named appropriately
dat$census <- factor(dat$CENSO)
dat$spp <- factor(dat$SPP)
dat$plot <- factor(dat$p_05)
dat$height <- as.numeric(dat$ALT_LAST_CENSUS)
dat$interval <- as.numeric(dat$TIME_SINCE_LAST_CENSUS)
dat$status <- as.numeric(dat$mort)
dat$con_dens <- as.numeric(dat$con_dens)
dat$tot_dens <- as.numeric(dat$tot_dens)
dat$BOSQUE <- factor(dat$BOSQUE)

# Pre-Filter

# How many species total?
length(unique(dat$spp)) # 399 across 55930 observations

# how many species across all sites?

dat %>% group_by(BOSQUE) %>%
  summarise(length(unique(spp)))

##Panama Pacifico = 67 species
##Metropolitano = 54 species
##ANAM = 86 Species
##Charco = 89 Species
##Oleoducto = 160 Species
##Buena Vista = 70 Species
##Sherman = 111 Species
##Santa Rita = 179 Species

# What's the median species richness at a cite for a given census?

dat %>% group_by(BOSQUE, census) %>%
  summarize(length(unique(spp))) %>%
  ungroup() %>%
  rename(spp = 'length(unique(spp))') %>%
  group_by(BOSQUE) %>%
  summarize(median(spp),
            sd(spp))

## Panama Pacifico = 54 (3.4) spp
## Metropolitano = 39 (1.6) spp
## ANAM = 57 (7.3) spp
## Charco = 67 (6)spp
## Oleoducto = 114 (20.7) spp
## Buena Vista = 53 (4.8) spp
## Sherman = 78 (9) spp
## Santa Rita = 152 (5.5) spp

# median number of observations?

dat %>% group_by(BOSQUE, census) %>%
  summarize(length(ID)) %>%
  ungroup() %>%
  rename(obs = 'length(ID)') %>%
  group_by(BOSQUE) %>%
  summarize(median(obs),
            sd(obs))
## Panama Pacifico = 477 (54.2)
## Metropolitano = 616 (164)
## ANAM = 518 (56)
## Charco = 402 (19.6)
## Oleoducto = 919 (248)
## Buena Vista = 584 (97)
## Sherman = 498 (112)
## Santa Rita = 1938 (243)


# Post - Filter

dat1 <- dat %>%
  # can't model mortality on a single census
  filter(FIRST_CENSUS != 1) %>%
  # height needs value for model to run
  filter(!is.na(height)) %>%
  # must be >= 200 mm to be within threshold
  filter(height >= 200) %>%
  # must be <= 1200 mm to be within threshold
  filter(height <= 1200) %>%
  mutate(b.spp = paste(BOSQUE, SPP, sep = "_"))

outliers <- dat1 %>% group_by(ID) %>% slice(1) %>% ungroup() %>%
  group_by(BOSQUE, SPP) %>% tally() %>%
  filter(n > 76) %>% #more than 76 obs is an outlier ( > 0.975 quantile)
  mutate(b.spp = paste(BOSQUE, SPP, sep = "_"))

outs <- outliers$b.spp

dat2 <- dat1 %>% filter(!b.spp %in% outs) %>%
  group_by(BOSQUE) %>%
  # add two to the high quantile to preserve a little bit more data for the K-estimation
  filter(con_dens <= quantile(con_dens, 0.975) + 2) %>%
  filter(tot_dens <= quantile(tot_dens, 0.975) + 2) %>%
  ungroup()

# with the filtering we have removed 34293 observations

# How many species total?
length(unique(dat2$spp)) # 379 spp across 21637 observations
# we only lose 20 species due to filtering

# how many species across all sites?

dat2 %>% group_by(BOSQUE) %>%
  summarise(length(unique(spp)))
## PP = 61 spp
## MET = 47 spp
## ANAM = 78 spp
## Charco = 83 spp
## Oleo = 149 spp
## BV = 65 spp
## Sherman = 104 spp
## SR = 167 spp

# What's the median species richness at a site across censuses?

dat2 %>% group_by(BOSQUE, census) %>%
  summarize(length(unique(spp))) %>%
  ungroup() %>%
  rename(spp = 'length(unique(spp))') %>%
  group_by(BOSQUE) %>%
  summarize(median(spp),
            sd(spp))
## PP = 42 (3.5)
## MET = 29 (4.2)
## ANAM = 43 (5)
## Charco = 52 (4.6)
## Oleo = 88 (16)
## BV = 44 (3)
## Sherman = 65 (7)
## SR = 135 (4.4)

# median number of observations?

dat2 %>% group_by(BOSQUE, census) %>%
  summarize(length(ID)) %>%
  ungroup() %>%
  rename(obs = 'length(ID)') %>%
  group_by(BOSQUE) %>%
  summarize(median(obs),
            sd(obs))

## PP = 138 (5.5)
## Met = 117 (13)
## ANAM = 168 (37.4)
## Charco = 197 (10.4)
## Oleo = 444 (103)
## BV = 358 (58)
## Sherman = 266 (51.5)
## SR = 826 (33)

#median total and conspecific densities?
dat2 %>% group_by(BOSQUE) %>%
  summarize(quantile(con_dens, na.rm = T, 0.5),
            quantile(con_dens, na.rm = T, 0.1),
            quantile(con_dens, na.rm = T, 0.9),
            quantile(tot_dens, na.rm = T, 0.5),
            quantile(tot_dens, na.rm = T, 0.1),
            quantile(tot_dens, na.rm = T, 0.9)) 

