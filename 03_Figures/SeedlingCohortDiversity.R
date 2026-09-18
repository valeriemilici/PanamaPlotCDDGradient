#Objective: Extract data on cohort survival over time, convert to measures of 
# diversity, assess whether cohort diversity increases over time (and at different
# rates across the sites). 

#Initialize Workspace ----------------------------------------------
rm(list =ls())

library(tidyverse)
library(ggplot2)
library(patchwork)
library(lme4)
library(lmerTest)
library(performance)
library(vegan)
library(ggpubr)
library(ggeffects)
library(mgcv) # For fitting gams

# Load Data 

dat <- read.csv("2023data/2023censusdat.csv")
sites <- read.csv("2023data/SiteData_Gradient_Jeffs.csv")
sites[3,1] <- "ANAM"

# Step 1: Extract data on cohort survival over time -----------------
dat$census <- factor(dat$CENSO)
dat$spp <- factor(dat$SPP)
dat$plot <- factor(dat$p_05)
dat$height <- as.numeric(dat$ALT_LAST_CENSUS)
dat$interval <- as.numeric(dat$TIME_SINCE_LAST_CENSUS)
dat$status <- as.numeric(dat$mort)
dat$con_dens <- as.numeric(dat$con_dens)
dat$tot_dens <- as.numeric(dat$tot_dens)
dat$BOSQUE <- factor(dat$BOSQUE)

#create a df of when and where each individual was first recorded. 
# Use to create cohorts
dat1 <- dat %>% 
  filter(FIRST_CENSUS == 1) |> #first observation for each seedling
  filter(CENSO <= 8) |> #want at least five years of change
  group_by(BOSQUE, census) %>%
  summarise(ID,
            SPP,
            GENUS_SPECIES,
            p_20,
            p_05) %>%
  ungroup()

dat1$census <- as.numeric(as.character(dat1$census))

#get each cohort time-series (there is no census 4)
C1 <- dat1 %>% filter(census == 1) 
C2 <- dat1 %>% filter(census == 2)
C3 <- dat1 %>% filter(census == 3)

C5 <- dat1 %>% filter(census == 5)
C6 <- dat1 %>% filter(census == 6)
C7 <- dat1 %>% filter(census == 7)
C8 <- dat1 %>% filter(census == 8)

ins1 <- C1$ID #a list of C1 IDs to filter main data
ins2 <- C2$ID
ins3 <- C3$ID

ins5 <- C5$ID
ins6 <- C6$ID
ins7 <- C7$ID
ins8 <- C8$ID

cohort1 <- dat %>% filter(ID %in% ins1)
cohort2 <- dat %>% filter(ID %in% ins2)
cohort3 <- dat %>% filter(ID %in% ins3)

cohort5 <- dat %>% filter(ID %in% ins5)
cohort6 <- dat %>% filter(ID %in% ins6)
cohort7 <- dat %>% filter(ID %in% ins7)
cohort8 <- dat %>% filter(ID %in% ins8)

# Step 2: Convert to diversity over time ---------------------------------------
#create a site.census column so that we can separate out each census
div1 <- cohort1 %>%
  group_by(BOSQUE, census, SPP) %>%
  tally() %>%
  pivot_wider(names_from = SPP,
              values_from = n) 
div1[is.na(div1)] <- 0

div2 <- cohort2 %>%
  group_by(BOSQUE, census, SPP) %>%
  tally() %>%
  pivot_wider(names_from = SPP,
              values_from = n) 
div2[is.na(div2)] <- 0

div3 <- cohort3 %>%
  group_by(BOSQUE, census, SPP) %>%
  tally() %>%
  pivot_wider(names_from = SPP,
              values_from = n) 
div3[is.na(div3)] <- 0

div5 <- cohort5 %>%
  group_by(BOSQUE, census, SPP) %>%
  tally() %>%
  pivot_wider(names_from = SPP,
              values_from = n) 
div5[is.na(div5)] <- 0

div6 <- cohort6 %>%
  group_by(BOSQUE, census, SPP) %>%
  tally() %>%
  pivot_wider(names_from = SPP,
              values_from = n) 
div6[is.na(div6)] <- 0

div7 <- cohort7 %>%
  group_by(BOSQUE, census, SPP) %>%
  tally() %>%
  pivot_wider(names_from = SPP,
              values_from = n) 
div7[is.na(div7)] <- 0

div8 <- cohort8 %>%
  group_by(BOSQUE, census, SPP) %>%
  tally() %>%
  pivot_wider(names_from = SPP,
              values_from = n) 
div8[is.na(div8)] <- 0

#prepare new data frame for diversity metrics

div.ind1 <- div1[,c(1:2)]
div.ind2 <- div2[,c(1:2)]
div.ind3 <- div3[,c(1:2)]
div.ind5 <- div5[,c(1:2)]
div.ind6 <- div6[,c(1:2)]
div.ind7 <- div7[,c(1:2)]
div.ind8 <- div8[,c(1:2)]

#remove first column, convert the rest to numeric
div1red <- data.frame(div1 %>% ungroup() %>% dplyr::select(-c(1:2)) %>%
                     sapply(as.numeric))
div2red <- data.frame(div2 %>% ungroup() %>% dplyr::select(-c(1:2)) %>%
                        sapply(as.numeric))
div3red <- data.frame(div3 %>% ungroup() %>% dplyr::select(-c(1:2)) %>%
                        sapply(as.numeric))
div5red <- data.frame(div5 %>% ungroup() %>% dplyr::select(-c(1:2)) %>%
                        sapply(as.numeric))
div6red <- data.frame(div6 %>% ungroup() %>% dplyr::select(-c(1:2)) %>%
                        sapply(as.numeric))
div7red <- data.frame(div7 %>% ungroup() %>% dplyr::select(-c(1:2)) %>%
                        sapply(as.numeric))
div8red <- data.frame(div8 %>% ungroup() %>% dplyr::select(-c(1:2)) %>%
                        sapply(as.numeric))

#calculate species richness
div.ind1$SR <- specnumber(div1red)
div.ind2$SR <- specnumber(div2red)
div.ind3$SR <- specnumber(div3red)
div.ind5$SR <- specnumber(div5red)
div.ind6$SR <- specnumber(div6red)
div.ind7$SR <- specnumber(div7red)
div.ind8$SR <- specnumber(div8red)
#shannon's
div.ind1$Shannon <- diversity(div1red)
div.ind2$Shannon <- diversity(div2red)
div.ind3$Shannon <- diversity(div3red)
div.ind5$Shannon <- diversity(div5red)
div.ind6$Shannon <- diversity(div6red)
div.ind7$Shannon <- diversity(div7red)
div.ind8$Shannon <- diversity(div8red)
#inverse simpson
div.ind1$inv.simpson <- diversity(div1red, "invsimpson")
div.ind2$inv.simpson <- diversity(div2red, "invsimpson")
div.ind3$inv.simpson <- diversity(div3red, "invsimpson")
div.ind5$inv.simpson <- diversity(div5red, "invsimpson")
div.ind6$inv.simpson <- diversity(div6red, "invsimpson")
div.ind7$inv.simpson <- diversity(div7red, "invsimpson")
div.ind8$inv.simpson <- diversity(div8red, "invsimpson")
#simpson
div.ind1$Simpson <- diversity(div1red, "simpson")
div.ind2$Simpson <- diversity(div2red, "simpson")
div.ind3$Simpson <- diversity(div3red, "simpson")
div.ind5$Simpson <- diversity(div5red, "simpson")
div.ind6$Simpson <- diversity(div6red, "simpson")
div.ind7$Simpson <- diversity(div7red, "simpson")
div.ind8$Simpson <- diversity(div8red, "simpson")
#Pielou's evenness (J')
div.ind1$evenness <- div.ind1$Shannon/log(div.ind1$SR)
div.ind2$evenness <- div.ind2$Shannon/log(div.ind2$SR)
div.ind3$evenness <- div.ind3$Shannon/log(div.ind3$SR)
div.ind5$evenness <- div.ind5$Shannon/log(div.ind5$SR)
div.ind6$evenness <- div.ind6$Shannon/log(div.ind6$SR)
div.ind7$evenness <- div.ind7$Shannon/log(div.ind7$SR)
div.ind8$evenness <- div.ind8$Shannon/log(div.ind8$SR)
#add cohort column
div.ind1$cohort <- 1
div.ind2$cohort <- 2
div.ind3$cohort <- 3
div.ind5$cohort <- 5
div.ind6$cohort <- 6
div.ind7$cohort <- 7
div.ind8$cohort <- 8

div.indices <- rbind( div.ind2, div.ind3, div.ind5,
                     div.ind6, div.ind7, div.ind8) #div.ind1,
#remove the first census because it seems suspiciously high. Prob includes way
#more obs than just 1 cohort of seedlings. 

div.indices <- left_join(div.indices, sites, by = "BOSQUE") #add precip info
div.indices <- div.indices[,c(1:9, 12, 15:17)] #reduce

div.indices$census <- as.numeric(as.character(div.indices$census))
div.indices <- div.indices %>% mutate(Census = (census + 1)- cohort)

# Part 3: Does Diversity Change over Time? ------------------------------------

div.indices$BOSQUE <- factor(div.indices$BOSQUE, levels = c("Panama Pacifico",
                                                            "Metropolitano",
                                                            "ANAM",
                                                            "Charco",
                                                            "Oleoducto",
                                                            "Buena Vista",
                                                            "Sherman",
                                                            "Santa Rita"))

#initial plots
Richness <-  ggplot(div.indices, aes(x = Census, y = SR, group = BOSQUE)) +
  geom_smooth(method = "lm", aes(color = BOSQUE,  fill = BOSQUE)) +
  geom_point(aes(color = BOSQUE),
             position = position_dodge2(width = 0.4)) +
  scale_fill_viridis_d(option = "cividis", direction = -1, name = "Site") +
  scale_color_viridis_d(option = "cividis", direction = -1, name = "Site") +
  ylab("Species Richness") +
  scale_x_continuous(breaks = seq(1,13, 1)) +
  theme_classic(12)

Shannon <- ggplot(div.indices, aes(x = Census, y = Shannon, group = BOSQUE)) +
  geom_smooth(method = "lm", aes(color = BOSQUE,  fill = BOSQUE)) +
  geom_point(aes(color = BOSQUE),
             position = position_dodge2(width = 0.4)) +
  scale_fill_viridis_d(option = "cividis", direction = -1, name = "Site") +
  scale_color_viridis_d(option = "cividis", direction = -1, name = "Site") +
  ylab("Shannon's Diversity") +
  scale_x_continuous(breaks = seq(1,13, 1)) +
  theme_classic(12)

Simpson <- ggplot(div.indices, aes(x = Census, y = Simpson, group = BOSQUE)) +
  geom_smooth(method = "lm", aes(color = BOSQUE,  fill = BOSQUE)) +
  geom_point(aes(color = BOSQUE),
             position = position_dodge2(width = 0.4)) +
  scale_fill_viridis_d(option = "cividis", direction = -1, name = "Site") +
  scale_color_viridis_d(option = "cividis", direction = -1, name = "Site") +
  ylab("Simpson's Diversity") +
  scale_x_continuous(breaks = seq(1,13, 1)) +
  theme_classic(12)

InvSimpson <- ggplot(div.indices, aes(x = Census, y = inv.simpson, group = BOSQUE)) +
  geom_smooth(method = "lm", aes(color = BOSQUE,  fill = BOSQUE)) +
  geom_point(aes(color = BOSQUE),
             position = position_dodge2(width = 0.4)) +
  scale_fill_viridis_d(option = "cividis", direction = -1, name = "Site") +
  scale_color_viridis_d(option = "cividis", direction = -1, name = "Site") +
  ylab("Inverse Simpson's Diversity") +
  scale_x_continuous(breaks = seq(1,13, 1)) +
  theme_classic(12)

InvSimpson

divplots <- Richness/ Shannon / InvSimpson +
  plot_layout(axis_titles = "collect_x", guides = "collect") +
  plot_annotation(tag_levels = "A")

divplots

ggsave(divplots, filename ="03_Figures/MS_Figs/seedlingcensusdiversity.png",
       height = 8, width = 6, units = "in")

m1 <- lmer(log(SR) ~ scale(MAP) + census + (1|BOSQUE) + (1|cohort),
           data = div.indices)


m2 <- lmer(Shannon ~ scale(MAP) + census + (1|cohort) + (1|BOSQUE),
           data = div.indices)

m3 <- lmer(inv.simpson ~ scale(MAP) + census + (1|cohort) + (1|BOSQUE),
           data = div.indices)

simp_mod <-lmer(inv.simpson ~ scale(MAP) * census + BOSQUE+  (1|cohort) ,
                data = div.indices)
summary(m1)

summary(m2)

summary(m3)

library(performance)
check_model(m1)

check_model(m3)

## BONUS: How does Evenness (J') change over time? --------

#plot the data --> expectations
Evenplot <- ggplot(div.indices, aes(x = Census, y = evenness, group = BOSQUE)) +
  geom_smooth(method = "lm", aes(color = BOSQUE,  fill = BOSQUE)) +
  geom_point(aes(color = BOSQUE),
             position = position_dodge2(width = 0.4)) +
  scale_fill_viridis_d(option = "cividis", direction = -1, name = "Site") +
  scale_color_viridis_d(option = "cividis", direction = -1, name = "Site") +
  ylab("Pielou's evenness (J')") +
  scale_x_continuous(breaks = seq(1,13, 1)) +
  theme_classic(12)

Evenplot
#according to the data, we can expect that there is a positive relationship
#where evenness increases over time. This relationship will be strongest at
#Metropolitano, presumably because it has the lowest initial evenness and 
#therefore is capable of the greatest change in evenness. Otherwise, I don't
#expect any mind-blowing trends to pop out of this. 

#simple model:
even_mod <- lmer(evenness ~ scale(MAP) + census + (1|BOSQUE) + (1|cohort),
                 data = div.indices)
#FYI- whether or not BOSQUE is included as a FE or RE does not change the value
#or significance of the FE estimates for MAP and census. The fact that MAP and
#BOSQUE are not independent causes some slight issues, so I'm keeping BOSQUE
#as a RE.

check_model(even_mod)
#the low values in the dataset (Metropolitano) are causing some non-ideal 
#deviations from the expectation of linearity. The "log" cheat does not help.
#could be that there is a non-linear response for Metropolitano. Consider GAM

summary(even_mod)

censusin <- c(2,13)
evenpreds <- predict_response(even_mod, c("MAP", "census"), type = "fixed") %>%
 filter(group %in% censusin) %>% #first and last census
  group_by(x) %>%
  mutate(difference = predicted - lag(predicted))


ggplot(data = evenpreds, mapping = aes(x = x, y = predicted)) +
         geom_rect(aes(xmin = as.numeric(x) + 10, 
                       xmax = as.numeric(x) - 10,
                       ymin = predicted - difference,
                       ymax = predicted)) +
  theme_bw()

simppreds <- predict_response(simp_mod, c("MAP", "census"), type = "fixed") %>%
  filter(group %in% censusin) %>% #first and last census
  group_by(x) %>%
  mutate(difference = predicted - lag(predicted))

ggplot(data = simppreds, mapping = aes(x = x, y = predicted)) +
  geom_rect(aes(xmin = as.numeric(x) + 10, 
                xmax = as.numeric(x) - 10,
                ymin = predicted - difference,
                ymax = predicted)) +
  theme_bw()
