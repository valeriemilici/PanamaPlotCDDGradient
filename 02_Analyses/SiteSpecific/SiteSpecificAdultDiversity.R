### Does the adult tree diversity vary with moisture? (Why are we working on this gradient?)
## Note: I wanted to get the census data from Panama Pacifico, Anam, and Santa Rita, but
## I couldn't find it on the CTFS portal. 

rm(list =ls())

library(tidyverse)
library(vegan)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(lme4)
library(lmerTest)
library(performance)

#load data
dat <- read.csv("2023data/siteadults.csv")
sites <- read.csv("2023data/SiteData_Gradient_Jeffs.csv")
sites[2,1] <- "Met"
sites[6,1] <- "BV"
sites[7,1] <- "SLZ"

#convert data into a format that is suitable for vegan package

#create a site.census column so that we can separate out each census
dat1 <- dat %>%
  group_by(Site, Mnemonic) %>%
  tally() %>%
  pivot_wider(names_from = Mnemonic,
              values_from = n) 
#replace nas with zeroes
dat1[is.na(dat1)] <- 0

#prepare new data frame for diversity metrics
div.indices <- data.frame(dat1$Site)
colnames(div.indices) <- "Site"

sites <- rename(sites, Site = BOSQUE)

div.indices <- left_join(div.indices, sites, by = "Site") #add precip info

#remove first column, convert the rest to numeric
dat2 <- data.frame(dat1 %>% ungroup() %>% dplyr::select(-c("Site")) %>%
                     sapply(as.numeric))

#calculate species richness
div.indices$SR <- specnumber(dat2)
#shannon's
div.indices$Shannon <- diversity(dat2)
#inverse simpson
div.indices$inv.simpson <- diversity(dat2, "invsimpson")
#simpson
div.indices$Simpson <- diversity(dat2, "simpson")

#make plots
Richness <- ggplot(div.indices, aes(x = MAP, y = SR)) +
  geom_smooth(method = "lm") +
  geom_point() +
  stat_regline_equation(label.x = 2000, label.y = 210) +
  stat_cor(aes(label = ..rr.label..), label.x = 2000, label.y = 200) +
  ylab("Species Richness") + 
  theme_bw()

Shannon <- ggplot(div.indices, aes(x = MAP, y = Shannon)) +
  geom_smooth(method = "lm") +
  geom_point() +
  stat_regline_equation(label.x = 2000, label.y = 4.25) +
  stat_cor(aes(label = ..rr.label..), label.x = 2000, label.y = 4.1) +
  ylab("Shannon's Diversity") +
  theme_bw()

Simpson <- ggplot(div.indices, aes(x = MAP, y = Simpson)) +
  geom_smooth(method = "lm") +
  geom_point() +
  stat_regline_equation(label.x = 2000, label.y = 1) +
  stat_cor(aes(label = ..rr.label..), label.x = 2000, label.y = 0.99) +
  ylab("Simpson's Diversity") +
  theme_bw()

divplots <- Richness/ Shannon / Simpson +
  plot_layout(axis_titles = "collect_x") +
  plot_annotation(tag_levels = "A")

divplots

ggsave(plot= divplots, filename = "figures/sitewidediversity_adults.jpg")


m1 <- lm(SR ~ scale(MAP),
           data = div.indices)

summary(m1) # no relationship

m2 <- lm(log(Shannon) ~ scale(MAP),
           data = div.indices)

summary(m2) # no relationship

m3 <- lm(Simpson ~ scale(MAP),
           data = div.indices)

summary(m3) # no relationship

m4 <- lm((SR) ~ scale(DRY_DEFICIT),
           data = div.indices)

summary(m4) # no relationship

m5 <- lm(Shannon ~ scale(DRY_DEFICIT),
           data = div.indices)

summary(m5) #no relationship 

m6 <- lm((Simpson) ~ scale(DRY_DEFICIT),
           data = div.indices)

summary(m6) # no relationship
