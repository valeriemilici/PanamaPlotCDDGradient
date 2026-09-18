### Diversity Measures Per Site
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
dat <- read.csv("2023data/2023censusdat.csv")
sites <- read.csv("2023data/SiteData_Gradient_Jeffs.csv")
sites[3,1] <- "ANAM"

#convert data into a format that is suitable for vegan package

#create a site.census column so that we can separate out each census
dat1 <- dat %>% mutate(site.census = paste(BOSQUE, CENSO, sep = "_")) %>%
  group_by(SPP, site.census) %>%
  tally() %>%
  pivot_wider(names_from = SPP,
              values_from = n) 
#replace nas with zeroes
dat1[is.na(dat1)] <- 0

#prepare new data frame for diversity metrics
div.indices <- data.frame(dat1$site.census)
div.indices[c("site", "census")] <- str_split_fixed(div.indices$dat1.site.census,
                                                    "_", n = 2)
sites <- rename(sites, site = BOSQUE)

div.indices <- merge(div.indices, sites, by = "site") #add precip info

#remove first column, convert the rest to numeric
dat2 <- data.frame(dat1 %>% dplyr::select(2:400) %>%
  sapply(as.numeric))
 
#calculate species richness
div.indices$SR <- specnumber(dat2)
#shannon's
div.indices$Shannon <- diversity(dat2)
#inverse simpson
div.indices$inv.simpson <- diversity(dat2, "invsimpson")
#simpson
div.indices$Simpson <- diversity(dat2, "simpson")


#summarize data for plotting
dat3 <- div.indices %>% group_by(site) %>%
  summarize( mu_SR = mean(as.numeric(SR)),
             sd_SR = sd(as.numeric(SR)),
             mu_H = mean(Shannon),
            sd_H = sd(Shannon),
             mu_d = mean(Simpson),
           sd_d = sd(Simpson),
             mu_id = mean(inv.simpson),
            sd_id = sd(inv.simpson),
             MAP = mean(MAP),
             DRY_DEFICIT = mean(DRY_DEFICIT))

write.csv(div.indices, "2023data/seedlingdiversityindices.csv")

#make plots
Richness <- ggplot(dat3, aes(x = MAP, y = mu_SR)) +
  geom_smooth(method = "lm") +
  geom_pointrange(aes(ymin = mu_SR - sd_SR, ymax = mu_SR + sd_SR)) +
  stat_regline_equation(label.x = 2800, label.y = 150) +
  stat_cor(aes(label = ..rr.label..), label.x = 2800, label.y = 140) +
  ylab("Species Richness") + 
  theme_bw()

Shannon <- ggplot(dat3, aes(x = MAP, y = mu_H)) +
  geom_smooth(method = "lm") +
  geom_pointrange(aes(ymin = mu_H - sd_H, ymax = mu_H + sd_H)) +
  stat_regline_equation(label.x = 2800, label.y = 4.25) +
  stat_cor(aes(label = ..rr.label..), label.x = 2800, label.y = 4.1) +
  ylab("Shannon's Diversity") +
  theme_bw()

Simpson <- ggplot(dat3, aes(x = MAP, y = mu_d)) +
  geom_smooth(method = "lm") +
  geom_pointrange(aes(ymin = mu_d - sd_d, ymax = mu_d + sd_d)) +
  stat_regline_equation(label.x = 2800, label.y = 0.8) +
  stat_cor(aes(label = ..rr.label..), label.x = 2800, label.y = 0.78) +
  ylab("Simpson's Diversity") +
  theme_bw()

divplots <- Richness/ Shannon / Simpson +
  plot_layout(axis_titles = "collect_x") +
  plot_annotation(tag_levels = "A")

divplots

ggsave(plot= divplots, filename = "figures/sitewidediversity.jpg")


## Finally some quick linear models to show there is no relationship
#all models show that there is no relationship. I think heteroscedasticity is 
#causing some issues with these models, but the diagnostics are decent (not perfect)
#and overall agree with what we see in the simple plots. There is no significant
#relationship between moisture and diversity. 

m1 <- lmer(SR ~ scale(MAP) + (1|census) + (1|site),
           data = div.indices)

summary(m1)


m2 <- lmer(log(Shannon) ~ scale(MAP) + (1|census) + (1|site),
           data = div.indices)

summary(m2)

check_model(m2)

qqnorm(resid(m2))
qqline(resid(m2))

plot(div.indices$MAP, div.indices$Shannon)
#No relationship between Shannon's diversity and MAP

m3 <- lmer(Simpson ~ scale(MAP) + (1|census) + (1|site),
           data = div.indices)

summary(m3)
#Marginal positive relationship between Simpson's and MAP

m4 <- lmer(log(SR) ~ scale(DRY_DEFICIT) + (1|census) + (1|site),
           data = filter(div.indices, SR <= 130))

summary(m4)


m5 <- lmer(Shannon ~ scale(DRY_DEFICIT) + (1|census) + (1|site),
           data = div.indices)

summary(m5)
#no relationship between Shannon and Deficit

m6 <- lmer((Simpson) ~ scale(DRY_DEFICIT) + (1|census) + (1|site),
           data = filter(div.indices, Simpson >= 0.75))

summary(m6)

#same as before a marginal positive relationship. It's really close to the 
#common cut-off though so it's essentially significant. Simpson is the best
#measure of diversity if we want to claim there is a moisture-diversity gradient,
#the slope is significant but not very strong. 

#There were some outliers on the low end of the data that were driving the relationship
#and causing the poor model fit. Thinning the data to remove the outliers greatly
#improves model fit, but now there is absolutely no relationship between Simpson's
#diversity and dry deficit