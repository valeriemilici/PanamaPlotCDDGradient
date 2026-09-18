#### Comparing log likelihood values for the gams at each site, and considering
#### density differently each time. (Null, linear, nonlinear)

###Initialize Workspace --------------------------------------------------------
rm(list =ls())
# Load libraries
library(tidyverse)
library(mgcv) # For fitting gams
library(gratia)
library(mgcViz)

options(scipen = 999)
#load data

load("02_Analyses/SpeciesSpecific/res_mod") #nonlinear results
load("02_Analyses/SpeciesSpecific/res_red_mod") #null model results
load("02_Analyses/SpeciesSpecific/lin_mod") #linear model results

# Loglikelihood estimates -------------------------------------

 #nonlinear
nl_list <- list()

for(i in 1:length(res_mod)){
 nll <- signif(logLik.gam(res_mod[[i]]), digits = 3)
  
  nl_list[[i]] <- nll
}

#null model

n_list <- list()

for(i in 1:length(res_red_mod)){
  nll <- signif(logLik.gam(res_red_mod[[i]]), digits = 3)
  
  n_list[[i]] <- nll
}

#linear
l_list <- list()

for(i in 1:length(lin_mod)){
  nll <- signif(logLik.gam(lin_mod[[i]]), digits = 3)
  
  l_list[[i]] <- nll
}

# extract p-values from nl model

pv_list <- list()

for(i in 1:length(res_mod)){
  s <- summary(res_mod[[i]])
 pv_list[[i]] <- s$s.table[3,4]
}

# extract p-values from l model
pv_list_l <- list()

for(i in 1:length(lin_mod)){
  s <- summary(lin_mod[[i]])
  pv_list_l[[i]] <- s$p.table[2,4]
}

#Combine Results
list <- data.frame(cbind(mod = names(res_mod),
              null = n_list, linear = l_list, nonlinear = nl_list))

list2 <- list %>% separate_wider_delim(cols = mod,
                                       delim = "_",
                                       names = c("Site", "Species"),
                                       too_many = "merge", 
                                       cols_remove = F )

#remove col 1 from list
list <- list[2:4]
list$null <- as.numeric(list$null)
list$linear <- as.numeric(list$linear)
list$nonlinear <- as.numeric(list$nonlinear)
#list$p.value <- as.numeric(list$p.value)

#identify max log likelihood
list3 <- cbind(list2, maxval = apply(list, 1, FUN = max))

#identify best model
list4 <- list3 %>% mutate(bestmod = case_when(maxval == null ~ "null",
                                              maxval == linear ~ "linear",
                                              maxval == nonlinear ~ "nonlinear")) %>%
  dplyr::select(-maxval)



p.value <- unlist(pv_list)
p.value.l <- unlist(pv_list_l)

list4 <- cbind(list4, p.value, p.value.l)

list4$null <- as.character(list4$null)
list4$linear <- as.character(list4$linear)
list4$nonlinear <- as.character(list4$nonlinear)
list4$mod <- as.character(list4$mod)
list4$p.value <- as.character(list4$p.value)
list4$p.value.l <- as.character(list4$p.value.l)

# the models that are removed from the meta regression -
MR_fails <- c("ANAM_MARIPA LEWISII", "ANAM_PIPECO", "ANAM_PSYCB2",
              "Buena Vista_PYSCHO", "Charco_MICOAF", "Metropolitano_ANACEX",
              "Oleoducto_HIRTRA", "Oleoducto_HYBAPR", "Oleoducto_PROTCO",
              "Panama Pacifico_PIPAE", "Panama Pacifico_PSYCHO", 
              "Panama Pacifico_RINOSY", "Santa Rita_CUPASC", "Santa Rita_PIPEA4",
              "Santa Rita_POSOLA", "Santa Rita_PSYCGL", "Santa Rita_RINODA",
              "Santa Rita_RINOSQ", "Santa Rita_TAB1ST", "Sherman_HIRTRA",
              "Sherman_MICOSI", "Sherman_TACHVE")

list4 <- list4 %>% filter(!mod %in% MR_fails)


#until I figure out what's up with the p-values, remove the columns
list5 <- list4 %>% dplyr::select(1:7)


#save output
write.csv(list5, file = "03_Figures/MS_Figs/modelcomparisontable.csv")
