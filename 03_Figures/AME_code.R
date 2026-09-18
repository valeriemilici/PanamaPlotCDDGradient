

# AME ------------

AME_eq <- AME %>%
  filter(change == "equilibrium") %>%
  rename(BOSQUE = b.sp)

AME_in <- AME %>%
  filter(change == "invasion")%>%
  rename(BOSQUE = b.sp)

# Join marginal effects and MAP data ------
AME_eq <- left_join(AME_eq, sites, by = "BOSQUE")
AME_in <- left_join(AME_in, sites, by = "BOSQUE")

# Reformat data for model fitting -------
# AME ----------
dat_meta_eq_A = metafor::escalc(measure = "GEN", # Set measure to generic, which passes the observed effect sizes or outcomes via the yi argument and the corresponding sampling variances via the vi argument (or the standard errors via the sei argument) to the function.
                                yi = estimate, # observed outcomes
                                sei = std.error, # standard errors
                                slab = BOSQUE, # label for site
                                data = AME_eq)

dat_meta_in_A = metafor::escalc(measure = "GEN", # Set measure to generic, which passes the observed effect sizes or outcomes via the yi argument and the corresponding sampling variances via the vi argument (or the standard errors via the sei argument) to the function.
                                yi = estimate, # observed outcomes
                                sei = std.error, # standard errors
                                slab = BOSQUE, # label for site
                                data = AME_in)

# Fit model (MAP) ---------
# AME --------
#Equilibrium
metamod_eq = metafor::rma(yi = yi,
                          vi = vi,
                          mods = ~ MAP,
                          method = "REML",
                          data = dat_meta_eq)

summary(metamod_eq)

forest(metamod_eq, 
       header = "Site", 
       xlab = "rAME",
       order = MAP)

# Invasion

metamod_in_A = metafor::rma(yi = yi,
                            vi = vi,
                            mods = ~ MAP,
                            method = "REML",
                            data = dat_meta_in_A)

summary(metamod_in_A)

forest(metamod_in_A, 
       header = "Site", 
       xlab = "AME",
       order = MAP)

# Fit model (Dry Deficit) ---------
# AME ----------

# Equilibrium
metamod2_eq_A = metafor::rma(yi = yi,
                             vi = vi,
                             mods = ~ DRY_DEFICIT,
                             method = "REML",
                             data = dat_meta_eq_A)

summary(metamod2_eq_A)

forest(metamod2_eq_A, 
       header = "Site", 
       xlab = "AME",
       order = DRY_DEFICIT)

# Invasion
metamod2_in_A = metafor::rma(yi = yi,
                             vi = vi,
                             mods = ~ DRY_DEFICIT,
                             method = "REML",
                             data = dat_meta_in_A)

summary(metamod2_in_A)

forest(metamod2_in_A, 
       header = "Site", 
       xlab = "AME",
       order = DRY_DEFICIT)