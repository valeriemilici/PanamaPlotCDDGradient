###Initialize Workspace --------------------------------------------------------
rm(list =ls())
# Load libraries
library(readr)
library(skimr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(parallel)
library(pbapply) # Adds progress bar to apply functions
library(here)
library(spatstat.geom)
library(mgcv) # For fitting gams
library(lubridate) # For calculating census intervals
library(broom) # For processing fit models
library(gratia)
library(boot)
library(mgcViz)
library(MASS)
library(stringr)

#load data
dat <- read.csv("2025data/2025censusdat.csv")

dds_list <- c("ANAM_MARIPA LEWISII", "ANAM_PIPECO", "ANAM_PSYCB2",
              "Buena Vista_PYSCHO", "Charco_MICOAF", "Metropolitano_ANACEX",
              "Oleoducto_HIRTRA", "Oleoducto_HYBAPR", "Oleoducto_PROTCO",
              "Panama Pacifico_PIPAE", "Panama Pacifico_PSYCHO", 
              "Panama Pacifico_RINOSY", "Santa Rita_CUPASC", "Santa Rita_PIPEA4",
              "Santa Rita_POSOLA", "Santa Rita_PSYCGL", "Santa Rita_RINODA",
              "Santa Rita_RINOSQ", "Santa Rita_TAB1ST", "Sherman_HIRTRA",
              "Sherman_MICOSI", "Sherman_TACHVE")
#the list of models that have too much SE to be included in the metaregression.
# note- When they are added to the data deficient seedlings for each site,
# it causes the data deficient models to fail. They must be removed entirely.

#make sure the data are the appropriate type and named appropriately
dat$census <- factor(dat$CENSO)
dat$spp <- factor(dat$SPP)
dat$plot <- factor(dat$PLOT)
dat$height <- as.numeric(dat$ht.std)
dat$interval <- as.numeric(dat$TIME_SINCE_LAST_CENSUS)
dat$status <- as.numeric(if_else(dat$STATUS == 1, 0, 1)) #dead == 1
dat$b.spp <- str_c(dat$BOSQUE, "_", dat$spp)
dat$con_dens.p <- as.numeric(dat$con_dens.p)
dat$tot_dens.p <- as.numeric(dat$tot_dens.p)


#I am looking at species within a site, so I need to group by Site and then Spp

nval = 3  ## this is the number of unique conspecific values req to run model
minrange = 1    # minimum range for conspecific density
# Val note- minrange is redundant given nval.
# If there are four unique conspecific values, then the range req. has been met
# nval = 4 is more important here for spline fitting. And setting the k-values

nsp <- dat %>% 
  group_by(BOSQUE, spp) %>% 
  summarise(
    range_con_dens = max(con_dens) - min(con_dens),
    max_con_dens = max(con_dens),
    unique_con_dens = length(unique(con_dens)),
    unique_total_dens = length(unique(tot_dens)),
    unique_height = length(unique(ALT_LAST_CENSUS))
  ) %>% 
  
  
  mutate(issue_nval = unique_con_dens < nval,              # less than nval unique values in consp densities
         issue_range = range_con_dens < minrange,                    # range should at least be equal to minrange
         trymodel = !(issue_nval|issue_range),   # 
         rare = !trymodel # preliminary assignment of rare species
  ) %>%
  unite("b.spp", BOSQUE:spp, remove = F) #combine site and species 

table(nsp$trymodel)
#126 sitexspecies combos have enough data to be modeled nonlinearly

# Functions for fitting models -----------
#Null and Nonlinear Models
model_fit = function(data, speciesinfo, reduced = F) {
  
  # create new factor with correct factor levels per species 
  data$census = factor(data$census)
  
  
  # create model formula
  term_c = ifelse(length(unique(data$census)) > 1, "+ s(census, bs = 're')", "") 
  term_p = ifelse(length(unique(data$plot)) > 1, "+ s(plot, bs = 're')", "")
  term_s = ifelse(length(unique(data$spp)) >1, " + s(spp, bs = 're')", "")
  
  
  if (reduced) {
    form =  as.formula(paste0("status ~ s(height, k = k1) +
                              s(tot_dens.p, k = k2)"
                              , term_c, term_p, term_s)) # reduced model
  } else {
    form =  as.formula(paste0("status ~ s(height, k = k1) +
                              s(tot_dens.p, k = k2)  +
                              s(con_dens.p, k = k3)" 
                              , term_c, term_p, term_s)) # full model
  }
  
  # Choose penalty
  # set to default k=10 
  k1 = k2 = k3 = 10
  if (k1 > speciesinfo$unique_height) k1 = speciesinfo$unique_height - 2
  if (k2 > speciesinfo$unique_total_dens) k2 = speciesinfo$unique_total_dens - 2
  if (k3 > speciesinfo$unique_con_dens) k3 = speciesinfo$unique_con_dens - 1
  # Fit model
  mod = try(gam(form
                , family = binomial(link=cloglog)
                , offset = log(interval)
                , data = data
                , method = "ML"
  ) , silent = T
  )
  
  return(mod)
  
}

#Linear Model 
linear_model_fit = function(data, speciesinfo, reduced = F) {
  
  # create new factor with correct factor levels per species 
  data$census = factor(data$census)
  
  
  # create model formula
  term_c = ifelse(length(unique(data$census)) > 1, "+ s(census, bs = 're')", "") 
  term_p = ifelse(length(unique(data$plot)) > 1, "+ s(plot, bs = 're')", "")
  term_s = ifelse(length(unique(data$spp)) >1, " + s(spp, bs = 're')", "")
  
  
  form =  as.formula(paste0("status ~ s(height, k = k1) +  s(tot_dens.p, k = k2)  +
                              con_dens.p" 
                            , term_c, term_p, term_s))
  
  
  # Choose penalty
  # set to default k=10 
  k1 = k2 = k3 = 10
  if (k1 > speciesinfo$unique_height) k1 = speciesinfo$unique_height - 2
  if (k2 > speciesinfo$unique_total_dens) k2 = speciesinfo$unique_total_dens - 2
  if (k3 > speciesinfo$unique_con_dens) k3 = speciesinfo$unique_con_dens - 1
  # Fit model
  lin_mod = try(gam(form
                , family = binomial(link=cloglog)
                , offset = log(interval)
                , data = data
                , method = "ML"
  ) , silent = T
  )
  
  return(lin_mod)
  
}

# check model run

model_convergence = function(model) {
  
  # gam not available
  if (!any(class(model)=="gam")) {
    print(paste(b.spp, "gam failed"))
  } else {
    
    # gam not converged
    if (!model$converged) {
      print(paste(b.spp, "no convergence"))
    } else {
      
      
      # Explore warning "glm.fit: fitted probabilities numerically 0 or 1 occurred (complete separation)"
      eps <- 10 * .Machine$double.eps
      glm0.resids <- augment(x = model) %>%
        mutate(p = 1 / (1 + exp(-.fitted)),
               warning = p > 1-eps,
               influence = order(.hat, decreasing = T))
      infl_limit = round(nrow(glm0.resids)/10, 0)
      # check if none of the warnings is among the 10% most influential observations, than it is okay..
      num = any(glm0.resids$warning & glm0.resids$influence < infl_limit)
      
      # complete separation
      if (num) {
        print(paste(b.spp, "complete separation is likely"))
      } else {
        
        # missing Vc
        if (is.null(model$Vc)) {
          print(paste(b.spp, "Vc not available"))
        } else {
          
          # successful model
          return(model)
        }
      }
    }
  }
}

## Fit models each species that is not data deficient ----------

#table(data_deficient = nsp$rare, trymodel = nsp$trymodel) 

#we have 68 species x site combinations for which they are not considered
#data deficient

# Extract species that are flagged as 'data deficient' from the nsp 
# dataframe
data_deficient_species <- nsp$b.spp[nsp$rare == TRUE]

# the NL GAMMs that do not converge. Put this data in with dd seedlings
fails <- c("Charco_ANACEX",
           "Charco_RINOSY",
           "Oleoducto_ANACEX",
           "Oleoducto_APEIME",
           "Oleoducto_PIPECO",
           "Santa Rita_APHELI",
           "Santa Rita_CHR2CO",
           "Santa Rita_FICUP1",
           "Sherman_ALCHLA",
           "Sherman_MAQUCO")

# Replace species names with "data_deficient_seedling" for those 
# flagged as 'data_deficient'. Be sure to retain BOSQUE information with each 
# DDseedling so that we can account for site by site differences.

dat1 <- dat %>% 
  mutate(b.spp = ifelse(b.spp %in% data_deficient_species, 
                        str_c(BOSQUE, "_", "data_deficient_seedling"), b.spp),
         b.spp = ifelse(b.spp %in% fails,
                        str_c(BOSQUE, "_", "data_deficient_seedling"), b.spp))

#Can data deficient seedlings be assessed for each site?

dat1 %>%
  group_by(b.spp) %>%
  summarise(
    range_con_dens = max(con_dens, na.rm = TRUE) - min(con_dens, 
                                                       na.rm = TRUE),
    max_con_dens = max(con_dens, na.rm = TRUE),
    unique_con_dens = length(unique(con_dens)),
    unique_total_dens = length(unique(tot_dens)),
    unique_height = length(unique(height))
  ) %>%
  mutate(
    # less than nval unique values in consp densities
    issue_nval = unique_con_dens < nval, 
    # range should at least be equal to minrange
    issue_range = range_con_dens < minrange,              
    trymodel = !(issue_nval | issue_range),
    # preliminary assignment of data_deficient species
    data_deficient = !trymodel                                      
  ) -> nsp_data_deficient #no. PP doesn't have enough condens variation

res_mod = list()      # main model fits (results of the model)
res_red_mod = list()  # reduced model fits for Pseudo R2 (results of the reduced model)
lin_mod = list() #linear model fits

# Fit models for individual species (null and non-linear)
for (b.spp in nsp_data_deficient$b.spp[nsp_data_deficient$trymodel]) {
  
  # select data for individual species and sites
  dat_sp = dat1[dat1$b.spp == b.spp, ]
  
  # model fit and reduced fit for Pseudo R2
  mod = model_fit(data = dat_sp,
                  speciesinfo = nsp_data_deficient[nsp_data_deficient$b.spp 
                                                   == b.spp, ])
  
  mod_red = model_fit(data = dat_sp,
                      speciesinfo = nsp_data_deficient[nsp_data_deficient$b.spp 
                                                       == b.spp, ], reduced = T)
  
  # check model success
  res = model_convergence(model = mod)
  res_red = model_convergence(model = mod_red)
  
  # save result
  if (is.character(res)) {
    nsp$rare[nsp$b.spp == b.spp] = T  
  } else {
    res_mod[[b.spp]] = res
    res_red_mod[[b.spp]] = res_red
  }
}


# Fit model for individual species (linear)
for (b.spp in nsp_data_deficient$b.spp[nsp_data_deficient$trymodel]) {
  
  # select data for individual species and sites
  dat_sp = dat1[dat1$b.spp == b.spp, ]
  
  # model fit and reduced fit for Pseudo R2
  mod = linear_model_fit(data = dat_sp,
                  speciesinfo = nsp_data_deficient[nsp_data_deficient$b.spp 
                                                   == b.spp, ])
  
  # check model success
  res = model_convergence(model = mod)
  
  # save result
  if (is.character(res)) {
    nsp$rare[nsp$b.spp == b.spp] = T  
  } else {
    lin_mod[[b.spp]] = res
    
  }
}


save(res_mod, file = "02_Analyses/SpeciesSpecific/res_mod") #nonlinear
save(res_red_mod, file = "02_Analyses/SpeciesSpecific/res_red_mod") #null
save(lin_mod, file = "02_Analyses/SpeciesSpecific/lin_mod") #linear
summary(res_mod[[4]])

## Summarize model fits ---------------

#Regression table via broom::tidy()
coefs = lapply(res_mod, broom::tidy)
coefs = Map(cbind, coefs, b.sp = names(coefs))
coefs = do.call(rbind, coefs)

# df logLik AIC BIC deviance df.residuals nobs 
sums = lapply(res_mod, broom::glance)
sums = Map(cbind, sums, b.sp = names(sums))
sums = do.call(rbind, sums)
#head(sums)


# AUC
aucs = lapply(res_mod, function(x) {
  roc <- performance::performance_roc(x, new_data = x$model)
  bayestestR::area_under_curve(roc$Spec, roc$Sens)
})
sums$AUC = unlist(aucs)


# Pseudo R2
sums$pseudoR2 = 1 - (unlist(lapply(res_mod, function(x) x$deviance)) /
                       unlist(lapply(res_red_mod, function(x) x$deviance)))

# plot splines in pdf ----------------------------------------------

# # Specify the name of the PDF file where the plots will be saved
pdf_file <- "mortality.pdf"

# Open the PDF file for writing
pdf(pdf_file)

# Loop through all the model results stored in 'res_mod'
for (i in 1:length(res_mod)) {
  
  # Get the vizmod for the current species
  vizmod <- getViz(res_mod[[1]], post = T, unconditional = T)
  pl <- plot(vizmod, nsim = 20, allTerms = T) + 
    # Add confidence interval line/ fit line/ simulation line
    l_ciLine() + l_fitLine() + l_simLine() +
    #Add confidence interval bar/# Add fitted points
    l_ciBar() + l_fitPoints(size = 1) +  
    l_rug() +                             # Add rug plot
    # Add title with the name of the current species
    labs(title = names(res_mod)[i])       
  
  # Print the plot to the R console only for the first 2 species
  if (i <= 2) {
    print(pl, pages = 1)
  }
  
  # Save the plot to the PDF
  print(pl, pages = 1)
}

# Close the PDF file
dev.off()

#This continues to give me just one page of output instead of one page per species
