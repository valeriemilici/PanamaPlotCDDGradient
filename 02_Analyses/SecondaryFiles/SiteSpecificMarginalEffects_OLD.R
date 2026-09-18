
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

### NOTES ###
# The current version of this code runs the rel and abs AME through the df
# titled "change2". Change2 sequentially estimates ME across common condens vals
# to revert back to the original ME estimates (i.e., "change"), delete "2" from
# all instances of "change" in the functions for AME and rAME. 13 Nov. '24

#load data
dat <- read.csv("2023data/2023censusdat.csv")

#make sure the data are the appropriate type and named appropriately
dat$census <- factor(dat$CENSO)
dat$spp <- factor(dat$SPP)
dat$plot <- factor(dat$p_05)
dat$height <- as.numeric(dat$ALT_LAST_CENSUS)
dat$interval <- as.numeric(dat$TIME_SINCE_LAST_CENSUS)
dat$status <- dat$mort
dat$b.spp <- str_c(dat$spp, "_", dat$BOSQUE)

#Visualize data ----------------
ggplot(dat, aes(x = con_dens)) +
  geom_histogram(binwidth = 1, color = "black", fill = "steelblue2") +
  labs(x = "Conspecific density", y = "Count", title = "Conspecific densities") +
  theme_bw(12)

ggplot(dat, aes(x = tot_dens)) +
  geom_histogram(binwidth = 1, color = "black", fill = "steelblue2") +
  labs(x = "Total density", y = "Count", title = "Total densities") +
  theme_bw(12)

#Define Data Deficient Species -----------
#I am looking at species within a site, so I need to group by Site and then Spp

nval = 4  ## this is the number of unique conspecific values
minrange = 1    # minimum range for conspecific density

#nval must be at least 4 so that with k = n - 2, we can still allow for a 
#nonlinear relationship, min = k = 2, rather than forcing a linear.

nsp <- dat %>% 
  group_by(spp) %>% 
  summarise(
    range_con_dens = max(con_dens) - min(con_dens),
    max_con_dens = max(con_dens),
    min_con_dens = min(con_dens),
    num_obs = length(con_dens),
    unique_con_dens = length(unique(con_dens)),
    unique_total_dens = length(unique(tot_dens)),
    unique_height = length(unique(ALT_LAST_CENSUS))
  ) %>% 
  
  
  mutate(issue_nval = unique_con_dens < nval,              # less than nval unique values in consp densities
         issue_range = range_con_dens < minrange,                    # range should at least be equal to minrange
         trymodel = !(issue_range),    
         rare = !trymodel # preliminary assignment of rare species
        )

table(nsp$trymodel)
# With this level of filtering, we can model 190/399 species in the dataset

# Function for fitting models -----------
model_fit = function(data,  reduced = F) {
  
  # create new factor with correct factor levels per species 
  data$census = factor(data$census)
  
  
  # create model formula
  term_c = ifelse(length(unique(data$census)) > 1, "+ s(census, bs = 're')", "") 
  term_p = ifelse(length(unique(data$plot)) > 1, "+ s(plot, bs = 're')", "")
  term_s = ifelse(length(unique(data$spp)) > 1, " +s(spp, bs = 're')", "")
  
  if (reduced) {
    form =  as.formula(paste0("status ~ s(height, k = 4) +
                              s(tot_dens, k = 4)"
                              , term_c, term_p, term_s)) # reduced model 
  } else {
    form =  as.formula(paste0("status ~ s(height, k = 4) +
                              s(tot_dens, k = 4)  +
                              s(con_dens, k = 4)" 
                              , term_c, term_p, term_s)) # full model 
  }
  
  # Choose penalty
  # set to default k=10 
 # k1 = k2 = k3 = 10
 # if (k1 > speciesinfo$unique_height) k1 = speciesinfo$unique_height - 2
 # if (k2 > speciesinfo$unique_total_dens) k2 = speciesinfo$unique_total_dens - 2
 # if (k3 > speciesinfo$unique_con_dens) k3 = speciesinfo$unique_con_dens - 2
  # Fit model
  mod = try(gam(form
                , family = binomial(link=cloglog)
                , offset = log(interval)
                , data = data
                , method = "REML"
  ) , silent = T
  )
  
  return(mod)
  
} 

# check model run

model_convergence = function(model) {
  
  # gam not available
  if (!any(class(model)=="gam")) {
    print(paste(BOSQUE, "gam failed"))
  } else {
    
    # gam not converged
    if (!model$converged) {
      print(paste(BOSQUE, "no convergence"))
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
        print(paste(BOSQUE, "complete separation is likely"))
      } else {
        
        # missing Vc
        if (is.null(model$Vc)) {
          print(paste(BOSQUE, "Vc not available"))
        } else {
          
          # successful model
          return(model)
        }
      }
    }
  }
}

## Fit models each species that is not rare ----------
res_mod = list()      # main model fits (results of the model)
res_red_mod = list()  # reduced model fits for Pseudo R2 (results of the reduced model)

# Fit models for individual forest sites

nsp2 <- nsp %>% filter(trymodel == "TRUE") 
splist <- nsp2$spp 
dat_red <- dat %>% filter(spp %in% splist) #reduce data to only include non-rare species
sites <- unique(dat$BOSQUE)

# Fit models for individual sites
for (BOSQUE in sites) {
  
  # select data for individual sites
  dat_s = dat_red[dat_red$BOSQUE == BOSQUE, ]
  
  # model fit and reduced fit for Pseudo R2
  mod = model_fit(data = dat_s) #, speciesinfo = nsp2[nsp2$spp == spp, ]
  mod_red = model_fit(data = dat_s,  reduced = T) #speciesinfo = nsp2[nsp2$spp == spp, ],
  
  # check model success
  res = model_convergence(model = mod)
  res_red = model_convergence(model = mod_red)
  
  # save result
  if (is.character(res)) {
    dat_red$BOSQUE = T  
  } else {
    res_mod[[BOSQUE]] = res
    res_red_mod[[BOSQUE]] = res_red
  }
    
  
}

#Fit models for all species (including rare) (Skip until fixed) ------------
table(rare = nsp$rare, trymodel = nsp$trymodel)

# Convert spp to character
dat$spp <- as.character(dat$spp)

# Create a vector of rare species
rare_species <- nsp$spp[nsp$rare == TRUE]

# Replace spp values to "rare_seedling" where nsp$rare is TRUE
dat2 <- dat %>% 
  mutate(spp = ifelse(spp %in% rare_species, "rare_seedling", spp))

# Fit models for all individual species included the rare group species

#for now this isn't working. I don't understand why. Something about the function
#is causing the dat2 object to convert the "BOSQUE" column into T/F after running,
#but everything looks ok to me. 
for(BOSQUE in sites) {
  
  # select data for individual species
  dat_s = dat2[dat2$BOSQUE == BOSQUE, ] #this might be the breaking point
  
  # model fit and reduced fit for Pseudo R2
  mod = model_fit(data = dat_s)
  mod_red = model_fit(data = dat_s, reduced = T)
  
  # check model success
  res = model_convergence(model = mod)
  res_red = model_convergence(model = mod_red)
  
  # save result
  
  if (is.character(res)) {
    dat2$BOSQUE = T  
  } else {
    res_mod[[BOSQUE]] = res
    res_red_mod[[BOSQUE]] = res_red
  }
  
}

## Summarize model fits ---------------

#Regression table via broom::tidy()
coefs = lapply(res_mod, broom::tidy) 
coefs = Map(cbind, coefs, BOSQUE = names(coefs))
coefs = do.call(rbind, coefs)

# df logLik AIC BIC deviance df.residuals nobs 
sums = lapply(res_mod, broom::glance)
sums = Map(cbind, sums, BOSQUE = names(sums))
sums = do.call(rbind, sums)
#head(sums)

# AUC
aucs = lapply(res_mod, function(x) {
  roc <- performance::performance_roc(x, new_data = x$model)
  bayestestR::area_under_curve(roc$Spec, roc$Sens)
})
sums$AUC = unlist(aucs)


# Pseudo R2
sums$pseudoR2 = 1 - (unlist(lapply(res_mod,
                                   function(x) x$deviance)) /
                       unlist(lapply(res_red_mod, function(x) x$deviance)))

## AMEs Absolute and Relative -------------------

#AMEs relative are used to compare CDD across sites, plots,
#or any other unit of interest with a meta-analyses approach
#while absolute AMEs are reported as they are commonly used in other studies.

#For AMEs we need the file "res_model" that have all the models results.

### Settings for AMEs

#### chose predictors for local density

predictors <- c(con_dens = "con_dens", 
                tot_dens = "tot_dens")

# change in conspecific density for AME calculations----

additive=1     #One more neighbor 

# different change settings for con-specific densities

interval = 1

change = list(equilibrium = data.frame(con_dens = "paste('+', additive)")
              , invasion = data.frame(con_dens = "c(0, additive)")
              , iqr = data.frame(con_dens = "c(q1, q3)")
)

#I'm going to try a multi-metaregression approach that examines changes across
# the most-common densities of seedlings (lwr - upr condens = 0 - 15)

#currently this is a test, so I'm only try every other number 0 -10
change2 = list(zero = data.frame(con_dens = "c(0, additive)"),
               one = data.frame(con_dens = "c(1, additive)"),
               two = data.frame(con_dens = "c(2, additive)"),
              three = data.frame(con_dens = "c(3, additive)"),
               four = data.frame(con_dens = "c(4, additive)"),
               five = data.frame(con_dens = "c(5, additive)"),
              six =  data.frame(con_dens = "c(6, additive)"),
              seven = data.frame(con_dens = "c(7, additive)"))
iter = 500

### Functions to calculate AMEs manually ------
# set step on x for numerical derivative
setstep = function(x) {
  eps = .Machine$double.eps
  x + (max(abs(x), 1, na.rm = TRUE) * sqrt(eps)) - x
}


# Average marginal effects ------------------------------------------------

get_AME = function(mod, data, term
                   , change2 = NULL
                   , at = NULL
                   , offset = 1
                   , relative = F
                   , iterations = 1000
                   , seed = 10
                   , samples = F) {
  
  #  Creating two copies of the data frames one of which will shift in term of interest:
  d0 = d1 = data
  
  # Adjusting the term of interests with the different options of change
  
  # if change is NULL, use numerical derivative
  if (is.null(change2)) {
    
    d0[[term]] = d0[[term]] - setstep(d0[[term]])
    d1[[term]] = d1[[term]] + setstep(d1[[term]])
    
  } 
  
  # if change includes "+", use explicit additive change
  if (grepl("\\+", paste(change2, collapse = "_"))) {
    
    d1[[term]] = d1[[term]] + as.numeric(gsub("\\+", "", change2))
    
  } 
  
  # if change has two values, use explicit change
  if (length(change2) == 2) {
    
    d0[[term]] = as.numeric(change2[1])
    d1[[term]] = as.numeric(change2[2])
    
  }
  
  # Fix at values: (allows the function to calculate the marginal effects at the specified values)
  if (!is.null(at)) {
    for (i in names(at))
      d0[[i]] = at[[i]]
    d1[[i]] = at[[i]]
  }
  
  # Matrices for prediction, map coefs to fitted curves
  Xp0 <- predict(mod, newdata = d0, type="lpmatrix")
  Xp1 <- predict(mod, newdata = d1, type="lpmatrix")
  
  # Model settings
  ilink <- family(mod)$linkinv
  beta <- coef(mod)
  vc <- mod$Vc # covariance matrix 
  
  
  # marginal effects
  pred0   <- 1 - (1-ilink(Xp0 %*% beta))^offset
  pred1   <- 1 - (1-ilink(Xp1 %*% beta))^offset
  ME <- (pred1-pred0)
  
  # if change is NULL, use numerical derivative
  if (is.null(change2)) {
    ME <- ME/(d1[[term]] - d0[[term]])
  } 
  
  
  # convert to relative if requested
  if (relative == T) ME = ME/pred0
  
  # average marginal effect
  AME = mean(ME)
  
  
  # generate AME samples for calculating uncertainty in the AME estimates
  
  # variance of average marginal effect via "posterior" simulation
  # simulate from multivariate normal using model beta means and covariance matrix
  if (!is.null(seed)) set.seed(seed)
  coefmat = mvrnorm(n = iterations
                    , mu = beta
                    , Sigma = vc)
  
  # estimate AME from each simulated coefficient vector
  AMEs = apply(coefmat, 1, function(coefrow) {
    
    # marginal effects
    pred0   <- 1 - (1-ilink(Xp0 %*% coefrow))^offset
    pred1   <- 1 - (1-ilink(Xp1 %*% coefrow))^offset
    ME <- (pred1-pred0)
    
    # if change is NULL, use numerical derivative
    if (is.null(change2)) {
      ME <- ME/(d1[[term]] - d0[[term]])
    } 
    
    # convert to relative if requested
    if (relative == T) ME = ME/pred0
    
    # average marginal effect
    AME = mean(ME)
    return(AME)
  })
  
  # Combine results
  if (!samples) {
    res = data.frame(term
                     , estimate = AME
                     , std.error = sqrt(var(AMEs))  
                     , estimate.sim = mean(AMEs)    
                     , offset
                     , change.value = paste(change2, collapse = "_"))
    return(res) 
    
  } else {
    
    res_sums = data.frame(term
                          , estimate = AME
                          , std.error = sqrt(var(AMEs)) 
                          , offset
                          , change.value = paste(change2, collapse = "_"))
    
    res_samples = data.frame(term
                             , estimate = AMEs
                             , MLE = AME
                             , offset
                             , change.value = paste(change2, collapse = "_"))
    res = list(res_sums, res_samples)
    return(res)  
    
  }
}
# Absolute AMEs -----------------------------------------------------------

# Calculate absolute AMEs based on manual function get_AME
AME = data.frame()
AMEsamples = data.frame()
for (i in names(predictors)[grepl("con_", names(predictors))]) { 
  for (j in names(change2)) {
    temp = lapply(res_mod, function(x){
      if (j == "iqr") {
        q1 = quantile(x$model$con_dens, probs = 0.25)
        q3 = quantile(x$model$con_dens, probs = 0.75)
      }
      get_AME(x
              , data = x$model
              , offset = interval
              , term = i
              , change2 = eval(parse(text = change2[[j]][,i]))
              , iterations = iter
              , samples = T
      )
    }
    )
    
    # AME
    tempAME = lapply(temp, function(x) x[[1]])
    tempAME = Map(cbind, tempAME, change2 = j, BOSQUE = names(tempAME))
    tempAME = do.call(rbind, tempAME)
    AME = rbind(AME, tempAME)
    
    # AME samples
    tempSamples = lapply(temp, function(x) x[[2]])
    tempSamples = Map(cbind, tempSamples, change2 = j, BOSQUE = names(tempSamples), iter = iter)
    tempSamples = do.call(rbind, tempSamples)
    AMEsamples = rbind(AMEsamples, tempSamples)
  }
}
head(AME) 

# Relative AMEs -----------------------------------------------------------


# Calculate relative AMEs based on manual function get_AME
rAME = data.frame()
rAMEsamples = data.frame()
for (i in names(predictors)[grepl("con_", names(predictors))]) { # for all predictors
  for (j in names(change2)) {
    temp = lapply(res_mod, function(x){
      if (j == "iqr") {
        q1 = quantile(x$model$con_dens, probs = 0.25)
        q3 = quantile(x$model$con_dens, probs = 0.75)
      }
      get_AME(x
              , data = x$model
              , offset = interval
              , term = i
              , change2 = eval(parse(text = change2[[j]][, i]))
              , iterations = iter
              , relative = T
              , samples = T
      )
    }
    )
    
    # rAME
    tempAME = lapply(temp, function(x) x[[1]])
    tempAME = Map(cbind, tempAME, change2 = j, BOSQUE = names(tempAME))
    tempAME = do.call(rbind, tempAME)
    rAME = rbind(rAME, tempAME)
    
    # rAME samples
    tempSamples = lapply(temp, function(x) x[[2]])
    tempSamples = Map(cbind, tempSamples, change2 = j, BOSQUE = names(tempSamples), iter = iter)
    tempSamples = do.call(rbind, tempSamples)
    rAMEsamples = rbind(rAMEsamples, tempSamples)
  }
}
head(rAME)

# plot splines in pdf -----------------------------------------------------

# Set the PDF file name
pdf_file <- "mortality.pdf" # save as pdf all the species figures

for (i in 1:length(res_mod))  {
  
  # with gamViz
  vizmod <- getViz(res_mod[[i]], post = T, unconditional = T)
  pl = plot(vizmod, nsim = 20, allTerms = T) + 
    l_ciLine() + l_fitLine() + l_simLine() + 
    l_ciBar() + l_fitPoints(size = 1) + 
    l_rug() + 
    labs(title = names(res_mod)[i]) 
  
  # Print the plot to the R console only for the first 3 species
  if (i <= 4) {
    print(pl, pages = 1)
  }
  
  # Save the plot to the PDF
  pdf(pdf_file)
  print(pl, pages = 1)
  dev.off()
}

#the PDF function isn't working and is only generating a 1 page PDF. I haven't
#figured out how to fix it (I haven't tried very hard), but this seems
#non-essential so move onto the next step. 

# Save results ------------------------------------------------------------


save(list = c("AME", "AMEsamples", "rAME", "rAMEsamples", "nsp", "coefs", "sums") # "nsp_rare"
     , file = paste0( "(r)AMEOutput/mortality_sitespecific.Rdata"))

#load("(r)AMEOutput/mortality.Rdata")
write.csv(AME, "(r)AMEOutput/AME_sitespecific.csv")
write.csv(rAME, "(r)AMEOutput/rAME_sitespecific.csv")
