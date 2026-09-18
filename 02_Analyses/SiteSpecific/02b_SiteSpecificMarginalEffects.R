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
dat <- read.csv("2023data/2023censusdat.csv")

#make sure the data are the appropriate type and named appropriately
dat$census <- factor(dat$CENSO)
dat$spp <- factor(dat$SPP)
dat$plot <- factor(dat$p_05)
dat$height <- as.numeric(dat$ALT_LAST_CENSUS)
dat$interval <- as.numeric(dat$TIME_SINCE_LAST_CENSUS)
dat$status <- as.numeric(dat$mort)
dat$con_dens.p <- as.numeric(dat$con_dens.p)
dat$tot_dens.p <- as.numeric(dat$tot_dens.p)
dat$con_dens <- as.numeric(dat$con_dens)
dat$tot_dens <- as.numeric(dat$tot_dens)
dat$BOSQUE <- factor(dat$BOSQUE)

# We can't analyze data from the first census, and seedlings must be
# within size threshold for the category (200 - 1200 mm)

dat1 <- dat %>% filter(FIRST_CENSUS != 1) %>%
  filter(!is.na(height)) %>%
  filter(height >= 200) %>%
  filter(height <= 1200) %>%
  mutate(b.spp = paste(BOSQUE, SPP, sep = "_"))

# are there any species over-represented in the dataset?

histdat <- dat1 %>% group_by(ID) %>% slice(1) %>% ungroup() %>%
  group_by(BOSQUE, SPP) %>% tally() 

outliers <- histdat %>% ungroup() %>% group_by(BOSQUE) %>%
  summarize(outlier = quantile(n, 0.975))

outs1 <- histdat %>% 
  filter(
    case_when(BOSQUE == "ANAM" & n >= 102 ~ T,
           BOSQUE == "Buena Vista" & n >= 73 ~ T,
           BOSQUE == "Charco" & n >= 36 ~ T,
           BOSQUE == "Metropolitano" & n >= 286 ~ T,
           BOSQUE == "Oleoducto" & n >= 52 ~ T,
           BOSQUE == "Panama Pacifico" & n >= 121 ~ T,
           BOSQUE == "Santa Rita" & n >= 70 ~ T,
           BOSQUE == "Sherman" & n >= 70 ~ T)) %>%
  mutate(b.spp = paste(BOSQUE, SPP, sep = "_"))

outs <- outs1$b.spp

#remove conspecific and total density observations if there is a >2 n gap in observations.
#the GAM doesn't estimate well across gaps in observations
cdgaps <- dat1  %>% filter(!b.spp %in% outs) %>%
  group_by(BOSQUE, con_dens) %>% tally() %>%
  mutate(dif = con_dens - lag(con_dens)) %>%
  filter(dif > 2) %>% 
  ungroup() %>%
  group_by(BOSQUE) %>%
  slice(1) %>%
  dplyr::select(1:2)

tdgaps <- dat1  %>% filter(!b.spp %in% outs) %>%
  group_by(BOSQUE, tot_dens) %>% tally() %>%
  mutate(dif = tot_dens - lag(tot_dens)) %>%
  filter(dif > 2) %>% 
  ungroup() %>%
  group_by(BOSQUE) %>%
  slice(1) %>%
  dplyr::select(1:2) 

#remove all outliers (high site-specific species observations and large gaps
dat2 <- dat1 %>% filter(!b.spp %in% outs) %>%
  filter(
    case_when(
      BOSQUE == "ANAM" & con_dens < 10 & tot_dens < 36 ~ T,
      BOSQUE == "Buena Vista" & con_dens < 11 ~ T,
      BOSQUE == "Santa Rita" & con_dens < 15 & tot_dens < 40 ~ T,
      BOSQUE == "Sherman" & con_dens < 10 ~ T,
      BOSQUE == "Metropolitano" & tot_dens < 24 ~ T,
      BOSQUE == "Charco" ~ T,
      BOSQUE == "Oleoducto" ~ T,
      BOSQUE == "Panama Pacifico" ~ T) )

#confirm there are no gaps
tdgaps <- dat2   %>%
  group_by(BOSQUE, tot_dens) %>% tally() %>%
  mutate(dif = tot_dens - lag(tot_dens)) %>%
  filter(dif > 2) %>% 
  ungroup() %>%
  group_by(BOSQUE) %>%
  slice(1) %>%
  dplyr::select(1:2) 

cdgaps <- dat2 %>%
  group_by(BOSQUE, con_dens) %>% tally() %>%
  mutate(dif = con_dens - lag(con_dens)) %>%
  filter(dif > 2) %>% 
  ungroup() %>%
  group_by(BOSQUE) %>%
  slice(1) %>%
  dplyr::select(1:2) #good!

siteinfo <- 
  dat2 %>% 
  group_by(BOSQUE) %>% 
  summarise(
    range_con_dens = max(con_dens) - min(con_dens),
    cd_q975 = quantile(con_dens, 0.975),
    max_con_dens = max(con_dens),
    min_con_dens = min(con_dens),
    num_obs = length(con_dens),
    unique_con_dens = length(unique(con_dens.p)),
    unique_total_dens = length(unique(tot_dens.p)),
    td_q975 = quantile(tot_dens, 0.975),
    unique_height = length(unique(ALT_LAST_CENSUS))
  )

# Function for fitting models -----------
model_fit = function(data, bosqueinfo,  reduced = F) {
  
  # create new factor with correct factor levels per species 
  data$census = factor(data$census)
  
  
  # create model formula
  term_c = ifelse(length(unique(data$census)) > 1, "+ s(census, bs = 're')", "") 
  term_p = ifelse(length(unique(data$plot)) > 1, "+ s(plot, bs = 're')", "")
  term_s = ifelse(length(unique(data$spp)) > 1, " +s(spp, bs = 're')", "")
  
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
  if (k1 > bosqueinfo$unique_height) k1 = bosque$unique_height - 2
  if (k2 > bosqueinfo$unique_total_dens) k2 = bosqueinfo$unique_total_dens - 2
  if (k3 > bosqueinfo$unique_con_dens) k3 = bosqueinfo$unique_con_dens - 2
  
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

# Fit models for individual sites

sites <- unique(dat$BOSQUE)

for (BOSQUE in sites) {
  
  # select data for individual sites
  dat_s = dat[dat$BOSQUE == BOSQUE, ]
  
  # model fit and reduced fit for Pseudo R2
  mod = model_fit(data = dat_s,
                  bosqueinfo = siteinfo[siteinfo$BOSQUE == BOSQUE, ]) 
  mod_red = model_fit(data = dat_s,
                      bosqueinfo = siteinfo[siteinfo$BOSQUE == BOSQUE, ],
                      reduced = T) 
  
  # check model success
  res = model_convergence(model = mod)
  res_red = model_convergence(model = mod_red)
  
  # save result
  if (is.character(res)) {
    dat$BOSQUE = T  
  } else {
    res_mod[[BOSQUE]] = res
    res_red_mod[[BOSQUE]] = res_red
  }
  
  
}

## Check the GAM fits

gam.check(res_mod[[2]])
summary(res_mod[[8]])
#ANAM has significant edf, but edf is far from k' so I'm not sure what the move is
#BV is fine p-value wise. The plots all look messy. Not sure if binomial cloglog affects expectations
#Charco is fine
#metropolitano has significant edf but k' is far from estimate
#Oleoducto tot_dens may be too low (same as ANAM)
#PP is fine
#SR is fine
#Sherman may be too low

#save model output to plot the splines
save(res_mod, file = "02_Analyses/SiteSpecific/res_mod")
save(res_red_mod, file = "02_Analyses/SiteSpecific/res_red_mod")
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


## AMEs Absolute and Relative -------------------

#AMEs relative are used to compare CDD across sites, plots,
#or any other unit of interest with a meta-analyses approach
#while absolute AMEs are reported as they are commonly used in other studies.

#For AMEs we need the file "res_model" that have all the models results.

### Settings for AMEs

#### chose predictors for local density

predictors <- c(con_dens.p = "con_dens.p", 
                tot_dens.p = "tot_dens.p")

# change in conspecific density for AME calculations----

additive=1     #One more neighbor 

# different change settings for con-specific densities

interval = 1

change = list(equilibrium = data.frame(con_dens.p = "paste('+', additive)")
              , invasion = data.frame(con_dens.p = "c(0, additive)")
              , iqr = data.frame(con_dens.p = "c(q1, q3)")
)

iter = 500

### Functions to calculate AMEs manually ------
# set step on x for numerical derivative
setstep = function(x) {
  eps = .Machine$double.eps
  x + (max(abs(x), 1, na.rm = TRUE) * sqrt(eps)) - x
}


# Average marginal effects ------------------------------------------------

get_AME = function(mod, data, term
                   , change = NULL
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
  if (is.null(change)) {
    
    d0[[term]] = d0[[term]] - setstep(d0[[term]])
    d1[[term]] = d1[[term]] + setstep(d1[[term]])
    
  } 
  
  # if change includes "+", use explicit additive change
  if (grepl("\\+", paste(change, collapse = "_"))) {
    
    d1[[term]] = d1[[term]] + as.numeric(gsub("\\+", "", change))
    
  } 
  
  # if change has two values, use explicit change
  if (length(change) == 2) {
    
    d0[[term]] = as.numeric(change[1])
    d1[[term]] = as.numeric(change[2])
    
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
  if (is.null(change)) {
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
    if (is.null(change)) {
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
                     , change.value = paste(change, collapse = "_"))
    return(res) 
    
  } else {
    
    res_sums = data.frame(term
                          , estimate = AME
                          , std.error = sqrt(var(AMEs)) 
                          , offset
                          , change.value = paste(change, collapse = "_"))
    
    res_samples = data.frame(term
                             , estimate = AMEs
                             , MLE = AME
                             , offset
                             , change.value = paste(change, collapse = "_"))
    res = list(res_sums, res_samples)
    return(res)  
    
  }
}

# Absolute AMEs -----------------------------------------------------------

# Calculate absolute AMEs based on manual function get_AME
AME = data.frame()
AMEsamples = data.frame()
for (i in names(predictors)[grepl("con_", names(predictors))]) { 
  for (j in names(change)) {
    temp = lapply(res_mod, function(x){
      if (j == "iqr") {
        q1 = quantile(x$model$con_dens.p, probs = 0.25)
        q3 = quantile(x$model$con_dens.p, probs = 0.75)
      }
      get_AME(x
              , data = x$model
              , offset = interval
              , term = i
              , change = eval(parse(text = change[[j]][,i]))
              , iterations = iter
              , samples = T
      )
    }
    )
    
    # AME
    tempAME = lapply(temp, function(x) x[[1]])
    tempAME = Map(cbind, tempAME, change = j, b.sp = names(tempAME))
    tempAME = do.call(rbind, tempAME)
    AME = rbind(AME, tempAME)
    
    # AME samples
    tempSamples = lapply(temp, function(x) x[[2]])
    tempSamples = Map(cbind, tempSamples, change = j, b.sp = names(tempSamples), iter = iter)
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
  for (j in names(change)) {
    temp = lapply(res_mod, function(x){
      if (j == "iqr") {
        q1 = quantile(x$model$con_dens.p, probs = 0.25)
        q3 = quantile(x$model$con_dens.p, probs = 0.75)
      }
      get_AME(x
              , data = x$model
              , offset = interval
              , term = i
              , change = eval(parse(text = change[[j]][, i]))
              , iterations = iter
              , relative = T
              , samples = T
      )
    }
    )
    
    # rAME
    tempAME = lapply(temp, function(x) x[[1]])
    tempAME = Map(cbind, tempAME, change = j, b.sp = names(tempAME))
    tempAME = do.call(rbind, tempAME)
    rAME = rbind(rAME, tempAME)
    
    # rAME samples
    tempSamples = lapply(temp, function(x) x[[2]])
    tempSamples = Map(cbind, tempSamples, change = j, b.sp = names(tempSamples), iter = iter)
    tempSamples = do.call(rbind, tempSamples)
    rAMEsamples = rbind(rAMEsamples, tempSamples)
  }
}

head(rAME)


#save output
save(list = c("AME", "AMEsamples", "rAME", "rAMEsamples", "coefs", "sums") # "nsp_rare"
     , file = paste0( "02_Analyses/(r)AMEOutput/mortality_sitespecific.Rdata"))

write.csv(AME, "02_Analyses/(r)AMEOutput/AME_sitespecific.csv")
write.csv(rAME, "02_Analyses/(r)AMEOutput/rAME_sitespecific.csv")
