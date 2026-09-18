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

load("02_Analyses/SpeciesSpecific/res_mod") #nonlinear results
load("02_Analyses/SpeciesSpecific/lin_mod") #linear model results


## AMEs Absolute and Relative -------------------

#AMEs relative are used to compare CDD across sites, plots,
#or any other unit of interest with a meta-analyses approach
#while absolute AMEs are reported as they are commonly used in other studies.

#For AMEs we need the file "res_mod" and "lin_mod" that have all the models results.

### Settings for AMEs ----

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

# Nonlinear AMEs -------
AME_nl = data.frame()
AMEsamples_nl = data.frame()
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
    AME_nl = rbind(AME_nl, tempAME)
    
    # AME samples
    tempSamples = lapply(temp, function(x) x[[2]])
    tempSamples = Map(cbind, tempSamples, change = j, b.sp = names(tempSamples), iter = iter)
    tempSamples = do.call(rbind, tempSamples)
    AMEsamples_nl = rbind(AMEsamples_nl, tempSamples)
  }
}

head(AME_nl) 

# Linear AMEs -----------
# Note: This is redundant. If a linear form is the best form, the GAM will have
# modeled the relationship as linear. Bookmark to remove this step.
AME_l = data.frame()
AMEsamples_l = data.frame()
for (i in names(predictors)[grepl("con_", names(predictors))]) { 
  for (j in names(change)) {
    temp = lapply(lin_mod, function(x){
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
    AME_l = rbind(AME_l, tempAME)
    
    # AME samples
    tempSamples = lapply(temp, function(x) x[[2]])
    tempSamples = Map(cbind, tempSamples, change = j, b.sp = names(tempSamples), iter = iter)
    tempSamples = do.call(rbind, tempSamples)
    AMEsamples_l = rbind(AMEsamples_l, tempSamples)
  }
}

head(AME_l)


# Relative AMEs -----------------------------------------------------------

# Nonlinear rAMEs ---------------
rAME_nl = data.frame()
rAMEsamples_nl = data.frame()
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
    rAME_nl = rbind(rAME_nl, tempAME)
    
    # rAME samples
    tempSamples = lapply(temp, function(x) x[[2]])
    tempSamples = Map(cbind, tempSamples, change = j, b.sp = names(tempSamples), iter = iter)
    tempSamples = do.call(rbind, tempSamples)
    rAMEsamples_nl = rbind(rAMEsamples_nl, tempSamples)
  }
}
head(rAME_nl)

# Linear rAMEs --------------------
# same comment as with Linear AME. Bookmark to remove. 
rAME_l = data.frame()
rAMEsamples_l = data.frame()
for (i in names(predictors)[grepl("con_", names(predictors))]) { # for all predictors
  for (j in names(change)) {
    temp = lapply(lin_mod, function(x){
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
    rAME_l = rbind(rAME_l, tempAME)
    
    # rAME samples
    tempSamples = lapply(temp, function(x) x[[2]])
    tempSamples = Map(cbind, tempSamples, change = j, b.sp = names(tempSamples), iter = iter)
    tempSamples = do.call(rbind, tempSamples)
    rAMEsamples_l = rbind(rAMEsamples_l, tempSamples)
  }
}
head(rAME_l)

# Save results ------------------------------------------------------------


#save(list = c("AME_nl", "AME_l", "AMEsamples_nl", "AMEsamples_l", "rAME_nl",
  #            "rAME_l", "rAMEsamples_nl", "rAMEsamples_l") # "nsp_rare", "nsp", "coefs", "sums"
 #    , file = paste0( "02_Analyses/(r)AMEOutput/mortality.Rdata"))

#load("(r)AMEOutput/mortality.Rdata")
write.csv(AME_nl, "02_Analyses/(r)AMEOutput/AMEmort_nl.csv")
write.csv(rAME_nl, "02_Analyses/(r)AMEOutput/rAMEmortnl.csv")
#write.csv(AME_l, "02_Analyses/(r)AMEOutput/AMEmort_l.csv")
#write.csv(rAME_l, "02_Analyses/(r)AMEOutput/rAMEmortl.csv")
