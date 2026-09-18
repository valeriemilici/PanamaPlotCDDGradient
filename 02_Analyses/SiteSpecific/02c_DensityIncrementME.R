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

# load gam results
load("02_Analyses/SiteSpecific/res_mod")
load("02_Analyses/SiteSpecific/res_red_mod")

## AMEs Absolute and Relative -------------------

#### chose predictors for local density

predictors <- c(con_dens = "con_dens", 
                tot_dens = "tot_dens")

change = list(one = data.frame(con_dens = "c(0, 1)")
              , four = data.frame(con_dens = "c(3, 4)")
              , eight = data.frame(con_dens = "c(7, 8)")
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
        q1 = quantile(x$model$con_dens, probs = 0.25)
        q3 = quantile(x$model$con_dens, probs = 0.75)
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
        q1 = quantile(x$model$con_dens, probs = 0.25)
        q3 = quantile(x$model$con_dens, probs = 0.75)
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
save(list = c("AME", "rAME")
     , file = paste0( "02_Analyses/(r)AMEOutput/mortality_increment_sitespecific.Rdata"))

