
################################################################################
### Example cross validation
################################################################################

## This script steps through carrying out cluster cross validation to compare 
## the predictive performance of the different models, using the growth vital 
## rate as an example

rm(list = ls(all = TRUE))

set.seed(110189)

source("Configure.R")

################################################################################
# Load and format data
################################################################################

## Load demographic data
load(file.path(datawd, "DemographicData.rda"))

## Load observed climate data
load(file.path(datawd, "ObservedClimate.rda"))

################################################################################
# Useful functions
################################################################################

## Format demographic data
growdata <- allsites %>% dplyr::select(bald, obsY, Sizem1, Size, tsf) %>% 
  filter(!is.na(Size), !is.na(Sizem1), !is.na(tsf))

## Function to fit the different possible models 
crossval <- function(year, climvar, demodata, climdata, fire, fullmod=NULL, meth = "REML", sknots = 7, 
                     fknots = 7, climknots = 8, nit=1000){
  ## Base model with size only (no time since fire or climate effects)
  if(is.na(fire)){
    usedata <- filter(demodata, !obsY==year) %>% mutate(yrpop = factor(paste0(bald, obsY)), obsY=factor(obsY), bald = factor(bald))
    model <- gam(Size ~ s(Sizem1, bs = "cr", k=sknots) + s(bald, bs="re") + s(yrpop, bs="re"), 
                 data=usedata, family="gaussian", method=meth)
    alldata <- demodata
  } else {
    ## Model with linear fire effect
    if(grepl("lin", fire)){
      if(is.na(climvar)){
        usedata <- filter(demodata, !obsY==year) %>% mutate(yrpop = factor(paste0(bald, obsY)), obsY=factor(obsY), bald = factor(bald))
        model <- gam(Size ~ s(Sizem1, bs = "cr", k=sknots)  + tsf + s(bald, bs="re") + s(yrpop, bs="re"), 
                     data=usedata, family="gaussian", method=meth)
        alldata <- demodata
      } else {
        ## Model with fire spline
        climdata <- climdata %>% dplyr::select(EryYear, EryFort, get(paste0(climvar, "_c"))) %>% tidyr::spread_("EryFort", paste0(climvar, "_c"))
        names(climdata) <- c("EryYear", paste0("P", 1:26))
        alldata <- inner_join(demodata, climdata, by=c("obsY" = "EryYear")) %>% 
          mutate(yrpop = factor(paste0(bald, obsY)), obsY=factor(obsY), bald = factor(bald))
        alldata$climate <- as.matrix(alldata[,paste0("P", 1:26)])
        alldata$fort <- matrix(rep(1:26, each=nrow(alldata)), nrow=nrow(alldata))
        usedata <- filter(alldata, !obsY==year)
        model <- gam(Size ~ s(Sizem1, bs = "cr", k=sknots)  + tsf + s(fort, bs="cr", by=climate, k=climknots) + s(bald, bs="re") + s(yrpop, bs="re") , 
                     data=usedata, family="gaussian", method=meth)  
      }} else {
        if(is.na(climvar)){
          usedata <- filter(demodata, !obsY==year) %>% mutate(yrpop = factor(paste0(bald, obsY)), obsY=factor(obsY), bald = factor(bald))
          model <- gam(Size ~ s(Sizem1, bs = "cr", k=sknots)  + s(tsf, bs="cr", k=fknots) + s(bald, bs="re") + s(yrpop, bs="re"), 
                       data=usedata, family="gaussian", method=meth)
          alldata <- demodata
        } else{
          ## Model with climate effect
          climdata <- climdata %>% dplyr::select_("EryYear", "EryFort", paste0(climvar, "_c")) %>% tidyr::spread_("EryFort", paste0(climvar, "_c"))
          names(climdata) <- c("EryYear", paste0("P", 1:26))
          alldata <- inner_join(demodata, climdata, by=c("obsY" = "EryYear")) %>% 
            mutate(yrpop = factor(paste0(bald, obsY)), obsY=factor(obsY), bald = factor(bald))
          alldata$climate <- as.matrix(alldata[,paste0("P", 1:26)])
          alldata$fort <- matrix(rep(1:26, each=nrow(alldata)), nrow=nrow(alldata))
          usedata <- filter(alldata, !obsY==year)
          model <- gam(Size ~ s(Sizem1, bs = "cr", k=sknots) + s(tsf, bs="cr", k=fknots) + s(fort, bs="cr", by=climate, k=climknots) + s(bald, bs="re") + s(yrpop, bs="re") , 
                       data=usedata, family="gaussian", method=meth)
        }}}
  ## Use gam for predictions
  predata <- filter(alldata, obsY==year) %>% mutate(yrpop = usedata$yrpop[1]) 
  ## Have to make "yrpop" something that has been seen in the data used to fit the model here - this isn't actually used in the predictions
  pred <- apply(predict(model, predata, type="terms", exclude="s(yrpop)"), 1, sum) + model$coefficients["(Intercept)"] ## model$coefficients[1] is the intercept
  ## Randomly sample random year effects 
  pred <- rep(pred, each=nit) + rnorm(nit, mean=0, sd=gam.vcomp(model)["s(yrpop)", "std.dev"])
  ## Calculate RMSE
  ind <- rep(1:nrow(predata), each=nit)
  like <- data.table(nit =rep(1:nit, times=nrow(predata)), rmse = (predata$Size[ind] - pred)^2)
  like <- like[, .(rmse = mean(rmse)), by=c("nit")][,.(rmse = sqrt(mean(rmse))),][, clim:=climvar,][,fire:=fire,]
  print(c(year, climvar, fire))
  return(like)
}

## Faster version of dnorm function (but with no checking)
fastdnorm <- function(values, means, stdevs){
  dem1 <- 2 * stdevs * stdevs
  dem2 <- sqrt(2*pi) * stdevs
  out <- exp(-((values - means)^2/dem1))/dem2
  return(out)
}

## Function to summarise cross validation results for each model
cvresults <- function(cv, demdata){
  nobs <- (demdata %>% group_by(obsY) %>% arrange(obsY) %>% summarise(n = n()))$n
  climmod <- cv %>% group_by(clim, fire) %>% mutate(weight = nobs) %>% 
    summarise(totrmse = sum(rmse*weight))
  return(climmod)
}
# Weighting to allow for different number of observations in the different years of the study

###########################################################################################
# Comparing time since fire models - using GAM, linear effect or no time since fire effect
###########################################################################################

param <- expand.grid(fire = c("6knots", "linear", NA), yrs = 1991:2014)
# Set comparisons - for time since fire we compared the GAM to a linear TSF effect and a base model with no TSF effect

# Run the cross validation
growfirecv <- do.call(rbind, mapply(crossval, year= param$yrs, fire = param$fire, 
                                     MoreArgs = list(climdata = formclim, climvar = NA, demodata = growdata, nit = 10),
                                    SIMPLIFY = FALSE))
# NB: nit sets the number of samples from the random year effect distribution - set to 10 here for due to computation time, 
# but 1000 samples used in the manuscript

# Summary of results
cvresults(growfirecv, growdata)

###########################################################################################
# Comparing the climate models
# Compare base model with tsf spline, but no climate effect to those including FLM for each climate variable
###########################################################################################

fire <- rep("6knots", times = 24*4)
clim <- rep(c("drought", "mintemp", "maxtemp", "precip"), each = 24)
yrs <- rep(1991:2014, length = length(fire))
# Set comparisons - for climate we compared the base model (with the TSF GAM) and no climate effects to climate models 
# containing each of the four potential climatic variables

# Run the cross validation
growclimcv <- do.call(rbind, mapply(crossval, year=yrs, fire=fire, climvar = clim,
                                MoreArgs = list(climdata = formclim, demodata = growdata, nit = 10), SIMPLIFY = FALSE))
# NB: nit sets the number of samples from the random year effect distribution - set to 10 here for due to computation time, 
# but 1000 samples used in the manuscript

# Summary of results
cvresults(growclimcv, growdata)

## Save output for future use
# allcv <- bind_rows(growfirecv, growclimcv)
# write.csv(allcv, file.path(dataout, "GrowCV.csv"), row.names = FALSE)

