
################################################################################
### Example IPM simulation
################################################################################

rm(list=ls(all=T))

source("Configure.R")

######################################################################################################
### Source IPM code and load GAMs & climate data
######################################################################################################

### Source IPM functions
source(file.path(codewd, "IPMfunctions.R"))

## This code steps through simulating data under present climatic conditions
## Code for future climatic conditions is the same apart from the data used and that 
## here the climate data were taken each year in turn (rather than sampling randomly)
## See enveffect function in IPMfunctions.R

### Load climate splines
load(file.path(dataout, "ClimModels_splines_sqrt.rda"))
# This is a list containing the models for the four vital rates (growth, survival, fecundity and recruit size)
# For each vital rate the climate predictor with the best predictive performance, using cross validation, is included (see Table 1 in manuscript)

### Load observed climate data and filter to contain the relevant years
load(file.path(datawd, "ObservedClimate.rda"))
  formclim <- ungroup(formclim) %>% filter(EryYear %in% 1990:2014)

## Predicted climate data 
load(file.path(datawd, "PredictedClimate_Present.rda"))  

## Remove means of observed data from predicted data as observed data were centered in the vital rate models
centclim <- ungroup(formclim) %>% group_by(EryFort) %>% 
  summarise(mintemp = mean(mintemp), maxtemp = mean(maxtemp), precip = mean(precip), drought = mean(drought)) %>%
  right_join(ungroup(pclim)) %>% mutate(mintemp_c = mintds - mintemp, maxtemp_c = maxtds - maxtemp, 
                               precip_c = precds - precip, drought_c = dr - drought)

######################################################################################################
### Simulate population dynamics
######################################################################################################

# Set implementation parameters for IPM - m is number of meshpoints, L and U are the lower and upper limits respectively
i.par <- calc_i_par(m = 100, L =0.1, U =7)

# Load random effects - this is a list containing the populations and year/population combinations from the raw data
load(file.path(datawd, "CommonRand.rda"))

## Set the different scenarios to simulate

nsim <- 5 # Set number of replicate simulations
# NB: this is set to 5 here to decrease computation time
# In the manuscript we run 1000 simulations for each parameter combination

# Set fertility parameters
fertscen <- data.frame(fgerm = 0, sbgerm = 0.005, mort = 0.3, sim = 1:nsim) 
# fgerm is germination in first year, sbgerm is germination from the seedbank and mort is mortality within seedbank

# Set fire return interval parameters
firescen <- expand.grid(firemed = seq(3, 30, by = 3), fireshape = c(2, 8, 32, 64), sim=1:nsim)
# firemed is the median RFI, fireshape is the shape of the weibull distribution
# nsim is the number of simulations to run, set above

allscen <- inner_join(fertscen, firescen)

##################################################################################################
### Simulate in parallel
##################################################################################################

## Set up parallel environment - if running in parallel
# np <- mpi.universe.size()
# np
# cl <- makeMPIcluster(np-1) ## For running on hpc
# cl <- makeCluster(2, type="SOCK") ## For running on laptop
# registerDoParallel(cl)

## Run models 
#popsize <- foreach(fgerm = allscen$fgerm, sbgerm = allscen$sbgerm, mort = allscen$mort, sim = allscen$sim,
#                   firemed = allscen$firemed, fireshape = allscen$fireshape, .combine = "rbind", .packages=c("dplyr", "mgcv")) %dopar%
#  run.ipm(fgerm = fgerm, sbgerm = sbgerm, mort = mort, pop =NA, sim = sim, firemed = firemed, fireshape = fireshape, 
#          num.years=85, i.par = i.par, rates = vr, init = 7000, yrpoprand = randef, climate = centclim)

# num.years is the length of the simulation 
# init is the initial population size (number of seeds)
# pop can be set to simulate a specific population. Set to NA, this randomly samples from the possible populations

### stop cluster - if running in parallel
# stopCluster(cl); print("Cluster stopped.")
# mpi.exit() 

##################################################################################################
### Simulate not in parallel
##################################################################################################

popsize <- do.call(rbind, mapply(run.ipm, fgerm = allscen$fgerm, sbgerm = allscen$sbgerm, mort = allscen$mort, sim = allscen$sim, 
                                 firemed = allscen$firemed, fireshape = allscen$fireshape,
                                 MoreArgs = list(pop = NA, num.years = 85, i.par = i.par, rates = vr, init = 7000, yrpoprand = randef, 
                                                 climate = centclim), SIMPLIFY = FALSE))


##################################################################################################
### Calculate min pop sizes & extinction
##################################################################################################

## Save the simulation output for future use
# save(popsize, file = file.path(ipmout, "IPMoutput.rda"))

## Output from the code above is the number of individuals at each timestep in the population 
## We can use this to calculate the minimum population size and quasi-extinction probability 

### Function for bootstrapped confidence intervals ###
require(boot)
meanfunc <- function(x, i){mn <- mean(x[i])}
bootfunc <- function(dataf){
  btci <- boot.ci(boot(dataf$minpop, statistic = meanfunc, R = 1000), type="norm")$normal[2:3]
  return(data.frame(minci = btci[1], maxci = btci[2], minp = mean(dataf$minpop), qe1 = mean(dataf$qe1)))
}


minpop <- popsize %>% group_by(sim, fgerm, sbgerm, sdmort, firemed, fireshape) %>%
  mutate(ntotal = nseed + nrose) %>% 
  summarise(minpop = min(ntotal), qe1 = ifelse(minpop<1, 1, 0)) %>% 
  group_by(fgerm, sbgerm, sdmort, firemed, fireshape) %>%
  do(bootfunc(.)) %>% 
  mutate(fireshape = factor(fireshape))

##################################################################################################
### Plot output
##################################################################################################

### Plot quasi-extinction probability
ggplot(minpop, aes(x=firemed, y = qe1, shape = fireshape, linetype = fireshape)) + 
  geom_point() + 
  geom_line() + 
  labs(x="Median fire return interval", y = "Extinction probability") + 
  theme_classic()


### Plot min pop sizes
ggplot(minpop, aes(x=firemed, y = minp, shape = fireshape, linetype = fireshape)) + 
  geom_point(position = position_dodge(1.5)) + 
  geom_line(position = position_dodge(1.5)) + 
  geom_errorbar(aes(ymin = minci, ymax = maxci), width =0.3, position = position_dodge(1.5)) + 
  labs(x="Median fire return interval", y = "Extinction probability") + 
  theme_classic()  
# Points are jittered to avoid overplotting

