################# IPM (density independent recruitment) ###############

## Use a weibull distribution to simulate fire frequencies
firefunc <- function(weib.median, weib.shape, nyears){
  if(is.numeric(weib.median)){
  weib.scale <- weib.median/(log(2)^(1/weib.shape))
  fire <- vector(mode='logical', length = nyears)
  tt <- 0
  while (tt < nyears) {
    tt1 <- rweibull(1, weib.shape, weib.scale)
    tt <- tt + tt1
    if (tt <= nyears) {
      fire[tt] <- TRUE
    }}
  tsf <- rep(1, times = nyears)
  for(i in 2:nyears){
    if(fire[i]) { tsf[i] <- 0 } else {tsf[i] <- tsf[i-1]+1} 
  } } else {
    tsf <- 1:nyears
  }
  return(tsf)
}

enveffect <- function(nyr, possrand, tsf, climd, vrmodels, pop = NA){
  ## Year effects - select at random at each year in the simulation
  yr <- sample(possrand$allyear, nyr, replace = TRUE)
  ## Population random effects - select pop at random from possible 11 (NB: no recruit size data from pop 95)
  pop <- ifelse(is.na(pop), sample(possrand$allpop, 1), pop)
  ## Seedling survival 
  sdlng <- logistic(-1.953 + -0.028 * tsf + rnorm(nyr, 0, 1.09))
  ## Sample from climate effects
  clyr <- sample(1:24, nyr, replace=TRUE) 
  # For the past climate simulations - at each year of the simulation we sample from the past years
  # For the future climate simulations - instead we would modify this to take each year in turn here, rather than sampling
  ## Climate effects
  topred <- data.frame(EryYear = climd$EryYear)
  topred$growclim <- predict(vrmodels[["grow"]], type='terms', exclude = c("s(Sizem1)", "s(tsf)", "s(bald)", "s(yrpop)"),
                      newdata = data.frame(Sizem1 = 1, tsf = 1, fort = climd$EryFort, 
                                           climate = climd$maxtemp_c, bald = 16, yrpop = 161992))
  topred$survclim <- predict(vrmodels[["surv"]], type='terms', exclude = c("s(Sizem1)", "s(tsf)", "s(bald)", "s(yrpop)"),
                      newdata = data.frame(Sizem1 = 1, tsf = 1, fort = climd$EryFort, 
                                           climate = climd$drought_c, bald = 16, yrpop = 161992))
  topred$fecclim <- predict(vrmodels[["fec"]], type='terms', exclude = c("s(Sizem1)", "s(tsf)", "s(bald)", "s(yrpop)"),
                      newdata = data.frame(Sizem1 = 1, tsf = 1, fort = climd$EryFort, 
                                           climate = climd$mintemp_c, bald = 16, yrpop = 161992))
  topred$recclim <- predict(vrmodels[["rec"]], type='terms', exclude = c("s(tsf)", "s(bald)", "s(yrpop)"),
                      newdata = data.frame(Sizem1 = 1, tsf = 1, fort = climd$EryFort, 
                                           climate = climd$maxtemp_c, bald = 16, yrpop = 161992))
  topred <- topred %>% group_by(EryYear) %>% 
    summarise(grow = sum(growclim), surv = sum(survclim), fec = sum(fecclim), rec = sum(recclim))
  return(data.frame(pop = pop, yr = yr, sdlng = sdlng, cgrow = topred$grow[clyr], csurv = topred$surv[clyr], 
                    cfec = topred$fec[clyr], crec = topred$rec[clyr]))
}

######################## IPM functions ###########################

# Growth function
g_lst1lst <- function(lst1, lst, growmod, fire, env)
{
  mean <- predict(growmod, newdata = data.frame(tsf = fire, Sizem1 = lst, bald = env$pop, yrpop = env$yr, 
                                                           fort = 1, climate = 0)) + env$cgrow
  sd <- sqrt(growmod$sig2)                                                     
  p.den.grow <- dnorm(lst1, mean = mean, sd = sd)                                
  return(p.den.grow)
}

# Survival function
s_lst <- function(lst, survmod, fire, env)
{
  p <- logistic(predict(survmod, newdata = data.frame(tsf = fire, Sizem1 = lst, bald = env$pop, yrpop = env$yr, 
                                             fort = 1, climate = 0)) + env$csurv)
  return(p)
}

# Fecundity function (number of flowering branches)

fn_lst <- function(lst, fecmod, fire, env)
{
  mu <- exp(predict(fecmod, newdata = data.frame(tsf = fire, Sizem1 = lst, bald = env$pop, yrpop = env$yr, 
                                                       fort = 1, climate = 0)) + env$cfec)
  return(mu)
}

# Recruit size function (not dependent on parental size)

c_lst1 <- function(lst1, recmod, fire, env)
{
  mean <- predict(recmod, newdata = data.frame(tsf = fire, bald = env$pop, yrpop = env$yr, 
                                                         fort = 1, climate = 0)) + env$crec
  sd <- sqrt(recmod$sig2)
  p.den.rcsz <- dnorm(lst1, mean=mean, sd=sd)
  return(p.den.rcsz)
}

################### Functions to build IPM kernels #################

## Survival-growth kernel
P_lst1lst <- function (lst1, lst, vrmodels, tsf, enviro) {
  return( s_lst(lst, vrmodels[["surv"]], tsf, enviro) *  g_lst1lst(lst1, lst, vrmodels[["grow"]], tsf, enviro))
}

## Seed production kernel
F_lst <- function (lst, vrmodels, tsf, enviro) 
{
  return(fn_lst(lst, vrmodels[["fec"]], tsf, enviro) * 182.5)  # NB: 182.5 is mean number of seeds per flowering branch (from 2004 paper)
}

## Fecundity kernel
F_lst1lst <- function (lst1, lst, vrmodels, tsf, enviro) 
{
  return(fn_lst(lst, vrmodels[["fec"]], tsf, enviro) * 182.5 * c_lst1(lst1 = lst1, recmod = vrmodels[["rec"]], fire = tsf, enviro)) 
  # NB: 182.5 is mean number of seeds per flowering branch (from 2004 paper)
}

# Calculate the 'implementation' parameters
calc_i_par <- function(m, L, U) {
  h <- (U - L) / m
  meshpts <- L + ((1:m) - 1/2) * h
  iS <- m+1 # index for seeds
  iR <- seq_len(m) # index for rosettes
  return(list(m = m, h = h, iS = iS, iR = iR, meshpts = meshpts))
}

# calculate a kernel for a particular year
mk_K <- function(i.par, vital, fire, d, gf, gsb, rand, clim) {
  with(i.par, {
    P <- F <- matrix(0, nrow = m+1, ncol = m+1)  # +1 as seedbank
    # Does fire occur that year?
    if(fire == 0) {
      # All rosettes will be killed in the fire (includes adults surviving from last year and this years recruits)
      #P[iR,] <- F[iR,] <- 0     
      # Survival: seeds from seeds - assume seed survival unaffected by fire
      P[iS, iS] <- (1-d) * (1-gsb)
      # No reproduction as all rosettes killed in fire
      #F[iS, iR] <- 0
    } else {
      # survival-growth: seeds from seeds (1)
      P[iS, iS] <- (1-d) * (1-gsb)
      # survival-growth: rosettes from seeds (2)
      P[iR, iS] <- gsb * rand$sdlng * c_lst1(meshpts, vital[["rec"]], fire, rand)    
      # survival-growth: rosettes from rosettes (3)
      P[iR, iR] <- outer(meshpts, meshpts, P_lst1lst, vital, fire, rand) * h
      # reproduction: seeds from rosettes (4)
      F[iS, iR] <- (1-d) * (1-gf) * F_lst(meshpts, vital, fire, rand) * h
      # reproduction: rosettes from rosettes (5)
      F[iR, iR] <- gf * rand$sdlng * outer(meshpts, meshpts, F_lst1lst, vital, fire, rand) * h
    }
    # keep anything we may need in a list
    list(F = F, P = P, K = F+P)
  })
}

############ Run IPM ############## 
run.ipm <- function(num.years, i.par, rates, pop, firemed, fireshape, fgerm, sbgerm, mort, init, yrpoprand, sim, climate, scpop=F, grrate=F){
  ## Simulate fires and select climate year
  fire <- firefunc(weib.median = firemed, weib.shape = fireshape, nyears = num.years)
  ## Simulate population
  with(i.par, {
    state.vec <- numeric(m+1)
    state.vec[iS] <- init  # Menges et al 2004 use 7000 seeds as starting population
    ## Sample from random year and population effects
    randef <- enveffect(num.years, pop = pop, yrpoprand, fire, climate,rates)
    nrose <- nseed <- meanlst <- nt <- grate <- numeric()
    ## Iterate population 
    for (tt in 1:num.years) {
      state.vec <- mk_K(i.par, fire = fire[tt], vital =rates, d = mort, gsb = sbgerm, gf = fgerm, 
                        rand = randef[tt,])$K %*% state.vec
      nseed[tt] <- state.vec[iS]  
      nt <- state.vec[iR]
      nrose[tt] <- sum(nt * h)
      grate[tt] <- log(nrose[tt]/nrose[tt-1])
      meanlst[[tt]] <- sum(nt * meshpts)/sum(nt)
    } 
    if(scpop==T){
      return(nrose/max(nrose))
    } else{
      if(grrate==T){
        return(grate)
      } else{
        output <- data.frame(nseed = nseed, nrose = nrose, meanz = meanlst, it = 1:num.years, sim, pop = randef$pop[1],
                             fgerm = fgerm, sbgerm = sbgerm, sdmort = mort, firemed = firemed, fireshape = fireshape)
        return(output)
      }
    }
  })
}
