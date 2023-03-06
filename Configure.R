
### Set up directories ###

datawd <- "./Data"
codewd <- "./Rcode"
dataout <- "./DataOut"

### Useful function ###
logistic <- function(p.vec)
{                
  p <- 1/(1+exp(-p.vec))                                          
  return(p)
}

### Load required packages ###

requiredPackageNames <- c("dplyr", "tidyr", "ggplot2", "lme4", "mgcv", "data.table", "Rmpi", "doParallel", "snow", "foreach")
# install (if we don't have it) and load the package 
sapply(requiredPackageNames, function(name) {
  if(!require(name, quietly=TRUE, character.only=TRUE)) {
    install.packages(name, repo = 'http://cran.uk.r-project.org')
    library(name, quietly=TRUE, character.only=TRUE)
  }
})