# 05_model_diagnnostics.R

library(coda)
library(tidyverse)
library(stringr)
library(runjags)

diagn <- readRDS("mcmc_chains/cat_dc_nest_inform_effect5.Rds")

summary(diagn)

## Select random sample of mu nodes
varnames <- varnames(diagn)
non_mu <- varnames[str_sub(varnames,1,3) != "mu["]
mu <- varnames[str_sub(varnames,1,3) == "mu["]
mu <- sample (x = mu, size = 10)

# Plot nodes
diagn_sub <- diagn[ , c(non_mu, mu) ]
diagn_pairs <- diagn[, non_mu] 
diagn_pairs <- map(diagn_pairs, function (x) x[4000:4100, ]) 
diagn_pairs <- do.call(rbind, diagn_pairs)

gelman.plot(diagn_sub, ylim = c(0.8,1.2), ask = FALSE)
pairs(diagn_pairs[, 1:10])
coda::acfplot(diagn_sub, ask = FALSE)
