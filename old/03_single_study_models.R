# "03_single_study_models.R"
# Choice of analyses - make in 01_simulate_outcome

# Read in sim_study
sim_data <- readRDS(file = "scratch_data/simulated_data1.Rds")

library(tidyverse)
library(runjags)
# load.module('glm')
library(stringr)

## Use double colon for these rather than library
# library("matrixcalc") 
# library(Matrix)
#library("clusterGeneration")


# Run models ----

#### Attempt to recover coefficients for a single study, fix precision very high ----
simplemodelstring <- "
model{
  coef[] ~ dmnorm(mu[], vcov_prec[,])
  for(i in 1:Ncoef) {
    mu[i] ~ dnorm(0, 0.0001)
  }
}# model
"
#Write model to text file
writeLines(simplemodelstring, con= "jags/simplemodelstring.txt")

## Create vcov with very low variance, then take recoprocals for
## very high precision
a<- matrix(rep(1e-10,100), nrow = 10)  # covariance
diag(a) <- 1e-6 #  variance
a <- solve(a) # invert to get precision, which is what JAGS expects
# Check matrix
Matrix::isSymmetric(a)
matrixcalc::is.positive.definite(a %>%  round(3))

# Run model for a single study
study_choose <- 1

jags <- jags.model('jags/simplemodelstring.txt',
                   data = list (coef = sim_data$con[[study_choose]]$coef,
                                vcov_prec = a,
                                Ncoef = 10
                   ),
                   n.chains = 2,
                   n.adapt = 1000)
data(LINE)
LINE$recompile()
LINE.out <- coda.samples(jags,
                         c('mu'),
                         5000)
CmprCoef <- function (){
  cmpr <- tibble (original = sim_data$con[[study_choose]]$coef,
          res = summary(LINE.out)$statistics[,"Mean"])
  cmpr$diff <- cmpr$res - cmpr$original
  cmpr
}
CmprCoef()
## Recovers coefficients

#### Attempt to recover coefficients for a single study, use actual precision matrix ----
# Run model for actual vcov
res <- vector(length = nrow(sim_data), mode = "list")
for (study_choose in seq_along(sim_data$study_id)){
  print(study_choose)
  jags <- autorun.jags('jags/simplemodelstring.txt',
                     data = list (coef = sim_data$con[[study_choose]]$coef,
                                  # Note round precision matrix to make sure symmetrical
                                  vcov_prec = round(sim_data$con[[study_choose]]$prec_matrix,3),
                                  Ncoef = 10
                     ),
                      n.chains = 2,
                     adapt = 1000,
                     startsample = 5000,
                   monitor = "mu")
  summary(jags)

  ## Recovers coefficients
  
  ## Runmodel with coefficients and variance estimates independently rather than part of a nultivariate normal
  
  separatemodelstring <- "
  model{
    dep  ~ dnorm(mu_dep, dep_prec)
    pain ~ dnorm(mu_pain, pain_prec)
    mu_dep  ~ dnorm(0, 0.0001)
    mu_pain ~ dnorm(0, 0.0001)
  }# model
  "
  
  vcov_prec <- round(sim_data$con[[study_choose]]$prec_matrix,3)
  dep_prec  <- vcov_prec["dep:alloc", "dep:alloc"]
  pain_prec <- vcov_prec["pain:alloc", "pain:alloc"]
  dep  <- sim_data$con[[study_choose]]$coef["dep:alloc"]
  pain <- sim_data$con[[study_choose]]$coef["pain:alloc"]
  
  #Write model to text file
  writeLines(separatemodelstring, con= "jags/sepmodelstring.txt")
  
  jags_sep <- autorun.jags('jags/sepmodelstring.txt',
                     data = list (dep = dep,
                                  dep_prec = dep_prec,
                                  pain = pain,
                                  pain_prec = pain_prec),
                     n.chains = 2,
                     startsample = 5000,
                     monitor = c('mu_dep', 'mu_pain'),
                     modules = "glm")
   summary(jags_sep)
  
  ## Compare single versus separate multiple times
  res[[study_choose]] <- (summary(jags)[c("mu[9]", "mu[10]"),] - summary(jags_sep))
}
saveRDS(res, file = "model_summaries/compare_separate_multivariate_models.Rds")


res2 <- do.call(rbind, res)
res2 <- res2[, 1:6]
any(res2 > 0.001)
