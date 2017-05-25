# "03_single_study_models.R"
# Choice of analyses - make in 01_simulate_outcome

# Read in sim_study
sim_data <- readRDS(file = "scratch_data/simulated_data.Rds")

library(tidyverse)
library(rjags)
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
study_choose <- sample(seq_along(sim_data$study_id), 1)
print(study_choose)
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
jags <- jags.model('jags/simplemodelstring.txt',
                   data = list (coef = sim_data$con[[study_choose]]$coef,
                                # Note round precision matrix to make sure symmetrical
                                vcov_prec = round(sim_data$con[[study_choose]]$prec_matrix,3),
                                Ncoef = 10
                   ),
                   n.chains = 2,
                   n.adapt = 1000)
data(LINE)
LINE$recompile()
LINE.out <- coda.samples(jags,
                         c('mu'),
                         5000)
cbind(original = sim_data$con[[study_choose]]$coef,
      summary(LINE.out)$statistics)
CmprCoef()
gelman.diag(LINE.out)

gelman.plot(LINE.out, ylim = c(0.9, 1.1), ask = FALSE)
pairs(window(LINE.out, start = 1000, end = 1200) %>% as.matrix(),
       labels = names(sim_data$con[[study_choose]]$coef))

# Model converged by 5000
update(jags, 5000)
data(LINE)
LINE$recompile()
LINE.out <- coda.samples(jags,
                         c('mu'),
                         20000)
CmprCoef()
print(study_choose)
print(sim_data$enrollment[study_choose])
## Recovers coefficients
