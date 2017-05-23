# 03 Categorical data analysis
# Choice of analyses - make in 01_simulate_outcome

drugs_select <- "Airways"
source('02_arrange_data_arrays.R') 

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
study_choose <- 118
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
res <- summary(LINE.out)
cbind(original = sim_data$con[[study_choose]]$coef,
      res$statistics)
## Recovers coefficients to 5 dps
vcov <- sim_data$con[[study_choose]]$vcov

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
gelman.diag(LINE.out)

pdf("figures/Diagnostic_plot_single_study.pdf")
gelman.plot(LINE.out, ylim = c(0.9, 1.1))
dev.off()

# Model converged by 5000
update(jags, 5000)
data(LINE)
LINE$recompile()
LINE.out <- coda.samples(jags,
                         c('mu'),
                         20000)

cbind(original = sim_data$con[[study_choose]]$coef,
      summary(LINE.out)$statistics)

