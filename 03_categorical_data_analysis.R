# 03 Categorical data analysis
# Choice of analyses - make in 01_simulate_outcome

drugs_select <- "Airways"
source('02_arrange_data_arrays.R') 

library(rjags)
load.module('glm')
library(stringr)
## Use double colon for these rather than library
# library("matrixcalc") 
# library(Matrix)

# Write model ----
modelstring <- "
model{
  for(Ddrug in 1:Ndrug) {
    for(Dcondition in Sdrug[Ddrug]:(Sdrug[Ddrug+1]-1)){
      for(Dstudy in Scondition[Dcondition]:(Scondition[Dcondition+1]-1)){
        coef[, Dstudy] ~ dmnorm(mu[,Dstudy], vcov[, , Dstudy])
        for(i in 1:Ncoef) { # 8 independent priors for first 8 coefficients
        mu[i, Dstudy] ~ dnorm(0, 0.0001)
        }
        mu[9, Dstudy]  ~ dnorm(mu9_mu[Dcondition],  mu9_prec[Dcondition])
        mu[10, Dstudy] ~ dnorm(mu10_mu[Dcondition], mu10_prec[Dcondition])

      } # each study
     mu9_mu[Dcondition] ~ dnorm(0, 0.0001)
     mu10_mu[Dcondition]~ dnorm(0, 0.0001)
    }# each condition
  }# each drug class
}# model
"

function (Simple_model){
### Simple model
# modelstring <- "
# model{
#   coef[] ~ dmnorm(mu[], vcov[,])
#   for(i in 1:Ncoef) {
#     mu[i] ~ dnorm(0, 0.0001)
#   }
# }# model
# "
}

writeLines(modelstring, con= "jags/model.txt")

# Run model
# Turn list of coefficeints into a matrix of coefficients
coef <-  map(myrag$con, ~ .x$coef)
coef <- sapply(coef, identity, simplify="matrix")
dim(coef)
# Reduce to vector for analysis
# coef <- coef[,1]

# Turn list of matriced into array of matrices
vcov <-  map(myrag$con, ~ .x$vcov) 
vcov <- sapply(vcov, identity, simplify="array")
dim(vcov)

# reduce to matrix
# vcov <- vcov[,,1]
apply(vcov, 3, matrixcalc::is.positive.definite)
apply(vcov, 3, Matrix::isSymmetric)

jags <- jags.model('jags/model.txt',
                 data = list (Ndrug = myrag$Ndrug,
                              Sdrug = myrag$Sdrug,
                              Scondition = myrag$Scondition,
                              coef = coef,
                              vcov = vcov,
                              Ncoef = 8
                   ),
                   n.chains = 2,
                   n.adapt = 1000)
update(jags, 10000) # Burn in

data(LINE)
LINE$recompile()
LINE.out <- coda.samples(jags,
                         c('mu'),
                         20000)
res <- summary(LINE.out)
compare_coefs <- res$statistics[, "Mean"]
compare_coefs <- tibble(lab = names(compare_coefs), mu = compare_coefs) %>% 
  mutate(lab = str_sub(lab, 4, -2)) %>% 
  separate(lab, into = c("coef", "study"), sep = ",")

mdld <- by(compare_coefs, compare_coefs$study, function (x) x$mu)
mdld <- tibble(study = names(mdld), mu = mdld)
cmpr <- mdld %>% 
  inner_join(mydf %>% mutate(study = as.character(study)), by = "study")
coefs <- map(cmpr$con, function (x) x$coef)

pdf("figures/examine coefficients.pdf")
par(mfrow = c(3,3))
pmap(list(coefs, cmpr$mu,  cmpr$enrollment), 
     function (x, y, z)  plot(x, y, xlim = c(-1,1), ylim = c(-1,1), main = z))
dev.off()

######## Diagnostics of model fit
pdf (file = "figures/Diagnostic plots2.pdf")
par(mfrow = c(4,3))
autocorr.plot(LINE.out)
gelman.plot(LINE.out, ylim = c(0.95, 1.2))
dev.off()
