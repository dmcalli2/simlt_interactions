## 04_nested_models.R


drugs_select <- "Airways"
source('02_arrange_data_arrays.R') 

library(rjags)
# load.module('glm')
library(stringr)



#### Model with study nested inside condition and drug-class ----
modelstringnested <- "
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

} # end study
# Priors for each condition
mu9_mu[Dcondition]  ~ dnorm(0, 0.0001)           # mean dep:alloc
mu10_mu[Dcondition] ~ dnorm(0, 0.0001)           # mean pain:alloc
mu9_sd[Dcondition]  ~ dnorm(0,0.00001)T(0,)      # sd dep:alloc
mu10_sd[Dcondition] ~ dnorm(0,0.00001)T(0,)      # sd pain:alloc
mu9_prec[Dcondition] <- 1/mu9_sd[Dcondition]^2   # prec dep:alloc
mu10_prec[Dcondition] <- 1/mu10_sd[Dcondition]^2 # prec pain:alloc
}# end condition
}# end drug class
}# model
"
#Write model to text file
writeLines(modelstringnested, con= "jags/modelstringnested.txt")

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
                                Ncoef = 8 # note 8 of 10 coefficeints are independent for each study
                   ),
                   n.chains = 2,
                   n.adapt = 1000)
update(jags, 10000) # Burn in

data(LINE)
LINE$recompile()
LINE.out <- coda.samples(jags,
                         c('mu[1,1]', 'mu9_mu', 'mu9_prec', 'mu10_mu', 'mu10_prec'),
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