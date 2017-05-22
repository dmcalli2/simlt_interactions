# 03 Categorical data analysis
# Choice of analyses - make in 01_simulate_outcome

drugs_select <- "Airways"
source('02_arrange_data_arrays.R') 

library(rjags)
load.module('glm')


# Write model ----
modelstring <- "
model{
  for(Ddrug in 1:Ndrug) {
    for(Dcondition in Sdrug[Ddrug]:(Sdrug[Ddrug+1]-1)){
      for(Dstudy in Scondition[Dcondition]:(Scondition[Dcondition+1]-1)){
        coef[, Dstudy] ~ dmnorm(mu[,Dstudy], vcov[, , Dstudy])
        for(i in 1:Ncoef) {
          mu[i, Dstudy] ~ dnorm(0, 0.0001)
        }
        } # each study
   }# each condition
   }# each drug class
}# model
"
writeLines(modelstring, con= "jags/model.txt")

# Run model
# Turn list of coefficeints into a matrix of coefficients
coef <-  map(myrag$con, ~ .x$coef)
coef <- sapply(coef, identity, simplify="matrix")
dim(coef)

# Turn list of matriced into array of matrices
vcov <-  map(myrag$con, ~ .x$vcov) 
vcov <- sapply(vcov, identity, simplify="array")
# turn variance into precision
vcov <- 1/vcov
dim(vcov)

jags <- jags.model('jags/model.txt',
                 data = list (Ndrug = myrag$Ndrug,
                              Sdrug = myrag$Sdrug,
                              Scondition = myrag$Scondition,
                              coef = coef,
                              vcov = vcov,
                              Ncoef = 10
                   ),
                   n.chains = 2,
                   n.adapt = 1000)
update(jags, 1000) # Burn in

data(LINE)
LINE$recompile()
LINE.out <- coda.samples(jags,
                         c('mu_indic3OR', 'mu_indic3_new', 'mu_indic3_mu', 'mu_indic3_sd'),
                         200000, thin = 10)

# Save model
save(LINE.out, file = paste0("Categorical Recover Interaction of ", gsub(".", "_", interaction_effect, fixed = TRUE), ".Rdata"))
a <- summary(LINE.out)$q
a <- a ["mu_indic3_mu", c("2.5%", "50%", "97.5%")]
a <- round(exp(a),2)
paste0(a[2], " (", a[1], " to ", a[3], ")")


pdf (file = paste0("Diagnostic plots for categorical of ", gsub(".", "_", interaction_effect, fixed = TRUE), ".pdf"))
plot( LINE.out, density = TRUE )
autocorr.plot (LINE.out)
dev.off()
