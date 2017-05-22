# 03 Categorical data analysis
# Choice of analyses - make in 01_simulate_outcome

drugs_select <- "Haemostasis"
source('02_arrange_data_arrays.R') # this sources "01_simulate_outcome_data.r"

library(rjags)
load.module('glm')


# Write model ----
modelstring <- "
model{
  for(Dindication in 1:Nindication) {
    for(Dstudyno in Sindication[Dindication]:(Sindication[Dindication+1]-1)){
      for(Dobservations in Sstudyno[Dstudyno]:(Sstudyno[Dstudyno+1]-1)){
        events[Dobservations] ~ dbin(p[Dobservations],n[Dobservations])
        logit(p[Dobservations]) <- 
        beta0[Dstudyno] + # intercept for study
        beta1[Dstudyno]*alloc[Dobservations] + # treatment effect
        beta2[Dstudyno]*positive[Dobservations] + # stratifying variable effect
        beta3[Dstudyno]*alloc[Dobservations]*positive[Dobservations] # interaction effect
        }#observations
      #study-specific priors
      beta0[Dstudyno] ~ dnorm(0,0.00001)
      beta1[Dstudyno] ~ dnorm(0,0.00001)
      beta2[Dstudyno] ~ dnorm (0,0.0001) # stratifying variable, no need to share this as not estimating it
      beta3OR[Dstudyno] <- exp(beta3[Dstudyno])
      # Shared priors (hyperpriors)
      # beta3[Dstudyno] ~ dt(mu_indic3[Dindication], prec_indic3[Dindication], 3) # interaction stratifying variable and treatment
      beta3[Dstudyno] ~ dnorm(mu_indic3[Dindication], prec_indic3[Dindication]) # with normal rather than t-distribution for block-updating
    }#studies
    #indication-specific priors, sharing information on mean stratifying and interaction effect across studies within indications
    mu_indic3[Dindication] ~ dnorm (mu_indic3_mu, mu_indic3_prec)
    mu_sd3[Dindication] ~ dnorm (mu_sd3_mu, mu_sd3_prec)T(0,)
    prec_indic3[Dindication] <- 1 / mu_sd3[Dindication]^2
#     mu_sd3[Dindication] ~  dlnorm(-4.06,1.45^2) # informative prior from Turner et al 2012
#     prec_indic3[Dindication] <- 1 / mu_sd3[Dindication] 
    mu_indic3OR[Dindication] <- exp(mu_indic3[Dindication])
  }#indications
  # Priors on indication, want to share information across indications
  # Priors on interacting variable parameter across indications
  mu_indic3_mu ~ dnorm (0, 0.0001) # average stratifying variable effect could range from 0.31 to 3.2 aross indications
  mu_indic3_sd ~ dnorm (0, 0.0001)T(0,) # plausible stratifying variable effect could range from same across indications to a 4-fold difference
  mu_indic3_prec <- 1 / mu_indic3_sd^2
  mu_sd3_mu ~ dnorm(0,0.0001)
  mu_sd3_sd ~ dnorm (0, 0.0001)T(0,) 
  mu_sd3_prec <- 1 / mu_sd3_sd^2

  mu_indic3_new ~ dnorm(mu_indic3_mu, mu_indic3_prec)
}# model
"
writeLines(modelstring,con="hierlog.txt")


# Run model
jags <- jags.model('hierlog.txt',
                 data = list (Nindication = myrag_cat$Nindication,
                                Sindication = myrag_cat$Sindication,
                                Sstudyno = myrag_cat$Sstudy,
                                alloc = myrag_cat$Xobservations[,"alloc"],
                                positive = myrag_cat$Xobservations[,"positive"],
                                n = myrag_cat$Xobservations[,"n"],
                                events = myrag_cat$events
                   ),
                   n.chains = 3,
                   n.adapt = 1000)
update(jags, 100000) # Burn in

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
