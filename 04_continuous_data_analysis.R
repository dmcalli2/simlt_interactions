# 04 Continuous data analysis
drugs_select <- "Antidiabetic"

source('02_arrange_data_arrays.R')
library(rjags)

# diabetes drugs

# write model ----
modelstring <- "
model{
  for(Dindication in 1:Nindication) {
    for(Dstudyno in Sindication[Dindication]:(Sindication[Dindication+1]-1)){
      for(Dobservations in Sstudyno[Dstudyno]:(Sstudyno[Dstudyno+1]-1)){
        y[Dobservations] ~ dnorm(theta[Dobservations],se_prec[Dobservations])
        theta[Dobservations] <- 
          beta0[Dstudyno] + # intercept for study
          beta1[Dstudyno]*alloc[Dobservations] + # treatment effect
          beta2[Dstudyno]*positive[Dobservations] + # stratifying variable effect
          beta3[Dstudyno]*alloc[Dobservations]*positive[Dobservations] # interaction effect
      }#observations
      #study-specific priors
      beta0[Dstudyno] ~ dnorm(0,0.00001)
      beta1[Dstudyno] ~ dnorm(0,0.00001)
      beta2[Dstudyno] ~ dnorm (0,0.0001) # stratifying variable, no need to share this as not estimating it
      # Shared priors (hyperpriors)
      # beta3[Dstudyno] ~ dt(mu_indic3[Dindication], prec_indic3[Dindication], 3) # interaction stratifying variable and treatment
      beta3[Dstudyno] ~ dnorm(mu_indic3[Dindication], prec_indic3[Dindication]) # with normal rather than t-distribution for block-updating
    }#studies
    #indication-specific priors, sharing information on mean stratifying and interaction effect across studies within indications
    mu_indic3[Dindication] ~ dnorm (mu_indic3_mu, mu_indic3_prec)
    mu_sd3[Dindication] ~ dnorm (mu_sd3_mu, mu_sd3_prec)T(0,)
    prec_indic3[Dindication] <- 1 / mu_sd3[Dindication]^2
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

writeLines(modelstring,con="hiercon.txt")


# Run model
jags <- jags.model('hiercon.txt',
                   data = list (Nindication = myrag_con$Nindication,
                                Sindication = myrag_con$Sindication,
                                Sstudyno = myrag_con$Sstudy,
                                alloc = myrag_con$Xobservations[,"alloc"],
                                positive = myrag_con$Xobservations[,"positive"],
                                se_prec = myrag_con$Xobservations[,"se_prec"],
                                y = myrag_con$y
                   ),
                   n.chains = 3,
                   n.adapt = 1000)
update(jags, 10000) # Burn in

data(LINE)
LINE$recompile()
LINE.out <- coda.samples(jags,
                         c('mu_indic3', 'mu_indic3_mu', 'mu_indic3_sd', 'mu_indic3_new'),
                         200000, thin = 10)

# Save model
save(LINE.out, file = paste0("Continuous Recover Interaction of ", gsub(".", "_", mean_interaction_effect, fixed = TRUE), ".Rdata"))
summary(LINE.out)

# 
a <- summary(LINE.out)$q
a <- a ["mu_indic3_mu", c("2.5%", "50%", "97.5%")]
a <- round((a),2)
paste0(a[2], " (", a[1], " to ", a[3], ")")
paste0(a[2], " ", a[1], " to ", a[3], "")


# 
# pdf (file = paste0("Diagnostic plots for continuous of ", gsub(".", "_", mean_interaction_effect, fixed = TRUE), ".pdf"))
# plot( LINE.out, density = TRUE )
# autocorr.plot (LINE.out)
# dev.off()
