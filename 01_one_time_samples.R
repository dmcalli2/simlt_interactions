# 11_simulate_data_base_r
# Initially generates random samples from normal distribution. THese will be the same for each scenario
# Next Calculates overall group effect and variance for each component that will make-up the final scenarios
# eg (intercept, comorbidity effect, treatment allocation and interaction)

## Read in diabetes trials
load("../Trial_identify/clinical_trials_august_2017/scratch_data/data_for_simulation.Rdata")

# Remove single alpha-glucosidase inhibitor trial
diabetes_final <-  
  subset(diabetes_final, diabetes_final$atc_5 != "A10BF")


## Make same answer each set of classes (not each class)
set.seed(1234)



### Each simulation scenario, overall effects
main_scen <- expand.grid(
  cept = 0,
  como = c(-0.2,-0.1, 0, 0.1, 0.2),
  allc = c(-0.2,-0.1, 0, 0.1, 0.2),
  actn = c(-0.2,-0.1, 0, 0.1, 0.2)
)

varn_scen <- list(cept = c(0.25, 0.5),
                     como = c(0.25, 0.5),
                     allc = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     actn = c(0.05, 0.10, 0.15, 0.20, 0.25))

### Determine number of samples needed for each scenario, do random sampling once for all analyses
### for all drugs, classes and trial

# Simulate variation around wider group (A10B) effect at trial, drug and class level
SampleVarn <- function(item, smpl_lvls, n_iter = 1000){
  my_mtrx <- matrix(nrow = length(item) * n_iter,
                    ncol = length(smpl_lvls))
  
  for(j in seq_along(smpl_lvls)){
    my_mtrx[,j] <- rnorm(length(item) * n_iter, 0, smpl_lvls[j])
  }
  rownames(my_mtrx) <- rep(item, n_iter)
  colnames(my_mtrx) <- smpl_lvls
  my_mtrx
}
trial <- lapply(list(cept = c(0.25, 0.5),
                     como = c(0.25, 0.5),
                     allc = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     actn = c(0.05, 0.10, 0.15, 0.20, 0.25)),
                function(x) {
                  SampleVarn(diabetes_final$nct_id, x)
                }
)

drug <- lapply(list(cept = c(0.25, 0.5),
                     como = c(0.25, 0.5),
                     allc = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     actn = c(0.05, 0.10, 0.15, 0.20, 0.25)),
                function(x) {
                  SampleVarn(unique(diabetes_final$drug), x)
                }
)

atc5 <- lapply(list(cept = c(0.25, 0.5),
                     como = c(0.25, 0.5),
                     allc = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     actn = c(0.05, 0.10, 0.15, 0.20, 0.25)),
                function(x) {
                  SampleVarn(unique(diabetes_final$atc_5), x)
                }
)

## Combine variation around A10B effects into single estimate 
# from trial, drug and class (A10BA, A10BB etc)
diabetes_final_smpls <- diabetes_final[rep(seq_along(diabetes_final$nct_id),1000),
                                       c("nct_id", "drug", "atc_5") ]
CombineVarn <- function(component = "cept", choice = list(atc5 = "0.25",
                                                       drug = "0.25",
                                                       trial = "0.25")){
   atc5[[component]][diabetes_final_smpls$atc_5,  choice$atc5] +
   drug[[component]][diabetes_final_smpls$drug, choice$drug] +
  trial[[component]][diabetes_final_smpls$nct_id, choice$trial]
}
diabetes_final_smpls$cept <- CombineVarn()
diabetes_final_smpls$como <- CombineVarn("como", list(atc5 = "0.5",
                                                            drug = "0.5",
                                                            trial = "0.5"))
diabetes_final_smpls$allc <- CombineVarn("allc", list(atc5 = "0.15",
                                                            drug = "0.2",
                                                            trial = "0.25"))
diabetes_final_smpls$actn <- CombineVarn("actn", list(atc5 = "0.05",
                                                            drug = "0.1",
                                                            trial = "0.15"))

## Convert effects into marginal estimates of difference
## from A10 effect (eg comoribdity = intercept + comorbidity) 
diabetes_final_smpls_obs <- within(diabetes_final_smpls,{
  cept <- cept
  como <- cept + como
  allc <- cept + allc
  actn <- cept + como + actn
})

diabetes_final_smpls_obs2 <- within(diabetes_final_smpls, {
  cept = cept + main_scen$cept[1]
  como = como + main_scen$como[1]
  allc = allc + main_scen$allc[1]
  actn = actn + main_scen$actn[1]
})


main_scen_obs <- within(main_scen, {
  cept <- cept
  como <- cept+como
  allc <- cept+allc
  actn <- cept + como + actn
})

diabetes_final_smpls_obs2 <- within(diabetes_final_smpls2, {
  cept = cept + main_scen$cept[1]
  como = como + main_scen$como[1]
  allc = allc + main_scen$allc[1]
  actn = actn + main_scen$actn[1]
})


### Calcualate error terms based on SD
# Standardised sd is approx 1 for change in HbA1c
comorbidity_prev <- 0.2
sd <- 1
diabetes_final_se <- within(diabetes_final, {
  ncomo_se = sd/ ((1-comorbidity_prev) * diabetes_final$n_per_grp)^0.5
  ycomo_se = sd/ (   comorbidity_prev  * diabetes_final$n_per_grp)^0.5
})

diabetes_final_se$inter_se <- with(diabetes_final_se,{
  como_se^2 + ycomo_se^2
})



