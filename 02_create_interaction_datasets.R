# 02_create_interaction_datasets.R
# Create effect estimates interactions with variation at trial, drug and class level
atc5 <- atc5["actn"]
drug <- drug["actn"]
trial <- trial["actn"]

count <- 0
res <- matrix(nrow = nrow(diabetes_final_smpls), ncol = 5^3)
res_names <- vector(length = 5^3)
for (i in 1:5){
  for(j in 1:5){
    for(k in 1:5){
      count <- count + 1
    res[, count] <- 
       atc5$actn[diabetes_final_smpls$atc_5,  c("0.05", "0.1", "0.15", "0.2", "0.25")[i]] +
       drug$actn[diabetes_final_smpls$drug,   c("0.05", "0.1", "0.15", "0.2", "0.25")[j]] +
      trial$actn[diabetes_final_smpls$nct_id, c("0.05", "0.1", "0.15", "0.2", "0.25")[k]] 
    res_names[count] <- paste(c("0.05", "0.1", "0.15", "0.2", "0.25")[i],
                       c("0.05", "0.1", "0.15", "0.2", "0.25")[j],
                       c("0.05", "0.1", "0.15", "0.2", "0.25")[k],
                       sep = "_")
    }
  }
}

colnames(res) <- res_names

scen1 <- diabetes_final_smpls

# Add in main effect to a chosen variation scenario
scen1$res <- res[, "0.1_0.1_0.1"] + -0.2

# Choose an iteration
scen1_iter1 <- scen1[scen1$iter == 1,]

# Add in se term for that iteration
comorbidity_prev <- 0.2
sd <- 1
# calculate SE for comorbidity adn non-comorbidity group (same for placebo and treatment)
ncomo_se = sd/ ((1-comorbidity_prev) * diabetes_final$n_per_grp)^0.5
ycomo_se = sd/ (   comorbidity_prev  * diabetes_final$n_per_grp)^0.5
# Calculate SE for interaction
inter_prec <- 1/(2*ncomo_se^2 + 2*ycomo_se^2)

library(INLA)

## Compare INLA
my_drug_n <- as.numeric(as.factor(scen1_iter1$drug))
study_id_n <- as.numeric(as.factor(as.character(scen1_iter1$nct_id)))
atc5_n <- as.numeric(as.factor(as.character(scen1_iter1$atc_5)))

my_data <- data.frame(y = scen1_iter1$res, y_prec = inter_prec, 
                      trial = study_id_n,
                      myatc4 = 1,
                      myatc5 = atc5_n,
                      mydrug = my_drug_n)


##Model formula trial within ATC4
myform1 <- y ~ -1 + # no intercept
  myatc4 + # single coefficient for overall treatment effect (prior for this specified below)
  f(trial, # Treatment in each trial (difference from overall treatment effect)
    model = "iid", # Standard assumption as trials uncorrelated
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 0.25)))) 
      # Half-normal prior for difference between trials

## linear combinations
# lc_list <-lapply(1:161, function (x) {
#   vect <- rep(NA, 161)
#   vect[x] <- 1
#   inla.make.lincomb(myatc4 = 1, trial = vect) 
# })
# lc_flat_list <- do.call(c, lc_list)
# names(lc_flat_list) <- paste0("lc", 1:161)

## Run model and specify rest of priors, and that likelihood fixed
mod1 <- inla(myform1, 
             data = my_data,
              # Calculate linear combinations
 #            lincomb = lc_flat_list,
            # control.inla(list(lincomb.derived.only = FALSE)), # Make more acurate estimate still half
             # Likelihood distribution
             family = "gaussian",
             # Fix likelihood hyperpars as the data is fixed with known precision.
             control.family = list(hyper = list(prec = list(fixed = TRUE, initial = 0))),
             # Likelihood precisions
             scale = my_data$y_prec,
             # Prior distribution for "fixed" effects - really for mu_mu
             control.fixed = list(mean = 0, prec = 0.25),
             # Optionally compute DIC
             verbose = FALSE)
summary(mod1)

## Run more complex model, trial within drug within ATC4 class

myform_nested <- y ~ -1 + myatc4 + 
  f(trial, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1)))) +
  f(mydrug, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 0.01))))

mod1_nested <- inla(myform_nested, 
             data = my_data,
             # Add linear combinations to estimate drug-class
             # lincomb = dc_all,
             # Likelihood distribution
             family = "gaussian",
             # Fix likelihood hyperpars as the data is fixed with known precision.
             control.family = list(hyper = list(prec = list(fixed = TRUE, initial = 0))),
             # Likelihood precisions
             scale = my_data$y_prec,
             # Prior distribution for "fixed" effects - really for mu_mu
             control.fixed = list(mean = 0, prec = 0.25),
             # Optionally compute DIC
             verbose = FALSE)
summary(mod1_nested)


## Run more complex model still, trial within drug within ATC5 class within ATC4 class
myform_nested2 <- y ~ -1 + myatc4 + 
  f(trial, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1)))) +
  f(mydrug, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 0.01)))) +
  f(myatc5, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 0.01))))

mod1_nested2 <- inla(myform_nested2, 
             data = my_data,
             # Add linear combinations to estimate drug-class
             # lincomb = dc_all,
             # Likelihood distribution
             family = "gaussian",
             # Fix likelihood hyperpars as the data is fixed with known precision.
             control.family = list(hyper = list(prec = list(fixed = TRUE, initial = 0))),
             # Likelihood precisions
             scale = my_data$y_prec,
             # Prior distribution for "fixed" effects - really for mu_mu
             control.fixed = list(mean = 0, prec = 0.25),
             # Optionally compute DIC
             verbose = FALSE)
summary(mod1_nested2)
