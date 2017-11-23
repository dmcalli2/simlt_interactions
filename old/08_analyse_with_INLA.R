##08_analyse_with_INLA

# Install INLA package (not available on CRAN)
update.packages("INLA", repos="http://www.math.ntnu.no/inla/R/stable")


library(tidyverse)
library(INLA)
library(forcats)
library(runjags)

drugs_select <- "Airways"
source('02_arrange_data_arrays.R') 

## Select the sixth simulation
res <- res_list[[6]]
res <- res %>% 
  filter(drug_group_allocation == drugs_select) %>% 
  group_by(my_drug) %>% 
  sample_n(4) %>% 
  arrange(my_drug) %>% 
  ungroup()
rm(res_list)

dep <- map_dbl(res$con, ~ .x$coef["dep:alloc"])
dep_prec <- map_dbl(res$con, ~ .x$prec_matrix["dep:alloc", "dep:alloc"])
res$my_drug_n <- as.numeric(as_factor(res$my_drug))
res$study_id_n <- as.numeric(as_factor(as.character(res$study_id)))

my_data <- data.frame(y = dep, y_prec = dep_prec, group = res$study_id_n, 
                      myclass = 1, myclass_act = res$my_drug_n)


## JAGS model
modelcompareinla <- "
model{
  for(i in 1:12) {
        y[i] ~ dnorm(mu[i], y_prec[i])
        mu[i]  ~ dnorm(mu_mu,  between_trial_prec)
  }# 
  mu_mu  ~ dnorm(0, 0.25)
  between_trial_sd ~ dnorm(0, 0.25)T(0,)
  between_trial_prec <- 1/between_trial_sd^2
}# model
"
writeLines(modelcompareinla, con= "jags/compare_inla.txt")

jags <- autorun.jags('jags/compare_inla.txt',
                     data = list(y = my_data$y, y_prec = my_data$y_prec),
                      n.chains = 2,
                     adapt = 1000,
                     startsample = 10000,
                   monitor = c("mu", "mu_mu", "between_trial_sd"))
summary(jags)

## Compare INLA
##Model formula
myform1 <- y ~ -1 + # no intercept
  myclass + # single coefficient for overall treatment effect (prior for this specified below)
  f(group, # Treatment in each trial (difference from overall treatment effect)
    model = "iid", # Standard assumption as trials uncorrelated
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 0.25)))) 
      # Half-normal prior for difference between trials

## linear combinations
lc_list <-map(1:12, function (x) {
  vect <- rep(NA, 12)
  vect[x] <- 1
  inla.make.lincomb(myclass = 1, group = vect) 
})
lc_flat_list <- do.call(c, lc_list)
names(lc_flat_list) <- paste0("lc", 1:12)

## Run model and specify rest of priors, and that likelihood fixed
mod1 <- inla(myform1, 
             data = my_data,
              # Calculate linear combinations
             lincomb = lc_flat_list,
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
#mod1_inc_acc <- inla.hyperpar(mod1) # increase precision of hyperparameter doesn't do much

### Compare JAGS and INLA model
summary(mod1)
summary(jags)["mu_mu",]

# Compare trial-specific effects
trial_inla_mean_manual  <- mod1$summary.random[[1]][,"mean"] + mod1$summary.fixed[,"mean"]
trial_inla_mean_derived <- mod1$summary.lincomb.derived$mean
trial_inla_sd <- mod1$summary.lincomb.derived$sd

trial_jags_mean <- summary(jags)[1:12, c("Mean")]
trial_jags_sd <- summary(jags)[1:12, c("SD")]

round(cbind(trial_inla_mean_manual, trial_inla_mean_derived, trial_jags_mean),3)
round(cbind(trial_inla_sd, trial_jags_sd),4)

mod1$marginals.fixed %>%  as.data.frame () %>%  plot()
re_mrg <- mod1$marginals.random[[1]]
re_mrg <- map(re_mrg, as.data.frame)
map(re_mrg, function (df) plot(df$x, df$y))

## Add drug-class effect
modelcompareinla_nested <- "
model{
  for(j in 1:3) {
    for(i in 1:4) {
        y[i,j] ~ dnorm(trial_effect[i,j], y_prec[i,j])
        trial_effect[i,j]  ~ dnorm(class_effect[j],  between_trial_prec)
    }
  class_effect[j]  ~ dnorm(group_effect, between_class_prec) ## drug-class effect
  }
  group_effect  ~ dnorm(0, 0.25) # overall effect
  between_trial_sd ~ dnorm(0, 1)T(0,)
  between_trial_prec <- 1/between_trial_sd^2

  between_class_sd ~ dnorm(0, 0.01)T(0,)
  between_class_prec <- 1/between_class_sd^2

}# model
"

writeLines(modelcompareinla_nested, con= "jags/modelcompareinla_nested.txt")

## Reshape data
MakeMatrix <- function (var){
  mtrx <- tapply(var, my_data$myclass_act, identity)
  do.call(cbind, mtrx)
}
y_mtrx <- MakeMatrix(my_data$y)
y_prec_mtrx <- MakeMatrix(my_data$y_prec)

# Run model
system.time({
jags_nested <- autorun.jags('jags/modelcompareinla_nested.txt',
                     data = list(y = y_mtrx, y_prec = y_prec_mtrx),
                      n.chains = 2,
                     adapt = 1000,
                     startsample = 10000,
                   monitor = c("trial_effect", "class_effect", "group_effect",
                               "between_trial_sd", "between_class_sd"))
})
## Runs in 2.74 seconds

summary(jags_nested)

## INLA model
myform_nested <- y ~ -1 + myclass + 
  f(group, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1)))) +
  f(myclass_act, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 0.01))))

# Set-up liner combinations
dc1 <- inla.make.lincomb(myclass = 1, myclass_act = c(1, NA, NA, NA)) 
dc2 <- inla.make.lincomb(myclass = 1, myclass_act = c(NA, 1, NA, NA)) 
dc3 <- inla.make.lincomb(myclass = 1, myclass_act = c(NA, NA, 1, NA)) 
dc_all <- c(dc1, dc2, dc3)
names(dc_all) <- paste0("dc", 1:3)

system.time({
mod1_nested <- inla(myform_nested, 
             data = my_data,
             # Add linear combinations to estimate drug-class
             lincomb = dc_all,
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
})
# Runs in 0.18 seconds

summary(mod1_nested)
summary(jags_nested)["group_effect", c("Mean", "SD")]

# Compare drug-class specific effects
class_inla_mean_manual <- mod1_nested$summary.random$myclass_act[,"mean"] + 
  mod1_nested$summary.fixed[,"mean"]
class_inla_mean_derived <- mod1_nested$summary.lincomb.derived$mean
class_inla_sd <- mod1_nested$summary.lincomb.derived$sd

class_jags_mean <- summary(jags_nested)[13:15, c("Mean")]
class_jags_sd <- summary(jags_nested)[13:15, c("SD")]

round(cbind(class_inla_mean_manual, class_inla_mean_manual, class_jags_mean),3)
round(cbind(class_inla_sd, class_jags_sd),3)