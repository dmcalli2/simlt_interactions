##08_analyse_with_INLA

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

my_data <- data.frame(y = dep, y_prec = dep_prec, group = res$study_id_n, myclass = 1, myclass_act = res$my_drug_n)


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

myform1 <- y ~ -1 + myclass + f(group, model = "iid", 
                     hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 0.25))))

mod1 <- inla(myform1, 
             data = my_data,
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
# mod1_inc_acc <- inla.hyperpar(mod1) # increase precision of hyperparameter doesn't do much

### Compare JAGS and INLA model
summary(mod1)
summary(jags)["mu_mu",]

# Compare trial-specific effects
trial_inla_mean <- mod1$summary.random[[1]][,"mean"] + mod1$summary.fixed[,"mean"]
trial_inla_sd <- mod1$summary.random[[1]][, "sd"]
trial_jags_mean <- summary(jags)[1:12, c("Mean")]
trial_jags_sd <- summary(jags)[1:12, c("SD")]

round(cbind(trial_inla_mean, trial_jags_mean),3)
round(cbind(trial_inla_sd, trial_jags_sd),3)

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
  between_trial_sd ~ dnorm(0, 0.001)T(0,)
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
system.time({
myform_nested <- y ~ -1 + myclass + 
  f(group, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 0.001)))) +
  f(myclass_act, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 0.01))))

mod1_nested <- inla(myform_nested, 
             data = my_data,
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
class_inla_mean <- mod1_nested$summary.random$myclass_act[,"mean"] + 
  mod1_nested$summary.fixed[,"mean"]
class_inla_sd <- mod1_nested$summary.random$myclass_act[, "sd"]

class_jags_mean <- summary(jags_nested)[13:15, c("Mean")]
class_jags_sd <- summary(jags_nested)[13:15, c("SD")]

round(cbind(class_inla_mean, class_jags_mean),3)
round(cbind(class_inla_sd, class_jags_sd),3)


# Compare trial specific effects
trial_inla_mean <- mod1_nested$summary.random$group[,"mean"] + # trial effect
  rep(mod1_nested$summary.random$myclass_act[,"mean"], each = 4) + # class effect
  mod1_nested$summary.fixed[,"mean"] # group effect

trial_jags_mean <- summary(jags_nested)[1:12, c("Mean")]
trial_jags_sd <- summary(jags_nested)[1:12, c("SD")]

round(cbind(trial_inla_mean, trial_jags_mean),3)

