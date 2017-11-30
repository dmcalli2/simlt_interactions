# 06_Show_impact_simulations
library(metRology)
library(tidyverse)
library(INLA)
library(stringr)

mean_scenario <- readRDS("scratch_data/mean_scenario")

mean_scenario <- mean_scenario %>% 
  filter(scenario %in% c("scenario_atc5_0.05_trial_0.05_drug_0.15.rds",
                         "scenario_atc5_0.15_trial_0.15_drug_0.15.rds",
                         "scenario_atc5_0.25_trial_0.25_drug_0.15.rds" )) %>% 
  mutate(scenario = str_sub(scenario, 10, -5))
# Identify mean for each scenarios

#02b_run_inla_model
library(INLA)
# INLA:::inla.dynload.workaround() 
load(file = "data/for_inla.Rdata")

## Loop through mean for 3 scenarios
mdls <- map(1:3, function (i) {
  scenario <- mean_scenario$scenario[i]
  iter <- mean_scenario$iter[i]
  
  # Add in main effect to a chosen variation scenario and select relevant classes
  diabetes$res <- res[, scenario] + -0.1
  a10bx <- diabetes %>% 
    filter(atc_5 == "A10BX")
  no_a10bx <- diabetes
  no_a10bx$res[diabetes$atc_5 %in% c("A10BX")] <- NA
  
  # Make linear combination
  dc1 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, NA, 1)) 
  
  # Select iteration
  ## Add values for specific iteration
  my_data$y <-  no_a10bx$res[no_a10bx$iteration == iter]
  
  ## Run model, trial within drug within ATC5 class within ATC4 class
  mod1_nested2 <- inla(myform_nested2, 
                       data = my_data,
                       # Add linear combinations to estimate drug-class
                       lincomb = dc1,
                       # Likelihood distribution
                       family = "gaussian",
                       # Fix likelihood hyperpars as the data is fixed with known precision.
                       control.family = list(hyper = list(prec = list(fixed = TRUE, initial = 0))),
                       # Likelihood precisions
                       scale = my_data$y_prec,
                       # Prior distribution for "fixed" effects - really for mu_mu
                       control.fixed = list(mean = 0, prec = 0.25),
                       # Optionally compute DIC
                       verbose = FALSE,
                       control.compute = list(config=TRUE),
                       control.inla = list(lincomb.derived.only=FALSE))
  # Return glinides and model
  list(a10bx = a10bx, mdl = mod1_nested2)
})
mdls <- transpose(mdls)
a10bx <- mdls$a10bx
mdls <- mdls$mdl
# Storing iterations in a single list
# scenario[[iter]] <- summary(mod1_nested2)
map(mdls, summary)
mdl_res <- map(mdls, ~ .x$summary.lincomb)
map(mdls, ~  plot(.x$marginals.lincomb[[1]], ylim = c(0,5), xlim = c(-2,2)))

## T distribution fits data well (as per BUGS)
mdl_t <- map(mdls, function(mdl_each){
  # Dataframe of values and densities
  res <- mdl_each$marginals.lincomb[[1]] %>% 
    as.data.frame() %>% 
    as_tibble()
  
  # Extract results
  mdl_res_each <- mdl_each$summary.lincomb
  # Fit data, using mean and sd and df 3 for initial values
  a <- nls(formula = y ~ metRology::dt.scaled(x, df = d, mean = m, sd = s),
         data = res, start = list(d = 3, m = mdl_res_each$mean, s = mdl_res_each$sd/2),
         algorithm = "port", lower = list(d = 2, m = -5, s = 0))
  a <- summary(a)
  a <- a$coefficients[,"Estimate"]
  
  # Add estimates to dataframe for comparison
  res$y_new <- metRology::dt.scaled(res$x, df = a["d"], mean = a["m"], sd = a["s"])

  # Plot and return parameters and data
  plot(res$x, res$y, main = paste(c("df", "m", "s"), round(a,2), sep = " = ", collapse = ", "))
  lines(res$x, res$y_new, col = "red")
  
  # Make new dataset with more points
  res_new <- tibble(x = seq(-2, 2, 0.005))
  res_new <- res_new %>% 
    mutate(y = metRology::dt.scaled(x, df = a["d"], mean = a["m"], sd = a["s"]))
  # Parameters and data
  list(param = a, res = res_new)
})

names(mdl_t) <- mean_scenario$scenario
mdl_t <- transpose(mdl_t)
mdl_t_g <- bind_rows(mdl_t$res, .id = "scenario") %>% 
  mutate(variation = factor(scenario, levels = mean_scenario$scenario,
                         labels = c("Small", "Moderate", "Large")))

plot1 <- ggplot(mdl_t_g, aes(x = x, y = y, colour = variation)) +
  geom_line() +
  scale_x_continuous("Treatment-covariate interaction", limits = c(-1, 1)) +
  scale_y_continuous("Density")
plot1
saveRDS(mdl_t, file = "scratch_data/Priors_for_examining_class.Rds")

## Random 