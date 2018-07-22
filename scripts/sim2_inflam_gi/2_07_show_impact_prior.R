#07_show_impact_prior
library(tidyverse)
library(coda)
library(ggplot2)
library(rjags)
library(stringr)
load.module("glm")

my_priors <- readRDS(file = "scratch_data/Priors_for_examining_class_sim2_sd1.Rds")
load(file = "data/sim2_for_inla.Rdata")

moa_ch <- rheum%>% 
   arrange(iteration)
rm(res, rheum)
moa_ch <- cbind(my_data, moa_ch)

## Select data and simulate effect for MoA class, with small between drug/trial variation
cd20 <- moa_ch %>% 
  filter(moa == "CD20-directed Antibody Interactions") %>% 
  as_tibble() 
rm(my_data)

## Downweight precision of that one trial, just to demonstrate more realistic effect
# Unusual to ahve a trial with 9000 participants
#cd20$y_prec <- cd20$y_prec * c(1, (1/3^2), 1)

# Simulate class effect
class_effect <- tibble(iteration = 1:1000, class_effect = rnorm(1000, 0, 0.10))
cd20 <- cd20 %>%
  inner_join(class_effect)
set.seed(1234)
cd20 <- cd20 %>% 
  mutate(y = rnorm(8000, class_effect, 0.05) + rnorm(1000, 0, 0.15))

## examine weighted mean
## Find that quite large treatment effects evident 
cd20 <- cd20 %>% 
  group_by(iteration) %>% 
  mutate(y_crude = weighted.mean(y, 1/y_prec)) %>% 
  ungroup() %>% 
  arrange(y_crude, iteration, nct_id)
hist(cd20$y_crude)
quantile(cd20$y_crude, probs = c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95))
# 20 % of samples, 10% either side

# Add vague non-informative prior to existing priors, remove moderate prior
# and set both informative priors to mean = -0.1, for simplicity for plot
my_priors <- my_priors$param
my_priors <- map(my_priors, function(x) {
  x["m"] <- -0.1
  x})
my_priors <- my_priors[-2]
my_priors$vague <- c(m = 0, s = 1, d = 3)

modelstring <- "
model{
for (z in 1:3){
   # likelihood
   y[z] ~ dnorm(wd[z], y_prec[z])
   # link and linear predictor
   # Shared priors (random effects)
   wd[z] ~ dnorm(wd_delta, wd_tau)
} #end of trials
   # Meta-analysis level priors
 
   wd_delta ~ dt(m, 1/s^2, d)
   wd_tau <-  1/wd_var
   wd_var ~ dlnorm(-3.23, 1/1.88^2)
   wd_sd <- wd_var^0.5
}# end of model
"
writeLines(modelstring, con = "jags/random_effects_new.txt")


# Select data for single run
df_choose <- map(c(q05 = 0.05, q50 = 0.5, q95 = 0.95), function(cutpoint){
  cd20 %>%
    arrange(y_crude, iteration) %>% 
    group_by(nct_id) %>% 
    slice(round(cutpoint * length(iteration))) %>% 
    ungroup() %>% 
    mutate(mydrug = 1:8)  %>% 
    select(y_prec, y) 
})

# Add NAs to generate plot priors (ie no data)
df_choose$non <- df_choose[[1]] %>% 
  mutate(y = NA) %>% 
  as_tibble()

system.time({
  
res <- map(df_choose, function(df_chosen){
  # Loop through dataframe choices
  map(my_priors, function(prior_slct){
  jags <- jags.model(paste0("jags/random_effects_new.txt"),
                       data = c(df_chosen,
                                prior_slct),
                       n.chains = 2,
                       n.adapt = 1000)
  update(jags, 10000)
  data(LINE)
  LINE$recompile()
  coda.samples(jags,
                c('wd', 'wd_delta', 'wd_sd'),
                40000)
 })
})

})

res2 <- do.call(c, res)

## Summary statistics
smrys <- map(res2, ~ summary(.x))
stats <- map(smrys, ~ .x$statistics["wd_delta",, drop = FALSE] %>%  as_tibble())
stats <- stats %>%  bind_rows(.id = "prior")
quants <- map(smrys, ~ .x$quantiles["wd_delta",, drop = FALSE] %>%  as_tibble())
quants <- quants %>%  bind_rows(.id = "prior")

## Densities for plots
res_plt <- map(res2, ~ .x[,"wd_delta"] %>%  as.matrix(drop = TRUE) %>%  as.vector())
res_plt <- stack(res_plt)
res_plt <- res_plt %>% 
  mutate(iteration = str_sub(ind, 1, 3),
         prior = str_sub(ind, 5)) %>% 
  mutate(iteration = factor(iteration, levels = c("non","q05", "q50", "q95"),
                            labels = c("Prior only","Worse", "Neutral", "Better")),
         prior = factor(prior, levels = c("vague",
                                          "path_0.15_moa_0.25_trial_0.25_drug_0.15",
                                          "path_0.15_moa_0.05_trial_0.05_drug_0.15"
                                          ),
                        labels = c("Non-informative",
                                   "Weak","Strong" ))) %>% 
  select(-ind) %>% 
  as_tibble()

res_plt <- res_plt %>% 
  group_by(iteration, prior) %>% 
  sample_n(10000)

plot_impact <- ggplot(res_plt, aes(x = prior, y = values, fill = prior)) + 
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975), adjust = 1) +
  facet_wrap(~iteration) +
  scale_y_continuous("Treatment-covariate interaction") +
  scale_x_discrete("Prior used") +
  scale_fill_discrete(guide = FALSE) +
  coord_cartesian(ylim = c(-0.75, 0.75))
plot_impact

tiff("figures/Impact_of_priors5_sim2_sd1.tiff", res = 600, compression = "lzw", unit = "in",
     height = 8, width = 8)
plot_impact
dev.off()
