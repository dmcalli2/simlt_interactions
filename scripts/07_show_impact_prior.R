#07_show_impact_prior
library(tidyverse)
library(coda)
library(ggplot2)
library(rjags)
library(stringr)
load.module("glm")

my_priors <- readRDS(file = "scratch_data/Priors_for_examining_class_dc3.Rds")
load(file = "data/sim1/std/for_inla.Rdata")

a10bx <- diabetes%>% 
   arrange(iteration)
rm(res, diabetes)
a10bx <- cbind(my_data, a10bx)

## Select data and simulate effect for atc5 class, with small between drug/trial variation
a10bx <- a10bx %>% 
  filter(atc_5 == "A10BX") %>% 
  as_tibble() 
rm(my_data)

## Downweight precision of that one trial, just to demonstrate more realistic effect
# Unusual to ahve a trial with 9000 participants
a10bx$y_prec <- a10bx$y_prec * c(1, (1/3^2), 1)

# Simulate class effect
class_effect <- tibble(iteration = 1:1000, class_effect = rnorm(1000, 0, 0.10))
a10bx <- a10bx %>%
  inner_join(class_effect)
set.seed(1234)
a10bx <- a10bx %>% 
  mutate(y = rnorm(3000, class_effect, 0.05) + rnorm(1000, 0, 0.15))

## examine weighted mean
## Find that quite large treatment effects evident 
a10bx <- a10bx %>% 
  group_by(iteration) %>% 
  mutate(y_crude = weighted.mean(y, 1/y_prec)) %>% 
  ungroup() %>% 
  arrange(y_crude, iteration, nct_id)
hist(a10bx$y_crude)
quantile(a10bx$y_crude, probs = c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95))
# 20 % of samples, 10% either side

# Add vague non-informative prior to existing priors, remove moderate prior
# and set both informative priors to mean = -0.1, for simplicity for plot
my_priors <- my_priors$param
my_priors <- map(my_priors, function(x) {
  x["m"] <- -0.1
  x["d"] <- 3
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
  a10bx %>%
    arrange(y_crude, iteration) %>% 
    group_by(nct_id) %>% 
    slice(round(cutpoint * length(iteration))) %>% 
    ungroup() %>% 
    mutate(mydrug = 1:3)  %>% 
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
                                          "atc5_0.25_trial_0.25_drug_0.25",
                                          "atc5_0.05_trial_0.05_drug_0.05"
                                          ),
                        labels = c("Non-informative",
                                   "Weak","Strong" ))) %>% 
  select(-ind) %>% 
  as_tibble()

res_plt <- res_plt %>% 
  group_by(iteration, prior) %>% 
  sample_n(10000)

res_plt_sum <- res_plt %>% 
  group_by(iteration, prior) %>%
  summarise(mean = mean(values))

hline_dat <- data.frame(iteration=levels(res_plt$iteration), hl=c(NA,-0.21,0.07,0.18)) 

save(res_plt,hline_dat, file = "scratch_data/new_drug_new_class.Rdata")


plot_impact <- ggplot(res_plt, aes(x = prior, y = values, fill = prior)) + 
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975), adjust = 1) +
  facet_wrap(~iteration) +
  scale_y_continuous("Treatment-covariate interaction") +
  scale_x_discrete("Prior used") +
  geom_hline(data=hline_dat, aes(yintercept=hl), linetype = "dashed", color = "red")+
    theme(axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)),
        text =element_text(size = 14.5),
        panel.background = element_rect(fill = "white", colour = "grey80")) +
  scale_fill_discrete(guide = FALSE) +
  coord_cartesian(ylim = c(-0.75, 0.75))
plot_impact

tiff("figures/Impact_of_priors5.tiff", res = 600, compression = "lzw", unit = "in",
     height = 8, width = 8)
plot_impact
dev.off()
