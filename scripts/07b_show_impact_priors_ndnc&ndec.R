#07_show_impact_prior
library(tidyverse)
library(coda)
library(ggplot2)
library(rjags)
library(stringr)
load.module("glm")

my_priors_ndnc <- readRDS(file = "scratch_data/Priors_for_newdrug_newclass.Rds")
my_priors_ndec <- readRDS(file = "scratch_data/Priors_for_newdrug_extclass.Rds")

load(file = "data/sim1/std/for_inla.Rdata")
load("Data/metadata_for_simulation.Rdata")
# Remove single alpha-glucosidase inhibitor trial, already out of diabetes file
diabetes_final <-  
  subset(diabetes_final, diabetes_final$atc_5 != "A10BF")
diabetes_final <- diabetes_final[with(diabetes_final, order(atc_5, drug, nct_id)),]

diabetes_orig <- diabetes%>% 
   arrange(iteration)
diabetes_orig <- cbind(my_data, diabetes_orig)

#append 3 new trials of a new drug in a new class

new_d <- as.data.frame(matrix(nrow=3,data=c(mean(my_data$y_prec),162,1,8,25,
                                            mean(my_data$y_prec),163,1,8,25,
                                            mean(my_data$y_prec),164,1,8,25), byrow=TRUE))
colnames(new_d) <- names(my_data)
my_data1 <- my_data %>%
  bind_rows(new_d) 

ndnc <- as.data.frame(matrix(ncol=4,data=c(rep('A10BN', 3000),
                                           rep('newdrug',3000),
                                           c(rep('NCT99999991',1000),rep('NCT99999992',1000),rep('NCT99999993',1000)),
                                           rep(seq(1,1000,1),3))), stringsAsFactors=FALSE)
colnames(ndnc) <- names(diabetes)
ndnc$iteration <- as.numeric(ndnc$iteration)
diabetes_ndnc <- diabetes %>%
  bind_rows(ndnc)%>%
  arrange(iteration)
diabetes_ndnc <- cbind(my_data1,diabetes_ndnc)

#append 3 new trials of a new drug in an exisitng class

new_d2 <- as.data.frame(matrix(nrow=3, ncol=5,data=c(mean(my_data$y_prec),162,1,3,25,
                                            mean(my_data$y_prec),163,1,3,25,
                                            mean(my_data$y_prec),164,1,3,25), byrow=TRUE))
colnames(new_d2) <- names(my_data)
my_data2 <- my_data %>%
  bind_rows(new_d2) 
ndec <- as.data.frame(matrix(ncol=4,data=c(rep('A10BG', 3000),
                                           rep('newdrug',3000),
                                           c(rep('NCT99999991',1000),rep('NCT99999992',1000),rep('NCT99999993',1000)),
                                           rep(seq(1,1000,1),3))), stringsAsFactors=FALSE)
colnames(ndec) <- names(diabetes)
ndec$iteration <- as.numeric(ndec$iteration)
diabetes_ndec <- diabetes %>%
  bind_rows(ndec) %>%
  arrange(iteration)
diabetes_ndec <- cbind(my_data2,diabetes_ndec)

### New Drug New Class

## Select data and simulate effect for atc5 class, with small between drug/trial variation
a10bn <- diabetes_ndnc %>% 
  filter(atc_5 == "A10BN") %>% 
  as_tibble() 
# Simulate class effect
class_effect <- tibble(iteration = 1:1000, class_effect = rnorm(1000, 0, 0.10))
a10bn <- a10bn %>%
  inner_join(class_effect)
set.seed(1234)
a10bn <- a10bn %>% 
  mutate(y = rnorm(3000, class_effect, 0.05) + rnorm(1000, 0, 0.15))

## examine weighted mean
## Find that quite large treatment effects evident 
a10bn <- a10bn %>% 
  group_by(iteration) %>% 
  mutate(y_crude = weighted.mean(y, 1/y_prec)) %>% 
  ungroup() %>% 
  arrange(y_crude, iteration, nct_id)
hist(a10bn$y_crude)
quantile(a10bn$y_crude, probs = c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95))
# 20 % of samples, 10% either side

# Add vague non-informative prior to existing priors, remove moderate prior
# and set both informative priors to mean = -0.1, for simplicity for plot -- this may not be advisable for new drug existing class
my_wdg_priors <- lapply(my_priors_ndnc$param, `[`, 1:2) 
my_dc_priors <- lapply(my_priors_ndnc$param, `[`, 3:4) 
my_drug_priors <- lapply(my_priors_ndnc$param, `[`, 5:6) 
my_priors_all <- list(my_wdg_priors, my_dc_priors, my_drug_priors)
priors <- list()
for (i in c(1:3)){
  my_priors <- map(my_priors_all[[i]], function(x) {
  x["m"] <- -0.1
  x["d"] <- 3
  if (!is.na(x["s"])) { x}
  })
  priors[[i]] <- my_priors
  }
my_priors <- priors
names(my_priors) <- c('WDG', 'nDC', 'nDg')

## A helper function that tests whether an object is either NULL _or_ 
## a list of NULLs
is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))
## Recursively step down into list, removing all such objects 
rmNullObs <- function(x) {
  x <- Filter(Negate(is.NullOb), x)
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
}
my_priors <- rmNullObs(my_priors)
my_priors <- unlist(my_priors, recursive = FALSE)
my_priors$vague__Full <- c(m = 0, s = 1, d = 3)
my_priors_names <-names(my_priors)
# This leaves 19 priors across 4 list elements:
#  - 6 at the WDG level (3 scenarios [hi,Me,Lo vrn] in 2 models [Full, Full no DC])
#  - 6 at the DC level (3 scenarios [hi,Me,Lo vrn] in 2 models [although in the no DC model these are just WDG priors])
#  - 6 at the drug level (3 scenarios [hi,Me,Lo vrn] in 2 models [Full, Full no DC])
#  - 1 vague

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
  a10bn %>%
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
                       data = c(df_chosen, prior_slct),
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
         prior_type_level = paste0(str_sub(ind, 5,7),str_sub(ind, 12)),
         prior_type_vrn = str_sub(ind, 9,10))%>% 
  mutate(iteration = factor(iteration, levels = c("non","q05", "q50", "q95"),
                            labels = c("Prior only","Worse", "Neutral", "Better")),
         prior_type_level = factor(prior_type_level, levels = c("vagFull", "WDGFull_noDCinfo","WDGFull", "nDCFull", "nDgFull", "nDgFull_noDCinfo"),
                        labels = c("Non-inf","WDG level - standard model","WDG level - DC model","New DC level - DC model","New Drug level - DC model","drop" )),
         prior_type_vrn = factor(prior_type_vrn, levels = c("e_", "Lo", "Me", "Hi"),
                                   labels = c("Non-inf","Strong","Medium","Weak" ))) %>% 
  select(-ind) %>% 
  as_tibble()

res_plt <- res_plt %>% 
  filter( prior_type_level != 'drop')

res_plt <- res_plt %>% 
  group_by(iteration, prior_type_level,prior_type_vrn) %>% 
  sample_n(10000)

res_plt_sum <- res_plt %>% 
  group_by(iteration, prior_type_level,prior_type_vrn) %>%
  summarise(mean = mean(values))

effects <- map(df_choose, function(x){
   colMeans(x['y'])
  })

hline_dat <- data.frame(iteration=levels(res_plt$iteration), hl=c(effects$non,
                                                                  effects$q05,
                                                                  effects$q50,
                                                                  effects$q95)) 

hline_dat2 <- data.frame(iteration=levels(res_plt$iteration), hl=c(rep(-0.1,4))) 

save(res_plt,hline_dat, file = "scratch_data/new_drug_new_class.Rdata")

dodge <- position_dodge(width=0.6)
plot_impact <- ggplot(res_plt, aes(x = prior_type_vrn, y = values, fill = prior_type_level)) + 
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975), adjust = 1, position = dodge,width = 1.7) +
  facet_grid(.~ iteration  , scales = "free_y") +
  scale_y_continuous("Treatment-covariate interaction") +
  scale_x_discrete("Prior used",expand = c(0,0)) +
  geom_hline(data=hline_dat2, aes(yintercept=hl, linetype = "At WDG level"), color = "grey80")+
  geom_hline(data=hline_dat, aes(yintercept=hl, linetype = "In new drug class"), color = "red")+
      theme(axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)),
        text =element_text(size = 14.5),
        panel.background = element_rect(fill = "white", colour = "grey80")) +
  scale_fill_discrete() +
  scale_linetype_manual(name = "'True' effect \nobserved in data", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = c( "grey80","red")))) + 
  coord_cartesian(ylim = c(-0.5, 0.5))
plot_impact

tiff("figures/Impact_of_priors5.tiff", res = 600, compression = "lzw", unit = "in",
     height = 8, width = 8)
plot_impact
dev.off()
