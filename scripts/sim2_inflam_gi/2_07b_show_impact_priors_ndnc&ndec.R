#07_show_impact_prior
library(tidyverse)
library(coda)
library(ggplot2)
library(rjags)
library(stringr)
load.module("glm")

my_priors_ndnc <- readRDS(file = "scratch_data/Priors_for_newdrug_newclass_sim2.Rds")
my_priors_ndec <- readRDS(file = "scratch_data/Priors_for_newdrug_extclass_sim2.Rds")

load(file = "data/sim2/std/for_inla.Rdata")
load("Data/rheum_final.Rdata")
load("scratch_data/Priors_for_newdrug_extclass_sim2_scns_itrs.Rdata")
# Remove single alpha-glucosidase inhibitor trial, already out of rheum file

rheum_orig <- rheum%>% 
   arrange(iteration) 

rheum_orig <- cbind(my_data, rheum_orig)

#append 3 new trials of a new drug in an existing class

new_d <- as.data.frame(matrix(nrow=3,data=c(mean(my_data$y_prec),163,1,6,33,2,
                                            mean(my_data$y_prec),164,1,6,33,2,
                                            mean(my_data$y_prec),165,1,6,33,2), byrow=TRUE))
colnames(new_d) <- names(my_data)
my_data1 <- my_data %>%
  bind_rows(new_d) 

ndec <- as.data.frame(matrix(ncol=5,data=c(rep('Conventional', 3000),
                                           rep('Folic acid analogues', 3000),
                                           rep('newdrug',3000),
                                           c(rep('NCT99999991',1000),rep('NCT99999992',1000),rep('NCT99999993',1000)),
                                           rep(seq(1,1000,1),3))), stringsAsFactors=FALSE)
colnames(ndec) <- names(rheum)
ndec$iteration <- as.numeric(ndec$iteration)
rheum_ndec <- rheum %>%
  bind_rows(ndec)%>%
  arrange(iteration)
rheum_ndec <- cbind(my_data1,rheum_ndec)

### New Drug Existing Class

## Select data and simulate effect for atc5 class, with small between drug/trial variation
newdrgs <- rheum_ndec %>% 
  filter(drug == "newdrug") %>% 
  as_tibble() 
# Simulate class effect - because class is existing, use distribution from prior in each scn
class_effect <- tibble(iteration = 1:1000, 
                       class_effect_lo = rnorm(1000, my_priors_ndec$param$Lo_Full[5], my_priors_ndec$param$Lo_Full[6]),
                       class_effect_me = rnorm(1000, my_priors_ndec$param$Me_Full[5], my_priors_ndec$param$Me_Full[6]),
                       class_effect_hi = rnorm(1000, my_priors_ndec$param$Hi_Full[5], my_priors_ndec$param$Hi_Full[6]))
newdrgs <- newdrgs %>%
  inner_join(class_effect)
set.seed(1234)

newdrgs_lo <- newdrgs %>% 
  mutate(y = rnorm(3000, class_effect_lo, 0.05) + rnorm(1000, 0, 0.05))
newdrgs_me <- newdrgs %>% 
  mutate(y = rnorm(3000, class_effect_me, 0.15) + rnorm(1000, 0, 0.15))
newdrgs_hi <- newdrgs %>% 
  mutate(y = rnorm(3000, class_effect_hi, 0.25) + rnorm(1000, 0, 0.25))
                      

## examine weighted mean
## Find that quite large treatment effects evident 
newdrgs_lo <- newdrgs_lo %>% 
  group_by(iteration) %>% 
  mutate(y_crude = weighted.mean(y, 1/y_prec)) %>% 
  ungroup() %>% 
  arrange(y_crude, iteration, nct_id)
hist(newdrgs_lo$y_crude)
quantile(newdrgs_lo$y_crude, probs = c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95))
# 20 % of samples, 10% either side
newdrgs_me <- newdrgs_me %>% 
  group_by(iteration) %>% 
  mutate(y_crude = weighted.mean(y, 1/y_prec)) %>% 
  ungroup() %>% 
  arrange(y_crude, iteration, nct_id)
hist(newdrgs_me$y_crude)
quantile(newdrgs_me$y_crude, probs = c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95))
# 20 % of samples, 10% either side
newdrgs_hi <- newdrgs_hi %>% 
  group_by(iteration) %>% 
  mutate(y_crude = weighted.mean(y, 1/y_prec)) %>% 
  ungroup() %>% 
  arrange(y_crude, iteration, nct_id)
hist(newdrgs_hi$y_crude)
quantile(newdrgs_hi$y_crude, probs = c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95))
# 20 % of samples, 10% either side

# Add vague non-informative prior to existing priors, remove moderate prior

my_wdg_priors <- lapply(my_priors_ndec$param, `[`, 1:2) 
my_pth_priors <- lapply(my_priors_ndec$param, `[`, 3:4) 
my_dc_priors <- lapply(my_priors_ndec$param, `[`, 5:6) 
my_drg_priors <- lapply(my_priors_ndec$param, `[`, 7:8) 
my_priors_all <- list(my_wdg_priors,my_pth_priors, my_dc_priors, my_drg_priors)
priors <- list()
for (i in c(1:4)){
  my_priors <- map(my_priors_all[[i]], function(x) {
#  x["m"] <- -0.1
  x["d"] <- 3
  if (!is.na(x["s"])) { x}
  })
  priors[[i]] <- my_priors
  }
my_priors <- priors
names(my_priors) <- c('WDG','Pth','nDC', 'nDg')

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
newdrgslst <- list(newdrgs_lo,newdrgs_me,newdrgs_hi)
df_choose <- list(list())
for (i in c(1:3)){
df_choose1 <- map(c(q05 = 0.05, q50 = 0.5, q95 = 0.95), function(cutpoint){
  newdrgslst[[i]] %>%
    arrange(y_crude, iteration) %>% 
    group_by(nct_id) %>% 
    slice(round(cutpoint * length(iteration))) %>% 
    ungroup() %>% 
    mutate(mydrug = 1:3)  %>% 
    select(y_prec, y) 
})
df_choose[i] <- list(df_choose1) 
}
df_choose <- unlist(df_choose, recursive = FALSE)
names(df_choose) <- paste0(c('q05','q50','q95'), rep(c('.lo','.me','.hi'),each=3) )

# Add NAs to generate plot priors (ie no data)
df_choose$non.lo <- df_choose[[1]] %>% 
  mutate(y = NA) %>% 
  as_tibble()
df_choose$non.me <- df_choose[[1]] %>% 
  mutate(y = NA) %>% 
  as_tibble()
df_choose$non.hi <- df_choose[[1]] %>% 
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
res_plt2 <- res_plt %>% 
  mutate(iteration = str_sub(ind, 1, 3),
         scenario = str_sub(ind, 5, 6),
         prior_type_level = paste0(str_sub(ind, 8,10),str_sub(ind, 15)),
         prior_type_vrn = str_sub(ind, 12,13)) %>% 
  mutate(iteration = factor(iteration, levels = c("non","q05", "q50", "q95"),
                            labels = c("Prior only","Worse", "Neutral", "Better")),
         scenario = factor(scenario, levels = c("lo", "me", "hi"),
                            labels = c("Low variation", "Medium variation", "High variation")),
         prior_type_level = factor(prior_type_level, levels = 
                                     c("vagFull", 
                                       "WDGFull_noDC_Path_info","PthFull_noDC_Path_info", "nDCFull_noDC_Path_info", "nDgFull_noDC_Path_info",
                                       "WDGNewDrgDC_only","PthNewDrgDC_only", "nDCNewDrgDC_only", "nDgNewDrgDC_only",
                                       "WDGFull","PthFull", "nDCFull", "nDgFull"),
                        labels = c("Non-inf", 
                                   "WDG level - standard model","drop1", "drop2", "drop3",
                                   "drop4","drop5", "DC level - single DC only model", "Drug level - single DC only model",
                                   "WDG level - DC model","Path level - DC model", "DC level - DC model", "Drug level - DC model")),
         prior_type_vrn = factor(prior_type_vrn, levels = c("e_", "Lo", "Me", "Hi"),
                                   labels = c("Non-inf","Strong","Medium","Weak" ))) %>% 
  select(-ind) %>% 
  as_tibble()

res_plt <- res_plt2 %>% 
  filter(!prior_type_level %in% c(paste0(rep('drop',5),seq(1,5,1))))

res_plt <- res_plt %>% 
  group_by(iteration, scenario, prior_type_level,prior_type_vrn) %>% 
  sample_n(10000)

res_plt_sum <- res_plt %>% 
  group_by(iteration, scenario, prior_type_level,prior_type_vrn) %>%
  summarise(mean = mean(values))

effects <- map(df_choose, function(x){
  colMeans(x['y'])
  }) 
  
hline_dat <- data.frame(iteration=rep(levels(res_plt$iteration),each=3),
                        scenario =rep(levels(res_plt$scenario), 4)) 
hline_dat$iteration <- factor(hline_dat$iteration,levels(res_plt$iteration))
hline_dat$scenario <- factor(hline_dat$scenario,levels(res_plt$scenario))

hline_dat <- hline_dat %>%
  arrange(scenario, iteration)
hline_dat <- hline_dat %>%
  mutate( hl = as.vector(c(effects$non.lo,
                           effects$q05.lo,
                           effects$q50.lo,
                           effects$q95.lo,
                           effects$non.me,
                           effects$q05.me,
                           effects$q50.me,
                           effects$q95.me,
                           effects$non.hi,
                           effects$q05.hi,
                           effects$q50.hi,
                           effects$q95.hi)))


hline_dat2 <- hline_dat %>%
  mutate( hl=c(rep(-0.1,12)))

hline_dat3 <- hline_dat %>%
  mutate( hl=c(rep(my_priors_ndec$param$Lo_Full[5],4),
               rep(my_priors_ndec$param$Me_Full[5], 4),
               rep(my_priors_ndec$param$Hi_Full[5], 4)))

save(res_plt,hline_dat, file = "scratch_data/new_drug_ext_class_sim2.Rdata")

res_plt <- res_plt %>%
  filter( !( scenario =="Low variation" & prior_type_vrn %in% c("Medium","Weak") ),
          !( scenario =="Medium variation" & prior_type_vrn %in% c("Strong","Weak") ),
          !( scenario =="High variation" & prior_type_vrn %in% c("Medium","Strong") ))

dodge <- position_dodge(width=0.8)
plot_impact <- ggplot(res_plt, aes(x = prior_type_vrn, y = values, fill = prior_type_level)) + 
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975), adjust = 0.5, position = dodge,width = 2.7, trim = TRUE) +
  facet_grid(scenario~ iteration  , scales = "free") +
  scale_y_continuous("Treatment-covariate interaction") +
  scale_x_discrete("Prior used",expand = c(0.2,0.2)) +
  geom_hline(data=hline_dat2, aes(yintercept=hl, linetype = "In original sample, at WDG level"), color = "red")+
  geom_hline(data=hline_dat3, aes(yintercept=hl, linetype = "In original sample, at DC level"), color = "forestgreen")+
  geom_hline(data=hline_dat, aes(yintercept=hl, linetype = "In new data, at DC level"), color = "grey80")+
        theme(axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)),
        text =element_text(size = 14.5),
        panel.background = element_rect(fill = "white", colour = "grey80")) +
  scale_fill_discrete() +
  scale_linetype_manual(name = "'True' effect observed", values = c(2, 2, 2), 
                        guide = guide_legend(override.aes = list(color = c( "grey80","forestgreen","red")))) +
  coord_cartesian(ylim = c(-0.75, 0.75))
plot_impact

tiff("figures/Impact_of_priors5.tiff", res = 600, compression = "lzw", unit = "in",
     height = 8, width = 8)
plot_impact
dev.off()
