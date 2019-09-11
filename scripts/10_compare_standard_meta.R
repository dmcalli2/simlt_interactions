# 09_one_drug

# install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(tidyverse)
library(INLA)

## Read in diabetes trials
load("Data/metadata_for_simulation.Rdata")
mean_scenarios <- readRDS("scratch_data/mean_scenario")
key_scens <- readRDS("scratch_data/key_scens")
scenario_res <- readRDS("scratch_data/scenario_res") %>% 
  filter(como_prev == "std",
         sim ==1)
load(file = "data/sim1/std/for_inla.Rdata")

# Remove single alpha-glucosidase inhibitor trial, already out of diabetes file
# Remove single alpha-glucosidase inhibitor trial, already out of diabetes file
diabetes_final <-  
  subset(diabetes_final, diabetes_final$atc_5 != "A10BF")
diabetes_final <- diabetes_final[with(diabetes_final, order(atc_5, drug, nct_id)),]

# Ascert
# Create effect estimates interactions with variation at trial, drug and class level
# Read in dataset with different interactions
# Each variable is for a different component
diabetes <- readRDS("scratch_data/interactn_opts.Rds")
diabetes <- as.data.frame(diabetes)


## Select only one drug saxagliptin in class A10BH

whichdata <- cbind(diabetes_final,my_data)

sel_drug_mydata <- whichdata %>%
  filter(drug == "saxagliptin") %>%
  select(names(my_data))

set.seed(2345)

n_its <- 15# change based on how many iterations you want

scenarios <- c('atc5_0.05_trial_0.05_drug_0.05','atc5_0.05_trial_0.25_drug_0.05','atc5_0.05_trial_0.25_drug_0.25')

mean_scenario_s1 <- scenario_res %>%
  filter(sim==1)%>% 
  mutate(scenario = str_sub(scenario, 4, -1L)) %>%
  arrange(scenario)

scenarios_random_iterations <- mean_scenario_s1 %>% 
  mutate(row = row_number()) %>%
  filter(scenario %in% scenarios) %>%
  group_by(scenario) %>% 
  sample_n(n_its) %>%
  ungroup() %>%
  pull(row)

## Loop through mean for 2 scenarios, single iteration for each
mdls <- map(scenarios_random_iterations, function (i) {
  #mdls <- map(c(sample(1:1000,n_its),sample(16001:17000,n_its),sample(26001:27000,n_its)), function (i) {
  scenario <- mean_scenario_s1$scenario[i]
  iter <- mean_scenario_s1$iter[i]
  
  # Add in main effect to a chosen variation scenario and select relevant classes
  diabetes$res <- res[, scenario] + -0.1
 
  # Select iteration
  ## Add values for specific iteration
  my_data$y <-  diabetes$res[diabetes$iteration == iter]
  
  # Make linear combination
  
  dc1 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA)) # DC level
  
  d1 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,20),1,NA,NA,NA) ) # drug level
  
  tr1 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,20),1,NA,NA,NA),
                           trial =c(rep(NA,10),1,rep(NA,161-11))) # trial level
  tr2 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,20),1,NA,NA,NA),
                           trial =c(rep(NA,16),1,rep(NA,161-17))) # trial level
  tr3 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,20),1,NA,NA,NA),
                           trial =c(rep(NA,35),1,rep(NA,161-36))) # trial level
  tr4 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,20),1,NA,NA,NA),
                           trial =c(rep(NA,42),1,rep(NA,161-43))) # trial level
  tr5 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,20),1,NA,NA,NA),
                           trial =c(rep(NA,48),1,rep(NA,161-49))) # trial level
  tr6 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,20),1,NA,NA,NA),
                           trial =c(rep(NA,55),1,rep(NA,161-56))) # trial level
  tr7 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,20),1,NA,NA,NA),
                           trial =c(rep(NA,89),1,rep(NA,161-90))) # trial level
  tr8 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,20),1,NA,NA,NA),
                           trial =c(rep(NA,130),1,rep(NA,161-131))) # trial level
  tr9 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,20),1,NA,NA,NA),
                           trial =c(rep(NA,148),1,rep(NA,161-149))) # trial level
  tr10 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,20),1,NA,NA,NA),
                           trial =c(rep(NA,151),1,rep(NA,161-152))) # trial level
  tr11 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,20),1,NA,NA,NA),
                           trial =c(rep(NA,154),1,rep(NA,161-155))) # trial level
  
  # f(myatc5, model = "iid", 
  # hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1))))
  myform_nested2 <- y ~ -1 + myatc4 + 
    f(trial, model = "iid", 
      hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1)))) +
    f(mydrug, model = "iid", 
      hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1)))) +
    f(myatc5, model = "iid", 
      hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1))))
  
  ## Run model, trial within drug within ATC5 class within ATC4 class
  mod1_nested2 <- inla(myform_nested2, 
                       data = my_data,
                       #data = my_data_test,            #Use instead if testing 
                       # Add linear combinations to estimate drug-class
                       #lincomb = dc4,
                       lincomb = c(a = dc1,b = d1,
                                   c=tr1,d=tr2,e=tr3,f=tr4,g=tr5,
                                   h=tr6,i=tr7,j=tr8,k=tr9,l=tr10,m=tr11),
                       #lincom = do.call(c, dc_tests),  #Use instead if testing 
                       # Likelihood distribution
                       family = "gaussian",
                       # Fix likelihood hyperpars as the data is fixed with known precision.
                       control.family = list(hyper = list(prec = list(fixed = TRUE, initial = 0))),
                       # Likelihood precisions
                       # scale = my_data_test$y_prec,   #Use instead if testing 
                       scale = my_data$y_prec,
                       # Prior distribution for "fixed" effects - really for mu_mu
                       control.fixed = list(mean = 0, prec = 0.25),
                       # Optionally compute DIC
                       verbose = FALSE,
                       control.compute = list(config=TRUE),
                       control.inla = list(lincomb.derived.only=FALSE))

  
  sel_drug_diabetes <- diabetes %>% 
    filter(drug == "saxagliptin")
  
  # Select iteration
  ## Add values for specific iteration
  sel_drug_mydata$y <-  sel_drug_diabetes$res[sel_drug_diabetes$iteration == iter]
  sel_drug_mydata$mydrug <-  1
  
  # Make linear combination
  
  tr1 <- inla.make.lincomb(mydrug = 1, 
                           trial =c(1,rep(NA,10))) # trial level
  tr2 <- inla.make.lincomb(mydrug = 1,
                           trial =c(rep(NA,1),1,rep(NA,9))) # trial level
  tr3 <- inla.make.lincomb(mydrug = 1, 
                           trial =c(rep(NA,2),1,rep(NA,8))) # trial level
  tr4 <- inla.make.lincomb(mydrug = 1, 
                           trial =c(rep(NA,3),1,rep(NA,7))) # trial level
  tr5 <- inla.make.lincomb(mydrug = 1, 
                           trial =c(rep(NA,4),1,rep(NA,6))) # trial level
  tr6 <- inla.make.lincomb(mydrug = 1,
                           trial =c(rep(NA,5),1,rep(NA,5))) # trial level
  tr7 <- inla.make.lincomb(mydrug = 1,
                           trial =c(rep(NA,6),1,rep(NA,4))) # trial level
  tr8 <- inla.make.lincomb(mydrug = 1, 
                           trial =c(rep(NA,7),1,rep(NA,3))) # trial level
  tr9 <- inla.make.lincomb(mydrug = 1,
                           trial =c(rep(NA,8),1,rep(NA,2))) # trial level
  tr10 <- inla.make.lincomb(mydrug = 1, 
                            trial =c(rep(NA,9),1,rep(NA,1))) # trial level
  tr11 <- inla.make.lincomb(mydrug = 1, 
                            trial =c(rep(NA,10),1)) # trial level
  
  ## Write model
  myform_nested3 <- y ~ -1 + mydrug + 
    f(trial, model = "iid", 
      hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 0.1)))) 
  
  ## Run model, trial within drug 
  mod1_nested3 <- inla(myform_nested3, 
                       data = sel_drug_mydata,
                       #data = my_data_test,            #Use instead if testing 
                       # Add linear combinations to estimate drug-class
                       #lincomb = dc4,
                       lincomb = c(c=tr1,d=tr2,e=tr3,f=tr4,g=tr5,
                                   h=tr6,i=tr7,j=tr8,k=tr9,l=tr10,m=tr11),
                       #lincom = do.call(c, dc_tests),  #Use instead if testing 
                       # Likelihood distribution
                       family = "gaussian",
                       # Fix likelihood hyperpars as the data is fixed with known precision.
                       control.family = list(hyper = list(prec = list(fixed = TRUE, initial = 0))),
                       # Likelihood precisions
                       # scale = my_data_test$y_prec,   #Use instead if testing 
                       scale = sel_drug_mydata$y_prec,
                       # Prior distribution for "fixed" effects - really for mu_mu
                       control.fixed = list(mean = 0, prec = 0.25),
                       # Optionally compute DIC
                       verbose = FALSE,
                       control.compute = list(config=TRUE),
                       control.inla = list(lincomb.derived.only=FALSE))
  
  
    # Return glinides and model
  list(mdl = mod1_nested2, mdl2 = mod1_nested3, scenario= scenario, iter=iter)
})

mdls <- transpose(mdls)
#a10bx <- as.data.frame(mdls$a10bx)
ydens <- mdls$ydens
scns <- mdls$scenario
itrs <- rep(seq(1, n_its, 1), times= length(unique(unlist(scns))))
mdl1 <- mdls$mdl
mdl2 <- mdls$mdl2

# The aims are 1) to compare standard meta-analysis of a single and hierarchcal meta-analyses in terms of 
# the average trial- and drug-level precision across iterations, within-scenarios and 2) to show this 
# comparison graphically for a specific iteration of three scenarios

drug_std <- data.frame()
trial_std <-data.frame()
for( i in 1:length(itrs)){
  drug_std <- rbind(drug_std, data.frame(mdl2[[i]]$summar.fixed))
  # Gets us the drug level effect from the drug only model
  trial_std <- rbind(trial_std, data.frame(mdl2[[i]]$summary.lincomb))
  # Gets us the trial level effect from the drug only model
}

drug_full <- data.frame()
trial_full <-data.frame()
for( i in 1:length(itrs)){
  drug_full <- rbind(drug_full, data.frame(mdl1[[i]]$summary.lincomb))
  # Gets us the drug level effect from the drug only model
  trial_full <- rbind(trial_full, data.frame(mdl1[[i]]$summary.lincomb))
  # Gets us the trial level effect from the drug only model
}


map(mdl1, summary)
mdl_res <- map(mdl1, ~ .x$summary.lincomb)
fixeds
mdl_res

mdl1_t <- map(mdl1, function(mdl_each){
  # Dataframe of values and densities
  res <- mdl_each$marginals.lincomb %>% 
    as.data.frame() %>% 
    as_tibble()
  # mdl_res_each <- mdl_each$summary.lincomb
  #list(mdl_t,mdl_res_each)
}) # abbreviated this function here to view distributions of marginals

mdl2_t <- map(mdl2, function(mdl_each){
  # Dataframe of values and densities
  res <- mdl_each$marginals.fixed %>% 
    as.data.frame() %>% 
    as_tibble()
  # mdl_res_each <- mdl_each$summary.lincomb
  #list(mdl_t,mdl_res_each)
}) # abbreviated this function here to view distributions of marginals

mdl_t_tdy <- bind_rows(mdl1_t) %>%
  mutate(scenario= rep(unlist(scns), each=75),
         iter= rep(itrs,each = 75)) %>%
  cbind(mdl2_t %>% bind_rows()) %>%
  rename(class.x = a.lc.x, class.y = a.lc.y,
         drug.x = b.lc.x, drug.y = b.lc.y) 

mdl_t_tdy2 <- mdl_t_tdy %>%
  select(contains('x'),scenario,iter)%>%
  gather(key = xtype, value = x, -iter,-scenario) %>%
  mutate(group =  str_sub(xtype, end=-3)) %>%
  bind_cols(mdl_t_tdy %>%
              select(contains('.y'),scenario,iter)%>%
              gather(key = ytype, value = y,-iter,-scenario) %>%
              select(ytype, y)) %>%
  mutate(myalpha = ifelse(group == "class", 0.7, 0.9),
         Effect = factor(case_when(group == "class" ~ "Drug class \n(full model)",
                                   group == "drug" ~ "Drug \n(full model)",
                                   group == "mydrug" ~ "Drug \n(drug only model)"), levels = c("Drug class \n(full model)",
                                                                                                   "Drug \n(full model)",
                                                                                                   "Drug \n(drug only model)")),
         Scenario = factor(case_when(scenario == "atc5_0.05_trial_0.05_drug_0.05" ~ "Low variation",
                                     scenario == "atc5_0.05_trial_0.25_drug_0.05" ~ "High between-trial \nvariation",
                                     scenario == "atc5_0.05_trial_0.25_drug_0.25" ~ "High between-trial, \nbetween-drug variation"), 
                           levels = c("Low variation",
                                      "High between-trial \nvariation",
                                      "High between-trial, \nbetween-drug variation")))



gg <- ggplot(mdl_t_tdy2,aes(x=x, y=y)) +
  geom_area(aes(group= Effect, fill=Effect)) + 
  facet_grid( iter ~ Scenario, scales = "fixed") + 
  scale_fill_manual(name="Effect at level of:",values = alpha(c("purple","blue","forestgreen"), c(0.3,0.7,0.7)))+
  theme_bw() +
  theme(axis.text.x = element_text(angle=0, vjust = 0, size = 11),
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text =element_text(size = 14.5),
        panel.background = element_rect(fill = "white", colour = "white"),
        strip.background.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        strip.background.x= element_rect(fill = "white", colour = "white"),
        panel.border = element_blank()) +
  scale_x_continuous(limits=c(-0.5, 0.5), breaks = c(-0.3,0,0.3)) + 
  ylim(0,17) + 
  xlab("Estimate of interaction effect (SDs)") +
  ylab("Posterior density")

tiff("figures/show_effect_drug_class_drug_level.tiff", res = 600, compression = "lzw", unit = "in",
     height = 10, width =8)
gg  
dev.off()


