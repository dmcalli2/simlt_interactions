# 09_one_drug

# install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

#install.packages("ggstance")
library(tidyverse)
library(ggstance)
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
  group_split(mydrug, keep=TRUE) 

set.seed(2345)

n_its <- 100# change based on how many iterations you want

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
  

  d1 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, 1, NA, NA), mydrug = c(1,rep(NA,24-1))) # drug level
  d2 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,1),1,rep(NA,24-2))) # drug level
  d3 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, 1, NA), mydrug = c(rep(NA,2),1,rep(NA,24-3))) # drug level
  d4 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, 1, NA), mydrug = c(rep(NA,3),1,rep(NA,24-4))) # drug level
  d5 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, 1, NA, NA), mydrug = c(rep(NA,4),1,rep(NA,24-5))) # drug level
  d6 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, 1, NA), mydrug = c(rep(NA,5),1,rep(NA,24-6))) # drug level
  d7 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, 1, NA, NA), mydrug = c(rep(NA,6),1,rep(NA,24-7))) # drug level
  d8 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,7),1,rep(NA,24-8))) # drug level
  d9 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, 1, NA, NA, NA, NA, NA), mydrug = c(rep(NA,8),1,rep(NA,24-9))) # drug level
  d10 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, 1, NA, NA, NA, NA, NA), mydrug = c(rep(NA,9),1,rep(NA,24-10))) # drug level
  d11 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,10),1,rep(NA,24-11))) # drug level
  d12 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, 1, NA, NA), mydrug = c(rep(NA,11),1,rep(NA,24-12))) # drug level
  d13 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, 1, NA, NA), mydrug = c(rep(NA,12),1,rep(NA,24-13))) # drug level
  d14 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(1, NA, NA, NA, NA, NA, NA), mydrug = c(rep(NA,13),1,rep(NA,24-14))) # drug level
  d15 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, NA, 1), mydrug = c(rep(NA,14),1,rep(NA,24-15))) # drug level
  d16 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, NA, 1), mydrug = c(rep(NA,15),1,rep(NA,24-16))) # drug level
  d17 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, 1, NA, NA, NA, NA), mydrug = c(rep(NA,16),1,rep(NA,24-17))) # drug level
  d18 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, NA, 1), mydrug = c(rep(NA,17),1,rep(NA,24-18))) # drug level
  d19 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, 1, NA, NA, NA, NA), mydrug = c(rep(NA,18),1,rep(NA,24-19))) # drug level
  d20 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, 1, NA, NA, NA, NA), mydrug = c(rep(NA,19),1,rep(NA,24-20))) # drug level
  d21 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,20),1,rep(NA,24-21))) # drug level
  d22 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,21),1,rep(NA,24-22))) # drug level
  d23 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, 1, NA, NA), mydrug = c(rep(NA,22),1,rep(NA,24-23))) # drug level
  d24 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,23),1) )# drug level

  
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
                       lincomb = c(d1= d1,d2= d2, d3= d3, d4= d4, d5= d5, d6= d6,
                                   d7= d7,d8= d8,d9= d9,d10= d10, d11= d11,
                                   d12= d12, d13= d13, d14= d14, d15= d15,
                                   d16= d16, d17= d17, d18= d18, d19= d19,
                                   d20= d20, d21= d21, d22= d22, d23= d23, d24= d24),
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

  drug_models <- list()
  
 for(j in 1:length(unique(my_data$mydrug)))
  {
  
    # Select iteration
    ## Add values for specific iteration
    sel_drug_mydata <-my_data %>% filter(mydrug==j)
 

    ## Write model
    myform_nested3 <- y ~ -1 + myatc4 + 
      f(trial, model = "iid", 
        hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 0.1)))) 
    
    ## Run model, trial within drug 
    mod1_nested3 <- inla(myform_nested3, 
                         data = sel_drug_mydata,
                         #data = my_data_test,            #Use instead if testing 
                         # Add linear combinations to estimate drug-class
                         #lincomb = d1,
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
    
    drug_models[[j]] <- mod1_nested3
    
  }
    # Return glinides and model
  list(mdl = mod1_nested2, mdl2 = drug_models, scenario= scenario, iter=iter)
})

mdls <- transpose(mdls)
#a10bx <- as.data.frame(mdls$a10bx)
ydens <- mdls$ydens
scns <- mdls$scenario
itrs <- rep(seq(1, n_its, 1), times= length(unique(unlist(scns))))
mdl1 <- mdls$mdl
mdl2 <- mdls$mdl2

# The aims are 1) to compare standard meta-analysis of a single and hierarchcal meta-analyses in terms of 
# the average drug-level precision across iterations, within-scenarios and 2) to show this 
# comparison graphically for a specific iteration of three scenarios

drug_std <- data.frame()

for( i in 1:length(itrs)){
  for( j in 1:length(unique(my_data$mydrug))){
    print(c(i,j))
  drug_std <- rbind(drug_std, data.frame(c(mdl2[[i]][[j]]$summary.fixed)))
  # Gets us the drug level effect from the drug only model
  }
}

drug_std <- drug_std %>%
  mutate(scenario= rep(unlist(scns), each=length(unique(my_data$mydrug))),
         iteration = rep(unlist(itrs),each=length(unique(my_data$mydrug))),
         mydrug = rep(seq(1:length(unique(my_data$mydrug))),length(itrs)),
         model= "drug") %>%
  select(mean, sd,lo=X0.025quant, hi=X0.975quant, scenario, iteration, mydrug,model)

drug_full <- data.frame()
for( i in 1:length(itrs)){
  drug_full <- rbind(drug_full, data.frame(mdl1[[i]]$summary.lincomb))
  # Gets us the drug level effect from the full model
}

drug_full <- drug_full %>%
  mutate(scenario= rep(unlist(scns), each=length(unique(my_data$mydrug))),
    iteration = rep(seq(1:length(itrs)),each=length(unique(my_data$mydrug))),
         mydrug = rep(seq(1:length(unique(my_data$mydrug))),length(itrs)),
    model="full")%>%
  select(mean=mean, sd=sd,lo=X0.025quant, hi=X0.975quant, scenario, iteration, mydrug, model)


drug <- drug_full %>%
  bind_rows(drug_std) %>%
  group_by(scenario,mydrug,model) %>%
  summarise(avg_mean = mean(mean),
            avg_lo = mean(lo),
            avg_hi = mean(hi),
            avg_sd = mean(sd),
            sd_sd = sd(sd)) %>%
  left_join(whichdata %>%
              select(-trial, -y_prec,-n_per_grp, -nct_id)) %>%
  distinct() %>%
  mutate(atc_5 = factor(atc_5),
         Model = factor(case_when(model=="drug"~ "Single drug\nmeta-analysis",
                                  model=="full"~ "Full hierarchical\nmodel"), levels = c("Single drug\nmeta-analysis","Full hierarchical\nmodel"))) %>%
  filter(scenario == "atc5_0.05_trial_0.05_drug_0.05") 

a <- my_data %>%
  left_join(whichdata %>% select(trial, n_per_grp)) %>% 
  group_by(mydrug) %>% summarise(n= n(), med_prec = median(y_prec), sum_npg = sum(n_per_grp))

drug <- drug %>% left_join(a) %>% mutate(n_prec = (n/med_prec))

gg <- ggplot(drug,aes(x=avg_mean, y=(drug), colour=Model, alpha=0.9) )+
  geom_point(aes(size=sqrt(sum_npg)/10),position = ggstance::position_dodgev(0.9))+
  geom_errorbarh(data=drug,aes(xmax=avg_hi, xmin=avg_lo),position = position_dodgev(0.9), width=0.4, size=0.9)  +
  scale_color_manual(values=c("forestgreen","blue"))+
  facet_grid(atc_5~., scales = "free", space = "free") +
  coord_cartesian(xlim=c(-0.3,0.1))+
  theme_bw() +
  theme(axis.text.x = element_text(),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(),
        text =element_text(size = 14.5),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.minor = element_blank(),
        strip.background.x= element_rect(fill = "white", colour = "white"),
        strip.text = element_text(size = 9),
        panel.border = element_blank()) +
  guides(colour=guide_legend(reverse = TRUE), alpha=FALSE, size=FALSE)+ 
  xlab("Estimate of interaction effect (SDs)") +
  ylab("Drug")

tiff("figures/show_increase_prec_vs_drug_lo_vrn.tiff", res = 600, compression = "lzw", unit = "in",
     height = 10, width =6.5)
gg  
dev.off()



