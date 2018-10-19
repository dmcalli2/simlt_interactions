# 06_Show_impact_priors_dc_vs_no_dc_info
library(metRology)
library(tidyverse)
library(INLA)
library(stringr)
library(magrittr)

mean_scenarios <- readRDS("scratch_data/mean_scenario")
key_scens <- readRDS("scratch_data/key_scens")
scenario_res <- readRDS("scratch_data/scenario_res") %>% 
  filter(como_prev == "std")


key_scens_list <-key_scens %>%
  filter(como_prev=="std") %>%
  select(scenario)
dput(as.vector(key_scens_list$scenario))


mean_scenario <- mean_scenarios %>% 
  filter(como_prev=="std") %>%
#  filter(scenario %in% key_scens_list$scenario) %>% 
  mutate(scenario = str_sub(scenario, 4, -1L))
# Identify mean for each scenarios


#Separate for each sim from here

mean_scenario_s1 <- scenario_res %>%
  filter(sim==1)%>% 
  mutate(scenario = str_sub(scenario, 4, -1L)) %>%
  arrange(scenario)


load(file = "data/sim1/std/for_inla.Rdata")


diabetes_abc <- diabetes %>% 
 # distinct(atc_5, drug) %>%
  mutate(newdrug = paste0(row_number(),drug))

diabetes <- diabetes %>%
  left_join(diabetes_abc) %>%
  select(-drug) %>%
  rename(drug = newdrug)


set.seed(523)

n_its <- 1 # change based on how many iterations you want

scenarios <- c('atc5_0.05_trial_0.05_drug_0.05','atc5_0.15_trial_0.05_drug_0.05','atc5_0.25_trial_0.05_drug_0.05') 

slct_itrs <- mean_scenario_s1 %>%
  filter(scenario %in% scenarios) %>%
  filter(between(mean,-0.101, -0.099))

scenarios_selected_iterations <- mean_scenario_s1 %>%
   mutate(row = row_number()) %>%
   right_join(mean_scenario %>%
                filter(scenario %in% scenarios)) %>%
   pull(row)

scenarios_random_iterations <- mean_scenario_s1 %>% 
  mutate(row = row_number()) %>%
  right_join(slct_itrs) %>%
  group_by(scenario) %>% 
  sample_n(n_its) %>%
  ungroup() %>%
  pull(row)
## Loop through mean for 3 scenarios, single iteration for each
mdls <- map(scenarios_random_iterations, function (i) {
#mdls <- map(scenarios_selected_iterations, function (i) {
  scenario <- mean_scenario_s1$scenario[i]
  iter <- mean_scenario_s1$iter[i]
 
  # Add in main effect to a chosen variation scenario and select relevant classes
  diabetes$res <- res[, scenario] + -0.1
#   a10bx <- diabetes %>% 
#     filter(atc_5 == "A10BX")
#  # Select iteration
  ## Add values for specific iteration

  my_data$y <-  diabetes$res[diabetes$iteration == iter]
   
#   Add a new, NA drug in each class to get class level priors
  
  library(truncnorm)
  
  y_prec_new <- rtruncnorm(n=7, a=quantile(my_data$y_prec, probs=0.10), b=quantile(my_data$y_prec, probs=0.90), mean=mean(my_data$y_prec), sd=sd(my_data$y_prec))
  
  new_d <- as.data.frame(matrix(nrow=7,data=c(y_prec_new[1],162,1,1,25,NA,
                                              y_prec_new[2],163,1,2,26,NA,
                                              y_prec_new[3],164,1,3,27,NA,
                                              y_prec_new[4],165,1,4,28,NA,
                                              y_prec_new[5],166,1,5,29,NA,
                                              y_prec_new[6],167,1,6,30,NA,
                                              y_prec_new[7],168,1,7,31,NA), byrow=TRUE))
  colnames(new_d) <- names(my_data)
  my_data1 <- my_data %>%
    bind_rows(new_d) 
  
  ydens <- density(my_data$y, bw=0.25,n=75)

  # Make linear combination
  
  #NB INLA will order each variable used in the lincomb first, so values below should correspond to the myatc5 and mydrug numbers, not their
  #position in the my_data dataframe - e.g., to make a linear combination of drug class 2/4 and drug 5/8 :
  # myatc4 = c(NA,1,NA,NA), mydrug = c(Na,NA,NA,NA,1,NA,NA,NA) - irrespective of where in the input dataframe these are
  
  dc1 <- inla.make.lincomb(myatc4 = 1) # WDG, should be equivalent to fixed effect
  
  dc2 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(1, NA, NA, NA, NA, NA, NA), mydrug = c(rep(NA,24),1,rep(NA,6)) ) # Posterior, DC 1
  dc3 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, 1, NA, NA, NA, NA, NA), mydrug = c(rep(NA,25),1,rep(NA,5)) ) # Posterior, DC 2 etc
  dc4 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, 1, NA, NA, NA, NA), mydrug = c(rep(NA,26),1,rep(NA,4)) )
  dc5 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,27),1,rep(NA,3)) )
  dc6 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, 1, NA, NA), mydrug = c(rep(NA,28),1,rep(NA,2)) ) 
  dc7 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, 1, NA), mydrug = c(rep(NA,29),1,rep(NA,1)) ) 
  dc8 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, NA, 1), mydrug = c(rep(NA,30),1) ) 
 


  ## Run model, trial within drug within ATC5 class within ATC4 class, using full data
  mod1_nested <- inla(myform_nested2, 
                       data = my_data1,
                       # Add linear combinations to estimate drug-class
                       lincomb = c(a = dc1, b = dc2, c= dc3, d=dc4,e=dc5,f=dc6,g=dc7,h=dc8),
                       # Likelihood distribution
                       family = "gaussian",
                       # Fix likelihood hyperpars as the data is fixed with known precision.
                       control.family = list(hyper = list(prec = list(fixed = TRUE, initial = 0))),
                       # Likelihood precisions
                       # scale = my_data_test$y_prec,   #Use instead if testing 
                       scale = my_data1$y_prec,
                       # Prior distribution for "fixed" effects - really for mu_mu
                       control.fixed = list(mean = 0, prec = 0.25),
                       # Optionally compute DIC
                       verbose = FALSE,
                       control.compute = list(config=TRUE),
                       control.inla = list(lincomb.derived.only=FALSE))


  # Return glinides and model
  list( mdl1 = mod1_nested, scenario= scenario, iter=iter, ydens = ydens)
})
mdls <- transpose(mdls)
#a10bx <- as.data.frame(mdls$a10bx)
ydens <- mdls$ydens
scns <- mdls$scenario
itrs <- rep(seq(1, n_its, 1), times= length(unique(unlist(scns))))
mdl1 <- mdls$mdl1
mdls <- as.list(c(mdl1))


# Storing iterations in a single list
# scenario[[iter]] <- summary(mod1_nested2)
map(mdls, summary)
mdl_res <- map(mdls, ~ .x$summary.lincomb)
#map(mdls, ~  plot(.x$marginals.lincomb[[]], ylim = c(0,5), xlim = c(-2,2)))

## T distribution fits data well (as per BUGS)
mdl_t <- map(mdls, function(mdl_each){
  # Dataframe of values and densities
  res <- mdl_each$marginals.lincomb %>% 
    as.data.frame() %>% 
    as_tibble()
 # mdl_res_each <- mdl_each$summary.lincomb
  #list(mdl_t,mdl_res_each)
}) # abbreviated this function here to view distributions of marginals
mdl_ts <- list()

mdl_t_tdy <- mdl_t[1] %>%
  as.data.frame() %>%
  mutate(model = 'Full',
         scenario = scns[[1]]) %>% 
  bind_rows(mdl_t[2] %>%
              as.data.frame() %>%
              mutate(model = 'Full',
                     scenario = scns[[2]])) %>%
  bind_rows(mdl_t[3] %>%
              as.data.frame() %>%
              mutate(model = 'Full',
                     scenario = scns[[3]])) %>%
  as_tibble()

#Drop the non-applicable group(s)

#saveRDS(mdl_t_tdy, "scratch_data/mdl_t_tdy4")


load("Data/metadata_for_simulation.Rdata")
# Remove single alpha-glucosidase inhibitor trial, already out of diabetes file
diabetes_final <-  
  subset(diabetes_final, diabetes_final$atc_5 != "A10BF")
diabetes_final <- diabetes_final[with(diabetes_final, order(atc_5, drug, nct_id)),]
unique(diabetes_final$atc_5)
library(RColorBrewer)
dput(brewer.pal(n = 7, name = "Set2"))

drugs <- cbind(my_data, diabetes_final) %>%
  group_by(myatc5) %>%
  sample_n(1) %>%
  arrange(myatc5)

mdl_t_tdy2 <- mdl_t_tdy %>%
  select(a.lc.x,b.lc.x,c.lc.x,d.lc.x,e.lc.x,f.lc.x,g.lc.x,h.lc.x,model,scenario)%>%
  gather(key = xtype, value = x, -model,-scenario) %>%
  mutate(group = ifelse(xtype == 'aa.lc.x', 'z', str_sub(xtype, 1,1))) %>%
  bind_cols(mdl_t_tdy %>%
              select(a.lc.y,b.lc.y,c.lc.y,d.lc.y,e.lc.y,f.lc.y,g.lc.y,h.lc.y,model,scenario)%>%
              gather(key = ytype, value = y, -model,-scenario) %>%
              select(ytype, y))%>%
  mutate(myatc5 = match(group, letters)-1) %>%
  left_join(drugs, 'myatc5')  %>%
  select(model, scenario, x ,y, group, myatc5, atc_5,drug)


d <- as.tibble(matrix(ncol=8,nrow=75*length(scns),
                          data=c(rep('Null',75*length(scns)),
                                 c(rep(scns[[1]],75),rep(scns[[2]],75),rep(scns[[3]],75)),
                                 c(ydens[[1]]$x,ydens[[2]]$x,ydens[[3]]$x),
                                 c(ydens[[1]]$y,ydens[[2]]$y,ydens[[3]]$y),
                                 rep('z',75*length(scns)),
                                 rep(NA,75*length(scns)),
                                 rep(NA,75*length(scns)),
                                 rep(NA,75*length(scns))), byrow=FALSE))


colnames(d) <- names(mdl_t_tdy2)
d$x <- as.numeric(d$x)
d$y <- as.numeric(d$y)
mdl_t_tdy3 <- mdl_t_tdy2 %>%
  bind_rows(d)


mdl_t_tdy3$model <- factor(mdl_t_tdy3$model)
mdl_t_tdy3$model <- factor(mdl_t_tdy3$model,levels(mdl_t_tdy3$model)[c(2,1)]) 


ggplot(mdl_t_tdy3,aes(x=x, y=y, colour=group)) +
  geom_line(size=0.9) +
  scale_colour_manual("Distribution of interaction\n effect estimate:",
                      values = c(z="grey80",a = "black", b="#FC8D62",c="#8DA0CB",d= "#E78AC3", e= "#A6D854", f= "#FFD92F", g= "#E5C494",h= "#B3B3B3" ),
                      labels = c(z = "Overall (across all trials)", a ="WDG level", b="Class 1", c="Class 2", d="Class 3", 
                                 e="Class 4", f="Class 5", g="Class 6", 
                                 h="Class 7")) +
  facet_grid(   model ~ scenario  , scales = "free") +
  xlim(-1.2,1.2) +
  theme(axis.text.x = element_text(angle=0, vjust = 0, size = 11),
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text =element_text(size = 14.5),
        panel.background = element_rect(fill = "white", colour = NULL),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(-0.2,"lines"))

saveRDS(mdl_t_tdy2, "scratch_data/mdl_t_tdy2_sim1")


## Only 1 iteration per scenario should be taken forward  

#mdls <- list(mdls[[1]],mdls[[21]])

## T distribution fits data well (as per BUGS)
mdl_t <- map(mdls, function(mdl_each){
  # Dataframe of values and densities
  res <- mdl_each$marginals.lincomb %>% 
    as.data.frame() %>% 
    as_tibble() %>%
    rename(wdg_x = a.lc.x, wdg_y = a.lc.y,
           dc_x = b.lc.x, dc_y = b.lc.y,
           ndec_x = c.lc.x, ndec_y = c.lc.y)
  # Extract results
  mdl_res_each <- mdl_each$summary.lincomb
  list(mdl_t,mdl_res_each)
  
  # Fit data, using mean and sd and df 3 for initial values
  
  a <- nls(formula = wdg_y ~ metRology::dt.scaled(wdg_x, df = 3, mean = m, sd = s),
           data = res, start = list(m = mdl_res_each$mean[[1]], s = mdl_res_each$sd[[1]]/2),
           algorithm = "port", lower = list(m = c(-5), s = c(0)))
  a <- summary(a)
  a <- a$coefficients[,"Estimate"]
  
  dc_a <- nls(formula = dc_y ~ metRology::dt.scaled(dc_x, df = 3, mean = m, sd = s),
              data = res, start = list(m = mdl_res_each$mean[[2]], s = mdl_res_each$sd[[2]]/2),
              algorithm = "port", lower = list(m = c(-5), s = c(0)))
  dc_a <- summary(dc_a)
  dc_a <- dc_a$coefficients[,"Estimate"]
  
  ndec_a <- nls(formula = ndec_y ~ metRology::dt.scaled(ndec_x, df = 3, mean = m, sd = s),
                data = res, start = list(m = mdl_res_each$mean[[3]], s = mdl_res_each$sd[[3]]/2),
                algorithm = "port", lower = list(m = c(-5), s = c(0)))
  ndec_a <- summary(ndec_a)
  ndec_a <- ndec_a$coefficients[,"Estimate"]
  
  # Add estimates to dataframe for comparison
  res$y_new <- metRology::dt.scaled(res$wdg_x, df = 3, mean = a["m"], sd = a["s"])
  res$dc_y_new <- metRology::dt.scaled(res$dc_x, df = 3, mean = dc_a["m"], sd = dc_a["s"])
  res$ndec_y_new <- metRology::dt.scaled(res$ndec_x, df = 3, mean = ndec_a["m"], sd = ndec_a["s"])
  
  # Plot and return parameters and data
  plot(res$wdg_x, res$wdg_y, main = paste(c("df", "m", "s"), c(3, round(a,2)), sep = " = ", collapse = ", "))
  lines(res$wdg_x, res$y_new, col = "red")
  
  plot(res$dc_x, res$dc_y, main = paste(c("df", "m", "s"), c(3, round(dc_a,2)), sep = " = ", collapse = ", "))
  lines(res$dc_x, res$dc_y_new, col = "red")
  
  plot(res$ndec_x, res$ndec_y, main = paste(c("df", "m", "s"), c(3, round(ndec_a,2)), sep = " = ", collapse = ", "))
  lines(res$ndec_x, res$ndec_y_new, col = "red")
  
  # Make new dataset with more points
  res_new <- tibble(x = seq(-2, 2, 0.005))
  res_new <- res_new %>% 
    mutate(wdg_y = metRology::dt.scaled(x, df = 3, mean = a["m"], sd = a["s"])) %>%
    mutate(dc_y = metRology::dt.scaled(x, df = 3, mean = dc_a["m"], sd = dc_a["s"])) %>%
    mutate(ndec_y = metRology::dt.scaled(x, df = 3, mean = ndec_a["m"], sd = ndec_a["s"]))
  # Parameters and data
  list(param = c(a,dc_a,ndec_a), res = res_new)
})


names(mdl_t) <- paste(rep(c("Lo","Me","Hi"),  3), c( rep("Full", 3), rep("Full_noDCinfo",3),rep("NewDrgDC_only",3)), sep = "_")
mdl_t <- transpose(mdl_t)
#Drop the non-applicable group(s)
mdl_t$res$Lo_Full_noDCinfo$dc_y <- NA
mdl_t$res$Me_Full_noDCinfo$dc_y <- NA
mdl_t$res$Hi_Full_noDCinfo$dc_y <- NA
mdl_t$res$Lo_NewDrgDC_only$dc_y <- NA
mdl_t$res$Me_NewDrgDC_only$dc_y <- NA
mdl_t$res$Hi_NewDrgDC_only$dc_y <- NA
mdl_t_g <- bind_rows(mdl_t$res, .id = "scenario_model") %>% 
  mutate(variation = factor(str_sub(scenario_model,1,2 ) , levels = c("Lo","Me","Hi"),
                            labels = c("Minimum", "Median", "Maximum")),
         model = str_sub(scenario_model,3,-1L )) %>%
  gather(key = y_type, value = y, -scenario_model,-x,-variation, -model)


plot1 <- ggplot(mdl_t_g, aes( x = x, y = y,colour = variation)) +
  geom_line() +
  facet_grid(model~ y_type) +
  scale_x_continuous("Treatment-covariate interaction", limits = c(-1, 1)) +
  scale_y_continuous("Density")
plot1

saveRDS(mdl_t, file = "scratch_data/Priors_for_newdrug_extclass.Rds")
