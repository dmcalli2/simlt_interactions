# 06_Show_impact_priors_dc_vs_no_dc_info
library(metRology)
library(tidyverse)
library(INLA)
library(stringr)

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
  mutate(scenario = str_sub(scenario, 4, -1L)) %>%
  filter(sim ==2)
# Identify mean for each scenarios


#Separate for each sim from here

mean_scenario_s2 <- scenario_res %>%
  filter(sim==2)%>% 
  mutate(scenario = str_sub(scenario, 4, -1L)) %>%
  arrange(scenario)


load(file = "data/sim2/std/for_inla.Rdata")


rheum_abc <- rheum %>% 
  distinct(moa, drug) %>%
  mutate(newdrug = paste0(row_number(),drug))

rheum <- rheum %>%
  left_join(rheum_abc) %>%
  select(-drug) %>%
  rename(drug = newdrug)

set.seed(375)

n_its <- 1 # change based on how many iterations you want

scenarios <- c('path_0.05_moa_0.05_trial_0.05_drug_0.05','path_0.15_moa_0.15_trial_0.15_drug_0.15','path_0.25_moa_0.25_trial_0.25_drug_0.25') 

slct_itrs <- mean_scenario_s2 %>%
  filter(scenario %in% scenarios) %>%
  filter(between(mean,-0.101, -0.099))

scenarios_selected_iterations <- mean_scenario_s2 %>%
  mutate(row = row_number()) %>%
  right_join(mean_scenario %>%
               filter(scenario %in% scenarios)) %>%
  pull(row)

scenarios_random_iterations <- mean_scenario_s2 %>% 
  mutate(row = row_number()) %>%
  right_join(slct_itrs) %>%
  group_by(scenario) %>% 
  sample_n(n_its) %>%
  ungroup() %>%
  pull(row)
## Loop through mean for 3 scenarios, single iteration for each
#mdls <- map(scenarios_selected_iterations, function (i) {
mdls <- map(scenarios_random_iterations, function (i) {
  scenario <- mean_scenario_s2$scenario[i]
  iter <- mean_scenario_s2$iter[i]
 
  # Add in main effect to a chosen variation scenario and select relevant classes
  rheum$res <- res[, scenario] + -0.1

  # Select iteration
  ## Add values for specific iteration

  my_data$y <-  rheum$res[rheum$iteration == iter]

  # Add in 3 trials of a new drug in an existing class
  
  new_d <- as.data.frame(matrix(nrow=3,data=c(mean(my_data$y_prec),163,1,6,33,2,NA,
                                              mean(my_data$y_prec),164,1,6,33,2,NA,
                                              mean(my_data$y_prec),165,1,6,33,2,NA), byrow=TRUE))
  colnames(new_d) <- names(my_data)
  my_data1 <- my_data %>%
    bind_rows(new_d) %>%
    arrange(wdg, mymoa,mydrug, trial)

  # Make linear combination
  
  #NB INLA will order each variable used in the lincomb first, so values below should correspond to the myatc5 and mydrug numbers, not their
  #position in the my_data dataframe - e.g., to make a linear combination of drug class 2/4 and drug 5/8 :
  # myatc4 = c(NA,1,NA,NA), mydrug = c(Na,NA,NA,NA,1,NA,NA,NA) - irrespective of where in the input dataframe these are
  
  dc1 <- inla.make.lincomb(wdg = 1) # WDG, should be equivalent to fixed effect
  
  dc2 <- inla.make.lincomb(wdg = 1, mypath = c(NA,1)) # path, should be equivalent to fixed effect
  
  dc3 <- inla.make.lincomb(wdg = 1, mypath = c(NA,1), mymoa = c(rep(NA,5),1,rep(NA,11) )) # DC level, new drug in a ext class
  
  dc4 <- inla.make.lincomb(wdg = 1, mypath = c(NA,1), mymoa = c(rep(NA,5),1,rep(NA,11)),  # drug level,  new drug ext class  
                            mydrug = c(rep(NA, 32), 1)) 
  

  ## Run model, trial within drug within ATC5 class within ATC4 class, using full data
  mod1_nested <- inla(myform_nested2, 
                       data = my_data1,
                       # Add linear combinations to estimate drug-class
                       lincomb = c(a = dc1, b = dc2, c= dc3, d= dc4),
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

 # Re-write model without MoA and path
  
  myform_nested3 <- y ~ -1 + wdg +    ## wdg = wider (/widest) drug grouping (top of hierarchy)
    f(trial, model = "iid", 
      hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1)))) +
    f(mydrug, model = "iid", 
      hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1)))) 

  
  # Re-make linear combination for no DC model
  dc1 <- inla.make.lincomb(wdg = 1) # WDG, should be equivalent to fixed effect
  
  dc2 <- inla.make.lincomb(wdg = 1) # path, NOT APPLICABLE - this is just a dummy
  
  dc3 <- inla.make.lincomb(wdg = 1) # path, NOT APPLICABLE - this is just a dummy
  
  dc4 <- inla.make.lincomb(wdg = 1, mydrug = c(rep(NA, 32), 1))  # drug level,  new drug new class  
                           
  
  ## Run model, trial within drug within ATC4 class (intermediate dc level info not included, but whole WDG used)
  mod2_nested <- inla(myform_nested3, 
                       data = my_data1,
                       # Add linear combinations to estimate drug-class
                       lincomb = c(a = dc1, b = dc2, c= dc3, d= dc4),
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
 
  my_data2 <- my_data1 %>%
    filter(mymoa==6)   
  # Re-make linear combination for single DC model
  
  dc1 <- inla.make.lincomb(wdg = 1) # WDG, should be equivalent to fixed effect
  dc2 <- inla.make.lincomb(wdg = 1) # Path level   #Not applicable, this included as dummy
  dc3 <- inla.make.lincomb(wdg = 1) # DC level   #Not applicable, this included as dummy
  dc4 <- inla.make.lincomb(wdg = 1,mydrug = c(NA, 1))    # DC new drug new class 
                                ## Drug 33 is now the 2nd drug in this class
  ## Run model, trial within drug within ATC4 class (intermediate dc level info not included, but whole WDG used)
  mod3_nested <- inla(myform_nested3, 
                      data = my_data2,
                      # Add linear combinations to estimate drug-class
                      lincomb = c(a = dc1, b = dc2, c= dc3, d= dc4),
                      # Likelihood distribution
                      family = "gaussian",
                      # Fix likelihood hyperpars as the data is fixed with known precision.
                      control.family = list(hyper = list(prec = list(fixed = TRUE, initial = 0))),
                      # Likelihood precisions
                      # scale = my_data_test$y_prec,   #Use instead if testing 
                      scale = my_data2$y_prec,
                      # Prior distribution for "fixed" effects - really for mu_mu
                      control.fixed = list(mean = 0, prec = 0.25),
                      # Optionally compute DIC
                      verbose = FALSE,
                      control.compute = list(config=TRUE),
                      control.inla = list(lincomb.derived.only=FALSE))
  
  # Return glinides and model
  list( mdl1 = mod1_nested, mdl2 = mod2_nested ,mdl3 = mod3_nested, scenario= scenario, iter=iter)
})
mdls <- transpose(mdls)  
#a10bx <- as.data.frame(mdls$a10bx)
scns <- mdls$scenario
itrs <- rep(seq(1, n_its, 1), times= length(unique(unlist(scns))))
mdl1 <- mdls$mdl1
mdl2 <- mdls$mdl2
mdl3 <- mdls$mdl3
mdls <- as.list(c(mdl1,mdl2,mdl3))

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
  bind_rows(mdl_t[4] %>%
              as.data.frame() %>%
              mutate(model = 'Full, no DC/path info',
                     scenario = scns[[1]])) %>%
  bind_rows(mdl_t[5] %>%
              as.data.frame() %>%
              mutate(model = 'Full, no DC/path info',
                     scenario = scns[[2]])) %>%
  bind_rows(mdl_t[6] %>%
              as.data.frame() %>%
              mutate(model = 'Full, no DC/path info',
                     scenario = scns[[3]])) %>%
  bind_rows(mdl_t[7] %>%
              as.data.frame() %>%
              mutate(model = 'DC of new drug only',
                     scenario = scns[[1]])) %>% 
  bind_rows(mdl_t[8] %>%
              as.data.frame() %>%
              mutate(model = 'DC of new drug only',
                     scenario = scns[[2]])) %>% 
  bind_rows(mdl_t[9] %>%
              as.data.frame() %>%
              mutate(model = 'DC of new drug only',
                     scenario = scns[[3]])) %>% 
  as_tibble()


#Drop the non-applicable group(s)

mdl_t_tdy <- mdl_t_tdy %>%
  mutate( b.lc.x = ifelse(model == "Full, no DC/path info",NA,b.lc.x),
          b.lc.y = ifelse(model == "Full, no DC/path info",NA,b.lc.y),
          c.lc.x = ifelse(model == "Full, no DC/path info",NA,c.lc.x),
          c.lc.y = ifelse(model == "Full, no DC/path info",NA,c.lc.y),
          b.lc.x = ifelse(model == "DC of new drug only",NA,b.lc.x),
          b.lc.y = ifelse(model == "DC of new drug only",NA,b.lc.y),
          c.lc.x = ifelse(model == "DC of new drug only",NA,c.lc.x),
          c.lc.y = ifelse(model == "DC of new drug only",NA,c.lc.y))

mdl_t_tdy$model <- as.factor(mdl_t_tdy$model )
mdl_t_tdy$model <- factor(mdl_t_tdy$model,levels(mdl_t_tdy$model)[c(3,1,2)]) 

ggplot(mdl_t_tdy) +
  geom_line( aes(x=a.lc.x,y=a.lc.y, colour="WDG")) +
  geom_line( aes(x=b.lc.x,y=b.lc.y, colour= "NewDrug_NewDC_Pathlevel")) +
  geom_line( aes(x=c.lc.x,y=c.lc.y, colour = "NewDrug_NewDC_DClevel")) +
  geom_line( aes(x=d.lc.x,y=d.lc.y, colour= "NewDrug_NewDC_Druglevel")) +
  scale_colour_manual("", 
                      breaks = c( "WDG", "NewDrug_NewDC_Pathlevel","NewDrug_NewDC_DClevel","NewDrug_NewDC_Druglevel","NewDrug_NewDC_DClevel"),
                      values = c("purple","red", "darkgreen","orange")) +
  facet_grid( model ~ scenario) +
  xlim(-1,1) +
  geom_vline(xintercept = -0.1, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = -0.1, linetype = "dashed", color = "grey50") +
  theme_bw()

## Only 1 iteration per scenario should be taken forward  

#mdls <- list(mdls[[1]],mdls[[21]])

## T distribution fits data well (as per BUGS)
mdl_t <- map(mdls, function(mdl_each){
  # Dataframe of values and densities
  res <- mdl_each$marginals.lincomb %>% 
    as.data.frame() %>% 
    as_tibble() %>%
    rename(wdg_x = a.lc.x, wdg_y = a.lc.y,
           pth_x = b.lc.x, pth_y = b.lc.y,
           dc_x = c.lc.x, dc_y = c.lc.y,
           ndnc_x = d.lc.x, ndnc_y = d.lc.y)
  # Extract results
  mdl_res_each <- mdl_each$summary.lincomb
  list(mdl_t,mdl_res_each)
  
  # Fit data, using mean and sd and df 3 for initial values

  a <- nls(formula = wdg_y ~ metRology::dt.scaled(wdg_x, df = 3, mean = m, sd = s),
         data = res, start = list(m = mdl_res_each$mean[[1]], s = mdl_res_each$sd[[1]]/2),
         algorithm = "port", lower = list(m = c(-5), s = c(0)))
  a <- summary(a)
  a <- a$coefficients[,"Estimate"]
  
  pth_a <- nls(formula = pth_y ~ metRology::dt.scaled(pth_x, df = 3, mean = m, sd = s),
              data = res, start = list(m = mdl_res_each$mean[[2]], s = mdl_res_each$sd[[2]]/2),
              algorithm = "port", lower = list(m = c(-5), s = c(0)))
  pth_a <- summary(pth_a)
  pth_a <- pth_a$coefficients[,"Estimate"]
  
  dc_a <- nls(formula = dc_y ~ metRology::dt.scaled(dc_x, df = 3, mean = m, sd = s),
           data = res, start = list(m = mdl_res_each$mean[[3]], s = mdl_res_each$sd[[3]]/2),
           algorithm = "port", lower = list(m = c(-5), s = c(0)))
  dc_a <- summary(dc_a)
  dc_a <- dc_a$coefficients[,"Estimate"]
  
  ndnc_a <- nls(formula = ndnc_y ~ metRology::dt.scaled(ndnc_x, df = 3, mean = m, sd = s),
           data = res, start = list(m = mdl_res_each$mean[[4]], s = mdl_res_each$sd[[4]]/2),
           algorithm = "port", lower = list(m = c(-5), s = c(0)))
  ndnc_a <- summary(ndnc_a)
  ndnc_a <- ndnc_a$coefficients[,"Estimate"]
  
  # Add estimates to dataframe for comparison
  res$y_new <- metRology::dt.scaled(res$wdg_x, df = 3, mean = a["m"], sd = a["s"])
  res$pth_y_new <- metRology::dt.scaled(res$pth_x, df = 3, mean = pth_a["m"], sd = pth_a["s"])
  res$dc_y_new <- metRology::dt.scaled(res$dc_x, df = 3, mean = dc_a["m"], sd = dc_a["s"])
  res$ndnc_y_new <- metRology::dt.scaled(res$ndnc_x, df = 3, mean = ndnc_a["m"], sd = ndnc_a["s"])

  # Plot and return parameters and data
  plot(res$wdg_x, res$wdg_y, main = paste(c("df", "m", "s"), c(3, round(a,2)), sep = " = ", collapse = ", "))
  lines(res$wdg_x, res$y_new, col = "red")
  
  plot(res$pth_x, res$pth_y, main = paste(c("df", "m", "s"), c(3, round(pth_a,2)), sep = " = ", collapse = ", "))
  lines(res$pth_x, res$pth_y_new, col = "red")
  
  plot(res$dc_x, res$dc_y, main = paste(c("df", "m", "s"), c(3, round(dc_a,2)), sep = " = ", collapse = ", "))
  lines(res$dc_x, res$dc_y_new, col = "red")
  
  plot(res$ndnc_x, res$ndnc_y, main = paste(c("df", "m", "s"), c(3, round(ndnc_a,2)), sep = " = ", collapse = ", "))
  lines(res$ndnc_x, res$ndnc_y_new, col = "red")

  # Make new dataset with more points
  res_new <- tibble(x = seq(-2, 2, 0.005))
  res_new <- res_new %>% 
    mutate(wdg_y = metRology::dt.scaled(x, df = 3, mean = a["m"], sd = a["s"])) %>%
    mutate(pth_y = metRology::dt.scaled(x, df = 3, mean = pth_a["m"], sd = pth_a["s"])) %>%
    mutate(dc_y = metRology::dt.scaled(x, df = 3, mean = dc_a["m"], sd = dc_a["s"])) %>%
    mutate(ndnc_y = metRology::dt.scaled(x, df = 3, mean = ndnc_a["m"], sd = ndnc_a["s"]))
  # Parameters and data
  list(param = c(a,pth_a,dc_a,ndnc_a), res = res_new)
})


names(mdl_t) <- paste(rep(c("Lo","Me","Hi"),  3), c( rep("Full", 3), rep("Full_noDC_Path_info",3),rep("NewDrgDC_only",3)), sep = "_")
mdl_t <- transpose(mdl_t)
#Drop the non-applicable group(s)
mdl_t$res$Lo_Full_noDC_Path_info$pth_y <- NA
mdl_t$res$Me_Full_noDC_Path_info$pth_y <- NA
mdl_t$res$Hi_Full_noDC_Path_info$pth_y <- NA
mdl_t$res$Lo_Full_noDC_Path_info$dc_y <- NA
mdl_t$res$Me_Full_noDC_Path_info$dc_y <- NA
mdl_t$res$Hi_Full_noDC_Path_info$dc_y <- NA
mdl_t$res$Lo_NewDrgDC_only$pth_y <- NA
mdl_t$res$Me_NewDrgDC_only$pth_y <- NA
mdl_t$res$Hi_NewDrgDC_only$pth_y <- NA
mdl_t$res$Lo_NewDrgDC_only$dc_y <- NA
mdl_t$res$Me_NewDrgDC_only$dc_y <- NA
mdl_t$res$Hi_NewDrgDC_only$dc_y <- NA
mdl_t_g <- bind_rows(mdl_t$res, .id = "scenario_model") %>% 
  mutate(variation = factor(str_sub(scenario_model,1,2 ) , levels = c("Lo","Me","Hi"),
                         labels = c("Minimum", "Median", "Maximum")),
         model = str_sub(scenario_model,4,-1L )) %>%
  gather(key = y_type, value = y, -scenario_model,-x,-variation, -model)


mdl_t_g$y_type <- as.factor(mdl_t_g$y_type )
mdl_t_g$y_type <- factor(mdl_t_g$y_type,levels(mdl_t_g$y_type)[c(4,3,1,2)]) 


plot1 <- ggplot(mdl_t_g, aes( x = x, y = y,colour = variation)) +
  geom_line() +
  facet_grid(model~ y_type) +
  scale_x_continuous("Treatment-covariate interaction", limits = c(-1, 1)) +
  scale_y_continuous("Density")
plot1

saveRDS(mdl_t, file = "scratch_data/Priors_for_newdrug_extclass_sim2.Rds")

