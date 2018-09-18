# 06_Show_impact_simulations
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
  filter(scenario %in% key_scens_list$scenario) %>% 
  mutate(scenario = str_sub(scenario, 4, -1L))
# Identify mean for each scenarios

#02b_run_inla_model

# INLA:::inla.dynload.workaround() 
# load(file = "data/for_inla.Rdata")

#Separate for each sim from here

mean_scenario_s1 <- scenario_res %>%
  filter(sim==1)%>% 
  mutate(scenario = str_sub(scenario, 4, -1L)) %>%
  arrange(scenario)

# mean_scenario_s1 <- mean_scenario %>%
#   filter(sim==1)

load(file = "data/sim1/std/for_inla.Rdata")

## NEW CLASS (dc1) NEW DRUG IN NEW CLASS (dc2) NEW DRUG IN EXISTING CLASS (dc3) 

diabetes_abc <- diabetes %>% 
  distinct(atc_5, drug) %>%
  mutate(newdrug = paste0(row_number(),drug))

diabetes <- diabetes %>%
  left_join(diabetes_abc) %>%
  select(-drug) %>%
  rename(drug = newdrug)


set.seed(2345)

n_its <- 5 # change based on how many iterations you want

scenarios <- c('atc5_0.05_trial_0.05_drug_0.05','atc5_0.05_trial_0.05_drug_0.15','atc5_0.05_trial_0.05_drug_0.25')

scenarios_random_iterations <- mean_scenario_s1 %>% 
  mutate(row = row_number()) %>%
  filter(scenario %in% scenarios) %>%
  group_by(scenario) %>% 
  sample_n(n_its) %>%
  ungroup() %>%
  pull(row)

## Loop through mean for 3 scenarios, single iteration for each
mdls <- map(scenarios_random_iterations, function (i) {
#mdls <- map(c(sample(1:1000,n_its),sample(16001:17000,n_its),sample(26001:27000,n_its)), function (i) {
  scenario <- mean_scenario_s1$scenario[i]
  iter <- mean_scenario_s1$iter[i]
 
  # Add in main effect to a chosen variation scenario and select relevant classes
  diabetes$res <- res[, scenario] + -0.1
  a10bx <- diabetes %>% 
    filter(atc_5 == "A10BX")
  
  # Create version where one DC is missing
  no_a10bx <- diabetes
  no_a10bx$res[diabetes$atc_5 %in% c("A10BX")] <- NA
  
   no_a10bx %>% 
    distinct(atc_5, drug)
  
  ###DM##
  # The drug is set to NA becuase all of the data in that class is set to NA
  no_a10bx %>%  filter(atc_5 == "A10BX") %>%  distinct(atc_5, nct_id,  drug, res)
#     atc_5      nct_id          drug res
# 1 A10BX NCT00519142 22mitiglinide  NA
# 2 A10BX NCT00097786 23nateglinide  NA
# 3 A10BX NCT00568984 24repaglinide  NA
  
     
  # Select iteration
  ## Add values for specific iteration
  my_data$y <-  no_a10bx$res[no_a10bx$iteration == iter]
  
  # Add in 'new' drug in class
  
  new_d <- as.data.frame(matrix(nrow=1,data=c(mean(my_data$y_prec),162,1,6,25,NA)))
  colnames(new_d) <- names(my_data)
  my_data1 <- my_data %>%
    bind_rows(new_d) 
  

  # Make linear combination
  dc1 <- inla.make.lincomb(myatc4 = 1) # WDG, should be equivalent to fixed effect
  
  dc2 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, NA, 1)) # DC level
  
  dc3 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, NA, 1),  # DC new drug new class
                            mydrug = c(rep(NA, 21), 1,NA,NA, NA)) 
  
  dc4 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, 1, NA), # DC new drug ext class
                           mydrug = c(rep(NA, 24), 1)) 

 ### This section used to test INLAs handling of lincomb instructions relative to data.frame ordering   
  ## Bit of code to create list of vectors for different linear combinations
  mydrug_vect_inld <- map(1:6, function(x){
    mydrug_vect <- rep(NA, 6)
    myatc5_vect <- c(NA, NA)
    
    mydrug_vect[x] <- 1
    if(x <=3) myatc5_vect[1] <- 1
    if(x>=4) myatc5_vect[2] <- 1
    list(myatc5_vect, mydrug_vect)
    })
 names(mydrug_vect_inld) <- paste0("i", 1:6)
 mydrug_vect_inld <- transpose(mydrug_vect_inld)   
 names(mydrug_vect_inld) <- c("myatc5_vect", "mydrug_vect") 
#myatc5_vect_incld <- map()
  ## create linear combinations
 dc_tests <- map2(mydrug_vect_inld$myatc5_vect, mydrug_vect_inld$mydrug_vect, function(atc5, drug){
    inla.make.lincomb(myatc4 =1, myatc5 = atc5, mydrug = drug)
 })
  
 ## create fake data
  my_data_test <- my_data %>% 
    filter(myatc5 <=2) %>% 
    group_by(myatc4, myatc5) %>% 
    slice(1:3) %>% 
    ungroup() %>% 
    mutate(y_prec = rep(100, 6),
           # Can vary the following to test what happens with INLA when given a dataframe
           y = c(-1, -0.5, 0, 0.5, 1, 2),
           trial = c(1:6),
           mydrug = c(1:3,6,5,4))
###
  
  ## Run model, trial within drug within ATC5 class within ATC4 class
  mod1_nested2 <- inla(myform_nested2, 
                       data = my_data1,
                       #data = my_data_test,            #Use instead if testing 
                       # Add linear combinations to estimate drug-class
                       #lincomb = dc4,
                       lincomb = c(a = dc1, b = dc2, c= dc3, d=dc4),
                       #lincom = do.call(c, dc_tests),  #Use instead if testing 
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
  list(a10bx = a10bx, mdl = mod1_nested2, scenario= scenario, iter=iter)
})
mdls <- transpose(mdls)
a10bx <- mdls$a10bx
scns <- mdls$scenario
itrs <- rep(seq(1, n_its, 1), times= length(unique(unlist(scns))))
mdls <- mdls$mdl

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
 
}) # abbreviated this function here to view distributions of marginals
mdl_ts <- list()

for (i in 1:length(mdl_t)){
  mdl_t2 <- mdl_t[[i]] %>%
  mutate(scenario = scns[[i]],
         iteration = itrs[[i]])
  mdl_ts[[length(mdl_ts)+1]] <- mdl_t2
  }
mdl_t_tdy <- bind_rows(mdl_ts)


#saveRDS(mdl_t_tdy, "scratch_data/mdl_t_tdy4")



ggplot(mdl_t_tdy) +
  geom_line( aes(x=a.lc.x,y=a.lc.y, colour="WDG")) +
  geom_line( aes(x=b.lc.x,y=b.lc.y, colour = "DC")) +
  geom_line( aes(x=c.lc.x,y=c.lc.y, colour= "NewDrug_NewDC")) +
  geom_line( aes(x=d.lc.x,y=d.lc.y, colour= "NewDrug_ExistDC")) +
  scale_colour_manual("", 
                      breaks = c("WDG", "DC", "NewDrug_NewDC","NewDrug_ExistDC"),
                      values = c("red", "orange", "blue","darkgreen")) +
  facet_grid( scenario~iteration) +
  xlim(-1,1) +
  geom_vline(xintercept = -0.1, linetype = "dashed", color = "grey50") +theme_bw()


## Modified to here

## T distribution fits data well (as per BUGS)
mdl_t <- map(mdls, function(mdl_each){
  # Dataframe of values and densities
  res <- mdl_each$marginals.lincomb %>% 
    as.data.frame() %>% 
    as_tibble() %>%
    rename(x = lc.x, y = lc.y)
  # Extract results
  mdl_res_each <- mdl_each$summary.lincomb
  list(mdl_t,mdl_res_each)
  

  # Fit data, using mean and sd and df 3 for initial values
  a <- nls(formula = y ~ metRology::dt.scaled(x, df = 3, mean = m, sd = s),
         data = res, start = list(m = mdl_res_each$mean, s = mdl_res_each$sd/2),
         algorithm = "port", lower = list(m = c(-5), s = c(0)))
  a <- summary(a)
  a <- a$coefficients[,"Estimate"]
  
  # Add estimates to dataframe for comparison
  res$y_new <- metRology::dt.scaled(res$x, df = 3, mean = a["m"], sd = a["s"])

  # Plot and return parameters and data
  plot(res$x, res$y, main = paste(c("df", "m", "s"), c(3, round(a,2)), sep = " = ", collapse = ", "))
  lines(res$x, res$y_new, col = "red")
  
  # Make new dataset with more points
  res_new <- tibble(x = seq(-2, 2, 0.005))
  res_new <- res_new %>% 
    mutate(y = metRology::dt.scaled(x, df = 3, mean = a["m"], sd = a["s"]))
  # Parameters and data
  list(param = a, res = res_new)
})

names(mdl_t) <- mean_scenario_s1$scenario
mdl_t <- transpose(mdl_t)
mdl_t_g <- bind_rows(mdl_t$res, .id = "scenario") %>% 
  mutate(variation = factor(scenario, levels = mean_scenario_s1$scenario,
                         labels = c("Minimum", "Median", "Maximum")))

plot1 <- ggplot(mdl_t_g, aes(x = x, y = y, colour = variation)) +
  geom_line() +
  scale_x_continuous("Treatment-covariate interaction", limits = c(-1, 1)) +
  scale_y_continuous("Density")
plot1

saveRDS(mdl_t, file = "scratch_data/Priors_for_examining_class_dc4.Rds")

## Random 


#quantile version
qres <- data.frame(y = c(0.025, 0.5, 0.975),
                   x = c( -0.2867,  -0.0963,     0.0631))

res2 <- nls(formula = y ~ metRology::pt.scaled(x, df = 3, mean = m, sd = s),
       data = qres, start = list(m = -0.1, s = 0.2),
       algorithm = "port", lower = list(m = -5, s = 0))
summary(res2)
a <- summary(res2)
a <- a$coefficients[,"Estimate"]

qres_new <- data.frame(x = seq(-1, 1, 0.1))
qres_new <- qres_new %>% 
  mutate(y = metRology::pt.scaled(x, df = 3, mean = a["m"], sd = a["s"]))
plot(qres_new$x, qres_new$y, type = "l")
points(qres$x, qres$y, col = "red")

qres_samples = metRology::rt.scaled(1000, df = 3, mean = a["m"], sd = a["s"])
hist(qres_samples)
mean(qres_samples)
sd(qres_samples)
