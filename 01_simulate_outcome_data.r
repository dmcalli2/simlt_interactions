#01 Simulate outcome data
### Run through 187 trials and simulat data with interactions
load("Dataset_of_trials.Rdata")

# Packages
library(tidyverse)
library("simstudy")
library(readxl)
library(Hmisc)
library(purrr)

# Functions
MakeCsv <- function (mypath) {
  # Functions to let me use excel rather than CSV files in constructing data description
  my_name <- round(runif(1, 10000, 90000),0)
  read_excel(path = mypath) %>%
    write_csv(paste0("scratch_data/", my_name, ".csv"))
  paste0("scratch_data/", my_name, ".csv")
}

# Take spreadsheet and make random samples
SampleFx <- function(desired_n, means, sds, varname = paste0("q", seq_along(means))) {
  a <- map2(means, sds, ~ rnorm(desired_n, .x, .y))
  a <- do.call(cbind, a) %>%  as_tibble()
  names(a) <- varname
  a
}

# Precision matrix
PrecMtrx <- function(model) solve(vcov(model))
sim_data <- select (trials, drug_group_allocation, my_condition, my_drug, enrollment)
rm(trials)

## Read in choice of parameters for simulation
sim_effects <- read_excel(path = "data/sim_effects.xlsx")  
print(sim_effects)

## Read in baseline variable definitions, note centred age
def <- defRead(MakeCsv("data/sim_df1.xlsx"))
print(def)

# Create vectors of effects at different levels
sim_data$study_id <- row_number(sim_data$drug_group_allocation)
sim_list <- sim_data[, c("drug_group_allocation", "my_condition", "my_drug")] %>%  as.list()
sim_list <- map(sim_list, unique)
sim_list_smpls <- map2(sim_list, names(sim_list), function (x, y) {
  sim_effects_select <- sim_effects[sim_effects$level_effect == y,]
  with(sim_effects_select,  SampleFx(desired_n = length(x), means = means, 
                                     sds = sds, varname = var_effect))
})
sim_list_smpls <- map2(sim_list, sim_list_smpls, function(x,y) {
  y$eff_lvl <- x
  y}
  )

sim_data_dg <- inner_join(sim_data %>% select(study_id, drug_group_allocation),
                          sim_list_smpls$drug_group_allocation,
                       by = c(`drug_group_allocation` = "eff_lvl"))
sim_data_cond <- inner_join(sim_data %>% select(study_id, my_condition),
                          sim_list_smpls$my_condition,
                       by = c(`my_condition` = "eff_lvl"))
sim_data_d <- inner_join(sim_data %>% select(study_id, my_drug),
                          sim_list_smpls$my_drug,
                       by = c(`my_drug` = "eff_lvl"))
sim_data <- bind_rows(sim_data_dg, sim_data_cond, sim_data_d) %>% 
  select(-drug_group_allocation, -my_condition, -my_drug) %>% 
  group_by(study_id) %>% 
  summarise_all(sum) %>% 
  inner_join(sim_data, by = "study_id")
rm(sim_data_cond, sim_data_d, sim_data_dg, sim_effects, sim_list, sim_list_smpls)

## Make list to loop through for dataframe
sim_list <- as.list(sim_data[, c("dep_eff", "pain_eff", "alloc_eff", "alloc_dep_eff", 
"alloc_pain_eff", "enrollment")]) 

## Loop through each trial
res <- pmap(sim_list, .f = function (dep_eff, pain_eff, alloc_eff, alloc_dep_eff, alloc_pain_eff, 
                                     enrollment) {
  # Generate data for each trial
  dt1 <- simstudy::genData(enrollment, def)

  ## Put age into quintiles to stratify by it
  dt1$age5 <- cut2(dt1$age, g = 5)

  ## Add treatent allocation, stratified by age and sex
  dt1 <- trtAssign(dt1, n = 2, balanced = TRUE, strata = c("age5", "sex"), 
            grpName = "alloc")

    ## Calculate outcome based on covariates and treatment, with interactions
  ## Each coefficient is simulated as coming from a normal distribution, except for age
  ## and sex
  ## dep_eff, pain_eff, alloc_eff, alloc_dep_eff, alloc_pain_eff
  trt_eff_mean <- with(dt1,
                       qlogis(0.2) + # USe this to create an intercept so is 20%, rest of data centred
                       age*0.01 +
                         sex*0.05 +
                         dep*dep_eff +
                         pain*pain_eff +
                         alloc*alloc_eff +
                         alloc*dep*alloc_dep_eff +
                         alloc*pain*alloc_pain_eff
                       )
 
  ## For logistic regression model sample from binomial error
  dt1$outcome_cat <- rbinom(length(trt_eff_mean), size = 1, prob = plogis(trt_eff_mean))
  ## Add a tiny bit of error to linear predictor so continuous model runs even if no error in coefficients
  dt1$outcome_con <- rnorm(length(trt_eff_mean), mean = trt_eff_mean, sd = 0.01)

  ## Centre all dichotmous variables (age already centred)
  dt1[ , c("alloc", "sex", "dep", "pain")] <- lapply(dt1[ , c("alloc", "sex", "dep", "pain")],
                                                     function(x) x - mean(x))
  ## Run regression model on data and save coefficients, variance covariance matrix
  ## and precision matrix
  mod_con <- glm(outcome_con ~ age + sex + dep + pain +
                  alloc +
                  alloc:age + alloc:sex + alloc:dep + alloc:pain,
                family = "gaussian",
                data = dt1)
  # mod_old <- update(mod_con, data = dt_old)
  mod_cat <- update(mod_con, outcome_cat ~ . , family = "binomial")

  mdls <- list(con = mod_con, cat = mod_cat)

    mdls  <- map(mdls, function (x) {
    map(list(coef = coef, vcov = vcov, prec_matrix = PrecMtrx), function (y) y(x))
  })
  mdls
})

## Save coefficients and variance-covariance matrix for each model to data
sim_data$con <- map(res, function (x) x$con)
sim_data$cat <- map(res, function (x) x$cat)

# Check precision matrix
PrecCheck <- function (model_type = sim_data$con) {
  map_lgl(model_type, ~ matrixcalc::is.positive.semi.definite(.x$prec_matrix %>%  
                                                                              round(3)))
}
any(!PrecCheck(sim_data$con)) # should be FALSE
any(!PrecCheck(sim_data$cat)) # should be FALSE

saveRDS(sim_data, file = "scratch_data/simulated_data.Rds")

## Examine effect on glm of centring variables

