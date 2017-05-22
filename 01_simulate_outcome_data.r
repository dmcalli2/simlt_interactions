#01 Simulate outcome data
### Run through 187 trials and simulat data with interactions


load("Dataset_of_trials.Rdata")

# Packages
library(tidyverse)
library("simstudy")
library(readxl)
library(Hmisc)

# Functions
MakeCsv <- function (mypath = "data/sim_df1.xlsx") {
  # Functions to let me use excel rather than CSV files in constructing data
  my_name <- round(runif(1, 10000, 90000),0)
  read_excel(path = mypath) %>%
    write_csv(paste0("scratch_data/", my_name, ".csv"))
  paste0("scratch_data/", my_name, ".csv")
}

sim_data <- select (trials, drug_group_allocation, my_condition, my_drug, enrollment)
rm(trials)

## Create baseline variables
def <- defRead(MakeCsv())
print(def)

## Loop through each trial
res <- map(sim_data$enrollment, function (enrol) {
  dt1 <- simstudy::genData(enrol, def)
  
  ## Dichotomise age to stratify by it
  dt1$age5 <- cut2(dt1$age, g = 5)
  
  ## Add treatent allocation, stratified by age and sex
  dt1 <- trtAssign(dt1, n = 2, balanced = TRUE, strata = c("age5", "sex"), 
            grpName = "alloc")
  
  ## Calculate outcome based on covariates and treatment, with interactions
  trt_eff_mean <- with(dt1,
                       age*0.01 +
                         sex*0.05 +
                         dep*0.05 +
                         pain* 0.1 +
                         alloc * -0.20 +
                         alloc*dep*0.10 +
                         alloc*pain*0.02)
  dt1$outcome_con <- rnorm(nrow(dt1), mean = trt_eff_mean, sd = 1)
  
  ## Continuous outcome as normal distribution with mean = 0, and SD 1
  dt1$outcome_cat <- dt1$outcome_con > quantile(dt1$outcome, probs = 0.9)
  dt1$outcome_log <- log(dt1$outcome_con + 0.0001- min(dt1$outcome_con))
  dt1$outcome_con <- dt1$outcome_con + rnorm(nrow(dt1), 0, 0.2)
  
  ## Run regression model on data adn save coefficient and variance covariance matrix
  mod_con <- glm(outcome_con ~ age + sex + dep + pain +
                  alloc +
                  alloc:age + alloc:sex + alloc:dep + alloc:pain,
                family = "gaussian",
                data = dt1)
  mod_cat <- update(mod_con, outcome_cat ~ . , family = "binomial")
  mod_log <- update(mod_con, outcome_log ~ . , family = "gaussian")
  
  mdls <- list(con = mod_con, cat = mod_cat, log = mod_log)
  mdls  <- map(mdls, function (x) {
    map(list(coef = coef, vcov = vcov), function (y) y(x))
  })
  mdls
})
sim_data$con <- map(res, function (x) x$con)
sim_data$cat <- map(res, function (x) x$cat)
sim_data$log <- map(res, function (x) x$log)

saveRDS(sim_data, file = "scratch_data/simulated_data.Rds")
