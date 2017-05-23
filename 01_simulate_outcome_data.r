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
MakeCsv <- function (mypath = "data/sim_df1.xlsx") {
  # Functions to let me use excel rather than CSV files in constructing data description
  my_name <- round(runif(1, 10000, 90000),0)
  read_excel(path = mypath) %>%
    write_csv(paste0("scratch_data/", my_name, ".csv"))
  paste0("scratch_data/", my_name, ".csv")
}

# Precision matrix
PrecMtrx <- function(model) solve(vcov(model))


sim_data <- select (trials, drug_group_allocation, my_condition, my_drug, enrollment)
rm(trials)

## Create baseline variables, note centred age
def <- defRead(MakeCsv())
print(def)

# Create vector of effects
sim_data$dep_tx        <- rnorm(nrow(sim_data), mean =  0.05, sd = 0.01)
sim_data$pain_tx       <- rnorm(nrow(sim_data), mean =  0.10, sd = 0.05)
sim_data$alloc_tx      <- rnorm(nrow(sim_data), mean = -0.20, sd = 0.05)
sim_data$alloc_dep_tx  <- rnorm(nrow(sim_data), mean =  0.10, sd = 0.05)
sim_data$alloc_pain_tx <- rnorm(nrow(sim_data), mean =  0.02, sd = 0.01)

## Make list to loop through for dataframe
sim_list <- as.list(sim_data[, c("enrollment", "dep_tx", "pain_tx", "alloc_tx",
             "alloc_dep_tx", "alloc_pain_tx")])

# hist(map_dbl(sim_list, ~ .x["alloc_pain_tx"]))

## Loop through each trial
res <- pmap(sim_list, .f = function (enrollment, dep_tx, pain_tx, alloc_tx, alloc_dep_tx, alloc_pain_tx) {
  dt1 <- simstudy::genData(enrollment, def)

  ## Put age into quintiles to stratify by it
  dt1$age5 <- cut2(dt1$age, g = 5)

  ## Add treatent allocation, stratified by age and sex
  dt1 <- trtAssign(dt1, n = 2, balanced = TRUE, strata = c("age5", "sex"), 
            grpName = "alloc")
  
  ## Calculate outcome based on covariates and treatment, with interactions
  ## Each coefficient is simulated as coming from a normal distribution, except for age
  ## and sex
  trt_eff_mean <- with(dt1,
                       age*0.01 +
                         sex*0.05 +
                         dep*dep_tx +
                         pain*pain_tx +
                         alloc*alloc_tx +
                         alloc*dep*alloc_dep_tx +
                         alloc*pain*alloc_pain_tx
                       )
  dt1$outcome_con <- rnorm(nrow(dt1), mean = trt_eff_mean, sd = 1)
  
  ## Continuous outcome as normal distribution with mean = 0, and SD 1
  dt1$outcome_cat <- dt1$outcome_con > quantile(dt1$outcome, probs = 0.9)
  dt1$outcome_log <- log(dt1$outcome_con + 0.0001- min(dt1$outcome_con))
  #dt1$outcome_con <- dt1$outcome_con + rnorm(nrow(dt1), 0, 0.2)
  
  ## Centre all dichotmous variables (age already centred)
  dt1[ , c("alloc", "sex", "dep", "pain")] <- lapply(dt1[ , c("alloc", "sex", "dep", "pain")],
                                                     function(x) x - mean(x))
  
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
    map(list(coef = coef, vcov = vcov, prec_matrix = PrecMtrx), function (y) y(x))
  })
  mdls
})

## Save coefficients and variance-covariance matrix for each model to data
sim_data$con <- map(res, function (x) x$con)
sim_data$cat <- map(res, function (x) x$cat)
sim_data$log <- map(res, function (x) x$log)

## Check precision matrix
# PrecCheck <- function (model_type = sim_data$con) {
#   each_study <- map_lgl(model_type, ~ matrixcalc::is.positive.semi.definite(.x$prec_matrix %>%  round()))
#   (each_study)
# }
# which(!PrecCheck(sim_data$con))

saveRDS(sim_data, file = "scratch_data/simulated_data.Rds")
