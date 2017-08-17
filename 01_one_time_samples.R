# 11_simulate_data_base_r
# Initially generates random samples from normal distribution. THese will be the same for each scenario
# Next Calculates overall group effect and variance for each component that will make-up the final scenarios
# eg (intercept, comorbidity effect, treatment allocation and interaction)


############## 1 Make samples for all scenarios
## Make same answer each time
set.seed(234)

### Determine number of samples needed for each scenario, do random sampling once for all analyses
### for all drug_classes
datasets <- 1000
drug_classes <- 6 # maximum number of classes
trials <- 10 # maximum number of trials per drug class

# Loop through parameters to keep it manageable, even though sampling 4 times
params <-  4 # intercept, treatment, 1 x comorbidities, 1 x interactions

param_smpls <- lapply(1:params, function (x) {
  ### Create array for classes, note 10 identical slices
  dc_var <-
    matrix(rnorm(datasets * drug_classes), ncol = drug_classes)
  dc_var <- lapply(1:trials, function(x)
    dc_var)
  dc_var <- sapply(dc_var, identity, simplify = "array")
  
  ## Create array for trials
  trl_var <-
    array(rnorm(datasets * drug_classes * trials),
          dim = c(datasets, drug_classes, trials))
  
  list(dc_var = dc_var, trl_var = trl_var)
})
names(param_smpls) <- c("intrcpt", "alloc", "com1", "intrct1")

## A few seconds have all samples need for each trial, and it is only 4MB
## It is 72,000 samples from random normal distributions,
## some of which (the drug class ones) are repeated

### Each simulation scenario, overall effects
main_scen <- expand.grid(
  intrcpt = 0,
  com1    = c(-0.2,-0.1, 0, 0.1, 0.2),
  alloc   = c(-0.2,-0.1, 0, 0.1, 0.2),
  intrct1 = c(-0.2,-0.1, 0, 0.1, 0.2)
)

### Each simulation scenario, variation at class and trial level
varn_scen <- list(
  intrcpt = expand.grid(btwn_class = c(0.25, 0.5),
                      btwn_trial = c(0.25, 0.5)),
## Comorbidity effect 1
  com1    = expand.grid(btwn_class = c(0.25, 0.5),
                      btwn_trial = c(0.25, 0.5)),
## Allocation
  alloc   = expand.grid(btwn_class = c(0.05, 0.10, 0.15, 0.20, 0.25),
                      btwn_trial = c(0.05, 0.10, 0.15, 0.20, 0.25)),
## Interactions effect
  intrct1 = expand.grid(btwn_class = c(0.05, 0.10, 0.15, 0.20, 0.25),
                      btwn_trial = c(0.05, 0.10, 0.15, 0.20, 0.25)))

## Component of linear predictor for combinations of variance, 1000 samples for each scenario
if(!all(names(varn_scen) %in% names(param_smpls))) print ("Warning mismatch in names")
varn_res <- lapply(names(varn_scen), function (param){ ## Loop through each paramter
  lapply(seq_along(varn_scen[[param]]$btwn_class), function(i){ ## Loop through each combination
      varn_scen[[param]]$btwn_class[i] * param_smpls[[param]]$dc_var + # between class
      varn_scen[[param]]$btwn_trial[i] * param_smpls[[param]]$trl_var  # between trial
  })
})
names(varn_res) <- names(varn_scen)
cumprod(sapply(varn_res, length) ) # gives us 10,000 scenarios alone
# At this point using <50MB
rm(list = c("param_smpls", "params"))

 # Reshapes arrays to a vector, where each array contributes one element in turn
  # The vector goes a1, b1, c1, d1, a2, b2, c2, d2, etc
  # First each array is collapsed, it goes columns, rows, slices etc
  # Then a matrix is created, which is again collapsed column wise
  # Final ordering is dataset 1:1000 (rows), for class 1(cols), for trial 1(slices),
  # then dataset 1:1000 for class 2, trial 1
  # then dataset 1:1000 for class 3, trial 1  etc
varn_res <- lapply(varn_res, function (x) { # apply to each parameter
  lapply(x, c)}) # apply to each combination
