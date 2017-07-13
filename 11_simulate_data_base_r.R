# 11_simulate_data_base_r
# Initially generates random samples from normal distribution. THese will be the same for each scenario
# Next calcualtes overall group effect and variance for each component that will make-up the final scenarios
# eg (intercept, comorbidity effect, treatment allocation and interaction)

## Make same answer each time
set.seed(234)

### Determine number of samples needed for each scenario, do random sampling once for all analyses
### for all drug_classes
datasets <- 1000
drug_classes <- 6 # maximum number of classes
trials <- 10 # maximum number of trials per drug class

# Loop through parameters to keep it manageable, even though sampling 4 times
params <-  4 # intercept, treatment, 2 x comorbidities, 2 x interactions

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
# At this point using <50MB and got 16GB, just create them all

## Create first comparison of interest, fix all other components to the middle value
## Set all the other effects to null (overall group, treatment allocation and comorbidity)
## and vary the comorbidity-treatment interaction
which(varn_scen$intrcpt$btwn_class == 0.50 &
      varn_scen$intrcpt$btwn_trial == 0.50)
which(varn_scen$com1$btwn_class == 0.50 &
      varn_scen$com1$btwn_trial == 0.50)
which(varn_scen$alloc$btwn_class == 0.15 &
      varn_scen$alloc$btwn_trial == 0.15)

## Linear predictor for other characteristics, no intercept
t0_g0_no_intrcp <- varn_res$intrcpt[[4]] 
t1_g0_no_intrcp <- varn_res$intrcpt[[4]]                      + varn_res$alloc[[13]] 
t0_g1_no_intrcp <- varn_res$intrcpt[[4]] + varn_res$com1[[4]] 
# Note this one does not include the interaction
t1_g1_no_intrcp <- varn_res$intrcpt[[4]] + varn_res$com1[[4]] + varn_res$alloc[[13]]

## Linear predictor for interactions
t1_g1_each_no_intrcp <- lapply(varn_res$intrct1, function (x) x + t1_g1_no_intrcp)
rm(t1_g1_no_intrcp)

#### Linear predictor for each of the 4 groups (treated or not, allocated or not)
n_trial <- 1000
n_trial_arm <- n_trial * 0.5
cmrbd_prop <- 0.2
cmrbd_prop_arm    <- cmrbd_prop*0.5
cmrbd_prop_arm_no <- (1-cmrbd_prop)*0.5

cmrbd_arm_no <- round(n_trial_arm*(1-cmrbd_prop))
cmrbd_arm    <- round(n_trial_arm*cmrbd_prop)

## Can I calculate the overall linear predictor like this
tg <- lapply(t1_g1_each_no_intrcp, function (x) {
  t0_g0_no_intrcp*cmrbd_prop_arm_no + t1_g0_no_intrcp*cmrbd_prop_arm_no + 
  t0_g1_no_intrcp*cmrbd_prop_arm    + x*cmrbd_prop_arm
  })
## In order to force the overall outcome to be 20%
intrcpt <- lapply(tg, function(x) qlogis(0.2) - x)

## The intercepts are similar for each element in the list, so just take the mean
intrcpt <- sapply(intrcpt, identity, simplify = "array")
intrcpt <- apply(intrcpt, 1:3, mean)

t0_g0 <- t0_g0_no_intrcp + intrcpt
t1_g0 <- t0_g0_no_intrcp + intrcpt
t0_g1 <- t0_g1_no_intrcp + intrcpt
t1_g1_each <- lapply(t1_g1_each_no_intrcp, function(x) x + intrcpt)

## Binomial
## No treatment, no comorbidity
t0_g0_r   <- t0_g0
t0_g0_r[] <- rbinom(n = length(t0_g0), size = cmrbd_arm_no, prob = plogis(t0_g0))
## Treatment, no comorbidity
t1_g0_r   <- t1_g0
t1_g0_r[] <- rbinom(n = length(t1_g0), size = cmrbd_arm_no, prob = plogis(t1_g0))
## No treatment, comorbidity
t0_g1_r   <- t0_g1
t0_g1_r[] <- rbinom(n = length(t0_g1), size = cmrbd_arm, prob = plogis(t0_g1))
## Treatment, comorbidity
t1_g1_r_each   <- t1_g1_each
t1_g1_r_each <- lapply(t1_g1_each, function (x) {
  x_res <- x
  x_res[] <- rbinom(n = length(x), size = cmrbd_arm, prob = plogis(x))
  x_res
})

# Continuous, note lm and glm implement inverse variance (ie precision weights)
# All the studies should have the same standard error, ASSUMING that the 
# within study standard deviation is the same for all groups
# This is incompatible with the overall standard deviation being the same as in all groups
MakeSE <- function (n) 1/n^0.5
t0_g0_se <- MakeSE(n = cmrbd_arm_no)
t1_g0_se <- MakeSE(n = cmrbd_arm_no)
t0_g1_se <- MakeSE(n = cmrbd_arm)
t1_g1_se <- MakeSE(n = cmrbd_arm)

## Bind data into a dataframe to run the glms
new_array <- sapply(list(t0_g0_r, t0_g1_r, t1_g0_r, t1_g1_r_each[[1]]), identity,
                    simplify = "array")
new_df <- as.data.frame.table(new_array)
rm(new_array)
names(new_df) <- c("datasets", "dc", "trials", "grps", "r")
new_df <- new_df[order(new_df$datasets, new_df$dc, new_df$trials, new_df$grps),]

new_df <- cbind(new_df, com = c(0,1,0,1))
new_df <- cbind(new_df, alloc = c(0,0,1,1))
new_df <- cbind(new_df, n = c(cmrbd_arm_no, cmrbd_arm, cmrbd_arm_no, cmrbd_arm)  )
new_df$my_indx <- as.numeric(factor(paste(new_df$datasets, new_df$dc, new_df$trials, sep = "_")))

mdls <- by(data = new_df, INDICES = new_df[, "my_indx"],
   FUN = function (x){
     mdl <- glm(cbind(r, n - r) ~ com + alloc + com:alloc, data = x, family = "binomial")
     print(x$my_indx)
     list(coef(mdl), vcov(mdl))
     # list(coef = coef(mdl), varcov = vcov(mdl))}
   }, simplify = FALSE)
saveRDS(mdls, file = "scratch_data/vary_interaction_constant_others.Rds") # 5 MB
