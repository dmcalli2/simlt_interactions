# 11_simulate_data_base_r
# Initially generates random samples from normal distribution. THese will be the same for each scenario
# Next Calculates overall group effect and variance for each component that will make-up the final scenarios
# eg (intercept, comorbidity effect, treatment allocation and interaction)
library(tidyverse)
## Read in diabetes trials
load("../Trial_identify/clinical_trials_august_2017/scratch_data/data_for_simulation.Rdata")

# Remove single alpha-glucosidase inhibitor trial
diabetes_final <-  
  subset(diabetes_final, diabetes_final$atc_5 != "A10BF")


## Make same answer each set of classes (not each class)
set.seed(1234)

### Each simulation scenario, overall effects
main_scen <- expand.grid(
  cept = 0,
  como = c(-0.2,-0.1, 0, 0.1, 0.2),
  allc = c(-0.2,-0.1, 0, 0.1, 0.2),
  actn = c(-0.2,-0.1, 0, 0.1, 0.2)
)
### Each simulation scenario, variation around the overall effect
varn_scen <- list(cept = c(0.25, 0.5),
                     como = c(0.25, 0.5),
                     allc = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     actn = c(0.05, 0.10, 0.15, 0.20, 0.25))

### Determine number of samples needed for each scenario, do random sampling once for all analyses
### for all drugs, classes and trial

# Simulate variation around wider group (A10B) effect at trial, drug and class level
SampleVarn <- function(item, smpl_lvls, n_iter = 1000){

  res <- lapply(smpl_lvls, function (x) {
    my_mtrx <- matrix(nrow = n_iter, ncol = length(item))
    for(j in seq_along(item)){
      my_mtrx[, j] <- rnorm(n_iter, 0, x)
    }
    rownames(my_mtrx) <- 1:n_iter
    colnames(my_mtrx) <- item
    my_mtrx
  })
  names(res) <- smpl_lvls
  res <- lapply(res, function(x) stack(as.data.frame(x)))
  res_names <- lapply(res, function(x) x$ind)[[1]]
  res <- do.call(cbind, lapply(res, function(x) x[, "values"]))
  res <- as.data.frame(res)
  res$item <- res_names
  res
}
trial <- lapply(list(cept = c(0.25, 0.5),
                     como = c(0.25, 0.5),
                     allc = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     actn = c(0.05, 0.10, 0.15, 0.20, 0.25)),
                function(x) {
                  SampleVarn(diabetes_final$nct_id, x)
                }
)

drug <- lapply(list(cept = c(0.25, 0.5),
                     como = c(0.25, 0.5),
                     allc = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     actn = c(0.05, 0.10, 0.15, 0.20, 0.25)),
                function(x) {
                  SampleVarn(unique(diabetes_final$drug), x)
                }
)

atc5 <- lapply(list(cept = c(0.25, 0.5),
                     como = c(0.25, 0.5),
                     allc = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     actn = c(0.05, 0.10, 0.15, 0.20, 0.25)),
                function(x) {
                  SampleVarn(unique(diabetes_final$atc_5), x)
                }
)

## Combine variation around A10B effects for intercept, comorbidity, allocation and itneraction
## into single estimate from trial, drug and class (A10BA, A10BB etc)
ExtractFx <- function(filename, filename_label = filename){
  x <- get(filename)
  x <- x$actn
  names(x)[!names(x) == "item"] <- paste(filename, 
                                             names(x)[!names(x) == "item"], sep = "_")
  names(x)[names(x)=="item"] <- filename_label
  x
}
atc5 <- ExtractFx("atc5", "atc_5")
drug <- ExtractFx("drug")
trial <- ExtractFx("trial", "nct_id")


trial <- trial %>% 
  arrange(nct_id)

drug <- drug %>% 
  inner_join(diabetes_final %>%  select(nct_id, drug)) %>% 
  arrange(nct_id) %>% 
  select(-nct_id)

atc5 <- atc5 %>% 
  inner_join(diabetes_final %>%  select(nct_id, atc_5)) %>% 
  arrange(nct_id) %>% 
  select(-nct_id)

## Combine all interaction effects
actn <- bind_cols(atc5, trial, drug) %>% 
  arrange(atc_5, drug, nct_id) %>% 
  mutate(iteration = rep(1:1000, 161)) %>% 
    select(atc_5, drug, nct_id, iteration, everything())

## Examine interaction effects at drug-class level
actn %>%
  select(iteration, atc_5, starts_with("atc5")) %>% 
  distinct() %>%
  group_by(atc_5) %>% 
  summarise_at(vars(starts_with("atc5")), mean)

atc5 %>%
  mutate(iteration = rep(1:1000, 161)) %>% 
  select(iteration, atc_5, starts_with("atc5")) %>% 
  distinct() %>%
  (function(x) {
    print(nrow(x))
    x
  }) %>% 
  group_by(atc_5) %>% 
  summarise_at(vars(starts_with("atc5")), mean)

## Fix the mean for a single ATC5-level class, choose J, so that it has a KNOWN
## difference compared to the overall mean Make two versions, a no difference
## and a large difference to see effects fo shrinkage
## Note that he overall effect (at ATC4 level) which has not yet been added
## First take unique value for each class for each iteration

atc_5_names <- c("atc5_0.05", "atc5_0.1", "atc5_0.15", "atc5_0.2", "atc5_0.25")

actn_one_fixed <- actn %>% 
  filter(atc_5 == "A10BJ") %>% 
  distinct(iteration, atc5_0.05, atc5_0.1, atc5_0.15, atc5_0.2, atc5_0.25) 

## Rename to disitinguish from other classes
names(actn_one_fixed)[-1] <- paste0("fx", atc_5_names)

## Join BJ class to other dataframe
actn_fixed <- actn %>% 
  inner_join(actn_one_fixed)

## Add fixed class effects
actn_fixed[, paste(names(actn_one_fixed)[-1], "f_same")] <- map2(actn_fixed[, atc_5_names],
                     actn_fixed[, names(actn_one_fixed)[-1]],
                     ~ .x - .y)
actn_fixed[, paste(names(actn_one_fixed)[-1], "f_diff")] <- map2(actn_fixed[, atc_5_names],
                     actn_fixed[, names(actn_one_fixed)[-1]],
                     ~ .x - .y - 0.2)

saveRDS(actn, file = "scratch_data/interactn_opts.Rds")
save(diabetes_final, file = "scratch_Data/interaction_opts_ordering.Rds")
