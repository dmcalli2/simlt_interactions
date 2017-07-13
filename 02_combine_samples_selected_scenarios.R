# 02_combine_samples_selected_scenarios

source("01_one_time_samples.R")

### Proportion with comorbidity
cmrbd_prop <- 0.2
cmrbd_prop_arm    <- cmrbd_prop*0.5
cmrbd_prop_arm_no <- (1-cmrbd_prop)*0.5

## Linear predictor for other characteristics, no intercept as will calculate this for whole 
## sample not for baseline group

## Loop through each combination
cmbns <- c(sapply(varn_scen, nrow), main_eff = nrow(main_scen))
count <- 0
cmbns[] <- c(1,1,1,2,5)

system.time({
  for(cpt_i in 1:cmbns["intrcpt"]){
    for(com_i in 1:cmbns["com1"]){
        for(alloc_i in 1:cmbns["alloc"]){
            for(rct_i in 1:cmbns["intrct1"]){
              new_dfs <- vector(length = cmbns["main_eff"], mode = "list")
               for(me_i in 1:cmbns["main_eff"]){
               main_eff <- main_scen[me_i, , drop = FALSE]
t0_g0_no_intrcp <- main_eff$intrcpt + varn_res$intrcpt[[cpt_i]]  
t1_g0_no_intrcp <- t0_g0_no_intrcp + main_eff$alloc + varn_res$alloc[[alloc_i]] 
t0_g1_no_intrcp <- t0_g0_no_intrcp + main_eff$com1  + varn_res$com1[[com_i]] 
t1_g1_no_intrcp <- t0_g0_no_intrcp +
  main_eff$alloc   + varn_res$alloc[[alloc_i]] +
  main_eff$com1    + varn_res$com1[[com_i]] +
  main_eff$intrct1 + varn_res$intrct1[[rct_i]]

## Note above the same for each comorbidity proportion

## Can I calculate the overall linear predictor like this
tg <- t0_g0_no_intrcp*cmrbd_prop_arm_no + t1_g0_no_intrcp*cmrbd_prop_arm_no + 
      t0_g1_no_intrcp*cmrbd_prop_arm    + t1_g1_no_intrcp*cmrbd_prop_arm

## In order to force the overall outcome to be 20%
intrcpt <- qlogis(0.2) - tg

t0_g0 <- t0_g0_no_intrcp + intrcpt
t1_g0 <- t0_g0_no_intrcp + intrcpt
t0_g1 <- t0_g1_no_intrcp + intrcpt
t1_g1 <- t1_g1_no_intrcp + intrcpt

#### Linear predictor for each of the 4 groups (treated or not, allocated or not)
n_trial <- 1000
cmrbd_arm_no <- round(n_trial * cmrbd_prop_arm_no)
cmrbd_arm    <- round(n_trial * cmrbd_prop_arm)

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
new_array <- sapply(list(t0_g0, t0_g1, t1_g0, t1_g1), identity,
                    simplify = "array")
new_df <- as.data.frame.table(new_array)
rm(new_array)
names(new_df) <- c("datasets", "dc", "trials", "grps", "lp")

## Order and add linear predictor and proportions
new_df <- new_df[order(new_df$datasets, new_df$dc, new_df$trials, new_df$grps),]
new_df$indx <- 1:nrow(new_df)

## Save variable labels only once, same for every loop
count <- count + 1                        
# if(count ==1) saveRDS(new_df[ , c("datasets", "dc", "trials", "grps", "indx")],
#                       file = "scratch_data/dataset_0_labels.Rds")
print(count)
## Create a list of dataframes for each combination of main effects
new_dfs[[me_i]] <- cbind(new_df[, c("indx", "lp")],
                com = c(0,1,0,1),
                alloc = c(0,0,1,1),
                prop = c(cmrbd_prop_arm_no, cmrbd_prop_arm, cmrbd_prop_arm_no, cmrbd_prop_arm))
        }
saveRDS(new_dfs, file = paste("scratch_data/dataset", cpt_i, "_", com_i, "_", alloc_i, "_", rct_i, ".Rds"))
      }
    }
  }
  }
  print(paste0(cpt_i, "_", com_i, "_", alloc_i, "_", rct_i,"_", me_i))
})


                        
#                 n = c(cmrbd_arm_no, cmrbd_arm, cmrbd_arm_no, cmrbd_arm))
# new_df$r  <- rbinom(n = length(new_df$lp), size = new_df$n, prob = plogis(new_df$lp))
# new_df$se <- 1/new_df$n^0.5

# 
# # Run regression models
# new_df$my_indx <- as.numeric(factor(paste(new_df$datasets, new_df$dc, new_df$trials, sep = "_")))
# mdls <- by(data = new_df, INDICES = new_df[, "my_indx"],
#    FUN = function (x){
#      mdl <- glm(cbind(r, n - r) ~ com + alloc + com:alloc, data = x, family = "binomial")
#      print(x$my_indx)
#      list(coef(mdl), vcov(mdl))
#      # list(coef = coef(mdl), varcov = vcov(mdl))}
#    }, simplify = FALSE)
# saveRDS(mdls, file = "scratch_data/vary_interaction_constant_others.Rds") # 5 MB

  
## Combinations, 1.25 million, not including proportions
# scnds <- prod(sapply(varn_scen, nrow)) * 125
# scnds / (60*60*25)
