# 02_combine_samples_selected_scenarios
# setwd("fship_sim")

source("01_one_time_samples.R")

# ### Functions

## Linear predictor for other characteristics, not intercept as will take this from actual trials

## Place each combination of variances in a matrix (intrcpt, com1, alloc, intrct1)
## Where each matrix is saved as one of the 125 mean effect combinations 
cmbns <- c(sapply(varn_scen, nrow))

# Create empty matrix to hold values
mtrx_res <- array(NA, dim = c(cmbns, datasets*drug_classes*trials*4)) # 4 is for treatment arms

system.time({
for(me_i in 1:2){ # nrow(main_scen)
  main_eff <- main_scen[me_i, , drop = FALSE]
  for(cpt_i in 1:cmbns["intrcpt"]){
    for(com_i in 1:cmbns["com1"]){
        for(alloc_i in 1:cmbns["alloc"]){
            for(rct_i in 1:cmbns["intrct1"]){
              browser()
            t0_g0 <- varn_res$intrcpt[[cpt_i]]                            # No comorbidity, no tx
            t1_g0 <- t0_g0 + main_eff$alloc + varn_res$alloc[[alloc_i]]   # No comorbidity, yes tx
            t0_g1 <- t0_g0 + main_eff$com1  + varn_res$com1[[com_i]]      # Comorbidity, no tx
            t1_g1 <- t0_g0 +                                              # Comorbidity, tx
              main_eff$alloc   + varn_res$alloc[[alloc_i]] +
              main_eff$com1    + varn_res$com1[[com_i]] +
              main_eff$intrct1 + varn_res$intrct1[[rct_i]]

            ## Make single vector of linear predictors
            mtrx_res[cpt_i, com_i, alloc_i, rct_i, ] <- c(t0_g0, t1_g0, t0_g1, t1_g1)
            }
        }
    }
  }
saveRDS(mtrx_res, file = paste0("scratch_data/mean_effects", me_i, ".Rds"),
        compress = FALSE)
      }
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

  
# ## Combinations, 1.25 million, not including proportions
#  scnds <- prod(sapply(varn_scen, nrow)) * 125 * 0.3
#  scnds / (60*60*25)

