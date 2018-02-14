## Run glm models
library(metafor)
## Using rma fixed effects function from metafor package to check results

# INITAL code demonstrates for LM model summations of estimates and standard errors will give effects
# If reshape data to wide (4 treatment comparisons) these are vectorised and will be very fast
# Will need to do so AFTER applying N's and SDs as these differ by trial

# Simulate data
mydf <- expand.grid(alloc = 0:1, com = 0:1)
mydf$est <- rnorm(4)/4
mydf$sd <- abs(rnorm(4))
mydf$n <- c(50, 10, 50, 10)
mydf$se <- mydf$sd/mydf$n^0.5
mydf_sim <- mydf

## Actual data
reg1 <- diabetes_final_smpls_obs[1, c("cept", "como", "allc", "actn")]
se <- diabetes_final_smpls_obs[1, c("ncomo_se", "ycomo_se", "ncomo_se", "ycomo_se")]
prec <- diabetes_final_smpls_obs[1, c("ncomo_prec", "ycomo_prec", "ncomo_prec", "ycomo_prec")]
mydf <- data.frame(est = as.vector(t(reg1)),
                   se = as.vector(t(se)),
                   prec = as.vector(t(prec)))
mydf$com <- c(0, 1, 0, 1)
mydf$alloc <- c(0, 0, 1, 1)
mydf <- mydf[c(1,3,2,4),]

# RUn RMA model
mod1_rma <- rma(yi = est, vi = se^2, mods = ~ alloc*com, method = "FE", data = mydf)
summary(mod1_rma)

# Calculate contrasts
mydf$est_rma <- summary(mod1_rma)$b

mydf$est_manual <- c(mydf$est[1], 
                     mydf$est[2] - mydf$est[1],
                     mydf$est[3] - mydf$est[1],
                     mydf$est[4] - mydf$est[2] - mydf$est[3] + mydf$est[1])

# Calculate standard errors
mydf$se_rma <- summary(mod1_rma)$se
SEDiff <- function (x = c(1,1)) sum(x^2)^0.5

alloc_se <- SEDiff (mydf$se[1:2])
com_se   <- SEDiff (mydf$se[c(1,3)])
inter_se   <- SEDiff (mydf$se[1:4])
mydf$se_manual <- c(mydf$se[1], alloc_se, com_se, inter_se)

## Examine result
mydf

## Standard errors for all 161 trials * 1000 samples
diabetes_final$inter_calc <- (
  diabetes_final$ncomo_se^2 +
  diabetes_final$ycomo_se^2 +
  diabetes_final$ncomo_se^2 +
  diabetes_final$ycomo_se^2)^0.5

## Recover interaction
## Demonstrates that only depends on SE (function of n and sd) and interaction
## Don't need to worry about the rest of the variables so is a 
## 3 DIMENSIONAL SIMULATION (trial, drug and class variation)
diabetes_final_smpls_obs$inter_calc <- 
  diabetes_final_smpls_obs$actn -
  diabetes_final_smpls_obs$como -
  diabetes_final_smpls_obs$allc +
  diabetes_final_smpls_obs$cept
check <- cbind(calc = diabetes_final_smpls_obs$inter_calc,
               orig = diabetes_final_smpls$actn)



