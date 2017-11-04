## Run glm models
library(metafor)
## Using rma fixed effects function from metafor package to check results


# INITAL code demonstrates for LM model summations of estimates and standard errors will give effects
# If reshape data to wide (4 treatmetn comparisons) these are vectorised and will be very fast
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

#### For logistic regression, will need to fit the models
# But can re-use the same design matrix
mydf <- expand.grid(alloc = 0:1, com = 0:1)
dsgn <- model.matrix(~ mydf$alloc * mydf$com)

mydf$n <- c(450, 50, 450, 50)
mydf$prop <- runif(4,0.05, 0.2)
mydf$x <- round(mydf$prop * mydf$n)
mydf$prop <- mydf$x / mydf$n

a <- glm.fit(y = cbind(mydf$x, mydf$n), x = dsgn, family = binomial(link = logit))
b <- glm(cbind(mydf$x, mydf$n) ~ mydf$alloc* mydf$com, family = binomial(link = logit))

## QR factorization of model matrix
qr.X <- qr(dsgn)

## get estimated coefficient
b <- qr.qty(qr.X, y)
beta <- as.vector(backsolve(qr.X$qr, b))

## compute residuals (I don't use `qr.resid` and `qr.fitted`, though I can)
res <- as.vector(y - X %*% beta)

## residual standard error
se2 <- sum(res ^ 2) / (nrow(X) - qr.X$rank)

## full variance-covariance matrix
chol2inv(qr.X$qr) * se2

########### Calculate on PUTTY
# for linear model
res$n <- c(10, 50, 100,540)
res$se <- res$sd * res$n^-1

res2 <- as.data.frame(matrix(res$est, ncol = 4))
names(res2) <- paste("est", 1:4, sep = "_")
res3 <- as.data.frame(matrix(res$se, ncol = 4))
names(res3) <- paste("se", 1:4, sep = "_")
res <- cbind(res2, res3)
head(res)

SEDiff <- function (x, y) (x^2 + y^2)^0.5
res$alloc <- SEDiff(res$se_1, res$se_2)
res$cmrbd <- SEDiff(res$se_1, res$se_4)
res$inter <- (res$se_1^2 + res$se_2^2 + res$se_3^2 + res$se_4^2)^0.5
head(res)


###########
# for glm
# too slow. Use Wald approximation for odds ratio
system.time({
  # SImulate data
res$n <- c(10, 50, 100,540)
res$x <- rbinom(nrow(res), res$n, plogis(res$est))
# Odds ratio se for each component 
res$n_recip <- 1/res$n
res$x_recip <- 1/res$x
res$x_and_n_recip <- res$n_recip + res$x_recip


res <- as.data.frame(matrix(res$x_and_n_recip, ncol = 4))
names(res) <- paste("all_eff", 1:4, sep = "_")
SEDiff <- function(x, y) (x+y)^0.5
res$cept  <- res$all_eff_1^0.5
res$alloc <- SEDiff(res$all_eff_1, res$all_eff_2)
res$cmrbd <- SEDiff(res$all_eff_1, res$all_eff_4)
res$inter <- res$all_eff_1 + res$all_eff_2 + res$all_eff_3 + res$all_eff_4 
})
# user  system elapsed
#    0.15    0.00    0.15




