## Run glm models

# Note thet for LM model simple summations will give effects

mydf <- expand.grid(alloc = 0:1, com = 0:1)
mydf$est <- rnorm(4)/4
mydf$sd <- abs(rnorm(4))
mydf$n <- c(50, 10, 50, 10)
mydf$se <- mydf$sd/mydf$n^0.5


mod1 <- lm(est ~ alloc*com, data = mydf, weights = n)
summary(mod1)

library(metafor)

mod1_rma <- rma(yi = est, vi = se^2, mods = ~ alloc*com, method = "FE", data = mydf)
summary(mod1_rma)


mydf$se_rma <- summary(mod1_rma)$se


SEDiff <- function (x = c(1,1)) sum(x^2)^0.5
SEDiff(c(3,4)) # shoudl be 5

alloc_se <- SEDiff (mydf$se[1:2])
com_se   <- SEDiff (mydf$se[c(1,3)])
inter_se   <- SEDiff (mydf$se[1:4])
SEDiff(c(0.714923, 0.32773))
mydf$se_manual <- c(mydf$se[1], alloc_se, com_se, inter_se)

mydf
