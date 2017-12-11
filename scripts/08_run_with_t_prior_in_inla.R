library(INLA)
library(metRology)
library(rjags)


### First reproduce normal distribution
x <- seq(-5, 5, 0.01)
# x_dt <- dt.scaled(x, df = 3, mean = 1, sd = 1)
PriorFunction <- function(x) {
  res <- dnorm(x, 0, (1/50)^0.5)
  log(res)
}

prior.table <-  paste(c("table:", cbind(x, PriorFunction(x))),
                      sep = "", collapse = " ")
plot(x, exp(PriorFunction(x)))
n <- 100

x <- sample(0:1, n, replace = TRUE)
y <- rnorm(n) + x*rnorm(n, 100)
tapply(y, x, mean)

myform <- y ~ f(x, model = "iid", hyper = list(prec = list(prior = "gaussian", 
    param = c(mean = 0, prec = 50))))
rr <- inla(myform,
           data = data.frame(x = x, y = y))
rr$summary.random


myform2 <- y ~ f(x, model = "iid", hyper = list(prec = list(prior = prior.table)))
rr2 <- inla(myform2,
           data = data.frame(x = x, y = y))
rr2$summary.random


modelstring <- "
model{
for (z in 1:100){
   # likelihood
   y[z] ~ dnorm(wd[z], y_prec[z])
   # link and linear predictor
   # Shared priors (random effects)
   wd[z] ~ dnorm(wd_delta, wd_tau)
} #end of trials
   # Meta-analysis level priors
 
   wd_delta ~ dt(m, 1/s^2, d)
   wd_tau <-  1/wd_var
   wd_var ~ dlnorm(-3.23, 1/1.88^2)
   wd_sd <- wd_var^0.5
}# end of model
"
writeLines(modelstring, con = "jags/random_effects_new.txt")
