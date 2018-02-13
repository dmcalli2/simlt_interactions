library(INLA)
library(metRology)
library(rjags)


### First reproduce normal distribution
x <- seq(-5, 5, 0.01)
# x_dt <- dt.scaled(x, df = 3, mean = 1, sd = 1)
PriorFunctionT <- function(x) {
  res <- dt.scaled(x, df = 40, mean = 0, sd = (1/10)^0.5)
  log(res)
}

PriorFunctionN <- function(x) {
  res <- dt.scaled(x, df = 1000, mean = 0, sd = (1/10)^0.5)
  log(res)
}
prior.matrixN <- cbind(x, PriorFunctionN(x))
prior.matrixT <- cbind(x, PriorFunctionT(x))
plot(prior.matrixT, col = "red")
points(prior.matrixN)

prior.table <-  paste(c("table:", prior.matrixN),
                      sep = "", collapse = " ")
n <- 100

x <- sample(0:1, n, replace = TRUE)
y <- rnorm(n) + x*rnorm(n, 20)
tapply(y, x, mean)

my_data <- data.frame(x = x, y = y, y_prec = 10)

## Built-in prior
myform <- y ~ -1 + f(x, model = "iid", hyper = list(prec = list(prior = "gaussian", 
    param = c(mean = -5, prec = 10))))
dc1 <- inla.make.lincomb(x = 1) 

rr <- inla(myform,
           data = my_data,
           lincomb = dc1,
           control.family = list(hyper = list(prec = list(fixed = TRUE, initial = 0))),
           scale = my_data$y_prec)
rr$summary.random
# summary(rr)

# 
# ## Home made prior
# myform2 <- y ~ -1 + f(x, model = "iid", hyper = list(prec = list(prior = prior.table)))
# rr2 <- inla(myform2,
#            data = my_data,
#            control.family = list(hyper = list(prec = list(fixed = TRUE, initial = 0))),
#            scale = my_data$y_prec)
# rr2$summary.random
# summary(rr2)
# hist(rr2$marginals.random$x$index.2[,"x"])

modelstring <- "
model{
for (z in 1:100){
   # likelihood
   y[z] ~ dnorm(mu[z], y_prec[z])
   # link and linear predictor
   mu[z] <- beta1*x[z]
} #end of trials
   # Prior
   beta1 ~ dnorm(-5, 10)
}# end of model
"
writeLines(modelstring, con = "jags/t_dist.txt")

jags <- jags.model(paste0("jags/t_dist.txt"),
                       data = my_data,
                       n.chains = 2,
                       n.adapt = 1000)
  update(jags, 1000)
 data(LINE)
  LINE$recompile()
  res <- coda.samples(jags,
                c('beta1'),
                10000)
summary(res)
tapply(y, x, mean)
