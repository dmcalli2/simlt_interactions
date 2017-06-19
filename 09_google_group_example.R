#09 example for google group


mydata <- data.frame(y = c(0.24, 0.233, 0.222, 0.231,
                           0.215, 0.214, 0.213, 0.211,
                           0.154, 0.153, 0.151),
                     y_prec = c(444418, 366373, 764921, 722141,
                                295114, 591535, 497774, 1389837, 
                                956262, 1257805, 390475),
                     group = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3))

formula <- y ~ 1 + f(group, model = "iid")

mod1 <- inla(formula, family = "gaussian",
             data = mydata,
             scale = mydata$y_prec, )
summary(mod1)

## Extract overall mean
mod1$summary.fixed
## Extract mean for each group
mod1$summary.random
## Extract SD on mean
mod1$summary.hyperpar

