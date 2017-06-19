##08_analyse_with_INLA

library(tidyverse)
library(INLA)
library(forcats)

drugs_select <- "Airways"
source('02_arrange_data_arrays.R') 

## Select the first simulation
res <- res_list[[1]]
res <- res %>% 
  filter(drug_group_allocation == drugs_select) %>% 
  sample_frac(0.25)
rm(res_list, res_rag)

dep <- map_dbl(res$con, ~ .x$coef["dep:alloc"])
dep_prec <- map_dbl(res$con, ~ .x$prec_matrix["dep:alloc", "dep:alloc"])

res$my_drug_n <- as.numeric(as_factor(res$my_drug))
res$study_id_n <- as.numeric(as_factor(as.character(res$study_id)))


my_data <- data.frame(y = dep, y_prec = dep_prec, group = res$study_id_n)
dput(my_data)

mydata <- structure(list(y = c(0.212328690693095, 0.212879178793947, 0.211761674185828, 
0.211368272612108, 0.214514875041842, 0.214431264419052, 0.213332755925144, 
0.210925591511177, 0.183937577901477, 0.183304019573961, 0.210580996742981
), y_prec = c(444417.744720562, 366373.128269875, 764921.410454027, 
722140.673944748, 295114.45772035, 591534.531828539, 497774.057496792, 
1389836.67050187, 956261.725720684, 1257804.53027412, 390474.587406024
), group = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)), .Names = c("y", 
"y_prec", "group"), row.names = c(NA, -11L), class = "data.frame")

formula <- y ~ 1 + f(group, model = "iid") 

mod1 <- inla(myform1, family = "gaussian",
             data = res,
             scale = res$dep_prec)
myform2 <- dep ~ 1 + f(study_id_n, model = "iid") + 
                    f(my_drug_n, model = "iid")


summary(mod1)

# dnorm(dep, 0.25) 

# *** Code for Figure 5.13 top
global_intercept <- inla.rmarginal(1000,marg = mod1$marginals.fixed$`(Intercept)`)

MarginalsFx <- function(x = varname, model = mod1) {
  lngth <- length(unique(res[, x]))
  output <- matrix(NA, 1000, lngth)
      for(i in 1:lngth){
        output[,i] <- inla.rmarginal(1000, marg = model$marginals.random[[x]][[i]])
      }
  output
}

group_low <- MarginalsFx("group_low")
# group_high <- MarginalsFx("group_high")
  

group_low_effect <- global_intercept + group_low


drug_class_effect <- global_intercept + drug_class_specific_intercept


study_effect_quartiles <- t(apply(study_effect, MARGIN=2,
                function(x) quantile( x, probs= c(0.025,0.5,0.975))))

drug_class_effect_quartiles <- t(apply(drug_class_effect, MARGIN=2,
                function(x) quantile( x, probs= c(0.025,0.5,0.975))))
