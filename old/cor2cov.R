# compare model where model full variance-covariance matrix 


library(tidyverse)
library(forcats)

drugs_select <- "Airways"
source('02_arrange_data_arrays.R') 

## Select the first simulation
res <- res_list[[1]]
res <- res %>% 
  filter(drug_group_allocation == drugs_select)
rm(res_list, res_rag)

coef <- map(res$con, ~ .x$coef)
covar <- map(res$con, ~ .x$vcov)

covar1 <- covar[[1]]
sds <- diag(covar1)^0.5
cor1 <- cov2cor(covar1)

covar1_rec <- t(sds) %*% cor1

