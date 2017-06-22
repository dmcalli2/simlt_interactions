# compare model where model full variance-covariance matrix and separate effects


library(tidyverse)
library(coda)
library(purrr)
library(stringr)
library(ggplot2)

models_want <- c("cat_dc_nest_inform_effect1.Rds",
                 "cat_dc_nest_inform_effect2.Rds",
                 "cat_dc_nest_inform_effect3.Rds", 
                 "cat_dc_nest_inform_effect4.Rds", "cat_dc_nest_inform_effect5.Rds", 
                 "cat_dc_nest_inform_effect6.Rds", "cat_dc_nest_inform_sep_effect1.Rds", 
                 "cat_dc_nest_inform_sep_effect2.Rds", "cat_dc_nest_inform_sep_effect3.Rds", 
                 "cat_dc_nest_inform_sep_effect4.Rds", "cat_dc_nest_inform_sep_effect5.Rds", 
                 "cat_dc_nest_inform_sep_effect6.Rds", 
                 "con_dc_nest_inform_effect1.Rds", 
                 "con_dc_nest_inform_effect2.Rds",
                 "con_dc_nest_inform_effect3.Rds", 
                 "con_dc_nest_inform_effect4.Rds",
                 "con_dc_nest_inform_effect5.Rds", 
                 "con_dc_nest_inform_effect6.Rds",
                 "con_dc_nest_inform_sep_effect1.Rds", 
                 "con_dc_nest_inform_sep_effect2.Rds",
                 "con_dc_nest_inform_sep_effect3.Rds", 
                 "con_dc_nest_inform_sep_effect4.Rds",
                 "con_dc_nest_inform_sep_effect5.Rds", 
                 "con_dc_nest_inform_sep_effect6.Rds")


ExtractParam <- function(model_name){
  res  <- readRDS(paste0("model_summaries/",model_name))
  res <- do.call(cbind, res)
  res <- res[c("dep", "pain"),]
  row_res <- row.names(res)
  res <- as.tibble(res)
  res$param <- row_res
  res
}

res <- map(models_want, ExtractParam)
names(res) <- models_want
res <- bind_rows(res, .id = "model")
res$outcome <- str_detect(res$model, "cat") 
res$outcome <- ifelse(res$outcome, "categorical", "countinuous")
res$effect_num <- str_extract(res$model, "[0-9]") %>%  parse_integer()
res$model_type <- str_detect(res$model, "sep")
res$model_type <- ifelse(res$model_type, "Separate", "Full")
res <- res %>% 
  arrange(outcome, param, effect_num)


res_quant <- res %>% 
  select(model_type, outcome, effect_num, param, `2.5%`, `25%`, `50%`, `75%`, `97.5%`) %>% 
  gather(key = "var", value = "value", `2.5%`, `25%`, `50%`, `75%`, `97.5%`)

plot1 <- ggplot(res_quant, aes(x = var, y = value, colour = model_type)) + 
  geom_point(position= position_dodge(width = 0.5)) + 
  facet_grid(effect_num ~ outcome + param) +
  scale_x_discrete("Quantiles of posterior") +
  scale_y_continuous("Full MVN model and separate models, log-odds ratio")
plot1

res_quant_spread <- res_quant %>%
  spread(key = model_type, value = value) %>% 
  mutate(full_minus_sep = Full - Separate)
mean(res_quant_spread$full_minus_sep >0) # 54% greater than zero, so not apparently systematic

plot2 <- ggplot(res_quant_spread, aes(x = var, y = full_minus_sep)) + 
  geom_point() + 
  facet_grid(effect_num ~ outcome + param) +
  scale_x_discrete("Quantiles of posterior") +
  scale_y_continuous("Full MVN model versus separate models, log-odds ratio")
plot2

pdf("compare full mvn versus separate likelihoods.pdf")
plot1
plot2
dev.off()
