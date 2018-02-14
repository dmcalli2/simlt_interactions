# 04_se_odds

library(tidyverse)

# Simulate log odds ratios for comorbidity, treatment and interaction effects
log_ors <- expand.grid(
  como = c(-0.5, 0, 0.5),
  allc = c(-0.2,-0.1, 0, 0.1, 0.2),
  actn = c(-0.2,-0.1, 0, 0.1, 0.2)
)

# Calculate log odds ratio for interaction-group versus no comorbidity, no treatment group
log_ors$como_allc <- log_ors$como + log_ors$allc + log_ors$actn

# Simulate number of people in each group
ns <- tibble(cept = 0.4,
             como = 0.1,
             allc = 0.4,
             como_allc = 0.1)
ns <- ns %>% 
  slice(rep((1:4), nrow(log_ors)))
ns[] <- map(ns, ~ 300*.x)

log_odds <- log_ors
log_odds[] <- lapply(log_ors, function (x) -2.19 + x)
log_odds$cept <- -2.19
prob <- log_odds %>%
  select(-actn)
prob[] <- lapply(prob, plogis)

x <- map2(prob ,
          ns,
          ~ round(.x*.y)) %>% 
  bind_cols()

## Estimate standard error for interaction using delta method
p_m <- as.matrix(prob[, c("cept", "como", "allc", "como_allc")])
n_m <- as.matrix(ns[, c("cept", "como", "allc", "como_allc")])
x_m <- n_m * p_m

VarLogOdds <- function(x, n) {
  a <- x
  b <- n-x
  (1/a + 1/b)
}
log_var_odds <- VarLogOdds(x_m, n_m)
log_actn_se <- rowSums(log_var_odds)^0.5
log(4/96)
VarLogOdds(4, 96)^0.5
exp(-3.17 + 2*0.51)
log_actn_se[1]

## Reshape data to long to run logistic regression as check on calculations
x_lng <- x %>% 
  mutate(scenario = seq_along(allc)) %>% 
  gather(key = "key", value = "x", -scenario)
n_lng <- ns %>% 
  mutate(scenario = seq_along(allc)) %>% 
  gather(key = "key", value = "n", -scenario)

xn <- x_lng %>% 
  inner_join(n_lng) %>% 
  arrange(scenario, key) %>% 
  mutate(allc = if_else(key %in% c("allc", "como_allc"), 1, 0),
         como = if_else(key %in% c("como", "como_allc"), 1, 0))

xn<- xn %>% 
  arrange(scenario, como, allc) %>% 
  select(-key)

xn_mdl <- xn %>% 
  group_by(scenario) %>% 
  nest()

## RUn glm on each dataset
models <- map(xn_mdl$data, function(sset) (glm(cbind(x, n-x) ~ allc*como, data = sset, 
                                              family = "binomial")))
models_smry <- map(models, broom::tidy)
models_smry <- map(models_smry, function(x) {
  x %>% 
    mutate(term = c("cept", "allc", "como", "actn"))
})
names(models_smry) <- stringr::str_pad(string = seq_along(models_smry), width = 2, pad = "0")
models_smry <- bind_rows(models_smry, .id = "scenario")

## Summarise model results
est <- models_smry %>% 
  as_tibble() %>% 
  select(-statistic, -p.value, -std.error) %>%
  spread(key = term, value = estimate)

se <-  models_smry %>% 
  as_tibble() %>% 
  select(-statistic, -p.value, -estimate) %>%
  spread(key = term, value = std.error)

est <- bind_cols(log_ors, est)
est <- est[, sort(names(est))]
models_smry <- models_smry[, sort(names(models_smry))]


se <- bind_cols(delta_method = log_actn_se, se)
plot((se$delta_method + se$actn)/2, se$delta_method-se$actn, ylim = c(-0.5, 0.5))

rm(list = ls())

## Calculate variation in SE by numbers ----
n_study_group <- c(300, 500, 1000)/2
prev_com <- seq(0.1, 0.3, 0.05)
event_rate_base <- seq(0.05, 0.20, 0.05) #
log_or_como <- seq(-0.7, 0.7, 0.1)
log_or_allc <- seq(-0.5, 0.5, 0.1)
log_or_actn <- seq(-0.5, 0.5, 0.1)

scenarios <- expand.grid(n_study_group, prev_com, event_rate_base, log_or_como,
                     log_or_allc, log_or_actn)
names(scenarios) <- c("n_study_group", "prev_com", "event_rate_base", "log_or_como",
                  "log_or_allc", "log_or_actn")
EventOdds <- function(x_curr, n_curr, log_or, n_new){
  odds_curr <- x_curr/(n_curr - x_curr)
  log_odds_new <- log(odds_curr) + log_or
  x_new <- plogis(log_odds_new) * n_new
  x_new
}
se_df <- scenarios %>% 
  transmute(n_cept = n_study_group * prev_com,
            n_como = n_study_group * (1-prev_com),
            n_allc = n_cept,
            n_actn = n_como,
            x_cept = n_cept*event_rate_base,
            x_como = EventOdds(x_cept, n_cept, log_or_como, n_como),
            x_allc = EventOdds(x_cept, n_cept, log_or_allc, n_allc),
            x_actn = EventOdds(x_cept, n_cept, log_or_allc + log_or_actn, n_allc))
se_df <- se_df %>% 
  mutate_at(vars(starts_with("x")), round)
se_df <- as.matrix(se_df)
x <- se_df[ , c("x_cept", "x_como", "x_allc", "x_actn")] %>% rowSums()

se_df[se_df ==0] <- 0.5 
se_df <- 1/se_df

scenarios <- scenarios %>% 
  mutate(se_actn = rowSums(se_df)^0.5,
         x = x)
rm(se_df)

library(ggplot2)

# scenarios <- scenarios %>% 
#                      mutate(log_or_como = paste0("log OR comorbidity ", log_or_como),
#                             log_or_allc = paste0("log OR treatment ", log_or_allc))

plot_grey <-  ggplot(scenarios , 
                 aes(x = event_rate_base, y = se_actn,
                     group = interaction(log_or_actn, log_or_como, log_or_allc))) + 
  geom_point() +
  geom_line(alpha = 0.2) +
  facet_grid(n_study_group ~ prev_com) +
  scale_y_continuous("se for log OR interaction") 

plot_actn <- plot_grey + aes(colour = log_or_actn)
plot_como <- plot_grey + aes(colour = log_or_como)
plot_allc <- plot_grey + aes(colour = log_or_allc)

pdf("figures/Plots of standard errors given itneractions.pdf", width = 10, height = 8)
plot_actn
plot_como
plot_allc
dev.off()

## Simplify
plot_actn2 <-  ggplot(scenarios %>% 
                        filter(log_or_como %>%  round() ==0,
                               log_or_allc %in% c(-0.5)), 
                 aes(x = event_rate_base, y = se_actn*(n_study_group*prev_com)^0.5,
                     colour = log_or_actn,
                     group = interaction(log_or_actn, log_or_allc))) + 
  geom_point() +
  geom_line(alpha = 0.2) +
  facet_grid(n_study_group ~ prev_com) +
  scale_y_continuous("se for log OR interaction x sqrt of number in one arm wiht comorbidity") 


# summarise standard errors by n, event rate in baseline, prevalence of
# comorbidity and size of interaction This means need to do around a mean of 10
# different standard errors per trial per simulation to cover the 165 different
# combinations per trial (based on the range of odds ratios for the 
# main treatment effect and comorbidity effect ( What we are assuming here is that a
# difference in standradr error for an individual trial of 0.025 is trivial
se_smry <- scenarios %>% 
  group_by(n_study_group, prev_com, event_rate_base, log_or_actn) %>% 
  summarise_at(vars(se_actn), .funs = funs(mean, sd, min, max, length)) %>% 
  ungroup() %>% 
  mutate(unq_se = round((max-min)/0.025, 0))
                         
stem(se_smry$unq_se)
mean(se_smry$unq_se)

## To simplify further, for each interaction, we need to run a mean of 53
## simulations per trial (as opposed to one simulation per trial for the
## continuous analysis). We could start with the top and bottom scenarios and
## work intelligently through until no differences are apparent.
se_smry_interaction <-  scenarios %>% 
  group_by(n_study_group, log_or_actn) %>% 
  summarise_at(vars(se_actn), .funs = funs(mean, sd, min, max, length)) %>% 
  ungroup() %>% 
  mutate(unq_se = round((max-min)/0.025, 0))
stem(se_smry_interaction$unq_se)
mean(se_smry_interaction$unq_se)
sum(se_smry_interaction$unq_se) # total number of simulations to cover all possibilities for 3 trials
unique(scenarios$log_or_actn)


