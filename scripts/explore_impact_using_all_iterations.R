# Examine variation in scenarios

scenarios <- readRDS("scratch_data/scenario_res")

set.seed(1234)
scenario_pick <- scenarios %>% 
  filter(scenario == sample(scenario, 1))


i <- runif(1, 1, 1000) %>% round()
one_example <- rnorm(10000, scenario_pick$mean[i], scenario_pick$sd[i])
all_possible_examples <-  rnorm(10000, scenario_pick$mean, scenario_pick$sd)

hist(all_possible_examples, col = "blue")
hist(one_example, add = TRUE, col = "white")
