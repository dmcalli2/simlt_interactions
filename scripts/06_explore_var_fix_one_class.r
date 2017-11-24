# explore varaicnes while fixing one class
library(tidyverse)

# 5 drug classes, one of which has the group mean
## In order to implement, need to correct every drug class samples by subtracting x[1]
## I don't add the group mean until later
set.seed(1234)
x <- rnorm(5, 0, 0.2)
x
sd(x)
mean(x)
x <- x - x[1]
x
sd(x)
mean(x)

# 5 drug classes, one of which differs by a known amount from the group mean
## In order to implement, need to correct every drug class samples by subtracting x[1]
## I don't add the group mean until later
## ONe consequence of doing it this way is that all of the drug classes
## have this relationship to the group mean (on average), but I think that is OK
set.seed(1234)
x <- rnorm(5, 0, 0.2)
x
sd(x)
x <- x - x[1] -0.2
x
sd(x)
mean(x)