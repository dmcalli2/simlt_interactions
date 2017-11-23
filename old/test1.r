print ("Hello world")
setwd("fship_sim")
mydf <- data.frame(a = 1:3)
saveRDS(mydf, file = "test1.rds")