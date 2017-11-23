## Experiment with matrices

mydims <- c(2,3,4)

A <- array(paste0("a", 1:24), dim = mydims)
B <- array(paste0("b", 1:24), dim = mydims)
C <- array(paste0("c", 1:24), dim = mydims)
D <- array(paste0("d", 1:24), dim = mydims)

            
c(rbind(c(A), c(B), c(C), c(D)))

FlattenLong <- function(array1, array2, array3, array4){
  c(rbind(c(array1), c(array2), c(array3), c(array4)))
}
FlattenLong(A, B, C, D)

system.time({
  for(i in 1:1000){
new_array <- sapply(list(A, B, C, D), identity,
                    simplify = "array")
new_df <- as.data.frame.table(new_array)
}
})

system.time({
  for(i in 1:1000){
FlattenLong(A, B, C, D)
}
})


dt <- paste0("dt", 1:15)
cl <- paste0("cl", 1:6)
tr <- paste0("tr", 1:10)
my_array <- array(NA, dim = c(15, 6, 10))
for (i in 1:15){
  for(j in 1:6){
    for(k in 1:10){
      my_array[i,j,k] <- paste(dt[i], cl[j], tr[k], sep = "_")
    }
  }
}

my_vect <- c(my_array)
my_array_new <- array(my_vect, dim = c(15, 6, 10))
all(my_array_new == my_array)
