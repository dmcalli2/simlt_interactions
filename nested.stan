// saved as simplest.stan
data {
  int n_coef;
  int n_studies;
  vector [n_coef] y[n_studies] ; // estimated coefficients, for each study
  cov_matrix [n_coef] sigma[n_studies]; // covariance matrix for each study
}
transformed data{
}
parameters {
  vector [n_coef] mu[n_studies]; 
   
}
transformed parameters {
  print(n_coef)
  print(n_studies)
}
model {
 for(j in 1:n_studies){
      y [j,]~ multi_normal(mu[j,], sigma[j,,]);
 //   for (i in 1:8){
 //    mu[i,j] ~ normal(0,5);
 //   }
 //   mu[9, j] ~ normal(mu9, 6);
 //   mu[10, j] ~ normal(mu10, 6);
 // }
 // mu9 ~ normal(0,1)
 // mu10 ~ normal(0,1)
 }
}
