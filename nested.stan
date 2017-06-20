// saved as nested.stan
data {
  int <lower=1> n_coef;
  int <lower=1> n_studies;
  vector [n_coef] y[n_studies] ; // estimated coefficients, for each study
  cov_matrix [n_coef] sigma[n_studies]; // covariance matrix for each study
}
transformed data{
}
parameters {
  vector [n_coef] mu[n_studies]; 
  real mu9;
  real mu10;
  real <lower=0> mu9_var;
  real <lower=0> mu10_var;
}
transformed parameters {
  print(n_coef)
  print(n_studies)
}
model {
 for(j in 1:n_studies){
    y [j,]~ multi_normal(mu[j,], sigma[j,,]);
    mu[j, 1] ~ normal(0, 5^2);  // intercept - study specific
    for (i in 2:8){ // next 9 coefficients, not intercept - study specific
      mu[j,i] ~ normal(0, 2^2);
   }
    mu[j, 9] ~ normal(mu9, mu9_var);  // one shared estimate for every study - dep
    mu[j, 10] ~ normal(mu10, mu10_var); // one shared estimate for every study - pain
 }
 mu9 ~ normal(0, 5^2); // dep
 mu10 ~ normal(0, 5^2); // pain
 mu9_var ~ cauchy(0, 5); // between study variance in depression:tx interaction
 mu10_var ~ cauchy(0, 5); // as above for pain:tx interaction
}

