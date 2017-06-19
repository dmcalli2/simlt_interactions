// saved as simplest.stan
data {
  int n_coef;
  vector [n_coef] y; // estimated coefficients
  cov_matrix [n_coef] sigma; // covariance 
}
transformed data{
}
parameters {
  vector[n_coef] mu; 
}
transformed parameters {
}
model {
  y ~ multi_normal(mu, sigma);
  for (i in 1:8)
    mu ~ normal(0,5);
  mu[9] ~ normal(0,6);
  mu[10] ~ normal(0,6);
 }

