
// Dynamic Logistic Regression
data {
  int D;
  real<lower=0> duration[D];
  int n[D];
  int y[D];
}

parameters{
  
  vector[D] theta_star;
  
  // real<lower=0, upper=1> theta_mu;
  real<lower=0> sigma_theta;
  real<lower=0> sigma_e;
  
  // real theta_1_mu;
  // real theta_1_raw[D-1];
  
}

transformed parameters{

  vector[D] theta;

  theta = cumulative_sum(theta_star) * sigma_theta;
}



model {
  
  target += normal_lpdf(sigma_theta | 0, 5);
  target += normal_lpdf(sigma_e | 0, 5);

  target += normal_lpdf( theta_star | 0, 1);

  target += binomial_logit_lpmf( y |  n, theta);

}



















