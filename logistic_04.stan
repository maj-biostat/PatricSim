
// Dynamic Logistic Regression
data {
  int D;
  real<lower=0> duration[D];
  int n[D];
  int y[D];
  real<lower=0> tau0;
}

transformed data {
  real<lower=0> duration_diff[D-1];
  
  for(i in 2:D) { duration_diff[i-1] = duration[i] - duration[i-1]; }
}

parameters {
  vector[D] eta_star;
  
  real<lower=0> tau;
}

transformed parameters {
  vector[D] eta;
  
  // the eta are the log odds of recovery accumulated from the 
  // first eta_star and then the noise that is represented by
  // eta_star[2:D]
  // so we get log odds as eta_1 + sum_2^D eta_star
  eta[1] = eta_star[1];
  // here, cumulative_sum is returning an array of length D-1
  eta[2:D] = eta[1] + cumulative_sum(eta_star[2:D]);
  
  // print("");
  // for(i in 1:D) { print("eta[i]      : ", eta[i]); }
  // for(i in 1:D) { print("eta_star[i] : ", eta_star[i]); }
  
}

model {
  // The first eta_star represents the starting point
  // eta_star[1] ~ normal(0, tau0);
  // eta_star[2:D] ~ normal(0, tau);
  
  target += normal_lpdf( eta_star[1]| 0, tau0);
  for(i in 2:D) { 
    // .* is non-vector multiplication
    target += normal_lpdf( eta_star[i]| 0, duration_diff[i-1] .* tau);
  }
  
  
  target += student_t_lpdf(tau | 3, 0, 2.5);
  
  target += binomial_logit_lpmf(y | n, eta);
  
  // tau ~ student_t(3, 0, 2.5);
  // // log binomial probability of y success given logit scaled
  // // chance of success eta
  // y ~ binomial_logit(n, eta);
}

