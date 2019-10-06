
// Dynamic Logistic Regression
data {
  int D;
  real<lower=0> duration[D];
  real n[D];
  real y[D];
}

transformed data {

  real<lower=0, upper = 1> y_mu[D];
  
  for (t in 1:D) {  y_mu[t] = y[t]/n[t];}
}

parameters{
  
  real<lower=0, upper = 1> theta_1[D];
  real<lower=0> sigma_1;
  real<lower=0> sigma_2;
  
  // real theta_1_mu;
  // real theta_1_raw[D-1];
  
}

// transformed parameters{
//   
//   real<lower=0, upper = 1> theta_1[D];
//   
//   // for(i in 2:D) { 
//   //   theta_1[i] = theta_1_mu + sigma_1 * theta_1_raw[i];
//   // }
// }



model {
  
  target += normal_lpdf(sigma_1 | 0, 5);
  target += normal_lpdf(sigma_2 | 0, 5);
  
  target += normal_lpdf( theta_1[1]| 0.5, sigma_2);
  for (t in 2:D) {
    theta_1[t] ~ normal(theta_1[t-1], sigma_2);
  }

  target += normal_lpdf( y_mu | theta_1, sigma_1);

}

