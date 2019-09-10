// generated with brms 2.9.0
functions {
}
data {
  int<lower=1> N;  // number of observations
  int y[N];  // number of successes
  int n[N];  // number trials each group
  int<lower=1> K;  // number of population-level effects (not incl intercept) 
  matrix[N, K] X;  // population-level design matrix (not incl intercept) 
  int prior_only;  // should the likelihood be ignored?
}
parameters {
  vector[K] b;  // population-level effects
  real b0;  // intercept
}

model {
  // N x K    and    K x 1 gives N x 1
  vector[N] mu = X * b;

  target += normal_lpdf(b | 0, 10);
  target += student_t_lpdf(b0 | 3, 0, 10);
  
  // print("");
  // print("        b0 : ",  b0);
  // for(i in 1:N) { print("mu[i]      : ", mu[i]); }
  // for(i in 1:N) { print("mu[i] + b0 : ", mu[i] + b0); }
  
  // likelihood including all constants
  if (!prior_only) {
    target += binomial_logit_lpmf(y | n, mu + b0);
  }
}



