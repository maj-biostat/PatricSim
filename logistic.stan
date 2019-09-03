// generated with brms 2.9.0
functions {
}
data {
  int<lower=1> N;  // number of observations
  int y[N];  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  // int Kc = K - 1;
  // matrix[N, Kc] Xc;  // centered version of X
  // vector[Kc] means_X;  // column means of X before centering
  // for (i in 2:K) {
  //   means_X[i - 1] = mean(X[, i]);
  //   Xc[, i - 1] = X[, i] - means_X[i - 1];
  // }
}
parameters {
  vector[K] b;  // population-level effects
  // real temp_Intercept;  // temporary intercept
}
transformed parameters {
}
model {
  vector[N] mu = X * b;
  // priors including all constants
  // target += student_t_lpdf(temp_Intercept | 3, 0, 10);
  // likelihood including all constants
  if (!prior_only) {
    target += bernoulli_logit_lpmf(y | mu);
  }
}
// generated quantities {
//   // actual population-level intercept
//   real b_Intercept = temp_Intercept - dot_product(means_X, b);
// }


