
// Two parameter dose response model.
// The parameterisation here allows you to identify the ED50 and the slope. 
// The response is the proportion recovered at day 7, which is a continuous
// value on [0, 1].
// In the two parameter DRM, the lower value is fixed at zero and the upper
// is fixed at 1. The next iteration of this model will allow both the upper
// and lower bounds to be estimated from the data.
// The ED50 is the dose for which we see a response equal to 50%.


functions {
  // generates dose response data governed by two parameters
  // inputs:
  // dose is the untransformed dose
  // ed50 corresponds to the dose where the response is 50%
  // slope controls the steepness of the transition between 0 and 1 or 1 to 0.
  // positive values for the slope align with data increasing from 0 to 1
  // negative values for the slope align with data decreasing from 1 to 0
  // high values for the slope makes the transition very steep
  // low values for the slope makes the transition very slow
  // real sigmoid(real dose, real ed50, real slope){
  //   return inv_logit((dose-ed50)*slope);
  // }
  
}

data {
  // sample size
  int<lower=0> N;
  real<lower=0> dose[N];
  int<lower=0> y[N];
  int<lower=0> trials[N];
}

parameters {
  
  real slope;
  real<lower=0> ed50;
  real<lower=0,upper=1> lwr;
  real<lower=0,upper=1> upr;
}

transformed parameters {

  real<lower=0,upper=1> yhat[N];

  for (j in 1:N){
    
    if(dose[j] == 0){
      yhat[j] = lwr ;
    } else {
      
      real numer = upr - lwr;
      real denom = 1 + exp(slope * (log(dose[j]) - log(ed50)));
      yhat[j] = lwr + (numer / denom);
    }

  }
}

model {
  
  // note that the parameterisation for the Beta distribution is in terms of 
  // the mean mu and the sample size v. Wikipedia says that this can be done 
  // by setting alpha=mu*v and beta=(1-mu)*v. 

  target += normal_lpdf(slope | 0, 10);
  
  // ed50 has lower bound of 0 therefore the next is half-normal 
  target += normal_lpdf(ed50 | 0, 10);

  target += uniform_lpdf(lwr | 0, 1);
  target += uniform_lpdf(upr | 0, 1);

  for (j in 1:N){

    // target += beta_lpdf(y[j] | yhat[j]*trials[j], (1-yhat[j])*trials[j]);
    target += binomial_lpmf(y[j] | trials[j], yhat[j]);
  }
}

