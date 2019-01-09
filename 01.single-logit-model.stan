data {
  int<lower=0> N;               // number of individuals
  vector[N] g;    // genotype
  int<lower=0,upper=1> y[N];    // case/control status
  vector[N] ind;  // indicator for heterozygotes
}
parameters {
  real mu;        // baseline odds of being an athlete
  real beta;      // additive effect
  real gamma;     // dominance effect
  // hyperparameters:
  real beta0;
  real gamma0;
  real<lower=0> taubeta;
  real<lower=0> taugamma;
  real<lower=-1,upper=1> rho;
}
model {
  y ~ bernoulli_logit(mu + beta * g + gamma * ind);
  {
    vector[2] coeffs = [ beta, gamma ]';
    vector[2] coeffs0 = [ beta0, gamma0 ]';
    matrix[2,2] sigma = [ [ taubeta*taubeta, rho*taubeta*taugamma ], [ rho*taubeta*taugamma,  taugamma*taugamma ] ];
    coeffs ~ multi_normal(coeffs0, sigma);
  }
}
