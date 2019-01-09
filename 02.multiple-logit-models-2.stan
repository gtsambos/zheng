data {
  int<lower=0> N1;               // number of individuals in study 1
  int<lower=0> N2;               // number of individuals in study 2
  int<lower=0> N3;               // number of individuals in study 3
  vector[N1] g1;    // genotypes in study 1
  vector[N2] g2;    // genotypes in study 2
  vector[N3] g3;    // genotypes in study 3
  int<lower=0,upper=1> y1[N1];    // case/control status in study 1
  int<lower=0,upper=1> y2[N2];    // case/control status in study 2
  int<lower=0,upper=1> y3[N3];    // case/control status in study 3
  vector[N1] ind1;  // indicator for heterozygotes in study 1
  vector[N2] ind2;  // indicator for heterozygotes in study 2
  vector[N3] ind3;  // indicator for heterozygotes in study 3
}
parameters {
  real mu1;        // baseline odds of being an athlete in study 1
  real mu2;        // baseline odds of being an athlete in study 2
  real mu3;        // baseline odds of being an athlete in study 3
  real beta1;      // additive effect in study 1
  real beta2;      // additive effect in study 2
  real beta3;      // additive effect in study 3
  real gamma1;     // dominance effect in study 1
  real gamma2;     // dominance effect in study 2
  real gamma3;     // dominance effect in study 3
  // hyperparameters:
  real beta0;
  real gamma0;
  real<lower=0> taubeta;
  real<lower=0> taugamma;
  real<lower=-1,upper=1> rho;
}
model {
  y1 ~ bernoulli_logit(mu1 + beta1 * g1 + gamma1 * ind1); # likelihood model for study 1
  y2 ~ bernoulli_logit(mu2 + beta2 * g2 + gamma2 * ind2); # likelihood model for study 2
  y3 ~ bernoulli_logit(mu3 + beta3 * g3 + gamma3 * ind3); # likelihood model for study 3
  {
    vector[2] coeffs1 = [ beta1, gamma1 ]';
    vector[2] coeffs2 = [ beta2, gamma2 ]';
    vector[2] coeffs3 = [ beta3, gamma3 ]';
    vector[2] coeffs0 = [ beta0, gamma0 ]';
    matrix[2,2] sigma = [ [ taubeta*taubeta, rho*taubeta*taugamma ], [ rho*taubeta*taugamma,  taugamma*taugamma ] ];
    coeffs1 ~ multi_normal(coeffs0, sigma);
    coeffs2 ~ multi_normal(coeffs0, sigma);
    coeffs3 ~ multi_normal(coeffs0, sigma);
  }
}
