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
}
model {
  y ~ bernoulli_logit(mu + beta * g + gamma * ind);
}
