data {
  int<lower=0> N[2];              // sample sizes in each study
  vector[N[1]] G1;                // genotypes in study 1
  vector[N[2]] G2;                // genotypes in study 2
  int<lower=0,upper=1> Y1[N[1]];  // case-control status in study 1
  int<lower=0,upper=1> Y2[N[2]];  // case-control status in study 2
  vector[N[1]] I1;                // indicators in study 1
  vector[N[2]] I2;                // indicators in study 2
}
parameters{
  real mu[2];       
  real beta[2];
  real gamma[2];
}
model{
  Y1 ~ bernoulli_logit(mu[1] + beta[1] * G1 + gamma[1] * I1); // likelihood model for study 1
  Y2 ~ bernoulli_logit(mu[2] + beta[2] * G2 + gamma[2] * I2); // likelihood model for study 2
}
