data {
  int<lower=1> n; // number of observations
  int<lower = 1> K;
  int <lower=0,upper=n> y[K]; // data
  vector<lower = 0>[K - 1] a;
  vector<lower = 0>[K - 1] b;
}
parameters {
  vector<lower = 0, upper = 1>[K - 1] z;
}
transformed parameters {
  simplex[K] pi; // do the stick breaking thing
  pi[1] = z[1];
  for (k in 2:(K - 1)) pi[k] = z[k] / (1 - sum(pi[1:(k - 1)]));
  pi[K] = 1 - sum(pi[1:(K - 1)]);
}
model {
  y ~ multinomial(pi);
  target += beta_lpdf(z | a, b);
}
