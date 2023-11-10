data{
  int<lower=1> n; // number of observations
  int<lower=1> K; // number of categories 
  int <lower=1,upper=n> y[K]; // data
  vector<lower=0>[K] alpha;
}
parameters{
  simplex[K] pi; // category probabilties
}
model{
  pi ~ dirichlet(alpha);
  y ~ multinomial(pi);
}
