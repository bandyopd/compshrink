data{
  int<lower=1> n; // number of observations
  int<lower=1> K; // number of categories 
  int <lower=0,upper=n> y[K]; // data
  real<lower=0> alpha; // Dirichlet shape 1
  real<lower=0> beta; // Dirichlet shape 2
}
parameters{
  simplex[K] pi; // category probabilties
}
transformed parameters{
  // stick-breaking construction 
  
  real<lower=0,upper=1> Z[K-1];
  Z[1] = pi[1];
  for(j in 2:(K-1)){
    Z[j] = pi[j]/(1-sum(pi[1:j-1]));
  }
}
model{
  for(j in 1:(K-1)){
    Z[j] ~ beta(alpha/K, alpha*(1-j/K)); // Dirichlet
  }
  y ~ multinomial(pi);
}
