data{
  int<lower=1> n; // number of observations
  int<lower=1> K; // number of categories 
  int <lower=0,upper=n> y[K]; // data
  real<lower=0> alpha; // TPB shape
  real<lower=0> beta; // TPB shape
  real<lower=0> phi; // TPB shape
}
parameters{
  simplex[K] pi; // category probabilties
  real<lower=0> lambda[K-1];
  real<lower=0> tau[K-1];
  // real<lower=0, upper=1> phi;
}
transformed parameters{
  // stick-breaking construction 
  
  real<lower=0,upper=1> Z[K-1];
  Z[1] = pi[1];
  for(j in 2:(K-1)){
    Z[j] = pi[j]/(1-sum(pi[1:j-1]));
  }
for(j in 1:(K-1)){
   Z[j] = 1/(1+tau[j]);
}
}
model{
  for(j in 1:(K-1)){
    // Z[j] ~ beta(alpha/K, alpha*(1-j/K)); // Dirichlet
    lambda[j] ~ gamma(alpha, phi);
    tau[j] ~ gamma(beta, lambda[j]);
  }
  y ~ multinomial(pi);
}
