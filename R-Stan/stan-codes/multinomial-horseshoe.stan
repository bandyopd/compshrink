data {
  int<lower=1> N; // total number of observations
  int<lower=2> ncolY; // number of categories
  int<lower=2> ncolX; // number of predictor levels
  matrix[N,ncolX] X; // predictor design matrix
  int <lower=0> Y[N,ncolY]; // data // response variable
  real sd_prior; // Prior standard deviation
}
parameters {
  matrix[ncolX, ncolY] beta_raw; // coefficients (raw)
  matrix[ncolX,ncolY] lambda; // local shrinkage
  real<lower=0> tau[ncolY]; // global shrinkage
  simplex[ncolY] pi[N];

}
transformed parameters{

matrix[ncolX,ncolY] beta; // coefficients

for (k in 1:ncolX) {
  for (l in 1:ncolY) {
    beta[k,l] = beta_raw[k,l]*lambda[k,l];
    }
  }
}
model {
// prior:
    for (k in 1:ncolX) {
      for (l in 1:ncolY) {
        tau[k] ~ cauchy(0, 1); // flexible 
        lambda[k,l] ~ cauchy(0, tau[k]);
        beta_raw[k,l] ~ normal(0,sd_prior);
    }
  }
// likelihood
for (i in 1:N) {
    vector[ncolY] logits;
    for (j in 1:ncolY){
      logits[j] = X[i,] * beta[,j];
      }
     pi[i,] ~ dirichlet(softmax(logits)); 
     Y[i,] ~ multinomial(pi[i,]);
    }
}

