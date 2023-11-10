data {
  int<lower=1> N; // total number of observations
  int<lower=2> ncolY; // number of categories
  int<lower=2> ncolX; // number of predictor levels
  matrix[N,ncolX] X; // predictor design matrix
  matrix[N,ncolY] Y; // response variable
  real sd_prior; // Prior standard deviation
}
parameters {
  matrix[ncolY-1,ncolX] beta_raw; // coefficients (raw)
  real theta;
  matrix[(ncolY-1),ncolX] lambda; // local shrinkage
  vector[(ncolY-1)] tau; // global shrinkage
}
transformed parameters{
  real exptheta = exp(theta);
  
  matrix[ncolY,ncolX] beta; // coefficients

  for (l in 1:ncolX) {
      beta[ncolY,l] = 0.0;
  }
for (k in 1:(ncolY-1)) {
  for (l in 1:ncolX) {
    beta[k,l] = beta_raw[k,l]*lambda[k,l];
    }
  }
}
model {
// prior:
  theta ~ normal(0,sd_prior);
    for (k in 1:(ncolY-1)) {
      for (l in 1:ncolX) {
        tau[k] ~ cauchy(0, 1); // flexible 
        lambda[k,l] ~ cauchy(0, 1);
        beta_raw[k,l] ~ normal(0,sd_prior);
    }
  }
// likelihood
for (n in 1:N) {
    vector[ncolY] logits;
    for (m in 1:ncolY){
      logits[m] = X[n,] * transpose(beta[m,]);
      }
    transpose(Y[n,]) ~ dirichlet(softmax(logits)); 
    }
}
