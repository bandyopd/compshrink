//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions {
// for likelihood estimation
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    real alpha_plus = sum(alpha);
    return lgamma(alpha_plus) + lgamma(sum(y)+1) + sum(lgamma(alpha + to_vector(y)))
                - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha))-sum(lgamma(to_vector(y)+1));
  }
}
// The input data is a matrix 'Y' of size 'N' by 'ncolY'.
data {
  int<lower=1> N; // total number of observations
  int<lower=2> ncolY; // number of categories
  int<lower=2> ncolX; // number of predictor levels
  matrix[N,ncolX] X; // predictor design matrix
  int <lower=0> Y[N,ncolY]; // data // response variable
  real<lower=0> sd_prior;
  real<lower=0> psi;
  real<lower=0> d;
}

// The parameters accepted by the model. Our model
// accepts the following parameters 
parameters {
  matrix[ncolX, ncolY] beta_raw; // coefficients (raw)
  vector[N] beta0; // intercept
  matrix<lower=0>[ncolX,ncolY] lambda_tilde; // truncated local shrinkage
  vector<lower=0>[ncolY] tau; // global shrinkage
}
transformed parameters{
  matrix[ncolX,ncolY] beta; // coefficients
  matrix<lower=0>[ncolX,ncolY] lambda; // local shrinkage

  lambda = diag_post_multiply(lambda_tilde, tau);
  beta = beta_raw .* lambda;
}

// The model to be estimated. We model the output
// 'Y' to have a Dirichlet-Multinomial Distribution 
model {
// prior:
for(k in 1:N){
  beta0 ~ normal(0, 10);
}
    for (k in 1:ncolX) {
      for (l in 1:ncolY) {
        tau[l] ~ inv_gamma(psi/2,psi*d^2/2);; // flexible 
        lambda_tilde[k,l] ~ exponential(1/(2*tau[l]^2));
        beta_raw[k,l] ~ normal(0,sd_prior);
    }
  }
// likelihood
for (i in 1:N) {
    vector[ncolY] logits;
    for (j in 1:ncolY){
      logits[j] = beta0[i] + X[i,] * beta[,j];
      }
     Y[i,] ~ dirichlet_multinomial(softmax(logits)*(1-psi)/psi);
    }
}

