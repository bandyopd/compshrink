//
// This Stan program defines a Generalized Dirichlet Multinomial model, with a
// matrix of values 'Y' modeled as GDM distributed
// with horseshoe prior on beta coefficients 
// ** Old version - directly modeling GDM parameters is bad **
// ** Fixed GDM definition: alpha, beta are K-1 dimensional. Fixed it**
// ** TL;DR: do NOT use this, bad convergence/estimates **
//

// Dirichlet likelihood 
functions {
// for likelihood estimation
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    real alpha_plus = sum(alpha);
    return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
                - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }

// Generalized Dirichlet Likelihood - for two shape paramaters alpha and beta
real generalized_dirichlet_multinomial_lpmf(int[] y, vector alpha, vector beta_) {
    // y is num_categories dimensional, alpha and beta are num_categories-1 dimensional
    int D = dims(alpha)[1]; // D = num_categories-1
    vector[D+1] x = cumulative_sum(to_vector(y));
    vector[D] z = rep_vector(x[D+1],D) - x[1:D];
    vector[D] alpha_prime = alpha + to_vector(y)[1:D];
    vector[D] beta_prime;

    beta_prime = beta_ + z;

    return (lgamma(x[D+1]+1) - sum(lgamma(to_vector(y)+rep_vector(1,D+1)))
      + sum(lgamma(alpha_prime)) - sum(lgamma(alpha)) 
      + sum(lgamma(beta_prime)) - sum(lgamma(beta_))
      - sum(lgamma(alpha_prime + beta_prime)) + sum(lgamma(alpha+beta_))
      );

 }
}

data {
  int<lower=1> N; // total number of observations
  int<lower=2> ncolY; // number of categories
  int<lower=2> ncolX; // number of predictor levels
  matrix[N,ncolX] X; // predictor design matrix
  int <lower=0> Y[N,ncolY]; // data // response variable
  real<lower=0> scale_icept;
  // real<lower=0> sd_prior;
  real<lower=0> psi;
}
parameters {
  vector[N] beta0; // intercept
  matrix[ncolX, ncolY] beta_raw; // coefficients (raw)
  matrix<lower=0>[ncolX,ncolY] lambda_tilde; // truncated local shrinkage
  vector<lower=0>[ncolY] tau; // global shrinkage
  
  vector[N] alpha0; // intercept
  matrix[ncolX, ncolY] alpha_raw; // coefficients (raw)
  matrix<lower=0>[ncolX,ncolY] gamma_tilde; // truncated local shrinkage
  vector<lower=0>[ncolY] phi; // global shrinkage
  
  real<lower=0> sigma;
}
transformed parameters{
  matrix[ncolX,ncolY] beta; // coefficients
  matrix<lower=0>[ncolX,ncolY] lambda; // local shrinkage 1
  matrix[ncolX,ncolY] alpha; // coefficients
  matrix<lower=0>[ncolX,ncolY] gamma; // local shrinkage 2
  
  lambda = diag_post_multiply(lambda_tilde, tau);
  beta = beta_raw .* lambda*sigma;
  
  gamma = diag_post_multiply(gamma_tilde, phi);
  alpha = alpha_raw .* gamma*sigma;
}

model {
   // psi ~ uniform(0,1);
 sigma ~ inv_gamma(0.5, 0.5);
// prior:
for(k in 1:N){
  beta0[k] ~ cauchy(0, scale_icept);
  alpha0[k] ~ cauchy(0, scale_icept);
}
// prior:
    for (k in 1:ncolX) {
      for (l in 1:ncolY) {
        tau[l] ~ cauchy(0.1, 5); // flexible 
        lambda_tilde[k,l] ~ cauchy(0, 5);
        beta_raw[k,l] ~ normal(0,1);
        
        phi[l] ~ cauchy(0.1, 5); // flexible 
        gamma_tilde[k,l] ~ cauchy(0, 5);
        alpha_raw[k,l] ~ normal(0,1);
    }
  }
// likelihood
for (i in 1:N) {
    vector[ncolY] logits_1;
    vector[ncolY] logits_2;
    
    for (j in 1:ncolY){
      logits_1[j] = beta0[i]+X[i,] * beta[,j];
      logits_2[j] = alpha0[i]+X[i,] * alpha[,j];
      }
     Y[i,] ~ generalized_dirichlet_multinomial(softmax(logits_1)[1:(ncolY-1)],softmax(logits_2)[1:(ncolY-1)]);
    }
}

