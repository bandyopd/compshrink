functions {
// for likelihood estimation
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    real alpha_plus = sum(alpha);
    return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
                - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }
}


data {
  int<lower=1> N; // total number of observations
  int<lower=2> ncolY; // number of categories
  int<lower=2> ncolX; // number of predictor levels
  matrix[N,ncolX] X; // predictor design matrix
  int <lower=0> Y[N,ncolY]; // data // response variable
}
parameters {
  matrix[ncolX, ncolY] beta_raw; // coefficients (raw)
  matrix<lower=0>[ncolX,ncolY] lambda_tilde; // truncated local shrinkage
  vector<lower=0>[ncolY] tau; // global shrinkage
  simplex[ncolY] pi[N]; // category probabilties
}
transformed parameters{
  matrix[ncolX,ncolY] beta; // coefficients
  matrix<lower=0>[ncolX,ncolY] lambda; // local shrinkage
  matrix[N, (ncolY-1)] Z; // stick-breaking construction 
  
  lambda = diag_post_multiply(lambda_tilde, tau);
  beta = beta_raw .* lambda;
  
  for(i in 1:N){
  real sum_pi;
  sum_pi = 0;
      Z[i,1] = pi[i,1];
      for(j in 2:(ncolY-1)){
        sum_pi += pi[i, j-1];
        if(sum_pi > 0.99)
          Z[i,j] = 0;
        else 
          Z[i,j] = pi[i,j]/(1-sum_pi);
      }
  }
}
model {
// prior:
    for (k in 1:ncolX) {
      for (l in 1:ncolY) {
        tau[k] ~ cauchy(0, 1); // flexible 
        lambda_tilde[k,l] ~ cauchy(0, 1);
        beta_raw[k,l] ~ normal(0,1);
    }
  }
// likelihood
for (i in 1:N) {
    vector[ncolY] alphas;
    for (j in 1:(ncolY-1)){
        alphas[j] = exp(X[i,] * beta[,j]);
           // A general beta = general Dirichlet 
        Z[i,j] ~ beta(alphas[j]/ncolY, alphas[j]*(1-j/(ncolY*1.0))); 
        
      }
     Y[i,] ~ multinomial(pi[i,]);
    }
}

