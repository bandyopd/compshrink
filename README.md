# compshrink

This is the repository for implementing a Bayesian Dirichlet-Multinomial regression using global-local shrinkage priors, which includes the horseshoe, horseshoe+ and Bayesian Lasso, as described in the paper: Datta J and Bandyopadhyay D. (2024). Bayesian variable shrinkage and selection in compositional data regression: Application to oral microbiome, *Journal of the Indian Society for Probability and Statistics*, 25, 491-515

## Abstract

We propose a variable selection and estimation framework for Bayesian compositional regression model using state-of-the-art continuous shrinkage priors to identify the significant associations between available covariates and taxonomic abundance from microbiome data. We use a generalized Dirichlet and Dirichlet distribution for modeling the compositional component and compare the popular horseshoe (Carvalho et al., 2010) and horseshoe+ (Bhadra et al., 2017) priors along with the Bayesian Lasso as a benchmark. We use Hamiltonial Monte Carlo for posterior sampling and posterior credible intervals and pseudo posterior inclusion probabilities for variable selection. Our simulation studies show excellent recovery and estimation accuracy for sparse parameter regime, and we apply our method to human microbiome data from NYC-Hanes study.

## References

1. Carvalho CM, Polson NG, Scott JG (2010). The horseshoe estimator for sparse signals, *Biometrika*, 97(2), 465–480
2. Bhadra A, Datta J, Polson NG, Willard B. (2017). The Horseshoe+ Estimator of Ultra-Sparse Signals, *Bayesian Analysis*, 12(4), 1105-1131

## Details 

The `Stan` codes are provided in the `stan-codes` folder. 

```{stan}
//
// This Stan program defines a Dirichlet Multinomial model, with a
// matrix of values 'Y' modeled as GDM distributed
// with horseshoe prior on beta coefficients 
// ** Wadsworth's model - with HS instead of Spike-Slab**
//

functions {
// for likelihood estimation
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    real alpha_plus = sum(alpha);
    return lgamma(alpha_plus) + lgamma(sum(y)+1) + sum(lgamma(alpha + to_vector(y)))
                - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha))-sum(lgamma(to_vector(y)+1));
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
  matrix[ncolX, ncolY] beta_raw; // coefficients (raw)
  vector[N] beta0; // intercept
  matrix<lower=0>[ncolX,ncolY] lambda_tilde; // truncated local shrinkage
  vector<lower=0>[ncolY] tau; // global shrinkage
  real<lower=0> sigma;
  // real<lower = 0> psi;
}
transformed parameters{
  matrix[ncolX,ncolY] beta; // coefficients 
  matrix<lower=0>[ncolX,ncolY] lambda; // local shrinkage
  lambda = diag_post_multiply(lambda_tilde, tau);
  beta = beta_raw .* lambda*sigma;
}

model {
 // psi ~ uniform(0,1);
 sigma ~ inv_gamma(0.5, 0.5);
// prior:
for(k in 1:N){
  beta0[k] ~ cauchy(0, scale_icept);
}
for (k in 1:ncolX) {
      for (l in 1:ncolY) {
        tau[l] ~ cauchy(0, 1); // flexible 
        lambda_tilde[k,l] ~ cauchy(0, 1);
        beta_raw[k,l] ~ normal(0,1);
    }
  }
// likelihood
for (i in 1:N) {
    vector[ncolY] logits;
    for (j in 1:ncolY){
      logits[j] = beta0[i]+X[i,] * beta[,j];
      }
     Y[i,] ~ dirichlet_multinomial(softmax(logits)*(1-psi)/psi);
    }
}
```

```{r}
stan.hs.fit <- stan_model(file = 'multinomial-horseshoe-marg.stan', 
                          model_name = "Dirichlet Horseshoe")
```

For any of the three candidate priors, we can sample from the posterior using the `sampling` function in R-Stan. 

```{r}
n.iters = 1000
n.chains = 1
rng.seed = 12345

set.seed(rng.seed)
dirfit <- dirmult(Ymat)

NYC.data = list(N = nrow(Ymat), ncolY = ncol(Ymat), ncolX = ncol(Xmat),
                 X = Xmat, Y = Ymat, psi = dirfit$theta, scale_icept = 2, d=1) 

ptm = proc.time()
smpls.hs.res = sampling(stan.hs.fit, 
                        data = NYC.data, 
                        iter = n.iters,
                        init = 0,
                        seed = rng.seed,
                        cores = 2,
                        warmup = floor(n.iters/2),
                        chains = n.chains,
                        control = list(adapt_delta = 0.85),
                        refresh = 100)
proc.time()-ptm
# summarize results

beta.smpls.hs <- rstan::extract(smpls.hs.res, pars=c("beta"), permuted=TRUE)[[1]]
beta.mean.hs <- apply(beta.smpls.hs, c(2,3), mean)
beta.median.hs <- apply(beta.smpls.hs, c(2,3), median)
beta.mode.hs <- apply(beta.smpls.hs, c(2,3), Mode)
beta.sd.hs <- apply(beta.smpls.hs, c(2,3),sd)
beta.hs.ci <- apply(beta.smpls.hs, c(2,3), quantile, probs=c(0.025,0.5,0.975)) #the median line with 95% credible intervals
```

