# compshrink

This is the repository for Bayesian Dirichlet-Multinomial regression using global-local shrinkage priors: including horseshoe, horseshoe+ and Bayesian Lasso as described in the paper entitled "Shrinkage and Selection for Compositional Data" by Jyotishka Datta and Dipankar Bandopadhyay. 

## Abstract

We propose a variable selection and estimation framework for Bayesian compositional regression model using state-of-the-art continuous shrinkage priors to identify the significant associations between available covariates and taxonomic abundance from microbiome data. We use a generalized Dirichlet and Dirichlet distribution for modeling the compositional component and compare the popular horseshoe (Carvalho et al., 2010) and horseshoe+ (Bhadra et al., 2017) priors along with the Bayesian Lasso as a benchmark. We use Hamiltonial Monte Carlo for posterior sampling and posterior credible intervals and pseudo posterior inclusion probabilities for variable selection. Our simulation studies show excellent recovery and estimation accuracy for sparse parameter regime, and we apply our method to human microbiome data from NYC-Hanes study.


```{r}
stan.hs.fit <- stan_model(file = 'multinomial-horseshoe-marg.stan', 
                          model_name = "Dirichlet Horseshoe")
stan.hsplus.fit <- stan_model(file = 'multinomial-hsplus-marg.stan', 
                              model_name = "Dirichlet HS+")
stan.laplace.fit <- stan_model(file = 'multinomial-laplace-marg.stan', 
                               model_name = "Dirichlet Laplace")
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

