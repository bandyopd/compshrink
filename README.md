# compshrink

This is the repository for Bayesian Dirichlet-Multinomial regression using global-local shrinkage priors: including horseshoe, horseshoe+ and Bayesian Lasso as described in the paper entitled "Shrinkage and Selection for Compositional Data" by Jyotishka Datta and Dipankar Bandopadhyay. 

## Abstract

We propose a variable selection and estimation framework for Bayesian compositional regression model using state-of-the-art continuous shrinkage priors to identify the significant associations between available covariates and taxonomic abundance from microbiome data. We use a generalized Dirichlet and Dirichlet distribution for modeling the compositional component and compare the popular horseshoe (Carvalho et al., 2010) and horseshoe+ (Bhadra et al., 2017) priors along with the Bayesian Lasso as a benchmark. We use Hamiltonial Monte Carlo for posterior sampling and posterior credible intervals and pseudo posterior inclusion probabilities for variable selection. Our simulation studies show excellent recovery and estimation accuracy for sparse parameter regime, and we apply our method to human microbiome data from NYC-Hanes study.


