#' simulate data from a Dirichlet-Multinomial regression model
#'
#' @note Requires the dirmult and MASS packages
#'
#' @param n_obs: the number of samples
#' @param n_vars: number of covariates excluding the intercept
#' @param n_taxa: number of species
#' @param n_relevant_vars: number of relevant nutrients
#' @param n_relevant_taxa: number of relevant species
#' @param beta_min: minimum absolute value of the regression parameters
#' @param beta_max: maximum absolute value of the regression parameters
#' @param signoise: scalar multiplier on the regression parameters
#' @param n_reads_min: lower bound on uniform distribution for number of reads in each sample
#' @param n_reads_max: upper bound on uniform distribution for number of reads in each sample
#' @param theta0: the dispersion parameter
#' @param rho: the correlation between covariates
#'
#' @return XX: (design matrix) with intercept: n_obs * (n_vars + 1)
#' @return YY: (count matrix) rows: n_obs samples, columns: n_taxa species
#' @return alphas: simulated intercept vector
#' @return betas: simulated coefficient matrix n_taxa * (n_vars + 1)
#' @return n_reads_min, n_read_max: row sum parameters
#' @return theta0, phi, rho, signoise: simulation inputs
#'
#' @export
simulate_dirichlet_multinomial_regression = function(n_obs = 100,
                                                     n_vars = 100,
                                                     n_taxa = 40,
                                                     n_relevant_vars = 4,
                                                     n_relevant_taxa = 4,
                                                     beta_min = 0.5,
                                                     beta_max = 1.0,
                                                     signoise = 1.0,
                                                     n_reads_min = 1000,
                                                     n_reads_max = 2000,
                                                     theta0 = 0.01,
                                                     rho = 0.4){
  
  # check for required packages
  if(!require(dirmult)){
    stop("dirmult package required")
  }
  if(!require(MASS)){
    stop("MASS package required")
  }
  # covariance matrix for predictors
  Sigma = matrix(1, n_vars, n_vars)
  Sigma = rho^abs(row(Sigma) - col(Sigma))
  # include the intercept
  XX = cbind(rep(1, n_obs),
             scale(MASS::mvrnorm(n = n_obs, mu = rep(0, n_vars), Sigma = Sigma)))
  # empties
  YY = matrix(0, n_obs, n_taxa)
  betas = matrix(0, n_taxa, n_vars)
  phi = matrix(0, n_obs, n_taxa)
  # parameters with signs alternating
  st = 0
  low_side = beta_min
  high_side = beta_max
  if(n_relevant_taxa != 1){
    # warning if the lengths don't match
    coef = suppressWarnings(seq(low_side, high_side, len = n_relevant_taxa) * c(1, -1))
  }else{
    coef = (low_side + high_side) / 2
  }
  coef_g = rep(1.0, len = n_relevant_vars)
  for(ii in 1:n_relevant_vars){
    # overlap species
    betas[(st:(st + n_relevant_taxa - 1)) %% n_taxa + 1, 3 * ii - 2] = coef_g[ii] * sample(coef)[((ii - 1):(ii + n_relevant_taxa - 2)) %% n_relevant_taxa + 1]
    st = st + 1
  }
  # -2.3 and 2.3 so that the intercept varies over three orders of magnitude
  intercept = runif(n_taxa, -2.3, 2.3)
  Beta = cbind(intercept, signoise * betas)
  # row totals
  ct0 = sample(n_reads_min:n_reads_max, n_obs, rep = T)
  for(ii in 1:n_obs){
    thisrow = as.vector(exp(Beta %*% XX[ii, ]))
    phi[ii, ] = thisrow/sum(thisrow)
    YY[ii, ] = dirmult::simPop(J = 1, n = ct0[ii], pi = phi[ii, ], theta = theta0)$data[1, ]
  }
  
  return(list(XX = XX, YY = YY, alphas = intercept, betas = Beta,
              n_reads_min = n_reads_min, n_reads_max = n_reads_max,
              theta0 = theta0, phi = phi, rho = rho, signoise = signoise))
  
}
