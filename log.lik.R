#' Log-likelihood function
#' 
#' Log-likelihood of the data under a Poisson distribution with log-linear model for the mean
#' 
#' @param data matrix of counts, with rows equal to the number of observations and columns equal to the number of categories
#' @param params vector of the alpha, beta, and gamma parameters
#' @param event vector of labels with one entry per row of data, where each entry identifies the event to which the corresponding row of data belongs
#' 
#' @return total log-likelihood of the data under the current parameter values
#' 
log.lik <- function(data, params, event) {
  #Basic constants:
  n = nrow(data)
  p = ncol(data)
  d = length(unique(event))
  
  #Extract parameters from the parameter vector:
  alpha = params[1:(p*d)]
  beta = params[(p*d+1):(p*d+p)]
  gamma = params[(p*d+1+p):(p*d+p+n)]
  sigma = params[(p*d+1+p+n)]
  
  #Event-level alphas:
  alpha.local = matrix(0, n, p)
  
  for (k in 1:d) {
    indx = which(event==unique(event)[k])
    alpha.local[indx,] = matrix(alpha[p*(k-1) + 1:p], length(indx), p, byrow=TRUE)
    
  }
  sumLogLik = -n/2*log(sigma)-1/(2*sigma)*(sum(gamma^2))
  #Anonymous function to compute log likelihood:
  (function (x) sum(data*x - exp(x), na.rm=TRUE))(gamma %*% t(beta) + alpha.local) -n/2*log(sigma)-1/(2*sigma)*(sum(gamma^2))
}