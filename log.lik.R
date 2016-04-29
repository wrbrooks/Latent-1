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
  #Change the lengths of the vectors
  #Beta is 2 per column
  #Gamma is 2 per row
  #sigma is 2
  #1 lambda(lagrangian multiplier)
  alpha = params[1:(p*d)]
  beta = params[(p*d+1):(2*p+p*d)]
  beta1 = beta[1:p]
  beta2 = beta[(p+1):(2*p)]
  gamma = params[(2*p+p*d+1):(2*p+p*d+2*n)]
  gamma1 = gamma[1:n]
  gamma2 = gamma[(n+1):(2*n)]
  sigma1 = params[(2*p+p*d+2*n+1)]
  sigma2 = params[(2*p+p*d+2*n+2)]
  lambda = params[(2*p+p*d+2*n+3):(2*p+p*d+2*n+2+d)]
  
  
  #Event-level alphas:
  #Victor, do loglikelihood of bivariate normal distribution here
  #Add lagrangians here as well
  alpha.local = matrix(0, n, p)
  sumLogLikNormal = 0
  for (k in 1:d) {
    indx = which(event==unique(event)[k])
    alpha.local[indx,] = matrix(alpha[(p*(k-1)+1):(p*k)], length(indx), p, byrow=TRUE)
    
    sumLogLikNormal = sumLogLikNormal - lambda[k]*sum(gamma1[indx]*gamma2[indx])
  }
  
  sumLogLikNormal = sumLogLikNormal -n/2*log(sigma1*sigma2) - 1/(2*sigma1)*(sum(gamma1^2))-1/(2*sigma2)*(sum(gamma2^2))
  
  #Anonymous function to compute log likelihood:
  (function (x) sum(data*x - exp(x), na.rm=TRUE))(gamma1 %*% t(beta1) + gamma2 %*% t(beta2) + alpha.local)+sumLogLikNormal
}