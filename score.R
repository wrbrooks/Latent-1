#' Gradient of the log-likelihood
#' 
#' Evaluates the score functions (gradients of the log-likelihood) for a 
#' Poisson distribution with log-linear model for the mean.
#'
#' @param data matrix of observations, with rows equal to the number of 
#'     observations and columns equal to the number of categories
#' @param params vector of model parameters, consisting of 
#'     \eqn{\boldsymbol{\alpha}, \boldsymbol{\beta}, \boldsymbol{\gamma}}.
#' @param event vector of labels with one entry per row of data, where 
#'     each entry identifies the event to which the corresponding row belongs
#' @param specific vector of logical values, one entry per column of data, 
#'     with each entry indicating whether the corresponding column of the 
#'     data is a human-specific bacterium
#' 
#' @return The score functions, evaluated at the current parameter values.
#' 
score <- function(data, params, event, specific=NULL) {
  #Basic constants:
  n = nrow(data)
  p = ncol(data)
  d = length(unique(event))
  if (is.null(specific))
    specific = rep(FALSE, p)
  if (length(specific) != p)
    stop("Argument 'specific' must have length equal to the number of 
         pathogens.")
  if (!all(typeof(specific) == "logical"))
    stop("Each element of 'specific' must be logical (boolean).")
  
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
  
  #Event-level alphas:'
  #Sort out lambdas here too
  alpha.local = matrix(0, n, p)
  lambda.local = vector()
  for (k in 1:d) {
    indx = which(event==unique(event)[k])
    alpha.local[indx,] = matrix(alpha[p*(k-1) + 1:p], length(indx), 
                                p, byrow=TRUE)
    #Derivative of lambda(lagrangian multiplier) is simply the dot product of the gammas in the event
    lambda.local = c(lambda.local, -sum(gamma1[indx]*gamma2[indx]))
  }
  
  #Just compute this once to save time:
  eta = exp(gamma1 %*% t(beta1) + gamma2 %*% t(beta2) + alpha.local)
  #Gradient in the direction of event-specific alphas:
  grad = vector()
  for (k in 1:d) {
    indx = which(event==unique(event)[k])
    grad = c(grad, as.integer(!specific) * colSums(data[indx,] - 
                                                     eta[indx,], na.rm=TRUE))
  }
  
  #Gradient in the direction of beta:
  grad = c(grad, colSums(sweep(data - eta, 1, gamma1, '*'), na.rm=TRUE))
  grad = c(grad, colSums(sweep(data - eta, 1, gamma2, '*'), na.rm=TRUE)[1:3])
  grad = c(grad, rep(0,2))
  
  #Gradient in the direction of gamma:
  #Treat these as zero, optimize them for IRLS
  grad = c(grad, rep(0,n))
  grad = c(grad, rep(0,n))
  
  #Victor, gradient in the direction of sigmas:
  grad = c(grad, (-n/(2*sigma1) + 1/(2*sigma1^2)*sum(gamma1^2)))
  grad = c(grad, (-n/(2*sigma2) + 1/(2*sigma2^2)*sum(gamma2^2)))
  
  #Lambdas as well
  grad = c(grad, lambda.local)
  
  grad
}