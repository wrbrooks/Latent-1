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
  # basic constants:
  n <- nrow(data)
  p <- ncol(data)
  d <- length(unique(event))
  
  # check for predictable error conditions
  if (is.null(specific))
    specific <- rep(FALSE, p)
  if (length(specific) != p)
    stop("Argument 'specific' must have length equal to the number of 
         pathogens.")
  if (!all(typeof(specific) == "logical"))
    stop("Each element of 'specific' must be logical (boolean).")
  
  # Extract parameters from the parameter vector.
  # sigma is length 2
  # lambda is the lagrange multiplier
  alpha <- params[1:(p * d)]
  beta <- params[(p * d + 1):(p * (2 + d))]
  gamma <- params[(p * (2 + d) + 1):(p * (2 + d) + 2 * n)]
  sigma1 <- params[p * (2 + d) + 2 * n + 1]
  sigma2 <- params[p * (2 + d) + 2 * n + 2]
  lambda <- tail(params, d)

  # read the separate elements of beta and gamma
  # beta is 2 per column
  # gamma is 2 per row
  beta1 <- beta[1:p]
  beta2 <- tail(beta, p)
  gamma1 <- gamma[1:n]
  gamma2 <- tail(gamma, n)

  
  # event-level alphas
  # sort out lambdas here too
  alpha.local <- matrix(0, n, p)
  lambda.local <- vector()
  for (k in 1:d) {
    indx <- which(event==unique(event)[k])
    alpha.local[indx,] <- matrix(alpha[(p * (k - 1) + 1):(p * k)], length(indx), p, byrow=TRUE)
    
    #Derivative w.r.t. lambda(lagrange multiplier) is simply the dot product of the gammas in the event
    lambda.local <- c(lambda.local, -sum(gamma1[indx]*gamma2[indx]))
  }
  
  # compute this once to save time:
  mu <- exp(gamma1 %*% t(beta1) + gamma2 %*% t(beta2) + alpha.local)
  
  # gradient in the direction of alpha
  grad <- rep(0, p*d)
  
  # gradient in the direction of beta:
  #grad <- c(grad, colSums(sweep(data - eta, 1, gamma1, '*'), na.rm=TRUE))
  #grad <- c(grad, colSums(sweep(data - eta, 1, gamma2, '*'), na.rm=TRUE)[1:3])
  #grad <- c(grad, rep(0, 2))
  grad <- c(grad, rep(0, p))
  grad <- c(grad, rep(0, p))
  
  # gradient in the direction of gamma:
  # treat these as zero, optimize them for IRLS
  #grad <- c(grad, rep(0, n))
  #grad <- c(grad, rep(0, n))
  
  # gradient in the direction of gamma1:
  grad <- c(grad, rowSums(sweep(data - mu, 2, beta1, '*')) - gamma2 * unlist(sapply(1:length(lambda), function(k) rep(lambda[k], sum(event==unique(event)[k])))))
  
  # gradient in the direction of gamma2:
  grad <- c(grad, rowSums(sweep(data - mu, 2, beta2, '*')) - gamma1 * unlist(sapply(1:length(lambda), function(k) rep(lambda[k], sum(event==unique(event)[k])))))
  
  # Victor, gradient in the direction of sigmas:
  grad <- c(grad, -n/2/sigma1 + sum(gamma1^2)/2/sigma1^2)
  grad <- c(grad, -n/2/sigma2 + sum(gamma2^2)/2/sigma2^2)
  
  # lambdas as well
  grad <- c(grad, lambda.local)
  
  grad
}
