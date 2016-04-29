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
  n <- nrow(data)
  p <- ncol(data)
  d <- length(unique(event))
  if (is.null(specific))
    specific <- rep(FALSE, p)
  if (length(specific) != p)
    stop("Argument 'specific' must have length equal to the number of 
         pathogens.")
  if (!all(typeof(specific) == "logical"))
    stop("Each element of 'specific' must be logical (boolean).")
  
  #Extract parameters from the parameter vector:
  alpha <- params[1:(p*d)]
  beta <- params[p*d + (1:p)]
  gamma <- params[p*d+p + (1:n)]
  sigma <- tail(params, 1)
  
  #Event-level alphas:
  alpha.local <- matrix(0, n, p)
  gamma.partial <- -sum(gamma)/(sigma)
  for (k in 1:d) {
    indx <- which(event==unique(event)[k])
    alpha.local[indx,] <- matrix(alpha[p*(k-1) + 1:p], length(indx), p, byrow=TRUE)
  }
  
  #Just compute this once to save time:
  eta <- exp(gamma %*% t(beta) + alpha.local)
  
  #Gradient in the direction of event-specific alphas:
  grad <- vector()
  for (k in 1:d) {
    indx <- which(event==unique(event)[k])
    grad <- c(grad, as.integer(!specific) * colSums(data[indx,] - eta[indx,], na.rm=TRUE))
  }
  
  #Gradient in the direction of beta:
  grad <- c(grad, colSums(sweep(data - eta, 1, gamma, '*'), na.rm=TRUE))
  
  #Gradient in the direction of gamma, set these to zero and deal with them in IRLS:
  # grad = c(grad, rowSums(sweep(data - eta, 2, beta, '*'), na.rm=TRUE)+gamma.partial)
  grad <- c(grad, rep(0,n))
  
  #Gradient in the direction of sigma:
  grad <- c(grad, -n/2/sigma + sum(gamma^2)/2/sigma^2)
  
  grad
}
