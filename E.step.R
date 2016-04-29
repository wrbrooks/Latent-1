#' Expectation step
#'
#' Replace censored values by their expectations of censored FIB counts
#' under the current parameter values, conditional on the counts being
#' below the minimum level of detection.
#' 
#' This function imputes the expectations of censored data values and
#' replaces the censored observations by their expectations, conditional on
#' the censored values being below the minimum level of detection.
#' The expectation is with respect to a Poisson 
#' distribution. It is calculated via the estimated parameters 
#' \eqn{(\hat{\mathbf{\alpha}}, \hat{\mathbf{\beta}}, \hat{\mathbf{\gamma}})}
#' through the relationship
#' \deqn{\eta = \alpha + \beta \gamma}
#' \deqn{Y ~ Poisson(\exp(\eta)).}
#' 
#' Let \eqn{c_j} be the minimum level of detection for column \eqn{j} of the
#' data, and suppose that observation \eqn{Y_{ij}} is censored, meaning that
#' our only knowledge of \eqn{Y_{ij}} is that \eqn{Y_{ij} < c_j}. Then the 
#' conditional expectation \eqn{E (Y_{ij} | Y_{ij} < c_j)} is
#' \deqn{E (Y_{ij} | Y_{ij} < c_j) \sum_{k=1}^c_j \{ k \lambda_{ij}^k 
#'     \exp (-\lambda_{ij}) \Gamma(k+1)^{-1} \} / \{ \sum_{k=1}^n
#'     \lambda_{ij}^k \exp (-\lambda_{ij}) \Gamma(k+1)^{-1} \},}
#' where \eqn{\lambda_{ij} = \alpha_{ij} + \beta_j \gamma_i}, \eqn{\alpha_{ij}}
#' is the intercept for observation \eqn{i}, bacteria species \eqn{j},
#' \eqn{\beta_j} is the slope for bacteria species \eqn{j} (meaning that
#' \eqn{\beta_j} is the marginal increase in expected long-concentration of 
#' species \eqn{j} when the contamination increases by one unit), and
#' \eqn{\gamma_i} is the contamination index of observation \eqn{i}.
#'
#' @param alpha vector of intercepts where each entry is the intercept for one species/event combination
#' @param beta slope vector with one entry per FIB species, which indicates the marginal increase in that species' log-mean for a unit increase in the contamination index
#' @param gamma contamination vector with one entry per row of data, indicating the contamination index for the corresponding row
#' @param data matrix of FIB counts with one row per sample and one column per FIB species
#' @param min.detect minimum detectable level of the FIB assay - counts below this threshold are censored 
#' @param event vector with one entry per row of data where each entry indicates the event with which the data row is associated
#'
#' @return A matrix of FIB counts, where the censored counts are replaced by their expectations under the current parameters.
#'
E.step = function(alpha, beta, gamma, data, min.detect, event) {
  #Basic constants:
  n = nrow(data)
  p = ncol(data)
  d = length(unique(event))
  
  #Event-level alphas:
  alpha.local = matrix(0, n, p)
  for (k in 1:d) {
    indx = which(event==unique(event)[k])
    alpha.local[indx,] = matrix(alpha[p*(k-1) + 1:p], length(indx), p, byrow=TRUE)
  }
  
  for (j in 1:ncol(data))
    for (t in 1:nrow(data))
      if (!is.na(data[t,j]))
        if (data[t,j] <= min.detect[j]) {
          #Get the probability that the count is 0,1,...,MLD
          #Where MLD is the minimum level of detection (censoring threshold).
          cens.log.lik = dpois((0:min.detect[j]), exp(gamma[t]*beta[j] +gamma[t+n]*beta[j+p] + alpha.local[t,j]), log=TRUE)
          cens.log.lik = cens.log.lik - max(cens.log.lik, na.rm=TRUE)
          
          #Compute the expected count, conditional on count being no greater than censoring threshold
          #If the total probability below the censoring threshold is indistinguishable from zero,
          #then set the expectation to the censoring threshold.
          if (sum(exp(cens.log.lik))==0) {
            data[t,j] = min.detect[j]
          } else data[t,j] = round(sum((0:min.detect[j]) * exp(cens.log.lik), na.rm=TRUE) / sum(exp(cens.log.lik), na.rm=TRUE))
        }
  
  return(data)
}