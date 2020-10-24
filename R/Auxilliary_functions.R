library(rmutil)
library(logcondens)

x <- sort(rgamma(500, shape=0.5, 2))
#' Scores of log-concave MLE
#'
#'Calculates \deqn{\frac{\hat\phi_n(X_{i+1})-\hat\phi_n(X_i)}{X_{i+1}-X_i}}
#'for i=1,...,n., where \eqn{\hat\phi_n=\log\hat f_n}, the latter being the
#'unsmoothed log-concave MLE given by \code{\link{logcondens}}.
#'
#' @param x a vector of length n; the dataset
#' @return   scores evaluated in between the data points, a vector of size n-1.
#' @export        
mle_score <- function(x)
{
mlef <- logcondens::logConDens(x, smoothed=FALSE)

#Evaluating the scores
my.phi <- mlef$phi
score <- diff(my.phi)/diff(x)
score
}

