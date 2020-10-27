# This file contains functions for choosing the
# Tuning parameters for Beran's  estimators

########### Beran's estimator selector  ######################################
#' Tuning parameter selector for Beran (1974)'s location estimator 
#'
#' \code{\link{beran.est}} computes the location estimator in a location
#' shift model using Beran (1974)'s estimator. This function depends on
#' tuning parameters theta and M.  We estimate the MSE of the estimators and
#' to choose the best (theta, M) from a set of options.
#' 
#' To this end, we generate B Bootstrap samples and for each pair
#' of (theta, M), we estimate the MSE by computing
#' \deqn{\frac{\sum_{i=1}^B(\hat\theta_i(theta, M)-\theta)^2}{B}}
#' where \eqn{\theta_i(theta, M)} is the estimator based on the i-th bootstrap sample
#' and \eqn{\theta} is the true theta.
#' 
#' We take theta to be in c(0.01, seq(0.10, 0.80, by=0.05)) and take M to
#' be in (2, 4, 6, ..., 10).
#' 
#' @param x An array of length n; the dataset.
#' @param inth A number; the initial estimator. The default is the median.
#' @param Mvec Optional, an array of integers, the number of basis functions to use. See details.
#' @param t Optional, an array of real numbers, contains values of theta, parameter needed for 
#'          estimating the Fourier coefficient.
#' @param B Optional, the number of bootstrap samples. The default is 100.
#' @seealso \code{\link{beran.est}}
#' @return   A vector of foure numbers.
#' \itemize{
#' \item \code{param:} A 2 columns, the number of rows is same as the
#'                    length of inth. Each row gives the optimal
#'                    (dn, tn) for the corresponding inth value.
#' \item \code{estimte:}  Gives the estimators of \eqn{\theta}
#'       using the optimal (dn, tn)'s.
#' }
#'@references Laha, N. \emph{ Location estimation fr symmetric and 
#'           log-concave densities}. Submitted.
#'@references Stone, C. (1975). \emph{Adaptive maximum likelihood estimators
#'of a location parameter}, Ann. Statist., 3, 267-284.
#'@author \href{https://connects.catalyst.harvard.edu/Profiles/display/Person/184207}{Nilanjana Laha}
#' (maintainer), \email{nlaha@@hsph.harvard.edu}.
#' @examples {}
#'   @export
beran.select <- function(x, t, Mvec, inth, B)
{
  #Missing values
  my.grid <- seq(2, 10, by=2)
  
  if(missing(Mvec)) Mvec <- my.grid
  if(missing(t)) t <- c(0.01, seq(0.10, 0.80, by=0.05))
  if(missing(B)) B <- 100
  if(missing(inth)) {
    warning("the preliminary estimator inth is missing; using the sample
            median")
    inth <- median(x)
  }
  
  # All possible (D, t) pairs
  #truncation
  n.D <- length(Mvec)
  #smoothening
  n.t <- length(t)
  gridmat <- numeric()
  for(i in 1 : n.D)
  {
    gridmat <- rbind(gridmat, cbind(rep(Mvec[i], n.t), t))
  }
  thmat <- parallel::mclapply(1:B, beran.select.in, x=x, gridmat=gridmat, inth=inth)
  
  #Selecting the one with the smallest MSE
  mse <- 1:nrow(gridmat)
  for(j in 1:nrow(gridmat))
  {
    ave <- 1:B
    for(b in 1:B)
    {
      ave[b] <- thmat[[b]][j]
    }
    mse[j] <- mean(ave)
  }
  
  #Selecting the best one
  mydntn <- gridmat[which.min(mse), ]
  Mn <- mydntn[1]
  tn <- mydntn[2]
  #The optimal dn and tn estimator and the default dn and tn estimator
  list( param=mydntn, estimate=beran.est(x, theta=tn, M=Mn, inth, 0.05))
}

beran.select.in <- function(seed, x, gridmat, inth)
{
  n <- length(x)
  set.seed(seed)
  thvec <- matrix(0, nrow(gridmat), length(inth))
  x <- sample(x, replace = TRUE)
  for(i in 1: nrow(gridmat))
  {
    #smoothening parameter
    tn <- gridmat[i, 2]
    #truncation parameter
    Mn <- gridmat[i, 1]
    thvec[i,] <- (beran.est(x, theta=tn, M=Mn, inth, 0.05)$estimate - 0)^2
  }
  set.seed(NULL)
  thvec
}