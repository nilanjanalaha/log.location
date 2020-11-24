# This file contains functions for choosing the
# Tuning parameters for Beran's and Stone's estimators

########### Stone's estimator selector  ######################################
#' Tuning parameter selector for Stone (1975)'s location estimator 
#'
#' \code{\link{giveth}} computes the location estimator in a location
#' shift model using Stone (1975)'s estimator. This function depends on
#' tuning parameters dn and tn.  We estimate the MSE of the estimators and
#' to choose the best (dn,tn) from a set of options.
#' 
#' To this end, we generate B Bootstrap samples and for each pair
#' of (dn, tn), we estimate the MSE by computing
#' \deqn{\frac{\sum_{i=1}^B(\hat\theta_i(dn,tn)-inth)^2}{B}}
#' where \eqn{\theta_i(dn,tn)} is the estimator based on the i-th bootstrap sample.
#' 
#' For asymptotic efficiency of the estimators,
#' a) dn\eqn{\to\infty} and tn\eqn{\to 0} b)
#' \deqn{\frac{(dn)^2}{n^{1-\epsilon}(tn)^6}=O(1)}
#' for some \eqn{\epsilon>0}. We take dn to be the minimum of
#'  \eqn{\min(D n^{1/2}(tn)^3 } and max(|x|)+2 (sd(x)*tn). Since dn is
#'  the truncation parameter and Stone (1975) uses kernel smoothening,
#'   we use the bound max(|x|)+2 (sd(x)tn).
#' and we take tn to be the minimum between one and \eqn{t n^{-1/7}/100}. 
#' We generate a set of (dn, tn) by varying
#' D and t acrros a grid. The default choice of the grid for t is c(0.01, seq(0.10, 0.80, by=0.05))
#' and that for D is
#' (0.5, 1, 1.5, 2, 3, 4,....., 20).
#' 
#' @param x An array of length n; the dataset.
#' @param inth A number; the initial estimator. The default is the median.
#' @param D Optional, an array of real numbers, contains values of D, parameter needed for tuning dn. See details.
#' @param t Optional, an array of real numbers, contains values of t, parameter needed for tuning tn. See details.
#' @param B Optional, the number of bootstrap samples. The default is 100.
#' @seealso \code{\link{giveth}}
#' @return   A list is returned:
#' \itemize{
#' \item \code{param:} A 2 columns, the number of rows is same as the
#'                    length of inth. Each row gives the optimal
#'                    (dn, tn) for the corresponding inth value.
#' \item \code{estimte:}  Gives the estimators of \eqn{\theta}
#'       using the optimal (dn, tn)'s.
#' \item\code{default_estimate:}A
#' }
#'@references Laha, N. \emph{ Location estimation fr symmetric and 
#'           log-concave densities}. Submitted.
#'@references Stone, C. (1975). \emph{Adaptive maximum likelihood estimators
#'of a location parameter}, Ann. Statist., 3, 267-284.
#'@author \href{https://connects.catalyst.harvard.edu/Profiles/display/Person/184207}{Nilanjana Laha}
#' (maintainer), \email{nlaha@@hsph.harvard.edu}.
#' @seealso \code{\link{giveth}}
#' @examples {}
#'   @export
stone.select <- function(x, inth, D, t, B)
{
  #Missing values
  my.grid <- seq(1,30, by=2)
  
  if(missing(D)) D <- seq(10,80, by=10)
  if(missing(t)) t <- seq(0.1,0.8, by=0.1)
  if(missing(B)) B <- 100
  if(missing(inth)) {
    warning("the preliminary estimator inth is missing; using the sample
            median")
    inth <- median(x)
  }
  
  # All possible (D, t) pairs
  #truncation
  n.D <- length(D)
  #smoothening
  n.t <- length(t)
  gridmat <- numeric()
  for(i in 1: n.D)
  {
    gridmat <- rbind(gridmat, cbind(rep(D[i], n.t), t))
  }
  
  th <- matrix(0, nrow(gridmat), length(inth))
  for(i in 1: nrow(gridmat))
  {
    #smoothening parameter
    tn <- gridmat[i, 2]
    #truncation parameter
    dn <- gridmat[i, 1]
    print(dn)
    th[i,] <- giveth(x, inth, tn, dn)
  }
thmat <- parallel::mclapply(1:B, stone.select.in, x=x, gridmat=gridmat, inth=inth, th=th)

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
tn <- mydntn[1]
dn <- min(mydntn[2]*n^{1/2}* (tn)^3, max(abs(x))+2*sd(x)*tn)
#The optimal dn and tn estimator and the default dn and tn estimator
c(mydntn, giveth(x, inth, mydntn[1], dn), giveth(x, inth))
}

stone.select.in <- function(seed, x, gridmat, inth, th)
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
    dn <- gridmat[i, 1]
    thvec[i,] <- (giveth(x, inth, tn, dn)- inth)^2
   }
   set.seed(NULL)
   thvec
}
########### Stone's estimator selector (discarded) ######################################
#' Tuning parameter selector for Stone (1975)'s location estimator 
#'
#' \code{\link{giveth}} computes the location estimator in a location
#' shift model using Stone (1975)'s estimator. This function depends on
#' tuning parameters dn and tn.  We estimate the MSE of the estimators and
#' to choose the best (dn,tn) from a set of options.
#' 
#' To this end, we split the dataset in tan parts. Then we compute  estimators
#' \eqn{\hat \theta_i(dn,tn)} for each (dn,tn) under consideration from the i-th part
#' of the data, where i=1,...,10. We then estimate the MSE corresponding to each (dn, tn)
#' by computing
#' \deqn{\frac{\sum_{i=1}^10(\hat\theta_i(dn,tn)-\hat\theta(dn, tn))^2}{10}.}
#' We choose the (dn, tn) pairs which minimize the estimated MSE.
#' 
#' For asymptotic efficiency of the estimators,
#' a) dn\eqn{\to\infty} and tn\eqn{\to 0} b)
#' \deqn{\frac{(dn)^2}{n^{1-\epsilon}(tn)^6}=O(1)}
#' for some \eqn{\epsilon>0}. We take \eqn{dn=Dn^{1/2}(tn)^3}
#' and \eqn{tn=t n^{-1/7}}. We generate a set of (dn, tn) by varying
#' D and t in the set 
#' (0.5, 1, 1.5, 2, 3, 4,....., 20).
#' 
#' @param x An array of length n; the dataset.
#' @param inth A vector of initial estimators
#' @param D An array of real numbers, contains values of D, parameter needed for tuning dn. See details.
#' @param t An array of real numbers, contains values of t, parameter needed for tuning tn. See details.
#' @seealso \code{\link{giveth}}
#' @return   A list is returned:
#' \itemize{
#' \item \code{param:} A matrix with 2 columns, the number of rows is same as the
#'                    length of inth. Each row gives the optimal
#'                    (dn, tn) for the corresponding inth value.
#' \item \code{estimte:} An array of same length as inth. Gives the estimators of \eqn{\theta}
#'       using the optimal (dn, tn)'s.
#' }
#'@references Laha, N. \emph{ Location estimation fr symmetric and 
#'           log-concave densities}. Submitted.
#'@references Stone, C. (1975). \emph{Adaptive maximum likelihood estimators
#'of a location parameter}, Ann. Statist., 3, 267-284.
#'@author \href{https://connects.catalyst.harvard.edu/Profiles/display/Person/184207}{Nilanjana Laha}
#' (maintainer), \email{nlaha@@hsph.harvard.edu}.
#' @examples
#' x <- rnorm(100); inth <- mean(x);
#'  giveth(x, inth=inth)
#'  giveth(x)
stone.select.ww <- function(x, inth)
{
  for(th in inth)
    catch.th <- stone.select.in(x, th)
}

stone.select.in.ww <- function(x, th)
{
  #Missing values
  my.grid <- c(0.5, 1, 1.5, 2, seq(3, 20, by=1))
  if(missing(D)) D <- my.grid
  if(missing(t)) t <- my.grid
  
  
  n <- length(x)
  #Randomly shuffle the data
  x <- base::sample(x)
  
  # Partitioning the data in ten equal parts
  folds <- cut(seq(1,n),breaks=10,labels=FALSE)
  test <- list()
  for(i in 1:10){
    #Segement x
    testIndexes <- which(folds==i, arr.ind=TRUE)
    test[[i]] <- x[testIndexes, ]
      }
  
  
  
}
  