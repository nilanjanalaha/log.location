#################################### MLE of theta #############################################
library(logcondens.mode)
giveth4ml=function(inth,x)
{
  n <- length(x)
  x=x-inth
  ind=sign(x)
  sortx <- sort(abs(x))
  flag=0
  if(sortx[1]!=0) 
  { flag=1
  sortx=c(0,sortx)
  }
  mlef= logcondens.mode::activeSetLogCon.mode(sortx,print=FALSE)
  L= sum(mlef$phi)/2
  L
}

# This function is called from the workspace
#Function for calculating the location estimator
#' The  MLE estimator of Laha (2020).
#'
#' Suppose n univariate observations are sampled from a density \eqn{f(x-m)} where
#' m is the location parameter and f is an unknown symmetric density. This function computes
#' a one step estimator to estimte m. This estimator uses
#' the log-concave MLE estimator from the package \code{\link{logcondens.mode}} to estimate f.
#' 
#' 
#' @param x       An array of length n; represents the data.
#' @param  alpha  The confidence level for the confidence bands. An (1-alpha) percent
#'               confidence interval is constructed. Alpha should lie in the interval (0, 0.50). The default
#'               value is 0.05.
#' 
#'                        
#'@return   A list of length two.   
#'\itemize{
#'\item \code{estimate:}  A scalar.
#'\item \code{CI:} A vector of length two.
#'}
#' @seealso \code{\link{p.mle}}, \code{\link{s.sym}}
#'@references Laha N. (2020). \emph{Location estimation for symmetric and
#' log-concave densities}. submitted.
#'@author \href{https://connects.catalyst.harvard.edu/Profiles/display/Person/184207}{Nilanjana Laha}
#' (maintainer), \email{nlaha@@hsph.harvard.edu}.
#' @examples
#' x <- rnorm(100); MLE_location(x)
#' @export 
MLE_location <- function(x, alpha = 0.05)
{
  x <- sort(x)
  inth <- mean(x)
  
  inx <- x
  n <- length(x)
  x=x-inth
  ind=sign(x)
  sortx <- sort(abs(x))
  flag=0
  if(sortx[1]!=0) 
  { flag=1
  sortx=c(0,sortx)
  }
  mlef= logcondens.mode::activeSetLogCon.mode(sortx,print=FALSE)
  phip <- diff(mlef$phi)/diff(mlef$x)
  Fhi <- mlef$Fhat
  kn<- sortx[mlef$IsKnot==1]
  kn <- sort(c(-kn,kn))
  
  rankv <- rank(abs(x)) # getting the ranks of x-th
  
  #calculating the L4
  L4 <- x
  L4[sign(x)<0] <- phip[(rankv[sign(x)<0])]
  L4[sign(x)>0] <- -phip[rankv[sign(x)>0]]
  Fh <-x
  Fh[sign(x)<0] <- 1/2*(1-Fhi[rankv[sign(x)<0]+1])
  Fh[sign(x)>0] <- 1/2*Fhi[rankv[sign(x)>0]+1]+1/2
  th4full = giveone_t(L4,inth,n)
  
  # Grids for MLE
  #Calculating I 
  I <- sum(L4^2)/n
  # The grid which depends on I
  evar <- 3/sqrt(I*n)
  thv <- seq(-evar+th4full,th4full+evar,length.out=500)
  Lvec <- sapply(thv,giveth4ml,x=inx)
  ml <- thv[which.max(Lvec)]
  CI <- conf.int(x, ml, alpha)
list( estimate=ml, CI=CI)
}

# Function for finding the confidence interval
conf.int <- function(x, ml, alpha)
{
  x <- x-ml
  x <- abs(x)
  x <- c(0,x)
 mlef <- logcondens.mode::activeSetLogCon.mode(x, mode=x[1])
 L4 <-   diff(mlef$phi)/diff(mlef$x)
 x <- mlef$Fhat
 len <- diff(x)
 I <- sum(L4^2*len)
 sd1 <- 1/sqrt(I)
 n <- length(x)-1
 #Gives almost similar result
 #I <- sum(L4^2)/length(L4)
# sd2 <-  1/sqrt(I)
 c(ml+qnorm(alpha/2)*sd1/sqrt(n),ml-qnorm(alpha/2)*sd1/sqrt(n))
}

