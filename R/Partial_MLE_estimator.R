#Function for calculating the location estimator
#' The partial MLE estimator of Laha (2020).
#'
#' Suppose n univariate observations are sampled from a density \eqn{f(x-m)} where
#' m is the location parameter and f is an unknown symmetric density. This function computes
#' a one step estimator to estimte m. This estimator uses
#' the log-concave MLE estimator from the package \code{\link{logcondens.mode}} to estimate f, and  is root-n consistent
#' for m provided f is log-concave. 
#' 
#' 
#' @param x       An array of length n; represents the data.
#' @param q A fraction between 0 and 1/2. Corresponds to the truncation parameter.
#'              The default is 0, which indicates no truncation.
#' @param  init Optional. An initial estimator of m. The default value is the sample median.
#' @param  alpha  The confidence level for the confidence bands. An (1-alpha) percent
#'               confidence interval is constructed. Alpha should lie in the interval (0, 0.50). The default
#'               value is 0.05.
#' 
#' @details \code{q:} If q is positive, the function
#'                      (1-2*q) percent observations from both tails while
#'                      computing the one step estimator. The scores are
#'                      estimated using the full data. See Laha et al. for more details.
#' @details \code{init:}  The default is mean. If init is the median, some jitter
#'                        (0.0001) is added.
#'                        
#'@return   A list of length two.   
#'\itemize{
#'\item \code{estimate:}  A matrix of two columns, and the rownumber equals the length 
#'                        of q. The firs column is q, and the second column is estimates corresponding to q.
#'\item \code{CI:} A matrix of three columns, the first column is q,
#'the second column is the left point and the third column is the right point
#'of the confidence intervals corresponding to q.
#'}
#'  
#'@references Laha N. (2020). \emph{Location estimation for symmetric and
#' log-concave densities}. submitted.
#'@author \href{https://connects.catalyst.harvard.edu/Profiles/display/Person/184207}{Nilanjana Laha}
#' (maintainer), \email{nlaha@@hsph.harvard.edu}.
#' @examples
#' x <- rlogis(100); p.mle(x, q=c(0, 0.001, 0.01))
#' @export 
 p.mle  <- function(x, init, q=0, alpha=0.05)
{
   if(missing(init)) init <- mean(x)
   if(missing(q)) q <- 0
   #Some jitter is added if init is the median
   if(init == median(x)) init <- init+0.0001
     
   
   #computation starts here
   inth <- init
  n <- length(x)
  x <- sort(x)
  x=x-inth
  ind=sign(x)
  sortx <- sort(abs(x))
  flag=0
  if(sortx[1]!=0) 
  { flag=1
  sortx=c(0,sortx)
  }
  mlef= logcondens.mode::activeSetLogCon.mode(sortx,print=FALSE)
  
  #computing the derivative of logf
  phip <- diff(mlef$phi)/diff(mlef$x)
  Fhi <- mlef$Fhat
  
  # what if mlef works on some x , less size than sortx
 
  
  rankv <- rank(abs(x)) # getting the ranks of abs(x)
  
  #calculating the L4
  L4 <- x
  L4[sign(x)<0] <- phip[(rankv[sign(x)<0])]
  L4[sign(x)>0] <- -phip[rankv[sign(x)>0]]
  Fh <-x
  Fh[sign(x)<0] <- 1/2*(1-Fhi[rankv[sign(x)<0]+1])
  Fh[sign(x)>0] <- 1/2*Fhi[rankv[sign(x)>0]+1]+1/2
  
 #Calculating the estimate and the confidence intervals 
 ans <- sapply( q, th4_truncate, L4=L4, inth=inth, n=n, alpha=alpha, mlef=mlef, x=x)
 list( estimate = data.frame(q=q, est=ans[1,]), CI = data.frame(q=q, lb=ans[2, ], ub=ans[3,]))
  
 }
 
 th4_truncate <- function(q, L4, inth, n, alpha, mlef, x)
 {
   if(q==0)
   {
     #estimate
     th4 = giveone_t(L4,inth,n)
     #CI
     I = mean(L4^2)
     lb <- th4+qnorm(alpha/2)/(sqrt(I)*n^(1/2))
     ub <- th4-qnorm(alpha/2)/(sqrt(I)*n^(1/2))
     return(c(th4, lb, ub))
   } else
   {
    q <- 1-q
     up <- quantilesLogConDens(2*q-1, mlef)[2]
     #Note: though not documented, it actually works for logcondens.mode
     #See the manual of the latter package, p. 11 for instance
     low <- -up
     upl <- findup(up,x)
     vecr <- c(findlow(low,x):upl)
     #mean
     th4=giveone_t(L4[vecr],inth,n)
     #CI
     I <- sum(L4[vecr]^2)/n
     lb <- th4+qnorm(alpha/2)/(sqrt(I)*n^(1/2))
     ub <- th4-qnorm(alpha/2)/(sqrt(I)*n^(1/2))
     return(c(th4, lb, ub))
   }
 }
 
 