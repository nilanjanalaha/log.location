#Function for calculating the location estimator
#' Smoothed symmetrized estimator of location from Laha (2020).
#'
#' Suppose n univariate observations are sampled from a density \eqn{f(x-m)} where
#' m is the location parameter and f is an unknown symmetric density. This function computes
#' a one step estimator to estimte m. This estimator uses
#' the smoothed log-concave MLE estimator from the package \code{\link{logcondens}} to estimate f, and  is root-n consistent
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
 s.sym  <- function(x, q=0, inth, alpha=0.05) 
{
  x <- sort(x)
  if(missing(inth)) inth <- mean(x)
  n <- length(x)
  mlef <- logcondens::logConDens(x,smoothed=FALSE,print=FALSE)
  J=asign(x,mlef)
  phi=mlef$phi
  b=-(gam(x,J,phi,inth)-var(x))
  b <- sqrt(b)
  
  L3=sapply(x,L3in,J=J,x=x,b=b,inth=inth,phi=phi)
  ans <- sapply( q, th3_truncate, L3=L3, inth=inth, n=n, alpha=alpha, mlef=mlef, x=x)
  list( est <- data.frame(q=q, est=ans[1,]), CI <- data.frame(q=q, lb=ans[2, ], ub=ans[3,]))
  
 }
 
 th3_truncate <- function(q, L3, inth, n, alpha, mlef, x)
 {
   if(q==0)
   {
     #estimate
     th3 = giveone_t(L3,inth,n)
     #CI
     I = mean(L3^2)
     lb <- th3+qnorm(alpha/2)/sqrt(I*n)
     ub <- th3-qnorm(alpha/2)/sqrt(I*n)
     return(c(th3, lb, ub))
   } else
   {
     q <- 1-q
     qt <- findq2(q, mlef, inth)
     if(qt==-1)
     { vecr <- 1:n } else { vecr <- c(findlow(-qt+inth,x): findup(inth+qt,x))}
     th3=giveone_t(L3[vecr],inth,n)
     #CI
     I <- sum(L3[vecr]^2)/n
     lb <- th3 + qnorm(alpha/2)/sqrt(n*I)
     ub <- th3 - qnorm(alpha/2)/sqrt(n*I)
     return(c(th3, lb, ub))
   }
 }
################################################################
################# Auxiliary functions ##########################
################################################################
 
g3in=function(y,J,x,b,phi)
{
  xp=x
  n=length(x)
  xp[1:(n-1)]=x[2:n]
  xp=xp[1:(n-1)]
  x=x[-n]
  phi=phi[-n]
  temp=exp(((J*b)^2)/2+y*J-x*J+phi)
  t2=pnorm((xp-J*b^2-y)/b)
  t3=pnorm((x-J*b^2-y)/b)
  l2=1/b*dnorm((xp-J*b^2-y)/b)
  l3=1/b*dnorm((x-J*b^2-y)/b)
  f=(temp*(t2-t3))
  ans=sum(f)
  t4=temp*(l3-l2)
  fhat=t4+J*f
  ans=c(ans,sum(fhat))
  ans
}
L3in=function(y,J,x,b,inth,phi)
{
  f=g3in(y,J,x,b,phi)[1]+g3in(2*inth-y,J,x,b,phi)[1]
  f=f/2
  fhat=g3in(y,J,x,b,phi)[2]-g3in(2*inth-y,J,x,b,phi)[2]
  fhat=fhat/2
  -fhat/f
}

L3inf=function(y,J,x,b,inth,phi)
{
  f=g3in(y,J,x,b,phi)[1]+g3in(2*inth-y,J,x,b,phi)[1]
  f/2
}



###################################################################################################################
# giving the density function smoothed

giveth3tf=function(x,mlef,inth) 
{
  n <- length(x)
  J=asign(x,mlef)
  phi=mlef$phi
  b=-(gam(x,J,phi,inth)-var(x))
  b <- sqrt(b)
  f <- sapply(x,L3inf,J=J,x=x,b=b,inth=inth,phi=phi)
  
  x2=2*inth-x
  cbind(c(x2,x),c(f,f))
  
}


gam=function(x,J,phi,inth)
{
  n=length(x)
  phi=phi[-n]
  x2=x[2:n]
  mea=mean(x)
  x=x[-n]
  sum(exp(phi-x*J)*gamv(x,x2,J))-mea^2
}

gamvin=function(a,b,c)
{
  if(c==0)
    return((b^3-a^3)/3)
  t1=(b^2*exp(b*c)-a^2*exp(a*c))/c
  t2=(b*exp(b*c)-a*exp(a*c))/c^2
  t3=(exp(b*c)-exp(a*c))/c^3
  t1-2*(t2-t3)
}
gamv=function(a,b,c) mapply(gamvin,a,b,c)