#----------calculating fourier coefficient for the score function------
#' Fourier coefficients of the score in a location model
#'
#' Estimates the Fourier coefficients of the scores in a location shift model,
#' with is the class of densities obtained by a location shift of
#' a fixed unknown density \eqn{f}. The score,
#'  given by 
#' \deqn{\phi'(x)=-f'\circ F^{-1}(t)/f\circ F^{-1}(x),}
#' does not depend on the unknown shift. The Fourier coefficients of the score
#' is estimated from the data, whose density is assumed to be a location
#' shift of the unknown density \eqn{f}. We use Beran (1974)'s nonparametric
#' estimator here.
#' 
#' @param x       An array of length n; represents the data from a location model whose score the Fourier coefficients
#'                seek to estimate. 
#' @param theta A small number, should be of order \eqn{O_p(n^{-1/2})}. The default
#'              is \eqn{4n^{-1/2}}
#' @param indices An array of positive integers, for each integer j in "indices", 
#'                Fourier coefficient  corresponding to the
#'                basis function \eqn{t\mapsto\exp(i2\pi jt)} is computed.
#' @param  which Optional. Takes value 1 or 2. If "which" is 1, only the real
#'                part of the Fourier coefficient is computed. If "which" is 2,
#'                only the imaginary part of the k-th coefficient is calculated. 
#'                The default is to calculate both real and imaginary
#'                parts.   
#' @details \code{theta:}   theta do not need to depend on the range of the data because
#'                         the estimators depend only on the rank of the dataponts.  If \eqn{theta=z_n} which 
#'                         equals \eqn{z n^{-1/2}}, then the estimated coefficients are 
#'                         root-n consistent by Theorem 2.1 of Beran (1974).
#' @details \code{which:}   If it is known that the density \eqn{f} is symmetric,
#'                          then the Fourier coefficients are real. Therefore, there is
#'                          no need to calculate the imaginary part. Similarly, if the density 
#'                          is odd, then the Fourier coefficinets of the scores are purely  
#'                          imaginary. Otherwise, one generally requires both the real and the imaginary parts, 
#'                          and "which" should be left unspecified in those cases.                            
#'@return       
#' \itemize{
#' \item If "which=1", an array (with the same length as "indices"), give the real parts of the Fourier coefficients.
#' \item If "which=2, an array (with the same length as "indices"), give the  imaginary parts of the Fourier coefficients.
#' \item The default is to return a matrix whose first and second column give the
#'       real and imaginary parts of the Fourier coefficients, respectively. The number
#'       of rows equal the length of the array "indices". 
#' }    
#'@references Beran, R. (1974). \emph{Asymptotically efficient adaptive 
#'rank estimates in location models}. Ann. Statist., 2, 63-74.
#'@author \href{https://connects.catalyst.harvard.edu/Profiles/display/Person/184207}{Nilanjana Laha}
#' (maintainer), \email{nlaha@@hsph.harvard.edu}.
#' @examples
#' x <- rnorm(100); score.coeff(x, length(x)^(-1/2), 1)
#' @export          
score.coeff <- function(x,  theta, indices, which)
{
  n <- length(x) 
  # Missing values
  if(missing(which)) which <- 0
  if(missing(theta)) theta <- 4/sqrt(n)
  if(min(indices) <= 0) stop("indices must be an array of positive integers")
  
  # sorting and computing ECDF
   x<- sort(x)
  
   nk <- length(indices)
   im.part <-  real.part <- 1: nk
   
   # Finding the real part
   if(which !=2)
   {
   for (k in 1: nk)
   {
     phi <- function(x) cos(2*pi*k*x)
     real.part[k] <- give.coeff(x, phi, theta, n)
   }
   }
   
   if(which==1) return(real.part)
   
   #Finding the imaginary part
   for (k in 1: nk)
   {
     phi <- function(x) sin(2*pi*k*x)
     im.part[k] <- give.coeff(x, phi, theta, n)
   }
   if(which==2) return(im.part)
   
   # returning the output
   return(cbind(real.part, im.part))
}

#Function for calculating T_N(Z,\phi) for any phi, as in Beran's paper (2.12)
give.coeff <- function(x, phi, theta, n)
{
  e <- ecdf(x)
  G2 <- function(x, n) n*e(x)/(n-1)
  G1 <- function(x, n) G2(x, n)-1/(n-1)
  
  mean(phi(G1(x+theta, n)) - phi(G2(x-theta, n)))/(2*theta)
}

#Function for calculating the location estimator
#' The ocation estimator of Beran (1974)
#'
#' Suppose the observations are sampled from a density \eqn{f(x-m)} where
#' m is the location parameter and f is an unknown density. This function, based on
#' Beran (1974), gives a nonparametric estimator of m. This estimator uses Fourier basis
#' expansion to estimate the  score function corresponding to the location model, which has the form
#' \deqn{\phi_F(t)=f'(F^{-1}(t))/f(F^{-1}(t)).}
#'   The Fourier
#' coefficients of the score function  are computed 
#' by \code{\link{score.coeff}}.
#' 
#' @param x       An array of length n; represents the data.
#' @param theta A small number, should be of order \eqn{O_p(n^{-1/2})}. The default
#'              is \eqn{4n^{-1/2}}. This is the tuning parameter for Fourier coefficient
#'              estimation.
#' @param M     The number of Fourier basis to be used to approximate \eqn{\phi_F}.
#' @param  init Optional. A vector of initial estimators of m. The default value is the sample median.
#' @param  alpha  The confidence level for the confidence bands. An (1-alpha) percent
#'               confidence interval is constructed. Alpha should lie in the interval (0, 0.50). The default
#'               value is 0.05.
#' 
#' @details \code{theta:}   theta do not need to depend on the range of the data because
#'                         the estimators depend only on the rank of the dataponts.  If \eqn{theta=z_n} which 
#'                         equals \eqn{z n^{-1/2}}, then the estimated coefficients are 
#'                         root-n consistent by Theorem 2.1 of Beran (1974).
#' @details \code{M:}   A higher value of M decreases bias, but increases the variance. 
#'                      Theorem 4.1 of Beran (1974) states that if \eqn{M\to\infty} as
#'                      the sample size n grows   
#'                      and \eqn{\lim_{n\to \infty} M^6/n=0}, then this estimator is asymptotically efficient.
#'                      A larger M always gives a more conservative confidence interval.
#' @details \code{init:}  Beran (1974) recommends using sample median and warns that
#'                        the method will be sensitive to the choice of the preliminary estimator.
#'                        This should be a root-n consistent estimator of m.
#'                          Some other choices are the mean, or the trimmed mean. 
#'@return   A list of length two.   
#'\itemize{
#'\item \code{estimate:} An array of same length as init , giving the estimators of m
#'            based on the corresponding initial estimators in init. If init is missing,
#'              only one estimate is produced, which uses the sample median as the initial
#'              estimator.
#'\item \code{CI:} A matrix giving the (1-alpha) percent confidence intervals (CI). Each row corresponds to an initial estimator
#'                 in init (if missing, the sample median is used). The first column
#'                 corresponds to the lower CI, and the second column corresponds to the 
#'                 upper CI.
#'}
#'  
#'@references Beran, R. (1974). \emph{Asymptotically efficient adaptive 
#'rank estimates in location models}. Ann. Statist., 2, 63-74.
#'@author \href{https://connects.catalyst.harvard.edu/Profiles/display/Person/184207}{Nilanjana Laha}
#' (maintainer), \email{nlaha@@hsph.harvard.edu}.
#' @examples
#' x <- rlogis(100); beran.est(x, M=10)
#' beran.est(x, M=20)
#' @export 
beran.est <- function(x, theta, M, init, alpha)
{
  #Missing elements
  n <- length(x)
  if(missing(init)) init <- median(x)
  if(missing(theta)) theta <- 4*n^{-1/2}
  if(missing(M)) {
    warning("M is missing; default value 10 is being used")
    M <- 10
  }
  if(missing(alpha)) alpha <- 0.05
  
  x <- sort(x)
  #The fourier coefficients, we only need the imaginary part
  coeff <- score.coeff(x, theta, 1:M, which=2)
  #The number of basis indices
  basis <- 1:M
  phi <- function (t) 2*sum(coeff*sin(2*pi*basis*t))
  phi <- Vectorize(phi)
  hatvar <- sum(4*coeff^2)
 temp <- sapply(init, beran.est.in, f=phi, x=x, hatvar= hatvar, alpha=alpha)
 list(estimate=temp[1,], CI=t(temp[2:3,]))
}

beran.est.in <- function(x, th, f,  hatvar, alpha)
{
  my.x <- abs(x-th)
  #The ranks
  rank.v <- rank(my.x)
  n <- length(x)
  
  temp <- f(rank.v/(n+1))*sign(x-th)
  #The estimate
  est <- th+ mean(temp)/hatvar
  #CI
  lb <- est+qnorm(alpha/2)*sqrt(hatvar)
  ub <- est-qnorm(alpha/2)*sqrt(hatvar)
  c(est, lb, ub)
}

