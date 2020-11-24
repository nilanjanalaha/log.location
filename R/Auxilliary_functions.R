

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

findup <- function(x,vec)
{
  i <- length(vec)
  while (vec[i]>x)
  {
    i=i-1
  }
  return(i)
}

findlow <- function(x,vec)
{
  i <- 1
  while (vec[i]<x)
  {
    i=i+1
  }
  return(i)
}

asign=function(x,mlef)
{
  n=length(x)
  phi=mlef$phi
  sup=mlef$x
  m=length(sup)
  if(m!=n)
    print("Warning: there are probably ties in data")
  diff(phi)/diff(x)
}

giveone_t=function(L,inth,n)
{
  I=sum(L^2)/n
  inf=L/I
  inth+mean(inf)
}
# find if x is in the interval given by y
Ininterval <- function(x, y)
{
  if(x <= y[1] | x>= y[2]) return(0)
  1
}

#######################################################################
############## Symmetrized beta density and information ###############
#######################################################################

deng <- function(x,r)
{
  f <- gamma((3+r)/2)/(sqrt(pi*r)*gamma(1+r/2))*(1-x^2/r)^(r/2) # The density of the symmetrized beta
  phip <- -x/(1-x^2/r)
  f*phip^2
}
nlg <- function(eta,r)
{
  L <- sqrt(r)*(2*qbeta(eta,1+r/2,1+r/2)-1)
  U <-sqrt(r)*(2*qbeta(1-eta,1+r/2,1+r/2)-1)
  integrate(deng, lower=L,upper=U,r=r)$value 
}

lg <- function(r,eta) 
{
  sapply(eta,nlg,r=r)/integrate(deng, lower=-sqrt(r),upper=sqrt(r),r=r)$value
}

varif <- function(r)
{
  r/(r+3)
}

# The information of symmetrized beta with r
information_of_r_beta <- function(r) integrate(deng, lower=-sqrt(r),upper=sqrt(r),r=r,rel.tol=10^(-8))$value

# density of sym(r)
dsbeta <- function(x,r) 
{
  ifelse(abs(x)>sqrt(r),0, gamma((3+r)/2)/(sqrt(pi*r)*gamma(1+r/2))*(1-x^2/r)^(r/2))
}

