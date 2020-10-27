
################## STONE ESTIMATOR ####################################
init <- function(x) c(mean(x), median(x), mean(x, trim=0.15))
s=function(x,inth) median(abs(x-inth))
m= function(x) median(x)
g= function(x) 2*dunif(x,-1,1)
fhin=function(t,rex,x,inth)
{
  ans= mean(dnorm(t+inth-x,sd=rex)+dnorm(-t+inth-x,sd=rex))/2
  ans
}
fh=function(t,rex,x,inth)
{
  sapply(t,fhin,rex=rex,x=x,inth=inth)
}
fhd=function(t,rex,x,inth)
{
  ans=-(t+inth-x)/rex^2*dnorm(t+inth-x,sd=rex)+(-t+inth-x)/rex^2*dnorm(-t+inth-x,sd=rex)
  mean(ans)/2
}
fhdvec=function(t,rex,x,inth) sapply(t,fhd,rex=rex,x=x,inth=inth)
Lh=function(t,rex,x,inth) {
  
  fhdvec(t,rex,x,inth)/fh(t,rex,x,inth)}

Ahin=function(t,rex,c,x,inth)  Lh(t,rex,x,inth)^2*g(t/c)*fh(t,rex,x,inth)
Ah=function(r,c,x,inth) pracma::integral(Ahin,xmin=-c,xmax=c,c=c,rex=r,x=x,inth=inth)
thfunc=function(y,rex,c,t,x,inth)  
{
  Lh(y+t-inth,rex,x,inth)*g((y+t-inth)/c)*dnorm(t,0,rex)
}

hthfuncin=function(t,rex,c,x,inth)
{
  yx=x
  vec=sapply(yx,thfunc,t=t,x=x,rex=rex,c=c,inth=inth)
  #print(mean(vec))
  mean(vec)
}
hthfunc=function(t,rex,c,xx,inth) 
{
  sapply(t,hthfuncin,rex=rex,c=c,x=xx,inth=inth)
}
pfnin=function(y,x,r,inth)
{
  n=length(x)
  mean( dnorm(y+inth-x,0,r))
}
pfn=function(y,x,r,inth) sapply(y,pfnin,x=x,r=r,inth=inth)
######################################MAIN FUNCTION#################################
givethin2=function(x,inth, tn, dn)
{
  x=sort(x)
  
  # Alert : the tuning parameters
  r=s(x,inth)*tn
  cc=s(x,inth)*dn
  xmin=-cc
  xmax=cc
  grid1=seq(xmin,xmax,length.out=61)
  v=Ah(r,cc,x,inth)
  grid2=grid1[2:61]
  grid1=grid1[1:60]
  gridv=grid2
  v1=pfn(gridv,x,r,inth)
  v2=Lh(gridv,r,x,inth)
  rim=(grid2-grid1)*v1*v2/v
  # rim=integral(hthfunc,xmin=-cc+inth-min(x),xmax=cc+inth-max(x),rex=r,c=cc,xx=x,method="Richardson",inth=inth)/v
  temp=inth-sum(rim)
  temp
}

# If no preliminary value is provided, it is chosen
init <- function(x) c(mean(x), median(x), mean(x, trim=0.05))

########### Stone's estimator ######################################
#' Stone (1975)'s estimator
#'
#'Gives a truncated one step estimator for  the location parameter in a symmetric
#' location family. This estimator is constructed using Stone (1975)'s
#' methods. This estimator uses the kernel density estimator (KDE) to estimate
#' the scores. Similar to all other one-step estimators, this estimator 
#' requires a preliminary estimator.
#'
#' @param x A vector of length n; the dataset
#' @param dn A parameter for the truncation. The default value is 20.
#' @param tn A parameter associated with the bandwidth of the kernel
#'           density estimator. The default value is 0.60.
#' @param inth An array giving some preliminary estimators for \eqn{\theta}. 
#'             The function calculates different one-step estimators
#'             for each of the preliminary estimators. The default value is
#'             an array containing the mean, the median, and a trimmed mean (.05 from both tails).
#' 
#' @details   Stone (1975) uses \eqn{r_n=\hat\sigma t_n} 
#' as the bandwidth of the Gaussian kernel, where \eqn{\hat\sigma} is the median of 
#' the \eqn{X_i-inth}'s.
#' For asymptotic efficiency of the estimators,
#' a) dn\eqn{\to\infty} and tn\eqn{\to 0} b)
#' \deqn{\frac{(dn)^2}{n^{1-\epsilon}(tn)^6}=O(1)}
#' for some \eqn{\epsilon>0}.  A rule of thumb plug-in estimate
#' for the optimal kernel width is 1.059\eqn{\hat\sigma n^{-1/5}} where \eqn{\hat\sigma}
#' is an estimator of the standard deviation, which does not satisfy the 
#' above relation. Stone(1975) takes inth to be the median in his simulations. The default choices of dn and tn 
#' are taken from Stone (1975), who 
#' considered a sample of size 40. Stone(1975)'s estimator uses Gaussian 
#' Kernels to estimate the unknown symmetric density. 
#' 
#' 
#' @return   An array of the same length as inth. Each element is the one-step
#'           estimator based on the corresponding element in inth.
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
#' @export 
giveth=function(x, inth, tn, dn)
{
  if(missing(inth)) inth <- init(x)
  if(missing(tn)) tn <- 0.60
  if(missing(dn)) dn <- 20
  
  r=length(inth)
  ans=numeric(r)
  for (i in 1:r)
    ans[i]=givethin2(x,inth[i], tn, dn)
  ans
}



#Gives stone's estimator using init(s) as initial estimators
#givethv=function(x) sapply(init(x),giveth,x=x)
#givethl=function(x) sapply(initl(x),giveth,x=x)
#givethr=function(x) sapply(initr(x),giveth,x=x)

