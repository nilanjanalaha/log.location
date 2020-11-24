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
  mlef=activeSetLogCon.mode(sortx,print=FALSE)
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
  mlef=activeSetLogCon.mode(sortx,print=FALSE)
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
 mlef <- activeSetLogCon.mode(x, mode=x[1])
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

# k is 1 for mean, 2 for median, 3 for trimmed mean, 4 folr logistic MLE
giveMLE <- function(x,k,q,leg,distname)
{
    x <- sort(x)
    inth=initl(x)[k]
  mlef=logConDens(x,smoothed=FALSE,print=FALSE)
  th2 <- giveLt(x,mlef,inth,t)
  th5 <- giveth5t(x,inth,q)
  
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
  mlef=activeSetLogCon.mode(sortx,print=FALSE)
  phip <- diff(mlef$phi)/diff(mlef$x)
  Fhi <- mlef$Fhat
  kn<- sortx[mlef$IsKnot==1]
  kn <- sort(c(-kn,kn))
  
  # what if mlef works on some x , less size than sortx
  if(length(mlef$x)<length(sortx))
  {
    v <- setdiff(sortx,mlef$x)
    print(table(sortx))
    if(length(v)>1)
      return(rep(-100,8))
    m <- length(phip)
    j <- findInterval(v,mlef$x)
    Fj <- mlef$Fhat[j]+mlef$fhat[j]/phip[j]*(exp(v-mlef$x[j])-1)
    Fh <- c(Fh[1:j],Fj,Fh[(j+1):(m+1)])
    phip <- c(phip[1:j],phip[j],phip[(j+1):m])
  }
  
  rankv <- rank(abs(x)) # getting the ranks of x-th
  
  #calculating the L4
  L4 <- x
  L4[sign(x)<0] <- phip[(rankv[sign(x)<0])]
  L4[sign(x)>0] <- -phip[rankv[sign(x)>0]]
  Fh <-x
  Fh[sign(x)<0] <- 1/2*(1-Fhi[rankv[sign(x)<0]+1])
  Fh[sign(x)>0] <- 1/2*Fhi[rankv[sign(x)>0]+1]+1/2
  
  
  # calculating th4 and obtaining the truncation quantiles
  up <- quantilesLogConDens(2*q-1, mlef)[2]
  low <- -up
  upl <- findup(up,x)
  vecr <- c(findlow(low,x):upl)
  th4=giveone_t(L4[vecr],inth,n)
  th4full = giveone_t(L4,inth,n) # The full one discrete information
  
  thv <- seq(x[1],x[length(x)],length.out=1000)
  Lvec <- sapply(thv,giveth4ml,x=inx)
  ml <- thv[which.max(Lvec)]
  
  
  #Plotting truncated bigger plot
  mytitle=paste(c("MLEpointfar_",distname,"_n_",length(x),"trunc.jpeg"),collapse = '')
  jpeg(filename = mytitle,width = 900, height = 600)
  
  dist <- c(leg,expression(hat(theta)^MLE))
  plot(thv,Lvec,xlab=expression(theta),ylab="Log-likelihood",type='l',col='red',xlim=range(inx))
  segments(kn,min(Lvec),kn,min(Lvec)+diff(range(Lvec))*.2,lty=2)

  tht <- c(inth,th2[1],th2[3],th4,th5[3],th5[1],ml)
  thf <- c(inth,th2[2],th2[4],th4full,th5[4],th5[2],ml)
  print(cbind(tht,thf))
  Lv1 <- sapply(tht,giveth4ml,x=inx)
  Lv2 <- sapply(thf,giveth4ml,x=inx)
  for (i in 1:length(dist))
  {
    points(tht[i],Lv1[i],col=(i+2),pch=16)
    if(i<length(dist))
    { abline(v=tht[i],col=i+2)} else {
      abline(v=tht[i],col=i+2,lty=3)
   }
  }
  abline(v=0,lty=2,col='gray')
  rug(x)
  legend("topright",legend= dist,lty=rep(1,length(dist)),lwd=rep(2.5,length(dist)),col=colv)
  dev.off()
  
  
  # Grids for close plot
  #Calculating I 
  Itrunc <- sum(L4[vecr]^2)/n
  # The grid
  evar <- 3/sqrt(Itrunc*n)
  thv <- seq(-evar+th4full,th4full+evar,length.out=500)
  Lvec <- sapply(thv,giveth4ml,x=inx)
  ml <- thv[which.max(Lvec)]
  
  # Closer plot
  mytitle=paste(c("MLEpoint_",distname,"_n_",length(x),"trunc.jpeg"),collapse = '')
  jpeg(filename = mytitle,width = 900, height = 600)
  
  dist <- c(leg,expression(hat(theta)^MLE))
  plot(thv,Lvec,xlab=expression(theta),ylab="Log-likelihood",type='l',col='red')
  segments(kn,min(Lvec),kn,min(Lvec)+diff(range(Lvec))*.2,lty=2)
  
  tht <- c(inth,th2[1],th2[3],th4,th5[3],th5[1],ml)
  thf <- c(inth,th2[2],th2[4],th4full,th5[4],th5[2],ml)
  print(cbind(tht,thf))
  Lv1 <- sapply(tht,giveth4ml,x=inx)
  Lv2 <- sapply(thf,giveth4ml,x=inx)
  for (i in 1:length(dist))
  {
    points(tht[i],Lv1[i],col=(i+2),pch=16)
    if(i<length(dist))
    { abline(v=tht[i],col=i+2)} else {
      abline(v=tht[i],col=i+2,lty=3)
    }
  }
  abline(v=0,lty=2,col='gray')
  rug(x)
  legend("topright",legend= dist,lty=c(rep(1,length(dist)-1),3),lwd=rep(2,length(dist)),col=c(1:length(dist))+2,cex=1.2)
  dev.off()
  
  
  #Plotting full
  mytitle=paste(c("MLEpoint_",distname,"_n_",length(x),"full.jpeg"),collapse = '')
 jpeg(filename=mytitle,width = 900, height = 600)
  plot(thv,Lvec,xlab=expression(theta),ylab="Log-likelihood",type='l',col='red')
  for (i in 1:length(dist))
  {
    points(thf[i],Lv2[i],col=(i+2),pch=16)
    if(i<length(dist))
    { abline(v=thf[i],col=i+2)} else {
      abline(v=thf[i],col=i+2,lty=3)
    }
  }
  abline(v=0,lty=2,col='gray')
  rug(x)
  legend("topright",legend= dist,lty=c(rep(1,length(dist)-1),3),lwd=rep(2,length(dist)),col=c(1:length(dist))+2,cex=1.2)
 dev.off()
}

#k is 1 for mean, 2 for median, 3 for trimmed mean
comparef <- function(x,k,dist)
{
  
  x <- sort(x)
  inth=initl(x)[k]
  mlef <- logConDens(x,smoothed=FALSE)
  kn <- x[mlef$IsKnot==1]
  x1 <- sort(c(-kn,kn)) #the knots for plotting
  
  f2 <- giveLtf(x,mlef,inth)
  f3 <- giveth3tf(x,mlef,inth)
  f4 <- giveth4tf(x,inth)
  f5 <- giveth5tf(x,mlef,inth)
  flist <- list(f2,f3,f4,f5)
  low <- -max(abs(x))-.1
  up <- -low
  vec <- seq(low,up,.01)
  ylen <- c(min(log(f2[,2]))-1,max(log(f4[,2]))+.2)
  #Plot#############################################################################################
  colv <- c('blue','red','green','orange','black')
  colm <- colv[-1]
  mytitle=paste(c("MLE_",dist,"_n_",length(x),".jpeg"),collapse = '')
    
  jpeg(file=mytitle, width = 900, height = 600)
  plot(vec,log(dlogis(vec)),ylab="",xlab="x",type='l',col='blue',ylim=ylen,lwd=2)
  rug(x)
  segments(x1,ylen[1],x1,ylen[1]+diff(ylen)*.2,lty=2)
   for (i in 1:4)
  {
    fv <- sort_line(flist[[i]])
    lines(fv[,1],fv[,2],col=colm[i],lwd=2)  
   
  }
  
  dist <-c(expression("true log-density "*phi[0]),"symmetrized estimator","smoothed symmetrized estimator","symmetric estimator","symmetrized: geometric mean estimator")
  legend("center",legend= dist,lty=rep(1,length(dist)),lwd=rep(2.5,length(dist)),col=colv)
  dev.off()
  
  
}

comparec <- function(x,k,dist)
{
  
  x <- sort(x)
  inth=initl(x)[k]
  mlef <- logConDens(x,smoothed=FALSE)
  kn <- x[mlef$IsKnot==1]
  x1 <- sort(c(-kn,kn)) #the knots for plotting
  
  f2 <- giveLtf(x,mlef,inth)
  f3 <- giveth3tf(x,mlef,inth)
  f4 <- giveth4tf(x,inth)
  f5 <- giveth5tf(x,mlef,inth)
  flist <- list(f2,f3,f4,f5)
  low <- -max(abs(x))-.1
  up <- -low
  vec <- seq(low,up,.01)
  ylen <- c(min(log(f2[,2]))-1,max(log(f4[,2]))+.2)
  #Plot#############################################################################################
  colv <- c('blue','red','green','orange','black')
  colm <- colv[-1]
  mytitle=paste(c("cMLE_",dist,"_n_",length(x),".jpeg"),collapse = '')
  ylen <- c(min(log(f5[,2]))-.2,max(log(f4[,2]))+.2)
  xlen <- range(f5[,1])
  
  pdf(file=mytitle, dpi=700)
  plot(vec,log(dlaplace(vec)),ylab="",xlab="x",type='l',col='blue',ylim=ylen,xlim=xlen,lwd=2)
  rug(x)
  segments(x1,ylen[1],x1,ylen[1]+diff(ylen)*.2,lty=2)
  for (i in 1:4)
  {
    fv <- sort_line(flist[[i]])
    lines(fv[,1],fv[,2],col=colm[i],lwd=2)  
    
  }
  
  dist <-c(expression("true log-density "*phi[0]),"symmetrized estimator","smoothed symmetrized estimator","symmetric estimator","symmetrized: geometric mean estimator")
  legend("topleft",legend= dist,lty=rep(1,length(dist)),lwd=rep(2.5,length(dist)),col=colv)
  dev.off()
  
  
}

sort_line <- function(m)
{
m<- m[order(m[,1]),]
m[,2] <- log(m[,2])
m  
}