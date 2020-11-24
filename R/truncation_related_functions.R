q

#Calculates distribution function from a given density function
#Can not be used directly into find_quantile beacuse it has an additional argument
givedist <- function(x,density){
  integrate(f=density,
            lower=-Inf,
            upper=x)$value
  }

#calculates quantile function from a given distribution function

objective <- function(x, quantile, dist){
  (dist(x) - quantile)^2
}
find_quantile <- function(dist, quantile, a,b){
  result = nlminb(start=(a+b)/2, objective=objective,
                  quantile = quantile,
                  dist = dist, lower=a, upper=b)$par
  return(result)
}


####################################################################################################################################
# Finding quantile for symmetrized log-concave density (1st)
# q>o :IMPORTANT##################################################################
#quantilesLogConDens(q, mlef)[2]

# non-symmetrized, mlef is a result of call to logConDens

findq1 <- function(q, mlef, bth)#q is for quantile
{
  
  sym_Fn <- function(x) (evaluateLogConDens(bth+x,mlef,which=3)[4]+1-evaluateLogConDens(bth-x,mlef,which=3)[4])/2
  kn <- sort(abs(mlef$xn-bth))
  Fn <- sapply(kn,sym_Fn)
    j <- findup(q,Fn)
  find_quantile(sym_Fn,q,kn[j],kn[j+1])
}

# Finding quantile for smoothed symmetrized log-concave density (2nd)
#mlef is smoothed dlc object
findq2 <- function(q, mlef,bth)#q is for quantile
{
  
  sym_Fn <- function(x) (evaluateLogConDens(bth+x,mlef,which=5)[6]+1-evaluateLogConDens(bth-x,mlef,which=5)[6])/2
  kn <- sort(abs(mlef$xn-bth))
  Fn <- sapply(kn,sym_Fn)
  if(Fn[length(Fn)]>q)
  {
  j <- findInterval(q,Fn)
 return( find_quantile(sym_Fn,q,kn[j],kn[j+1]))
  }
  else
  {
    return(-1) # Indicator that q is attained outside the range of x
  }
}

# Finding quantile for  symmetric log-concave density (3rd))
#mlef is logConDens.mode object
#quantilesLogConDens(2*q-1, mlef)[2]#here q must be positive

# Finding quantile for smoothed  symmetric log-concave density (4th)

# Finding quantile for  GM (5th)

findq5 <- function(q,mlef,bth)
{
  xt <- mlef$xn-bth
  kn <- sort(c(0,abs(xt)))
  phi.c <- function(x) (evaluateLogConDens(bth+x,mlef,which=1)[2]/2+evaluateLogConDens(bth-x,mlef,which=1)[2]/2)
  l <- sapply(kn,phi.c)
  kn <- kn[-which(l==-Inf)]
  l=l[-which(l==-Inf)]
  f <- exp(l)
  
  deriv <- diff(l)/diff(kn)
 zero <- which(deriv==0)
 dist.part <- deriv
 posi <- 1:length(deriv)
 if(length(zero>0))
 {
   zerop <- max(zero)
   dist.part[1:zerop] <- f[1]*(diff(kn)[1:zerop])
   posi <- setdiff(posi,1:zerop)
 }
  
  
  dist.part[posi] <- diff(f)[posi]/deriv[posi]
  Fn <- cumsum(dist.part)
  C <- 2*max(Fn)
  Fn <- 1/2+ Fn/C
  
  j <- findup(q,Fn)
  print(j)
  sym_Fn <- function(x) Fn[j]+f[j]*(exp(deriv[j]*(x-kn[j]))-1)/deriv[j]
  find_quantile(sym_Fn,q,kn[j],kn[j+1])
  
}

giveth5tf <- function(x,mlef,inth)
{
  xt <- mlef$xn-bth
  kn <- sort(c(0,abs(xt)))
  phi.c <- function(x) (evaluateLogConDens(bth+x,mlef,which=1)[2]/2+evaluateLogConDens(bth-x,mlef,which=1)[2]/2)
  l <- sapply(kn,phi.c)
  kn <- kn[-which(l==-Inf)]
  l=l[-which(l==-Inf)]
  f <- exp(l)
  
  deriv <- diff(l)/diff(kn)
  zero <- which(deriv==0)
  dist.part <- deriv
  posi <- 1:length(deriv)
  if(length(zero>0))
  {
    zerop <- max(zero)
    dist.part[1:zerop] <- f[1]*(diff(kn)[1:zerop])
    posi <- setdiff(posi,1:zerop)
  }
  
  
  dist.part[posi] <- diff(f)[posi]/deriv[posi]
  Fn <- cumsum(dist.part)
  C <- 2*max(Fn)
  
  kn2 <- c(sort(-kn),kn[-1])
  f <- c(sort(f),f[-1])
  cbind(kn2+inth,f/C)
}











###################################################################################################
