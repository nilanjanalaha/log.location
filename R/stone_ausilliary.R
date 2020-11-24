
################## STONE ESTIMATOR ####################################

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
Ah=function(r,c,x,inth) integral(Ahin,xmin=-c,xmax=c,c=c,rex=r,x=x,inth=inth)
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
givethin2=function(x,inth)
{
  x=sort(x)
  
  # Alert : the tuning parameters
  r=s(x,inth)*.6
  cc=s(x,inth)*20
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
giveth=function(x,inthv)
{
  r=length(inthv)
  ans=numeric(r)
  for (i in 1:r)
    ans[i]=givethin2(x,inthv[i])
  ans
}

init <- function(x) c(mean(x), median(x), mean(x, trim=0.05))

#Gives stone's estimator using init(s) as initial estimators
givethv=function(x) sapply(init(x),giveth,x=x)



