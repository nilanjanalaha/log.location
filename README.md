
<!-- README.md is generated from README.Rmd. Please edit that file -->

# log.location

<!-- badges: start -->

<!-- badges: end -->

The goal of log.location is to …

## Installation

You can install the released version of log.location from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("log.location")
```

And the development version from [GitHub](https://github.com/) with:

``` r
#Uninstalling previous versions, if installed
devtools::install_github("nilanjanalaha/log.location", force=TRUE)
#> Downloading GitHub repo nilanjanalaha/log.location@master
#> 
#>      checking for file ‘/private/var/folders/vb/498bk3fd1x38tk9jpy4760c00000gn/T/RtmpTgMXaA/remotescf087619cc5d/nilanjanalaha-log.location-776ae33/DESCRIPTION’ ...  ✓  checking for file ‘/private/var/folders/vb/498bk3fd1x38tk9jpy4760c00000gn/T/RtmpTgMXaA/remotescf087619cc5d/nilanjanalaha-log.location-776ae33/DESCRIPTION’ (458ms)
#>   ─  preparing ‘log.location’:
#>      checking DESCRIPTION meta-information ...  ✓  checking DESCRIPTION meta-information
#>   ─  excluding invalid files
#>      Subdirectory 'R' contains invalid file names:
#>      ‘n100N500.RData’
#>   ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#> ─  building ‘log.location_0.0.0.9000.tar.gz’
#>      
#> 
```

Let us simulate some Toy data

``` r
library(log.location)
x <- rgamma(1000, 1, 1)
```

## Beran (1974) ’s estimator : the fourier coefficients

The following function computes the Fourier coefficients of the scores.
Here theta is a tuning parameter which controls the precision of the
estimation and “indices” is the vector of indexes corresponding to the
basis functions.

``` r
score.coeff(x, theta=0.01, indices= c(1,2,4))
#>       real.part    im.part
#> [1,] -1.0564957 -0.0307393
#> [2,] -1.1536626 -0.1317294
#> [3,] -0.9337857 -0.1525421
```

If we want only real parts or only imaginary parts, we can indicate so
by setting “which” to be one and two, respectively.

``` r
# To output only the real part of the cefficients
score.coeff(x, theta=0.01, indices= c(1,2,4), which=1)
#> [1] -1.0564957 -1.1536626 -0.9337857
```

We plot the real parts of the estimated Fourier coefficients as a
function of theta.
<img src="man/figures/README-plotgamma-1.png" width="100%" />

We repeat the same procedure for normal data. The true score for the
normal location odel is \(\Phi^{-1}(x)\) where \(\Phi\) is the standard
normal distribution function. We use the stats functon integrate to
calculate the imaginary part of its first Fourier coefficient (k=1). The
real part will be zero because standard normal density is symmetric
about 0.

``` r
f <- function(x) qnorm(x)*sin(2*pi*x)
val <- stats::integrate(f, 0, 1)$value
val
#> [1] -0.4638467
```

Overlaying the true value on the plot:

``` r
set.seed(5)
x <- rnorm(1000)
grid.theta <- seq(0.01, 0.5, length.out= 500)
reals <- sapply(grid.theta, score.coeff, x=x, indices=1, which=2)
plot(grid.theta, reals, xlab="theta", ylab="coefficients", type='l')
abline(h=val, col='blue')
```

<img src="man/figures/README-overlay value.-1.png" width="100%" />

The suggested value of theta in Beran (1974) is of order \(n^{-1/2}\).
We set \(\theta=cn^{-1/2}\). Let us plot the real part of the first
Fourier coefficient as a function of \(c\).

<img src="man/figures/README-overlay value 2-1.png" width="100%" />

We observe that for Gaussian samples, \(cn^{-1/2}\) is not a very good
choice for \(\theta\). We consider k=3 now.

``` r
set.seed(8)
grid.theta <- seq(1, 100, by= 1)
reals <- sapply(grid.theta/sqrt(1000), score.coeff, x=x, indices=3, which=2)
plot(grid.theta, reals, xlab="c", ylab="coefficients", type='l')

f <- function(x) qnorm(x)*sin(2*3*pi*x)
val <- stats::integrate(f, 0, 1)$value
print(val)
#> [1] -0.2052047
abline(h=val, col='blue')
```

<img src="man/figures/README-overlay value 3-1.png" width="100%" />

We run the same calculations by changing the seed.
<img src="man/figures/README-overlay value 4-1.png" width="100%" />

Note that the previous c no longer works.

``` r
set.seed(1)
grid.theta <- seq(1, 100, by= 1)
reals <- sapply(grid.theta/sqrt(1000), score.coeff, x=x, indices=3, which=2)
plot(grid.theta, reals, xlab="c", ylab="coefficients", type='l')

f <- function(x) qnorm(x)*sin(2*3*pi*x)
val <- stats::integrate(f, 0, 1)$value
print(val)
#> [1] -0.2052047
abline(h=val, col='blue')
```

<img src="man/figures/README-overlay value 5-1.png" width="100%" /> It
seems that the estimates are more stable with larger values of k.

## Beran (1974) ’s estimator

Beran (1974)’s estimator is given by the function beran.est. This
function also computes \((1-\alpha)\%\) confidence intervals. If
\(\alpha\) is missing, 0.95\(%\) confidence intervals are calculated.
One important tuning parameter is the number of basis functions, which
is given by M. Below we plot the estimates given by Beran’s estimator
versus M for a standard logistic sample, whose location parameter
\(\theta\) is 0. The true value of \(\theta\) is given by the blue line.
The confidence intervals are in red.

``` r
set.seed(1)
grid.M <- seq(1, 20, by= 1)
#The data
x=rlogis(100)
est <- CI.lb <- CI.ub <- numeric(length(grid.M))

#The estimates
for(M in grid.M)
{
  est[M] <- beran.est(x, M=M)$estimate
  CI.lb[M] <- beran.est(x, M=M)$CI[,1]
  CI.ub[M] <- beran.est(x, M=M)$CI[,2]
}
  
plot(grid.M, est, xlab="M", ylab="The estimates", type='l', ylim=c(-2,2))

#Confidence intervals are in red
lines(grid.M, CI.lb, col='red')
lines(grid.M, CI.ub, col='red')
abline(h=0, col='blue')
```

<img src="man/figures/README-beran estimate-1.png" width="100%" />

We observe that the confidence intervals are too large for n=100. We try
n=1000, which still has high confidence intervals.
<img src="man/figures/README-beran estimate 1000-1.png" width="100%" />

## Stone’s estimator

Stone (1975)’s estimators also give the estimator of the location
\(\theta\). The function to use is **giveth**. You have to supply an
initial estimator of \(\theta\), which is **inth**. The tuning
parameters are dn (truncation) and tn (smoothness). If not sapplied byy
the user, default values are used.

``` r
x <- rnorm(100)
giveth(x, inth = mean(x))
#> $estimate
#> [1] -0.04902123
#> 
#> $CI
#> [1] -0.11892722  0.02088475
giveth(x, inth = mean(x), dn=30, tn=0.50)
#> $estimate
#> [1] -0.04641573
#> 
#> $CI
#> [1] -0.11392271  0.02109125
```

\#Estimators in Laha et al. (2020)

## Partial MLE estimator

Use the function **p.mle**. The parameter q can be used to give the
truncation level. The default is 0. Please load the package `logcondens"
before using this. Although log.location adds logcondens as a
dependency, for some reason, the funvtion`PreProcess does not get
loaded".

``` r
x <- rnorm(100)
library(logcondens)
p.mle(x, init=mean(x))
#> $estimate
#>   q        est
#> 1 0 -0.2460157
#> 
#> $CI
#>   q         lb          ub
#> 1 0 -0.4663034 -0.02572806
p.mle(x, q=c(0, 0.01, 0.02))
#> $estimate
#>      q        est
#> 1 0.00 -0.2460157
#> 2 0.01 -0.2492603
#> 3 0.02 -0.2528789
#> 
#> $CI
#>      q         lb          ub
#> 1 0.00 -0.4663034 -0.02572806
#> 2 0.01 -0.4748238 -0.02369670
#> 3 0.02 -0.4841165 -0.02164123
```

## MLE estimator

The function **MLE\_location** gives the MLE under the log-concave
location model of laha et al. (2020). Please load the package
`logcondens" before using this. Although log.location adds logcondens as
a dependency, for some reason, the funvtion`PreProcess does not get
loaded".

``` r
library(logcondens)
MLE_location(x)
#> $estimate
#> [1] -0.1966886
#> 
#> $CI
#> [1] -0.41164290  0.01826566
```

## The smoothed symmetrized estimator

The function **s.sym** calculates this estimator. Pass the initial
estimator as inth and use q to set the truncation level.

``` r
library(logcondens)
s.sym(x)
#> $estimate
#>   q       est
#> 1 0 -0.199468
#> 
#> $CI
#>   q         lb           ub
#> 1 0 -0.3967695 -0.002166426
s.sym(x, q=c(0.01, 0.02))
#> $estimate
#>      q        est
#> 1 0.01 -0.2248547
#> 2 0.02 -0.2295700
#> 
#> $CI
#>      q         lb         ub
#> 1 0.01 -0.4730244 0.02331501
#> 2 0.02 -0.4875726 0.02843262
```
