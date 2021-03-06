% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Beran_estimator.R
\name{score.coeff}
\alias{score.coeff}
\title{Fourier coefficients of the score in a location model}
\usage{
score.coeff(x, theta, indices, which)
}
\arguments{
\item{x}{An array of length n; represents the data from a location model whose score the Fourier coefficients
seek to estimate.}

\item{theta}{A small number, should be of order \eqn{O_p(n^{-1/2})}. The default
is \eqn{4n^{-1/2}}}

\item{indices}{An array of positive integers, for each integer j in "indices", 
Fourier coefficient  corresponding to the
basis function \eqn{t\mapsto\exp(i2\pi jt)} is computed.}

\item{which}{Optional. Takes value 1 or 2. If "which" is 1, only the real
part of the Fourier coefficient is computed. If "which" is 2,
only the imaginary part of the k-th coefficient is calculated. 
The default is to calculate both real and imaginary
parts.}
}
\value{
\itemize{
\item If "which=1", an array (with the same length as "indices"), give the real parts of the Fourier coefficients.
\item If "which=2, an array (with the same length as "indices"), give the  imaginary parts of the Fourier coefficients.
\item The default is to return a matrix whose first and second column give the
      real and imaginary parts of the Fourier coefficients, respectively. The number
      of rows equal the length of the array "indices". 
}
}
\description{
Estimates the Fourier coefficients of the scores in a location shift model,
with is the class of densities obtained by a location shift of
a fixed unknown density \eqn{f}. The score,
 given by 
\deqn{\phi'(x)=-f'\circ F^{-1}(t)/f\circ F^{-1}(x),}
does not depend on the unknown shift. The Fourier coefficients of the score
is estimated from the data, whose density is assumed to be a location
shift of the unknown density \eqn{f}. We use Beran (1974)'s nonparametric
estimator here.
}
\details{
\code{theta:}   theta do not need to depend on the range of the data because
                        the estimators depend only on the rank of the dataponts.  If \eqn{theta=z_n} which 
                        equals \eqn{z n^{-1/2}}, then the estimated coefficients are 
                        root-n consistent by Theorem 2.1 of Beran (1974).

\code{which:}   If it is known that the density \eqn{f} is symmetric,
                         then the Fourier coefficients are real. Therefore, there is
                         no need to calculate the imaginary part. Similarly, if the density 
                         is odd, then the Fourier coefficinets of the scores are purely  
                         imaginary. Otherwise, one generally requires both the real and the imaginary parts, 
                         and "which" should be left unspecified in those cases.
}
\examples{
x <- rnorm(100); score.coeff(x, length(x)^(-1/2), 1)
}
\references{
Beran, R. (1974). \emph{Asymptotically efficient adaptive 
rank estimates in location models}. Ann. Statist., 2, 63-74.
}
\author{
\href{https://connects.catalyst.harvard.edu/Profiles/display/Person/184207}{Nilanjana Laha}
(maintainer), \email{nlaha@hsph.harvard.edu}.
}
