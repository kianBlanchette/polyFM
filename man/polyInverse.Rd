% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extraFunctions.R
\name{polyInverse}
\alias{polyInverse}
\title{polyInverse}
\usage{
polyInverse(V, Vdot, maxdeg, alpha)
}
\arguments{
\item{V}{An N x K matrix of the Slepian sequences (or another set of tapers)}

\item{Vdot}{An N x K matrix of the derivatives of the Slepian sequences}

\item{maxdeg}{The highest degree polynomial at which to compute the associated polynomials}

\item{alpha}{A parameter of the Gegenbauer polynomials used as a starting point. The Gegenbauer polynomials are orthogonal on [-1,1] with the       respect to the weight function (1-x^2)^(alpha - 1/2)}
}
\value{
A The N x K matrix used to produce the polynomial inverse

Adot The N x K matrix used to produce the derivative of the polynomial inverse
}
\description{
Derives the matrices used to compute the inverse based on the associated polynomials
}
\details{
The associated polynomials and their derivatives are derived, then combined with the Slepians and their derivatives to produce the matrices used in computing the polynomial inverse.
}
\examples{
 N <- 1024
 NW <- 4
 K <- 2*NW - 1
 DW <- dpss(N,K,NW,returnEigenvalues=FALSE)
 Vdot <- dpssDeriv(DW,NW)

 P <- 3

 polyInverse(DW$v, Vdot, maxdeg = P, alpha = 0.75)
}
\references{
D. J. Thomson, "Polynomial phase demodulation in multitaper analysis," 2009 IEEE/SP 15th Workshop on Statistical Signal Processing, 2009, pp. 401-404, doi: 10.1109/SSP.2009.5278553.

D. J. Thomson, "Inverse-constrained projection filters," Proc. SPIE 4478, Wavelets: Applications in Signal and Image Processing IX, (5 December 2001); https://doi.org/10.1117/12.449708

Blanchette, Kian. Multitaper statistical tests for the detection of frequency-modulated signals. MSc. Thesis. Queen's University (Canada), 2020
}
