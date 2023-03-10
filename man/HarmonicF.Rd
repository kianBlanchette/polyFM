% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extraFunctions.R
\name{HarmonicF}
\alias{HarmonicF}
\title{HarmonicF}
\usage{
HarmonicF(yk,dw)
}
\arguments{
\item{yk}{The matrix of eigencoefficients. M x K where M = number of Fourier frequencies, K = number of tapers.}

\item{dw}{N x K matrix of the Slepian sequences. N = the length of the series.}
}
\value{
HF the vector containing the Harmonic F statistic

cmv the vector containing the complex mean values
}
\description{
Computes the Harmonic F test statistic
}
\details{
The Harmonic F statistic is computed at all the Fourier frequencies for which the eigencoefficients were computed.
}
\examples{
N <- 1024
nFFT <- 2^ceiling(log2(2*N))
NW <- 4
K <- 2*NW - 1
DW <- dpss(N,K,NW,returnEigenvalues=FALSE)

X <- rnorm(N)
xzp <- c(X, rep(0,nFFT-N))
V <- rbind(DW$v, matrix(0, nrow = nFFT - N, ncol = K))
yk <- mvfft(V * X)

HF <- HarmonicF(yk, DW$v)
}
\references{
Thomson, David J. "Spectrum estimation and harmonic analysis." Proceedings of the IEEE 70.9 (1982): 1055-1096.
}
