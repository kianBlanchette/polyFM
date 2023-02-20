
.is.dpss <- function(obj){
  class(obj) == "dpss"
}

.is.deriv <- function(obj){
  class(obj) == "dpssDeriv"
}

.is.ap <- function(obj){
  class(obj) == "ass.poly"
}

#' @title polyInverse
#' @description Derives the matrices used to compute the inverse based on the associated polynomials
#' @usage  polyInverse(V, Vdot, maxdeg, alpha)
#' @param V An N x K matrix of the Slepian sequences (or another set of tapers)
#' @param Vdot An N x K matrix of the derivatives of the Slepian sequences
#' @param maxdeg The highest degree polynomial at which to compute the associated polynomials
#' @param alpha A parameter of the Gegenbauer polynomials used as a starting point. The Gegenbauer polynomials are orthogonal on [-1,1] with the       respect to the weight function (1-x^2)^(alpha - 1/2)
#' @details The associated polynomials and their derivatives are derived, then combined with the Slepians and their derivatives to produce the matrices used in computing the polynomial inverse.
#' @returns A The N x K matrix used to produce the polynomial inverse
#' @returns Adot The N x K matrix used to produce the derivative of the polynomial inverse
#' @references D. J. Thomson, "Polynomial phase demodulation in multitaper analysis," 2009 IEEE/SP 15th Workshop on Statistical Signal Processing, 2009, pp. 401-404, doi: 10.1109/SSP.2009.5278553.
#' @references D. J. Thomson, "Inverse-constrained projection filters," Proc. SPIE 4478, Wavelets: Applications in Signal and Image Processing IX, (5 December 2001); https://doi.org/10.1117/12.449708
#' @references Blanchette, Kian. Multitaper statistical tests for the detection of frequency-modulated signals. MSc. Thesis. Queen's University (Canada), 2020
#' @export
#' @examples
#'  N <- 1024
#'  NW <- 4
#'  K <- 2*NW - 1
#'  DW <- dpss(N,K,NW,returnEigenvalues=FALSE)
#'  Vdot <- dpssDeriv(DW,NW)
#'
#'  P <- 3
#'
#'  polyInverse(DW$v, Vdot, maxdeg = P, alpha = 0.75)

polyInverse <- function(V, Vdot, maxdeg = 0, alpha = 0.75){

  AP <- apDeriv(V=V, maxdeg=maxdeg, alpha=alpha)

  H <- AP$H
  G <- AP$G
  D <- AP$D
  K <- nrow(H)

  Hproj <- diag(x=1, nrow = K) - tcrossprod(H)
  Gproj <- tcrossprod(G, H)
  A <- tcrossprod(V, Hproj) + Gproj

  Dproj <- tcrossprod(D, H)
  Adot <- tcrossprod(Vdot, Hproj) + Dproj

  out = list(A=A, Adot=Adot)
  return(out)
}

#' @title apDeriv
#' @description Derives the associated polynomials and their derivatives.
#' @usage  apDeriv(V, maxdeg, alpha)
#' @param V An N x K matrix of the Slepian sequences (or another set of tapers)
#' @param maxdeg The highest degree polynomial at which to compute the associated polynomials
#' @param alpha A parameter of the Gegenbauer polynomials used as a starting point. The Gegenbauer polynomials are orthogonal on [-1,1] with the       respect to the weight function (1-x^2)^(alpha - 1/2)
#' @details The associated polynomials are created through a Gram-Schmidt type process involving the Slepian sequences and the Gegenbauer polynomials.
#' @returns H The K x P matrix of eigencoefficients of the associated polynomials
#' @returns G The N x P matrix of associated polynomials
#' @returns D The N x P matrix of the derivatives of the associated polynomials
#' @returns Hn The column sums of squares of R
#' @references D. J. Thomson, "Polynomial phase demodulation in multitaper analysis," 2009 IEEE/SP 15th Workshop on Statistical Signal Processing, 2009, pp. 401-404, doi: 10.1109/SSP.2009.5278553.
#' @references D. J. Thomson, "Inverse-constrained projection filters," Proc. SPIE 4478, Wavelets: Applications in Signal and Image Processing IX, (5 December 2001); https://doi.org/10.1117/12.449708
#' @export
#' @examples
#'  N <- 1024
#'  NW <- 4
#'  K <- 2*NW - 1
#'  DW <- dpss(N,K,NW,returnEigenvalues=FALSE)
#'
#'  P <- 3
#'
#'  apDeriv(DW$v,P,alpha=0.75)

apDeriv <- function(V, maxdeg, alpha=0.75) {

  # Sanity checks
  stopifnot(is.matrix(V), is.numeric(maxdeg), maxdeg>=0)
  N <- length(V[, 1])
  K <- length(V[1, ])
  P <- maxdeg + 1
  timeArr <- 1:N

  G <- D <- matrix(data=0, nrow=N, ncol=P)
  H <- matrix(data=0, nrow=K, ncol=P)


  # Setup centered time index
  midTime <- (1+N) / 2
  tscl <- 2/(N-1)
  timeArrC <- (timeArr - midTime) * tscl

  # Start with Gegenbauer polynomials; convergence is faster
  #alpha <- 0.75
  G[, 1] <- 1.0
  if(maxdeg > 0) {
    G[, 2] <- 2 * alpha * timeArrC
    D[, 2] <- tscl
    if(maxdeg > 1) {
      for(j in 2:maxdeg) {
        A1 <- 2 * ( (j-1) + alpha ) / j
        A2 <- ( (j-2) + 2 * alpha ) / j
        G[, (j+1)] <- A1 * timeArrC * G[, j] - A2 * G[, (j-1)]
        D[, (j+1)] <- tscl * j * G[, j]
      } # end of loop on higher orders
    } # end of maxdeg > 1
  } # end of maxdeg > 0

  # Inner Products of R and V
  for(L in 1:P) {
    Kmin <- ( (L-1) %% 2 ) + 1
    for(k in seq(Kmin, K, 2)) {  # loop on non-zero Slepians
      H[k, L] <- t(V[, k]) %*% G[, L]
    }
  }


  # Degree 0, 1 (manual) -- L = degree+1
  for(L in 1:min(2,P)) {
    scl <- 1 / sqrt( sum(H[, L]^2) )
    H[, L] <- H[, L] * scl # orthonormalize
    G[, L] <- G[, L] * scl
    D[, L] <- D[, L] * scl
  }



  # loop on higher degrees, applying Gram-Schmidt only on similar
  # parity functions (as even/odd are already orthogonal in U)
  if( P > 2 ) {
    for(L in 3:P) {
      if(L %% 2 == 0) {
        Kmin <- 2
      } else {
        Kmin <- 1
      }
      for(j in seq(Kmin, L-1, 2)) {
        scl <- sum( H[, L] * H[, j] )
        H[, L] <- H[, L] - scl * H[, j] # Gram-Schmidt
        G[, L] <- G[, L] - scl * G[, j]
        D[, L] <- D[, L] - scl * D[, j]
      }
      scl <- 1 / sqrt(sum(H[, L]^2))
      H[, L] <- H[, L] * scl  # orthonormalize
      G[, L] <- G[, L] * scl
      D[, L] <- D[, L] * scl
    }
  }


  Hn <- colSums(G^2)
  ap <- list(H=H,G=G,D=D,Hn=Hn)
  class(ap) <- "ass.poly"
  return(ap)
}

#' @title HarmonicF
#' @description Computes the Harmonic F test statistic
#' @usage HarmonicF(yk,dw)
#' @param yk The matrix of eigencoefficients. M x K where M = number of Fourier frequencies, K = number of tapers.
#' @param dw N x K matrix of the Slepian sequences. N = the length of the series.
#' @details The Harmonic F statistic is computed at all the Fourier frequencies for which the eigencoefficients were computed.
#' @returns HF the vector containing the Harmonic F statistic
#' @returns cmv the vector containing the complex mean values
#' @references Thomson, David J. "Spectrum estimation and harmonic analysis." Proceedings of the IEEE 70.9 (1982): 1055-1096.
#' @export
#' @examples
#' N <- 1024
#' nFFT <- 2^ceiling(log2(2*N))
#' NW <- 4
#' K <- 2*NW - 1
#' DW <- dpss(N,K,NW,returnEigenvalues=FALSE)
#'
#' X <- rnorm(N)
#' xzp <- c(X, rep(0,nFFT-N))
#' V <- rbind(DW$v, matrix(0, nrow = nFFT - N, ncol = K))
#' yk <- mvfft(V * X)
#'
#' HF <- HarmonicF(yk, DW$v)

HarmonicF <- function(yk, dw){

  k <- ncol(yk)
  Ukz <- colSums(dw)
  ssqUkz <- sum(Ukz^2)
  cmv <- (yk %*% Ukz) / ssqUkz
  ssqave <- ssqUkz * Mod(cmv)^2
  Ukz <- as.matrix(Ukz)

  ssqres <- apply( Mod(yk - (cmv %*% t(Ukz)))^2, MARGIN = 1, FUN = sum)
  HF <- (k-1) * ssqave / ssqres
  class(HF) <- "Ftest"

  out <- list(HF = HF, cmv = cmv)

  return(out)
}

.modF <- function(mxdeg, nord, FPcoef, Fcoef, nfreqs){
  output <- .Fortran("modF",
                     mxdeg = as.integer(mxdeg),
                     nord = as.integer(nord),
                     FPcoef = as.double(FPcoef),
                     Fcoef = as.double(Fcoef),
                     Fp = double( mxdeg*nfreqs ),
                     nfreqs = as.integer(nfreqs)
  )
  Fp <- matrix(data = output$Fp, nrow = mxdeg, ncol = nfreqs)
}

.jkFreq <- function(Phi, H, freq){
  K <- nrow(H)
  P <- ncol(H)
  stopifnot(K > P)

  mF <- matrix(0, nrow = P, ncol = length(freq))
  fkp <- matrix(0, nrow = K, ncol = P)
  for(k in 1:K){
    ssq1 <- colSums( Phi[-k,]^2 )
    PhiP <- crossprod(H[-k,], Phi[-k,])
    ssq2 <- 0
    for(p in 1:P){
      ssq2 <- ssq2 + PhiP[p,]^2
      mF[p,] <- (K-1-p) * PhiP[p,]^2 / (ssq1 - ssq2)
      fkp[k,p] <- freq[which.max(mF[p,])]
    }
  }
  jkMean <- apply(fkp, MARGIN = 2, FUN = function(x) mean(x))
  jkVar <- numeric(P)
  for(p in 1:P){
    jkVar[p] <- (1-1/K) * sum( (fkp[,p] - jkMean[p])^2 )
  }

  jk <- list(jkMean = jkMean, jkVar = jkVar)
  return(jk)
}

.jackknifeModulatedF <- function(yk, derivIN = NULL, apIN = NULL, dpssIN = NULL, freq){

  #speed version of this function. everything must be passed in.
  #just computes modified F3. note that H matrix is missing first column,
  #the zero degree column.
  # also does jackknife estimates of mean and variance of frequency
  K <- ncol(yk)
  V <- dpssIN$v
  H <- apIN$H[,-1, drop = FALSE]
  P <- ncol(H)

  stopifnot(K > P)
  nfreqs <- nrow(yk)
  #############################################################################

  phi <- IFcompute(yk, V, derivIN)

  Phi <- crossprod(V,phi)
  PhiP <- crossprod(H,Phi)

  mF <- .modF(P,K,PhiP,Phi,nfreqs)

  jk <- .jkFreq(Phi,H,freq)

  ModF <- list(mF = mF, jk = jk)

  return(ModF)

}
