#' @title dpssap
#' @description Computes a set of polynomials associated with the Slepian sequences. The polynomials themselves are not orthogonal, but when expanded in their eigencoeffiecent form U = V^T R, they are.
#' @usage dpssap(V,maxdeg,alpha=0.75)
#' @param V An NxK matrix of Slepian sequences as created in the dpss function of the multitaper package.
#' @param maxdeg The highest degree of polynomial to be generated.
#' @param alpha A parameter of the Gegenbauer polynomials used as a starting point. The Gegenbauer polynomials are orthogonal on [-1,1] with the       respect to the weight function (1-x^2)^(alpha - 1/2)
#' @details The associated polynomials are created through a Gram-Schmidt type process involving the Slepian sequences and the Gegenbauer polynomials.
#' @returns U The K x P matrix of eigencoefficients of the associated polynomials
#' @returns R The N x P matrix of associated polynomials
#' @returns Hn The column sums of squares of R
#' @references D. J. Thomson, "Polynomial phase demodulation in multitaper analysis," 2009 IEEE/SP 15th Workshop on Statistical Signal Processing, 2009, pp. 401-404, doi: 10.1109/SSP.2009.5278553.
#' @references D. J. Thomson, "Inverse-constrained projection filters," Proc. SPIE 4478, Wavelets: Applications in Signal and Image Processing IX, (5 December 2001); https://doi.org/10.1117/12.449708
#' @export
#' @examples
#'  N <- 1000
#'  NW <- 4
#'  K <- 2*NW - 1
#'  DW <- dpss(N,K,NW,returnEigenvalues=FALSE)
#'
#'  P <- 3
#'
#'  dpssap(DW$v,P,alpha=0.75)

dpssap <- function(V, maxdeg, alpha=0.75) {

  # Sanity checks
  stopifnot(is.matrix(V), is.numeric(maxdeg), maxdeg>=0)
  N <- length(V[, 1])
  K <- length(V[1, ])
  P <- maxdeg + 1
  timeArr <- 1:N

  R <- matrix(data=0, nrow=N, ncol=P)
  U <- matrix(data=0, nrow=K, ncol=P)

  # Setup centered time index
  midTime <- (1+N) / 2
  scl <- 2/(N-1)
  timeArrC <- (timeArr - midTime) * scl

  # Start with Gegenbauer polynomials; convergence is faster
  R[, 1] <- 1.0
  if(maxdeg > 0) {
    R[, 2] <- 2 * alpha * timeArrC
    if(maxdeg > 1) {
      for(j in 2:maxdeg) {
        A1 <- 2 * ( (j-1) + alpha ) / j
        A2 <- ( (j-2) + 2 * alpha ) / j

        R[, (j+1)] <- A1 * timeArrC * R[, j] - A2 * R[, (j-1)]
      } # end of loop on higher orders
    } # end of maxdeg > 1
  } # end of maxdeg > 0

  # Inner Products of R and V
  for(L in 1:P) {
    Kmin <- ( (L-1) %% 2 ) + 1
    for(k in seq(Kmin, K, 2)) {  # loop on non-zero Slepians
      U[k, L] <- t(V[, k]) %*% R[, L]
    }
  }

  # Degree 0, 1 (manual) -- L = degree+1
  for(L in 1:min(2,P)) {
    scl <- 1 / sqrt( sum(U[, L]^2) )
    U[, L] <- U[, L] * scl # orthonormalize
    R[, L] <- R[, L] * scl
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
        scl <- sum( U[, L] * U[, j] )
        U[, L] <- U[, L] - scl * U[, j] # Gram-Schmidt
        R[, L] <- R[, L] - scl * R[, j]
      }
      scl <- 1 / sqrt(sum(U[, L]^2))
      U[, L] <- U[, L] * scl  # orthonormalize
      R[, L] <- R[, L] * scl
    }
  }

  Hn <- colSums(R^2)
  ap <- list(U=U,R=R,Hn=Hn)
  class(ap) <- "ass.poly"
  return(ap)
}

#' @title dpssDeriv
#' @description Computes "derivatives" of the Slepian sequences at the same points the Slepians are defined at.
#' @usage  dpssDeriv(DW,NW)
#' @param DW The DPSS object as output by the dpss function of the multitaper package. This is essentially a list containing both the matrix of the first K Slepian sequences of length N and the associated eigenvalues. Both are needed here.
#' @param NW The multitaper time-bandwidth product (N = time = number of points * sampling period, W = bandwidth). This should have been selected already when computing the Slepian sequences. Typically, K = 2NW - 1, so a common choice of NW is (K+1)/2, where K is the number of tapers. The NW chosen here MUST match the NW used to generate the DW object.
#' @details The derivatives of the Slepians are defined through a matrix equation involving a skew-symmetrix Toeplitz matrix and the Slepian sequences scaled by their eigenvalues. Thus, we only need one row of the Toeplitz matrix and we can use the Fast Fourier Transform and the Convolution Theorem to generate the derivatives as opposed to computing them directly through a lengthy matrix multiplication.
#' @returns efnDeriv An N x K matrix whose columns contain the derivatives of the first K Slepian sequences.
#' @references  Blanchette, Kian. Multitaper statistical tests for the detection of frequency-modulated signals. MSc. Thesis. Queen's University (Canada), 2020.
#' @export
#' @examples
#' N <- 1024
#' NW <- 4
#' K <- 2*NW - 1
#'
#' DW <- dpss(N,K,NW,returnEigenvalues = TRUE)
#'
#' Vdot <- dpssDeriv(DW,NW)

dpssDeriv <- function(DW,NW){
  stopifnot(.is.dpss(DW))
  efn <- DW$v
  ev <- DW$eigen
  ndata <- nrow(efn)
  K <- ncol(efn)
  W <- NW / ndata
  tn <- 1:(ndata-1)

  b <- c(0, sin( 2*pi*W*tn) / (pi * tn^2) )
  a <- 2 * W * c(0, cos( 2*pi*W*tn ) / tn )
  y <- a - b

  efnDeriv <- matrix(data = NA, nrow = ndata, ncol = K)
  x <- as.matrix(c(y,0,-y[ndata:2]))
  for(k in 1:K){
    p <- c(efn[,k],rep(0,ndata))
    h <- as.vector(fft(p)*fft(x))
    efnDeriv[,k] <- Re(fft(h, inverse = TRUE)[1:ndata] / length(h)) / ev[k]
  }

  class(efnDeriv) <- "dpssDeriv"
  return(efnDeriv)
}



#' @title computeEigencoefficients
#' @description Computes the eigencoefficients of a series given a set of tapers
#' @usage computeEigencoefficients(x, dw, nFFT)
#' @param x A time series
#' @param dw An N x K matrix of orthogonal tapers
#' @param nFFT The length of the series after zero padding
#' @details The series is first zero padded and tapered. Then, the eigencoefficients are computed by Fast Fourier transform on the tapered, padded series.
#' @returns yk The K x nFFT matrix of eigencoefficients
#' @references Thomson, David J. "Spectrum estimation and harmonic analysis." Proceedings of the IEEE 70.9 (1982): 1055-1096.
#' @export
#' @examples
#' N <- 1024
#' NW <- 4
#' K <- 2*NW - 1
#' DW <- dpss(N,K,NW,returnEigenvalues=FALSE)
#'
#' X <- rnorm(N)
#' nFFT <- 2^ceiling(log2(2*N))
#'
#' yk <- computeEigencoefficients(X,DW$v, nFFT)

computeEigencoefficients <- function(x, dw, nFFT){
  N <- length(x)
  stopifnot(nFFT >= N)
  padX <- c(x, rep(0,nFFT-N))
  padDW <- rbind(dw, matrix(0, nrow = nFFT - N, ncol = ncol(dw)))
  taperX <- padDW * padX

  yk <- mvfft(taperX)
  return(yk)
}

#' @title IFcompute
#' @description Computes the estimate of the instantaneous frequency of a time series based on the Slepian sequences.
#' @usage IFcompute(yk,V,Vdot)
#' @param yk The matrix of eigencoefficients. M x K where M = number of Fourier frequencies, K = number of tapers.
#' @param V The matrix of Slepian sequences. N x K where N = length of the series in the time domains.
#' @param Vdot The matrix of derivatives of the Slepian sequences. N x K matrix.
#' @details For each frequency at which the eigencoefficients are computed, an inverse series and its derivative series are also computed. These are used together to produce a matrix of instantaneous frequency series in the time domain.
#' @returns phi The N x M matrix of instantaneous frequency series, one column for each frequency at which the eigencoefficients are computed.
#' @references D. J. Thomson, "Inverse-constrained projection filters," Proc. SPIE 4478, Wavelets: Applications in Signal and Image Processing IX, (5 December 2001); https://doi.org/10.1117/12.449708
#' @references  Blanchette, Kian. Multitaper statistical tests for the detection of frequency-modulated signals. MSc. Thesis. Queen's University (Canada), 2020.
#' @export
#' @examples
#' N <- 1024
#' NW <- 4
#' K <- 2*NW - 1
#' DW <- dpss(N,K,NW,returnEigenvalues=FALSE)
#'
#' Vdot <- dpssDeriv(DW,NW)
#'
#' X <- rnorm(N)
#' sp <- spec.mtm(X, dpssIN = DW, returnInternals = TRUE)
#' yk <- sp$mtm$eigenCoefs
#'
#' instFreq <- IFcompute(yk, DW$v, Vdot)

IFcompute <- function(yk, V, Vdot){
  U <- tcrossprod(V, Re(yk))
  W <- tcrossprod(V, Im(yk))
  Udot <- tcrossprod(Vdot, Re(yk))
  Wdot <- tcrossprod(Vdot, Im(yk))

  num <- U*Wdot - Udot*W
  amp2 <- U^2 + W^2

  phi <- num / (2 * pi * amp2)

  return(phi)
}

#' @title ModulatedF
#' @description Computes the multitaper polynomial modulated test statistics for a given time series
#' @usage ModulatedF(x, K, NW, P, alpha, nFFT, dpssIN, apIN, derivIN, returnInternals, plotEase)
#' @param x A time series
#' @param K The number of tapers
#' @param NW The multitaper time-bandwidth product
#' @param P The highest degree at which to compute the test statistics
#' @param alpha Parameter of the Gegenbauer polynomials used to generate the set of polynomials associated to the Slepian sequences. 0.75 by default
#' @param nFFT The length of the series after zero padding. 2^ceiling(log2(2*length(x))) by default
#' @param dpssIN A dpss (Slepian sequence) object that can be passed in instead of being generated within the function.
#' @param apIN An associated polynomial object that can be passed in instead of being generated within the function.
#' @param derivIN A dpss derivative object that can be passed in instead of being generated within the function.
#' @param returnInternals A Boolean object. If TRUE, return the intermediate objects created inside this function. If FALSE (the default), only return the test statistics.
#' @param plotEase A Boolean object. If TRUE, sets any elements of ModF < 1 to 1 for clearer plotting on a log scale.
#' @details The estimates of the instantaneous frequency are computed for each Fourier frequency. Then the eigencoefficients and polynomial coefficients are derived for the instantaneous frequencies, which are then used to compute an F statistic.
#' @returns ModF The P x nFFT matrix containing the multitaper polynomial modulated test statistics.
#' @references  Blanchette, Kian. Multitaper statistical tests for the detection of frequency-modulated signals. MSc. Thesis. Queen's University (Canada), 2020
#' @importFrom multitaper dpss
#' @export
#' @examples
#' N <- 1024
#' NW <- 4
#' K <- 2*NW - 1
#' P <- 3
#' dpssIN <- dpss(N,K,NW,returnEigenvalues=FALSE)
#' apIN <- apDeriv(V = dpssIN$v, maxdeg = P, alpha = 0.75)
#' derivIN <- dpssDeriv(DW = dpssIN, NW = NW)
#'
#' x <- rnorm(N)
#'
#' mF <- ModulatedF(x, K, NW, P, alpha = 0.75, dpssIN = dpssIN, apIN = apIN, derivIN = derivIN)

ModulatedF <- function(x, K, NW, P, alpha = 0.75, nFFT = "default", dpssIN = NULL, apIN = NULL, derivIN = NULL, returnInternals = FALSE, plotEase = TRUE){

  stopifnot(K > P)

  N <- length(x)

  if(nFFT == "default") {
    nFFT <- 2^ceiling(log2(2*N))
  } else {
    stopifnot(is.numeric(nFFT))
  }
  stopifnot(nFFT >= N)

  if(!.is.dpss(dpssIN)){
    dpssIN <- dpss(n = N, k = K, nw = NW)
  }
  if(!.is.deriv(derivIN)){
    derivIN <- dpssDeriv(DW = dpssIN, NW = NW)
  }
  if(!.is.ap(apIN)){
    apIN <- apDeriv(V = dpssIN$v, maxdeg = P, alpha = alpha)
  }

  V <- dpssIN$v
  H <- apIN$H[,-1, drop = FALSE]

  xEigenCoefs <- computeEigencoefficients(x, dw = V, nFFT = nFFT)

  instFreq <- IFcompute(xEigenCoefs, V, derivIN)

  IFEigenCoefs <- crossprod(V,instFreq)
  IFPolyEigenCoefs <- crossprod(H,IFEigenCoefs)

  ModF <- .modF(P,K,IFPolyEigenCoefs,IFEigenCoefs,nFFT)

  if(plotEase){
    for(j in 1:P){
      idx <- which(ModF[j,] < 1)
      ModF[j,idx] <- 1
    }
  }

  if(!returnInternals){
    out <- list(ModF = ModF)
    return(out)
  }
  else{
    comp <- list(dpss = dpssIN, dpssDeriv = derivIN, assPoly = apIN, N = N, nFFT = nFFT)
    out <- list(ModF = ModF, instFreq = instFreq, IFEigenCoefs = IFEigenCoefs, IFPolyEigenCoefs = IFPolyEigenCoefs, xEigenCoefs = xEigenCoefs, comp = comp)
    return(out)
  }

}

#' @title findLocalFMaxM
#' @description Finds local maxima that exceed a given threshold, typically a percentile of an F distribution.
#' @usage  findLocalFMaxM(obj, k, cutoff)
#' @param obj The F statistic object, a matrix of size (number of Fourier frequencies) x (number of polynomial degrees)
#' @param k The number of tapers used to compute the polynomial F statistics.
#' @param cutoff The threshold above which we want to find local maxima.
#' @details The function finds all the indices that are above the provided cutoff. The local maxima are the turning points of this set of indices. That is, the local maxima are the points above the threshold that are greater than the points just before and just after them.
#' @returns MAXES A list of size P, where P is the number of polynomial degrees for which the F statistics were computed. Each element of the list contains the indices corresponding to the Fourier frequencies where local maxima can be found.
#' @references  Blanchette, Kian. Multitaper statistical tests for the detection of frequency-modulated signals. MSc. Thesis. Queen's University (Canada), 2020
#' @export
#' @examples
#' N <- 1024
#' nFFT <- 2^ceiling(log2(2*N))
#' NW <- 4
#' K <- 2*NW - 1
#' P <- 3
#'
#' dpssIN <- dpss(N,K,NW,returnEigenvalues=FALSE)
#' apIN <- apDeriv(V = dpssIN$v, maxdeg = P, alpha = 0.75)
#' derivIN <- dpssDeriv(DW = dpssIN, NW = NW)
#'
#' noiz <- rnorm(N)
#' n <- 0:(N-1)
#' A <- 0.5
#' freq <- seq(1/nFFT, 0.5, by = 1/nFFT)
#' f <- freq[nFFT/2]
#' B <- 0.0025
#' tt <- 2 * n / (N-1) - 1
#' FMpoly <- B * (1 - 2 * tt^2)
#' PhMod <- 2*pi*f*n + 2*pi*cumsum(FMpoly)
#' sig <- A * cos(PhMod)
#' x <- sig + noiz
#'
#' mF <- ModulatedF(x, K, NW, P, alpha = 0.75, dpssIN = dpssIN, apIN = apIN, derivIN = derivIN)
#' 
#' maxes <- findLocalFMaxM(mF$ModF,K,cutoff=1-1/N)

findLocalFMaxM <- function(obj, k, cutoff){

  M <- nrow(obj)
  stopifnot(k > M)
  MAXES <- list()
  for(m in 1:M){
    Fval <- obj[m,]
    fMaxInd <- which(Fval > qf(cutoff, 1, k-m))
    maxes <- c()

    if (length(fMaxInd) == 0){
      next
    }

    for (i in 1:length(fMaxInd)){
      if (fMaxInd[i] == 1 || fMaxInd[i] == length(Fval)){
        next
      }

      if (Fval[fMaxInd[i]] > Fval[fMaxInd[i]-1] &&
          Fval[fMaxInd[i]] > Fval[fMaxInd[i]+1]){
        maxes <- c(maxes, fMaxInd[i])
      }
    }
    MAXES[[m]] <- maxes
  }
  return(MAXES)
}