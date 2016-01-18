#' @export BaselineCorrection
BaselineCorrection <- function (RawSpect_data, ptw.bc=TRUE, maxIter = 42,
                                lambda.bc=1e7, p=0.05, eps=1e-8, returnBaseline=F) {
  begin_info <- beginTreatment("BaselineCorrection", RawSpect_data, force.real=T)
  RawSpect_data <- begin_info[["Signal_data"]]
  if (ptw.bc) {
#     require("ptw")
    asysm <- ptw::asysm
  } else {
    difsmw <- function (y, lambda.bc, w, d) {
      # Weighted smoothing with a finite difference penalty
      # y:      signal to be smoothed
      # lambda.bc: smoothing parameter 
      # w:      weights (use0 zeros for missing values)
      # d:      order of differences in penalty (generally 2)
      m <- length(y)
      W <- Matrix::Diagonal(x=w)
      E <- Matrix::Diagonal(m)
      D <- diff(E, differences = d)
      C <- Matrix::chol(W + lambda.bc * t(D) %*% D)
      x <- Matrix::solve(C, Matrix::solve(t(C), w * y))
      return(as.numeric(x))
    }
    asysm <- function (y, d=2, lambda.bc, p, eps) {
      # Baseline estimation with asymmetric least squares
      # y:      signal
      # lambda.bc: smoothing parameter (generally 1e5 to 1e8)
      # p:      asymmetry parameter (generally 0.001)
      # d:      order of differences in penalty (generally 2)
      # eps:    1e-8 in ptw package
      m = length(y)
      w = rep(1, m)
      i <- 1
      repeat {
        z = difsmw(y, lambda.bc, w, d)
        w0 = w
        w <- p * (y > z+eps | y < 0) + (1 - p) * (y <= z+eps)
        if (sum(abs(w - w0)) == 0) {
          break
        }
        i <- i + 1
        if (i > maxIter) {
          warning("cannot find baseline estimation in asysm")
          break
        }
      }
      return(z)
    }
  }
  n <- nrow(RawSpect_data)
  m <- ncol(RawSpect_data)
  baseline = matrix(NA,nrow = n, ncol = m)
  for (k in 1:n) {
    baseline[k,] <- asysm(RawSpect_data[k,], lambda.bc, p, eps)
    if (F & k == 1) {
      m = ncol(RawSpect_data)
      plot(1:m,RawSpect_data[k,],type="l",col="red")
      lines(1:m,baseline[k,],type="l", col="blue")
      lines(1:m,RawSpect_data[k,] - baseline[k,],type="l", col="green")
    }
    RawSpect_data[k,] <- RawSpect_data[k,] - baseline[k,]
  }

RawSpect_data = endTreatment("BaselineCorrection", begin_info, RawSpect_data) # FIXME create removeImaginary filter ??

if (returnBaseline) {
  return(list(RawSpect_data=RawSpect_data, baseline=baseline))
} else {
  return(RawSpect_data)
  
 }
}
