#' @export BaselineCorrection
#' @import ptw
#' @import Matrix
BaselineCorrection <- function(Spectrum_data, ptw.bc = TRUE, maxIter = 42, 
                               lambda.bc = 1e+07, p.bc = 0.05, eps = 1e-08, 
                               returnBaseline = F) {
  
  # Data initialisation ----------------------------------------------
  begin_info <- beginTreatment("BaselineCorrection", Spectrum_data, force.real = T)
  Spectrum_data <- begin_info[["Signal_data"]]
  p <- p.bc
  lambda <- lambda.bc
  n <- dim(Spectrum_data)[1]
  
  # Data check
  checkArg(ptw.bc, c("bool"))
  checkArg(maxIter, c("int", "pos"))
  checkArg(lambda, c("num", "pos0"))
  checkArg(p.bc, c("num", "pos0"))
  checkArg(eps, c("num", "pos0"))
  checkArg(returnBaseline, c("bool"))
  
  
  
  # Baseline Correction implementation definition ----------------------
  
  # 2 Ways: either use the function asysm from the ptw package or by 
  # built-in functions 
  if (ptw.bc) {
    asysm <- ptw::asysm
  } else {
    difsmw <- function(y, lambda, w, d) {
      # Weighted smoothing with a finite difference penalty cf Eilers, 2003.
      # (A perfect smoother) 
      # y: signal to be smoothed 
      # lambda: smoothing parameter 
      # w: weights (use0 zeros for missing values) 
      # d: order of differences in penalty (generally 2)
      m <- length(y)
      W <- Matrix::Diagonal(x=w)
      E <- Matrix::Diagonal(m)
      D <- Matrix::diff(E, differences = d)
      C <- Matrix::chol(W + lambda * t(D) %*% D)
      x <- Matrix::solve(C, Matrix::solve(t(C), w * y))
      return(as.numeric(x))
      
    }
    asysm <- function(y, lambda, p, eps) {
      # Baseline estimation with asymmetric least squares
      # y: signal
      # lambda: smoothing parameter (generally 1e5 to 1e8)
      # p: asymmetry parameter (generally 0.001)
      # d: order of differences in penalty (generally 2)
      # eps: 1e-8 in ptw package
      m <- length(y)
      w <- rep(1, m)
      i <- 1
      repeat {
        z <- difsmw(y, lambda, w, d = 2)
        w0 <- w
        w <- p * (y > z + eps | y < 0) + (1 - p) * (y <= z + eps)
        if (sum(abs(w - w0)) == 0) {
          break
        }
        i <- i + 1
        if (i > maxIter) {
          warning("cannot find Baseline estimation in asysm")
          break
        }
      }
      return(z)
    }
  }
  
  # Baseline estimation ----------------------------------------------
  Baseline <- matrix(NA, nrow = nrow(Spectrum_data), ncol = ncol(Spectrum_data))

  for (k in 1:n) {
    Baseline[k, ] <- asysm(y = Spectrum_data[k, ], lambda = lambda, p = p, eps = eps)
    if (F & k == 1) {
      m <- ncol(Spectrum_data)
      graphics::plot(1:mm, Spectrum_data[k, ], type = "l", col = "red")
      graphics::lines(1:mm, Baseline[k, ], type = "l", col = "blue")
      graphics::lines(1:m, Spectrum_data[k, ] - Baseline[k, ], type = "l",
        col = "green")
    }

    Spectrum_data[k, ] <- Spectrum_data[k, ] - Baseline[k, ]
  }
  
  # Data finalisation ----------------------------------------------
  Spectrum_data <- endTreatment("BaselineCorrection", begin_info, Spectrum_data)  # FIXME create removeImaginary filter ??
  
  if (returnBaseline) {
    return(list(Spectrum_data = Spectrum_data, Baseline = Baseline))
  } else {
    return(Spectrum_data)
  }
}
