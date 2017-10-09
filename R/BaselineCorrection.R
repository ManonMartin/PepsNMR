#' @export BaselineCorrection
#' @import ptw
#' @import Matrix
BaselineCorrection <- function(Spectrum_data, ptw.bc = TRUE, maxIter = 42, 
                               lambda.bc = 1e+07, p.bc = 0.05, eps = 1e-08, 
                               ppm.bc = TRUE, exclude.bc = list(c(5.1,4.5)),
                               returnBaseline = F) {
  
  # Data initialisation ----------------------------------------------
  begin_info <- beginTreatment("BaselineCorrection", Spectrum_data, force.real = T)
  Spectrum_data <- begin_info[["Signal_data"]]
  p <- p.bc
  lambda <- lambda.bc
  n <- dim(Spectrum_data)[1]
  m <- dim(Spectrum_data)[2]
  
  
  # Data check
  checkArg(ptw.bc, c("bool"))
  checkArg(maxIter, c("int", "pos"))
  checkArg(lambda, c("num", "pos0"))
  checkArg(p.bc, c("num", "pos0"))
  checkArg(eps, c("num", "pos0"))
  checkArg(returnBaseline, c("bool"))
  checkArg(ppm.bc, c("bool"))
  checkArg(unlist(exclude.bc), c("num"), can.be.null = TRUE)
  
  # Define the interval where to search for (by defining Data)
  if (is.null(exclude.bc)) {
    exclude_index <- NULL
  } else  {
  # if ppm.bc == TRUE, then exclude.bc is in the colnames values, else, in the column
  # index
    if (ppm.bc == TRUE)  {
      colindex <- as.numeric(colnames(Spectrum_data))
    } else  {
      colindex <- 1:m
    }
    
    Int <- vector("list", length(exclude.bc))
    for (i in 1:length(exclude.bc))  {
      Int[[i]] <- indexInterval(colindex, from = exclude.bc[[i]][1], 
                                to = exclude.bc[[i]][2], inclusive = TRUE)
    }
    exclude_index <- unlist(Int)
  }
  
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
    asysm <- function(y, lambda, p, eps, exclude_index) {
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
        p_vect <- rep((1-p), m) # if y <= z + eps
        p_vect[y > z + eps | y < 0] <- p  # if y > z + eps | y < 0
        if(!is.null(exclude_index)){
          p_vect[exclude_index] <- 0 # if exclude area
        }
        
        w <- p_vect  
        # w <- p * (y > z + eps | y < 0) + (1 - p) * (y <= z + eps)
        
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

  # for (k in 1:n) {
    # Baseline[k, ] <- asysm(y = Spectrum_data[k, ], lambda = lambda, p = p, eps = eps)
    
  if (ptw.bc ){
    Baseline <- apply(Spectrum_data,1, asysm, lambda = lambda, p = p, 
                      eps = eps)
  }else {
    Baseline <- apply(Spectrum_data,1, asysm, lambda = lambda, p = p, 
                      eps = eps, exclude_index = exclude_index)
  }
    
    
    Spectrum_data <- Spectrum_data - t(Baseline)
  # }
  
  # Data finalisation ----------------------------------------------
  Spectrum_data <- endTreatment("BaselineCorrection", begin_info, Spectrum_data)  # FIXME create removeImaginary filter ??
  
  if (returnBaseline) {
    return(list(Spectrum_data = Spectrum_data, Baseline = Baseline))
  } else {
    return(Spectrum_data)
  }
}
