#' @export SolventSuppression
SolventSuppression <- function(Fid_data, lambda.ss=1e6, ptw.ss=TRUE, plotSolvent=F, returnSolvent=F) {
  begin_info <- beginTreatment("SolventSuppression", Fid_data)
  Fid_data <- begin_info[["Signal_data"]]
  checkArg(ptw.ss, c("bool"))
  checkArg(lambda.ss, c("num", "pos0"))
  if (ptw.ss) {
    # Use of the function in ptw that smoothes signals with a finite difference penalty of order 2
    difsm <- ptw::difsm
  } else {
    # Or manual implementation based on sparse matrices for large data series (cf. Eilers, 2003. "A perfect smoother")
    difsm <- function(y, d=2, lambda) {

      m <- length(y)
      # Sparse identity matrix m x m
      E <- Matrix::Diagonal(m)
      D <- Matrix::diff(E, differences = d)
      A <- E + lambda.ss * Matrix::t(D) %*% D
      # base::chol does not take into account that A is sparse
      # and is extremely slow
      C <- Matrix::chol(A)
      x <- Matrix::solve(C, Matrix::solve(Matrix::t(C), y))
      return(as.numeric(x))
    }
  }
  n <- dim(Fid_data)[1]
  if (returnSolvent) {
    SolventRe <- Fid_data
    SolventIm <- Fid_data
  }
  for (i in 1:n) {
    FidRe <- Re(Fid_data[i,])
    FidIm <- Im(Fid_data[i,])
    solventRe <- difsm(FidRe, lambda=lambda.ss)
    solventIm <- difsm(FidIm, lambda=lambda.ss)
    if (plotSolvent) {
      m = length(FidRe)
      graphics::plot(1:m,FidRe,type="l",col="red")
      graphics::lines(1:m,solventRe,type="l", col="blue")
      graphics::plot(1:m,FidIm,type="l",col="red")
      graphics::lines(1:m,solventIm,type="l", col="blue")
    }
    FidRe <- FidRe - solventRe
    FidIm <- FidIm - solventIm
    Fid_data[i,] <- complex(real=FidRe, imaginary=FidIm)
    if (returnSolvent) {
      SolventRe[i,] = solventRe
      SolventIm[i,] = solventIm
    }
  }
  #if (plotSolvent) { # device independent
  #  dev.off()
  #}
  Fid_data <- endTreatment("SolventSuppression", begin_info, Fid_data)
  if (returnSolvent) {
    return(list(Fid_data=Fid_data, SolventRe=Re(SolventRe), SolventIm=Im(SolventIm)))
  } else {
    return(Fid_data)
  }
}
