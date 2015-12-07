SolventSuppression <- function(Fid_data, lambda=1e6, use.ptw=TRUE, plotSolvent=F, returnSolvent=F) {
  begin_info <- beginTreatment("SolventSuppression", Fid_data)
  Fid_data <- begin_info[["Signal_data"]]
  checkArg(use.ptw, c("bool"))
  checkArg(lambda, c("num", "pos0"))
  if (use.ptw) {
#     require("ptw")
    difsm <- ptw::difsm
  } else {
    difsm <- function(y, d=2, lambda) {
#       require('Matrix')
      m <- length(y)
      # Sparse identity matrix m x m
      E <- Matrix::Diagonal(m)
      D <- Matrix::diff(E, differences = d)
      A <- E + lambda * Matrix::t(D) %*% D
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
    solventRe <- difsm(FidRe, lambda=lambda)
    solventIm <- difsm(FidIm, lambda=lambda)
    if (plotSolvent) {
      m = length(FidRe)
      plot(1:m,FidRe,type="l",col="red")
      lines(1:m,solventRe,type="l", col="blue")
      plot(1:m,FidIm,type="l",col="red")
      lines(1:m,solventIm,type="l", col="blue")
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