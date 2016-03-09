#' @export Warping
Warping <- function(RawSpect_data,
                    normalization.type=c("median","mean","firstquartile","peak","none"), from.normW=3.05, to.normW=4.05,
                    Reference.choosing=c("fixed", "before", "after"), Reference=1,
                    optim.crit=c("RMS", "WCC"), ptw.wp=F, K=3, L=40,
                    lambda.smooth=0, deg=3, lambda.bspline=0.01, kappa=0.0001,
                    max_it_Bspline=10, returnReference=F, returnWarpingfunc=F) {
  meanSqrDiff <- function(m, row) {
    # for the row ref, x - m[ref,] is 0
    # to get, the mean, we divide by nrow(data)-1 because the row ref is ignored
    # we could actually just do the sum and not the mean, both choices are good
    return (sum(apply(m, 1, function(x) sum((x - m[row,])^2, na.rm=T))) / (nrow(m)-1))
  }

  begin_info <- beginTreatment("Warping", RawSpect_data, force.real=T)
  RawSpect_data <- begin_info[["Signal_data"]]
  normalization.type <- match.arg(normalization.type)
  Reference.choosing <- match.arg(Reference.choosing)
  optim.crit <- match.arg(optim.crit)
  checkArg(K, c("int", "pos"))
  checkArg(L, c("int", "pos0"))
  checkArg(lambda.smooth, c("num", "pos0"))
  checkArg(deg, c("int", "pos"))
  checkArg(lambda.bspline, c("num", "pos0"))
  checkArg(kappa, c("num", "pos0"))
  checkArg(max_it_Bspline, c("int", "pos"))

  if (is.null(Reference) || !(Reference %in% row.names(RawSpect_data))) {
      checkArg(Reference, c("int", "pos"))
      Reference <- row.names(RawSpect_data)[Reference]
  }

  if (L > 0 && L <= deg) {
    stop("L should be greater than deg because with 1 interval, there is already deg+1 Bsplines.")
  }
  n <- nrow(RawSpect_data)
  m <- ncol(RawSpect_data)
  if (normalization.type != "none") {
    Normalization(RawSpect_data, normalization.type, from.normW, to.normW)
  }
  rnames <- rownames(RawSpect_data)
  if (n > 1) {
    if (Reference.choosing == "fixed") {
      pool <- Reference
    } else if (Reference.choosing == "before") {
      argmin <- which.min(sapply(rnames, function (i) meanSqrDiff(RawSpect_data,i)))
      pool <- names(argmin)
    } else {
      pool <- rnames
    }
    best.meanSqrDiff <- NULL
    best.Warped_data <- NULL
    decreasing <- FALSE
    for (Reference in pool) {
      ref <- RawSpect_data[Reference, ]
      samp_rnames <- rnames[rnames != Reference]
      sample <- RawSpect_data[samp_rnames, , drop=F]
      beta <- rep(0, K+1)
      cur.Warped_data <- RawSpect_data
      if (K >= 1) {
        # This is beta_1 which should be approximately 1 because
        # if they are perfectly aligned w(v) = v = 0*1 + 1*v + 0*v^2 + 0*v^3
        beta[2] = 1
      }
      if (ptw.wp) {
#         require("ptw")
        ptw.output <- ptw::ptw(ref, sample, optim.crit=optim.crit, init.coef=beta,
                               smooth.param=lambda.smooth)
        cur.Warped_data[samp_rnames, ] <- ptw.output$warped.sample
        # We transpose because 'diff' is along the columns
        w <- t(ptw.output$warp.fun)
      } else {
        if (optim.crit == "WCC") {
          stop("WCC is only implemented in ptw, set ptw.wp=T to use WCC.")
        }
        for (samp_rname in samp_rnames) {
          sw.output <- SingleWarp(ref=ref, sample=sample[samp_rname, ], beta=beta, L=L,
                                  lambda.smooth=lambda.smooth, deg=deg,
                                  lambda.bspline=lambda.bspline,
                                  kappa=kappa, max_it_Bspline=max_it_Bspline)
          warped <- sw.output$warped
          w <- sw.output$w
          cur.Warped_data[samp_rname, ] = warped
        }
      }
      cur.meanSqrDiff <- meanSqrDiff(cur.Warped_data, Reference)
      if (is.null(best.meanSqrDiff) || cur.meanSqrDiff < best.meanSqrDiff) {
        best.meanSqrDiff <- cur.meanSqrDiff
        best.Warped_data <- cur.Warped_data
        decreasing <- (FALSE %in% (diff(w) > 0))
      }
    }
    if (decreasing) {
      warning("The warping function is not increasing for the sample ", samp_rname, ".")
    }
  } else {
    best.Warped_data <- RawSpect_data
  }
  RawSpect_data = endTreatment("Warping", begin_info, best.Warped_data)
  if (returnReference & ! returnWarpingfunc) {
    return(list(RawSpect_data=RawSpect_data, Reference=Reference))
  } else if (returnReference & returnWarpingfunc) {
    return(list(RawSpect_data=RawSpect_data, Reference=Reference, Warpingfunc=w))
  } else if (! returnReference & returnWarpingfunc) {
    return(list(RawSpect_data=RawSpect_data, Warpingfunc=w))
  } else {return(RawSpect_data) }
}
