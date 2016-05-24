#' @export ZeroOrderPhaseCorrection
ZeroOrderPhaseCorrection <- function (RawSpect_data, plot_rms=NULL, returnAngle = FALSE) {
  begin_info <- beginTreatment("ZeroOrderPhaseCorrection", RawSpect_data)
  RawSpect_data <- begin_info[["Signal_data"]]

  debug_plot <- F

  rms <- function(ang, y) {
    if (debug_plot) {
      abline(v=ang, col="gray60")
    }
    roty <- y * exp(complex(real=0, imaginary=ang))
    Rey <- Re(roty)
    negRey <- Rey[Rey < 0]
    a <- mean(negRey)
    b <- mean(Rey)
    ssa <- sqrt(sum((negRey - a)^2))
    ssb <- sum((Rey - b)^2)
    return(ssa/ssb)
  }

  n <- nrow(RawSpect_data)
  rnames <- rownames(RawSpect_data)
  Angle <- c()
  for (k in 1:n) {
    # The function is rms is periodic (period 2pi) and it seems that there is a phase x
    # such that rms is unimodal (i.e. decreasing then increasing) on the interval [x; x+2pi].
    # However, if we do the optimization for example on [x-pi; x+pi],
    # instead of being decreasing then increasing, it might be increasing then decreasing
    # in which case optimize, thinking it is a valley will have to choose between
    # the left or the right of this hill and if it chooses wrong, it will end up at
    # like x-pi while the minimum is close to x+pi.

    # Supposing that rms is unimodal, the classical 1D unimodal optimization will
    # work in either [-pi;pi] or [0;2pi] (this is not easy to be convinced by that I agree)
    # and we can check which one it is simply by the following trick
    f0  <- rms(0,  RawSpect_data[k,])
    fpi <- rms(pi, RawSpect_data[k,])
    if (f0 < fpi) {
      interval <- c(-pi, pi)
    } else {
      interval <- c(0, 2*pi)
    }

    debug_plot <- F # rms should not plot anything now, only when called by optimize
    if (!is.null(plot_rms) && rnames[k] %in% plot_rms) {
      x <- seq(min(interval),max(interval),length.out=100)
      y <- rep(1,100)
      for (K in (1:100)) {
        y[K] <- rms(x[K], RawSpect_data[k,])
      }
      plot(x, y)
      debug_plot <- T
    }

    best <- optimize(rms, interval=interval, maximum=FALSE, y=RawSpect_data[k,])
    ang <- best[["minimum"]]

    if (debug_plot) {
      abline(v=ang, col="black")
      dev.off()
    }

    RawSpect_data[k,] <- RawSpect_data[k,] * exp(complex(real=0, imaginary=ang))
    Angle = c(Angle, ang)
  }
  RawSpect_data <- endTreatment("ZeroOrderPhaseCorrection", begin_info, RawSpect_data)
    if (returnAngle) {
    return(list(RawSpect_data=RawSpect_data, Angle=Angle))
  } else {
    return(RawSpect_data)
  }
}
