#' @export ZeroOrderPhaseCorrection
#' @importFrom stats quantile sd
#' @importFrom graphics par plot
#' 
ZeroOrderPhaseCorrection <- function (RawSpect_data, plot_rms=NULL, returnAngle = FALSE, createWindow=TRUE, 
                                      Angle = NULL,   p.zo=0.8, plot_spectra = FALSE, rotation = TRUE) {
  begin_info <- beginTreatment("ZeroOrderPhaseCorrection", RawSpect_data)
  RawSpect_data <- begin_info[["Signal_data"]]
  n <- nrow(RawSpect_data)
  rnames <- rownames(RawSpect_data)
  
  debug_plot <- F
    
    if (is.null(Angle)) {

  rms <- function(ang, y) {
    # if (debug_plot) {
    #   graphics::abline(v=ang, col="gray60")
    # }
    roty <- y * exp(complex(real=0, imaginary=ang))
    Rey <- Re(roty)
    # negRey <- Rey[Rey < 0]
    # a <- mean(negRey)
    # b <- mean(Rey)
    # ssa <- sqrt(sum((negRey - a)^2))
    # ssb <- sum((Rey - b)^2)
    # return(ssa/ssb)
    
    Rey[abs(Rey)>=quantile(abs(Rey), .95)] = quantile(abs(Rey), .95) # trim
    Rey = abs(Rey)*sign(Rey)
    ReyPos <- Rey[Rey >= 0]
    
    POSss = sum((ReyPos-mean(ReyPos))^2)
    
    ss = sum((Rey - mean(Rey))^2)
    
    return(POSss/ss)
  }

  
  
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
      if (createWindow==TRUE) {
        grDevices::dev.new(noRStudioGD = FALSE) 
      }
      graphics:: plot(x, y, main =  rownames(RawSpect_data)[k])
      debug_plot <- T
    }

    best <- stats::optimize(rms, interval=interval, maximum=FALSE, y=RawSpect_data[k,])
    ang <- best[["minimum"]]

    if (debug_plot) {
      graphics::abline(v=ang, col="black")
      # grDevices::dev.off()
    }

    RawSpect_data[k,] <- RawSpect_data[k,] * exp(complex(real=0, imaginary=ang))
    Angle = c(Angle, ang)
  }
  
  
    } else {
      
    if (!is.vector(Angle)) {
      stop("Angle is not a vector")
    }
      
    if (!is.numeric(Angle)) {
      stop("Angle is not a numeric")
    }
      
    if (length(Angle) != n) {
      stop(paste("Angle has length", length(Angle), "and there are", n, "spectra to rotate."))
    }
      
      for (k in 1:n) {
        RawSpect_data[k,] <- RawSpect_data[k,] * exp(complex(real=0, imaginary=Angle[k]))
      }
  }
  
  
  #================== Detect a 180° rotation due to the water signal
  MEAN_Q = c()
  for (i in 1:nrow(RawSpect_data)) {
    data = Re(RawSpect_data[i,])
    data_p = data[data >= stats::quantile(data[data >=0 ], p.zo)]
    data_n = data[data <= stats::quantile(data[data <0 ], (1-p.zo))]
    
    mean_quant = (sum(data_p) + sum(data_n)) / (length(data_p) +length(data_n))
    # mean(p.zo% higher pos and neg values)
    MEAN_Q = c(MEAN_Q, mean_quant)
  }
  
  vect = which(MEAN_Q < 0)
  if (length(vect)!=0) {
    warning("The mean of", p.zo," positive and negative quantiles is negative for ", paste0(rownames(RawSpect_data)[vect],"; "))
    if(rotation == TRUE) {
      warning(" An automatic 180 degree rotation is applied to these spectra")
      Angle[vect] = Angle[vect] + pi       
    }

  }
  
  vect_risk = which(MEAN_Q<0.1*mean(MEAN_Q[MEAN_Q>0])) # is there any MEAN_Q with a very low value copared to mean of positive mean values?
  if (length(vect_risk)!=0) {
    warning("the rotation angle for spectra", paste0(rownames(RawSpect_data)[vect_risk],"; "), "might not be optimal, you need to check visually for those spectra")
  }
  
  # result of automatic rotation
  for (k in vect_risk) {
    RawSpect_data[k,] <- RawSpect_data[k,] * exp(complex(real=0, imaginary=Angle[k]))
  }
  #==================
  
  
  
  #========= Draw spectra
  if (plot_spectra==TRUE) {
    nn = ceiling(n/4)
    i = 1
    for (k in 1:nn) {
      graphics::par(mfrow=c(4,2))
      while (i <= n)
      {
        last = min(i + 4-1, n)
        graphics::plot(Re(RawSpect_data[i,]), type="l", ylab = "intensity", xlab="Index",main = paste0(rownames(RawSpect_data)[i], " - Real part"))
        graphics::plot(Im(RawSpect_data[i,]), type="l", ylab = "intensity", xlab="Index",main = paste0(rownames(RawSpect_data)[i], " - Imaginary part"))
        i=i+1
      }
      i = last + 1
    }
  }
  
  #=========
  
  RawSpect_data <- endTreatment("ZeroOrderPhaseCorrection", begin_info, RawSpect_data)
    if (returnAngle) {
    return(list(RawSpect_data=RawSpect_data, Angle=Angle)) 
  } else {
    return(RawSpect_data)
  }
}
