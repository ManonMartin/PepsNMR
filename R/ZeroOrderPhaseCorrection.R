#' @export ZeroOrderPhaseCorrection
#' @importFrom stats quantile sd
#' @importFrom graphics par plot
#' 
ZeroOrderPhaseCorrection <- function(Spectrum_data, type.zopc = c("rms", "manual", "max"), 
                                     plot_rms = NULL, returnAngle = FALSE, createWindow = TRUE, 
                                     angle = NULL, plot_spectra = FALSE,  
                                     ppm.zopc = TRUE, exclude.zopc = list(c(5.1,4.5))) {
  
  
  # Data initialisation and checks ----------------------------------------------
  
  # Entry arguments definition:
  # plot_rms : graph of rms criterion returnAngle : if TRUE, returns avector of
  # optimal angles createWindow : for plot_rms plots angle : If angle is not NULL,
  # spectra are rotated according to the angle vector values
  # plot_spectra : if TRUE, plot rotated spectra  
 
  
  
  begin_info <- beginTreatment("ZeroOrderPhaseCorrection", Spectrum_data)
  Spectrum_data <- begin_info[["Signal_data"]]
  n <- nrow(Spectrum_data)
  m <- ncol(Spectrum_data)
  
  rnames <- rownames(Spectrum_data)
  
  # Check input arguments
  type.zopc <- match.arg(type.zopc)
  checkArg(ppm.zopc, c("bool"))
  checkArg(unlist(exclude.zopc), c("num"), can.be.null = TRUE)
  
  
  # type.zopc in c("max", "rms") -----------------------------------------
  if (type.zopc %in% c("max", "rms")) {
    # angle is found by optimization
    
    # rms function to be optimised
    rms <- function(ang, y, meth = c("max", "rms"))  {
      # if (debug_plot) { graphics::abline(v=ang, col='gray60') }
      roty <- y * exp(complex(real = 0, imaginary = ang))  # spectrum rotation
      Rey <- Re(roty)
      
      if (meth == "rms")  {
        ReyPos <- Rey[Rey >= 0]  # select positive intensities
        POSss <- sum((ReyPos)^2, na.rm = TRUE)  # SS for positive intensities
        ss <- sum((Rey)^2, na.rm = TRUE)  #  SS for all intensities
        return(POSss/ss)  # criterion : SS for positive values / SS for all intensities 
      } else  {
        maxi <- max(Rey, na.rm = TRUE)
        return(maxi)
      }
    }
    
   
    # Define the interval where to search for (by defining Data)
    if (is.null(exclude.zopc)) {
      Data <- Spectrum_data
    } else  {
      
      # if ppm.zopc == TRUE, then fromto is in the colnames values, else, in the column
      # index
      if (ppm.zopc == TRUE)  {
        colindex <- as.numeric(colnames(Spectrum_data))
      } else  {
        colindex <- 1:m
      }
      
      # Second check for the argument exclude.zopc
        diff <- diff(unlist(exclude.zopc))[1:length(diff(unlist(exclude.zopc)))%%2 !=0]
        for (i in 1:length(diff)) {
          if (ppm.zopc == TRUE & diff[i] >= 0)  {
            stop(paste("Invalid region removal because from <= to in ppm.zopc"))
          } else if (ppm.zopc == FALSE & diff[i] <= 0) {stop(paste("Invalid region removal because from >= to in column index"))}
        }
      
      
      Int <- vector("list", length(exclude.zopc))
      for (i in 1:length(exclude.zopc))  {
        Int[[i]] <- indexInterval(colindex, from = exclude.zopc[[i]][1], 
                                  to = exclude.zopc[[i]][2], inclusive = TRUE)
      }
      
      vector <- rep(1, m)
      vector[unlist(Int)] <- 0
      if (n > 1)  {
        Data <- sweep(Spectrum_data, MARGIN = 2, FUN = "*", vector)  # Cropped_Spectrum
      } else   {
        Data <- Spectrum_data * vector
      }  # Cropped_Spectrum
    }
    
    
    # angles computation
    Angle <- c()
    for (k in 1:n)
    {
      # The function is rms is periodic (period 2pi) and it seems that there is a phase
      # x such that rms is unimodal (i.e. decreasing then increasing) on the interval
      # [x; x+2pi].  However, if we do the optimization for example on [x-pi; x+pi],
      # instead of being decreasing then increasing, it might be increasing then
      # decreasing in which case optimize, thinking it is a valley will have to choose
      # between the left or the right of this hill and if it chooses wrong, it will end
      # up at like x-pi while the minimum is close to x+pi.
      
      # Supposing that rms is unimodal, the classical 1D unimodal optimization will
      # work in either [-pi;pi] or [0;2pi] (this is not easy to be convinced by that I
      # agree) and we can check which one it is simply by the following trick
      
      f0 <- rms(0, Data[k, ],meth = type.zopc)
      fpi <- rms(pi, Data[k, ], meth = type.zopc)
      if (f0 < fpi) {
        interval <- c(-pi, pi)
      } else {
        interval <- c(0, 2 * pi)
      }
      
      # graphs of rms criteria
      debug_plot <- F  # rms should not plot anything now, only when called by optimize
      if (!is.null(plot_rms) && rnames[k] %in% plot_rms) {
        x <- seq(min(interval), max(interval), length.out = 100)
        y <- rep(1, 100)
        for (K in (1:100))   {
          y[K] <- rms(x[K], Data[k, ],  meth = type.zopc)
        }
        if (createWindow == TRUE)  {
          grDevices::dev.new(noRStudioGD = FALSE)
        }
        graphics::plot(x, y, main = paste("Criterion maximization \n", 
                                          rownames(Data)[k]), ylim = c(0, 1.1),
                       ylab = "positiveness criterion", xlab = "angle ")
        debug_plot <- T
      }
      
      # Best angle
      best <- stats::optimize(rms, interval = interval, maximum = TRUE, 
                              y = Data[k,],  meth = type.zopc)
      ang <- best[["maximum"]]
      
      
      if (debug_plot)  {
        graphics::abline(v = ang, col = "black")
        graphics::text(x = (ang+0.1*ang), y = (y[ang]-0.1*y[ang]), labels = round(ang, 3))
      }
      
      # Spectrum rotation
      Spectrum_data[k, ] <- Spectrum_data[k, ] * exp(complex(real = 0, imaginary = ang))
      Angle <- c(Angle, ang)
    }
    
    
    
    
  } else {
    # type.zopc is "manual" -------------------------------------------------------
    # if Angle is already specified and no optimisation is needed
    Angle <- angle
    
    if (!is.vector(angle)) {
      stop("angle is not a vector")
    }
    
    if (!is.numeric(angle))  {
      stop("angle is not a numeric")
    }
    
    if (length(angle) != n) {
      stop(paste("angle has length", length(angle), "and there are", n, "spectra to rotate."))
    }
    for (k in 1:n)  {
      Spectrum_data[k, ] <- Spectrum_data[k, ] * exp(complex(real = 0, imaginary = - angle[k]))
    }
  }
  
  
  #  Draw spectra
  if (plot_spectra == TRUE) {
    nn <- ceiling(n/4)
    i <- 1
    for (k in 1:nn)  {
      if (createWindow == TRUE)  {
        grDevices::dev.new(noRStudioGD = FALSE)
      }
      graphics::par(mfrow = c(4, 2))
      while (i <= n)   {
        last <- min(i + 4 - 1, n)
        graphics::plot(Re(Spectrum_data[i, ]), type = "l", ylab = "intensity", 
          xlab = "Index", main = paste0(rownames(Spectrum_data)[i], " - Real part"))
        graphics::plot(Im(Spectrum_data[i, ]), type = "l", ylab = "intensity", 
          xlab = "Index", main = paste0(rownames(Spectrum_data)[i], " - Imaginary part"))
        i <- i + 1
      }
      i <- last + 1
    }
  }
  
  
  # Data finalisation ----------------------------------------------
  
  Spectrum_data <- endTreatment("ZeroOrderPhaseCorrection", begin_info, Spectrum_data)
  if (returnAngle) {
    return(list(Spectrum_data = Spectrum_data, Angle = Angle))
  } else {
    return(Spectrum_data)
  }
}
