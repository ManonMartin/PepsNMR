#' @export ZeroOrderPhaseCorrection
#' @importFrom stats quantile sd
#' @importFrom graphics par plot
#' 
ZeroOrderPhaseCorrection <- function(Spectrum_data, method = c("rms", "manual", "max"), 
                                     plot_rms = NULL, returnAngle = FALSE, createWindow = TRUE, 
                                     Angle = NULL, p.zo = 0.8, plot_spectra = FALSE, quant = 0.95, 
                                     freq = TRUE, fromto.0OPC = NULL) {
  
  
  # Data initialisation and checks ----------------------------------------------
  
  # Entry arguments definition:
  # plot_rms : graph of rms criterion returnAngle : if TRUE, returns avector of
  # optimal angles createWindow : for plot_rms plots Angle : If Angle is not NULL,
  # spectra are rotated according to the Angle vector values p.zo: idem que
  # rotation: à retirer plot_spectra : if TRUE, plot rotated spectra  
  # quant: probability for sample quantile used to trim the spectral intensities
  
  
  begin_info <- beginTreatment("ZeroOrderPhaseCorrection", Spectrum_data)
  Spectrum_data <- begin_info[["Signal_data"]]
  n <- nrow(Spectrum_data)
  m <- ncol(Spectrum_data)
  
  rnames <- rownames(Spectrum_data)
  
  # Check input arguments
  method <- match.arg(method)
  checkArg(freq, c("bool"))
  checkArg(unlist(fromto.0OPC), c("num"), can.be.null = TRUE)
  
  
  # method in c("max", "rms") -----------------------------------------
  if (method %in% c("max", "rms")) {
    # Angle is found by optimization
    
    # rms function to be optimised
    rms <- function(ang, y, meth = c("max", "rms"), p = 0.95)  {
      # if (debug_plot) { graphics::abline(v=ang, col='gray60') }
      roty <- y * exp(complex(real = 0, imaginary = ang))  # spectrum rotation
      Rey <- Re(roty)
      
      if (meth == "rms")  {
        si <- sign(Rey)  # sign of intensities
        
        Rey[abs(Rey) >= quantile(abs(Rey), p)] <- quantile(abs(Rey), p)  # trim the values
        Rey <- abs(Rey) * si  # spectral trimmed values
        ReyPos <- Rey[Rey >= 0]  # select positive intensities
        
        # POSss = sum((ReyPos-mean(ReyPos))^2) # centred SS for positive intensities
        POSss <- sum((ReyPos)^2)  # SS for positive intensities
        
        # ss = sum((Rey - mean(Rey) )^2) # centred SS for all intensities
        ss <- sum((Rey)^2)  #  SS for all intensities
        
        return(POSss/ss)  # criterion : SS for positive values / SS for all intensities 
      } else  {
        maxi <- max(Rey)
        return(maxi)
      }
    }
    
   
    # Define the interval where to search for (by defining Data)
    if (is.null(fromto.0OPC)) {
      Data <- Spectrum_data
    } else  {
      
      # if freq == TRUE, then fromto is in the colnames values, else, in the column
      # index
      if (freq == TRUE)  {
        colindex <- as.numeric(colnames(Spectrum_data))
      } else  {
        colindex <- 1:m
      }
      
      # Second check for the argument fromto.0OPC
      diff <- diff(unlist(fromto.0OPC))[1:length(diff(unlist(fromto.0OPC)))%%2 != 0]
      for (i in 1:length(diff))   {
        if (diff[i] <= 0)   {
          stop(paste("Invalid region removal because from > to"))
        }
      }
      
      Int <- vector("list", length(fromto.0OPC))
      for (i in 1:length(fromto.0OPC))  {
        Int[[i]] <- indexInterval(colindex, from = fromto.0OPC[[i]][1], 
                                  to = fromto.0OPC[[i]][2], inclusive = TRUE)
      }
      
      vector <- rep(0, m)
      vector[unlist(Int)] <- 1
      if (n > 1)  {
        Data <- sweep(Spectrum_data, MARGIN = 2, FUN = "*", vector)  # Cropped_Spectrum
      } else   {
        Data <- Spectrum_data * vector
      }  # Cropped_Spectrum
    }
    
    
    # Angles computation
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
      
      f0 <- rms(0, Data[k, ], p = quant, meth = method)
      fpi <- rms(pi, Data[k, ], p = quant, meth = method)
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
          y[K] <- rms(x[K], Data[k, ], p = quant, meth = method)
        }
        if (createWindow == TRUE)  {
          grDevices::dev.new(noRStudioGD = FALSE)
        }
        graphics::plot(x, y, main = rownames(Data)[k])
        debug_plot <- T
      }
      
      # Best angle
      best <- stats::optimize(rms, interval = interval, maximum = TRUE, 
                              y = Data[k,], p = quant, meth = method)
      ang <- best[["maximum"]]
      
      
      if (debug_plot)  {
        graphics::abline(v = ang, col = "black")
        # grDevices::dev.off()
      }
      
      # Spectrum rotation
      Spectrum_data[k, ] <- Spectrum_data[k, ] * exp(complex(real = 0, imaginary = ang))
      Angle <- c(Angle, ang)
    }
    
    
    
    
  } else {
    # method is "manual" -------------------------------------------------------
    # if Angle is already specified and no optimisation is needed
    
    if (!is.vector(Angle)) {
      stop("Angle is not a vector")
    }
    
    if (!is.numeric(Angle))  {
      stop("Angle is not a numeric")
    }
    
    if (length(Angle) != n) {
      stop(paste("Angle has length", length(Angle), "and there are", n, "spectra to rotate."))
    }
    for (k in 1:n)  {
      Spectrum_data[k, ] <- Spectrum_data[k, ] * exp(complex(real = 0, imaginary = ang))
    }
  }
  
  
  
  # #================== Detect a 180° rotation due to the water signal MEAN_Q = c()
  # for (i in 1:nrow(Spectrum_data)) { data = Re(Spectrum_data[i,]) data_p =
  # data[data >= stats::quantile(data[data >=0 ], p.zo)] data_n = data[data <=
  # stats::quantile(data[data <0 ], (1-p.zo))] mean_quant = (sum(data_p) +
  # sum(data_n)) / (length(data_p) +length(data_n)) # mean(p.zo% higher pos and neg
  # values) MEAN_Q = c(MEAN_Q, mean_quant) } vect = which(MEAN_Q < 0) if
  # (length(vect)!=0) { warning('The mean of', p.zo,' positive and negative
  # quantiles is negative for ', paste0(rownames(Spectrum_data)[vect],'; '))
  # if(rotation == TRUE) { warning(' An automatic 180 degree rotation is applied to
  # these spectra') Angle[vect] = Angle[vect] + pi } } vect_risk =
  # which(MEAN_Q<0.1*mean(MEAN_Q[MEAN_Q>0])) # is there any MEAN_Q with a very low
  # value copared to mean of positive mean values?  if (length(vect_risk)!=0)
  # { warning('the rotation angle for spectra',
  # paste0(rownames(Spectrum_data)[vect_risk],'; '), 'might not be optimal, you
  # need to check visually for those spectra') } # result of automatic rotation for
  # (k in vect_risk) { Spectrum_data[k,] <- Spectrum_data[k,] * exp(complex(real=0,
  # imaginary=Angle[k])) } #==================
  
  
  
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
