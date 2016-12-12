#' @export FourierTransform
FourierTransform <- function(Fid_data, Fid_info = NULL, SW_h = NULL, SW = NULL, O1 = NULL, reverse.axis = TRUE) {
  
  # Data initialisation and checks ----------------------------------------------
  begin_info <- beginTreatment("FourierTransform", Fid_data, Fid_info)
  Fid_data <- begin_info[["Signal_data"]]
  Fid_info <- begin_info[["Signal_info"]]
  
  SW_h <- getArg(SW_h, Fid_info, "SW_h")
  SW <- getArg(SW, Fid_info, "SW")  # Sweep Width in ppm (semi frequency scale in ppm)
  O1 <- getArg(O1, Fid_info, "O1")
  
  m <- ncol(Fid_data)
  n <- nrow(Fid_data)
  
  checkArg(reverse.axis, c("bool"))
  
  # Fourier Transformation ----------------------------------------------
  # mvfft does the unnormalized fourier transform (see ?mvfft), so we need divide
  # by m.  It does not matter a lot in our case since the spectrum will be
  # normalized.
  
  # FT
  RawSpect_data <- fftshift1D2D(t(stats::mvfft(t(Fid_data))))/m
  # recover the frequencies values
  f <- ((0:(m - 1)) - floor(m/2)) * Fid_info[1, "SW_h"]/m
  
  if(reverse.axis == TRUE) {
    revind <- rev(1:m)
    RawSpect_data <- RawSpect_data[,revind] # reverse the spectrum
  }
  
  colnames(RawSpect_data) <- f
  
  
  # PPM conversion ----------------------------------------------
  
  # The Sweep Width has to be the same since the column names are the same
  
  ppmInterval <- SW/m  # FIXME divide by two ??
  
  O1index <- binarySearch(a = f, target = O1, lower = TRUE)
  
  end <- O1index - m
  start <- O1index -1
  ppmScale <- (start:end) * ppmInterval
  RawSpect_data <- matrix(RawSpect_data, nrow = n, ncol =  -(end - start) + 1, dimnames = 
                              list(rownames(RawSpect_data), ppmScale))
   
  
  # Data finalisation ----------------------------------------------
  return(endTreatment("FourierTransform", begin_info, RawSpect_data))
}
