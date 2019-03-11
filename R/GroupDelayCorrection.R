#' @export GroupDelayCorrection
GroupDelayCorrection <- function(Fid_data, Fid_info = NULL, group_delay = NULL,
                                 verbose = FALSE) {
  
  
  # Data initialisation and checks ----------------------------------------------
  checkArg(verbose, c("bool"))
  begin_info <- beginTreatment("GroupDelayCorrection", Fid_data, Fid_info,verbose = verbose)
  Fid_data <- begin_info[["Signal_data"]]
  dimension_names <- dimnames(Fid_data)
  Fid_info <- begin_info[["Signal_info"]]
  checkArg(group_delay, c("num", "pos0"), can.be.null = TRUE)

  
  # if Fid_info and group_delay are NULL, getArg will generate an error
  
  group_delay <- getArg(group_delay, Fid_info, "GRPDLY", can.be.absent = TRUE)
  
  if (is.null(group_delay)) {
    
    # See DetermineBrukerDigitalFilter.m in matNMR MATLAB library
    group_delay_matrix <- matrix(c(44.75, 46, 46.311, 33.5, 36.5, 36.53, 66.625, 
      48, 47.87, 59.0833, 50.1667, 50.229, 68.5625, 53.25, 53.289, 60.375, 
      69.5, 69.551, 69.5313, 72.25, 71.6, 61.0208, 70.1667, 70.184, 70.0156, 
      72.75, 72.138, 61.3438, 70.5, 70.528, 70.2578, 73, 72.348, 61.5052, 70.6667, 
      70.7, 70.3789, 72.5, 72.524, 61.5859, 71.3333, NA, 70.4395, 72.25, NA, 
      61.6263, 71.6667, NA, 70.4697, 72.125, NA, 61.6465, 71.8333, NA, 70.4849, 
      72.0625, NA, 61.6566, 71.9167, NA, 70.4924, 72.0313, NA), nrow = 21, 
      ncol = 3, byrow = TRUE, dimnames = list(c(2, 3, 4, 6, 8, 12, 16, 24, 
        32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048), 
        c(10, 11, 12)))
    decim <- Fid_info[1, "DECIM"]
    dspfvs <- Fid_info[1, "DSPFVS"]
    if (!(toString(decim) %in% rownames(group_delay_matrix)))  {
      stop(paste("Invalid DECIM", decim, "it should be one of", rownames(group_delay_matrix)))
    }
    if (!(toString(dspfvs) %in% colnames(group_delay_matrix)))  {
      stop(paste("Invalid DSPFVS", dspfvs, "it should be one of", colnames(group_delay_matrix)))
    }
    group_delay <- group_delay_matrix[toString(decim), toString(dspfvs)]
    if (is.na(group_delay))  {
      stop(paste("Invalid DECIM", decim, "for DSPFVS", dspfvs))
    }
  }
  m <- ncol(Fid_data)
  n <- nrow(Fid_data) 
  
  # GroupDelayCorrection ----------------------------------------------
  
  # We do the shifting in the Fourier domain because the shift can be non-integer.
  # That way we automatically have the circular behaviour of the shift and the
  # interpolation if it is non-integer.
  
  Spectrum <- t(stats::mvfft(t(Fid_data)))
  
  # Spectrum <- FourierTransform(Fid_data, Fid_info)
  p <- ceiling(m/2)
  new_index <- c((p + 1):m, 1:p)
  Spectrum <- Spectrum[,new_index]
  Spectrum <- matrix(data = Spectrum, ncol = m, nrow = n)
  
  Omega <- (0:(m - 1))/m
  i <- complex(real = 0, imaginary = 1)
  
  if (n>1) {
    Spectrum <- sweep(Spectrum, MARGIN = 2, exp(i * group_delay * 2 * pi * Omega), `*`)
    Spectrum <- Spectrum[,new_index]
  }else {
    Spectrum <- Spectrum* exp(i * group_delay * 2 * pi * Omega)
    Spectrum <- Spectrum[new_index]
    Spectrum <- matrix(data = Spectrum, ncol = m, nrow = n)
  }
  
  
  Fid_data <- t(stats::mvfft(t(Spectrum), inverse = TRUE))/m
  colnames(Fid_data) <- dimension_names[[2]]
  rownames(Fid_data) <- dimension_names[[1]]
  
  # Data finalisation ----------------------------------------------
  
  return(endTreatment("GroupDelayCorrection", begin_info, Fid_data, verbose = verbose))
}
