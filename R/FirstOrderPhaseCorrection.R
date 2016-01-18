#' @export FirstOrderPhaseCorrection
FirstOrderPhaseCorrection <- function(Fid_data, Fid_info=NULL, group_delay=NULL) {
  begin_info <- beginTreatment("FirstOrderPhaseCorrection", Fid_data, Fid_info)
  Fid_data <- begin_info[["Signal_data"]]
  Fid_info <- begin_info[["Signal_info"]]
  checkArg(group_delay, c("num", "pos0"), can.be.null=TRUE)
  # if Fid_info and group_delay are NULL, getArg will generate an error
  group_delay <- getArg(group_delay, Fid_info, "GRPDLY", can.be.absent=TRUE)
  if (is.null(group_delay)) {
    # See DetermineBrukerDigitalFilter.m in matNMR MATLAB library
    group_delay_matrix <- matrix(c(
      44.7500, 46.0000, 46.311,
      33.5000, 36.5000, 36.530,
      66.6250, 48.0000, 47.870,
      59.0833, 50.1667, 50.229,
      68.5625, 53.2500, 53.289,
      60.3750, 69.5000, 69.551,
      69.5313, 72.2500, 71.600,
      61.0208, 70.1667, 70.184,
      70.0156, 72.7500, 72.138,
      61.3438, 70.5000, 70.528,
      70.2578, 73.0000, 72.348,
      61.5052, 70.6667, 70.700,
      70.3789, 72.5000, 72.524,
      61.5859, 71.3333, NA,
      70.4395, 72.2500, NA,
      61.6263, 71.6667, NA,
      70.4697, 72.1250, NA,
      61.6465, 71.8333, NA,
      70.4849, 72.0625, NA,
      61.6566, 71.9167, NA,
      70.4924, 72.0313, NA),
      nrow = 21, ncol = 3, byrow = TRUE,
      dimnames = list(
        c(2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048),
        c(10, 11, 12)))
    decim <- Fid_info[1,"DECIM"]
    dspfvs <- Fid_info[1,"DSPFVS"]
    if (!(toString(decim) %in% rownames(group_delay_matrix))) {
      stop(paste("Invalid DECIM", decim, "it should be one of", rownames(group_delay_matrix)))
    }
    if (!(toString(dspfvs) %in% colnames(group_delay_matrix))) {
      stop(paste("Invalid DSPFVS", dspfvs, "it should be one of", colnames(group_delay_matrix)))
    }
    group_delay <- group_delay_matrix[toString(decim), toString(dspfvs)]
    if (is.na(group_delay)) {
      stop(paste("Invalid DECIM", decim, "for DSPFVS", dspfvs))
    }
  }
  m <- ncol(Fid_data)
  # We do the shifting in the Fourier domain because the shift can be non-integer.
  # That way we automatically have the circular behaviour of the shift and the interpolation
  # if it is non-integer.
  Spectrum <- t(mvfft(t(Fid_data))) / m
  Omega <- (0:(m-1)) / m
  i <- complex(real=0,imaginary=1)
  Spectrum <- sweep(Spectrum, MARGIN=2, exp(i*group_delay*2*pi*Omega), `*`)
  Fid_data <- t(mvfft(t(Spectrum), inverse = TRUE))
  return(endTreatment("FirstOrderPhaseCorrection", begin_info, Fid_data))
}
