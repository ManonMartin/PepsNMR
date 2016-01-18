#' @export FourierTransform
FourierTransform <- function(Fid_data, Fid_info=NULL, SW_h=NULL) {
  begin_info <- beginTreatment("FourierTransform", Fid_data, Fid_info)
  Fid_data <- begin_info[["Signal_data"]]
  Fid_info <- begin_info[["Signal_info"]]
  getArg(SW_h, Fid_info, "SW_h")
  m <- ncol(Fid_data)
  # mvfft does the unnormalized fourier transform (see ?mvfft),
  # so we need divide by m.
  # It does not matter a lot in our case since the spectrum will be normalized.
  RawSpect_data <- fftshift1D2D(t(mvfft(t(Fid_data)))) / m
  f <- ((0:(m-1)) - floor(m/2)) * Fid_info[1,"SW_h"] / m
  colnames(RawSpect_data) <- f
  return(endTreatment("FourierTransform", begin_info, RawSpect_data))
}
