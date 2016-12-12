#' @export NegativeValuesZeroing
NegativeValuesZeroing <- function(Spectrum_data) {
  # Data initialisation and checks ----------------------------------------------
  begin_info <- beginTreatment("NegativeValuesZeroing", Spectrum_data, force.real = T)
  Spectrum_data <- begin_info[["Signal_data"]]
  
  # NegativeValuesZeroing ----------------------------------------------
  Spectrum_data[Spectrum_data < 0] <- 0
  
  # Data finalisation ----------------------------------------------
  return(endTreatment("NegativeValuesZeroing", begin_info, Spectrum_data))
}
