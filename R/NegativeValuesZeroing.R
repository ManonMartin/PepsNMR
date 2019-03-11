#' @export NegativeValuesZeroing
NegativeValuesZeroing <- function(Spectrum_data, verbose = FALSE) {
  # Data initialisation and checks ----------------------------------------------
  checkArg(verbose, c("bool"))
  begin_info <- beginTreatment("NegativeValuesZeroing", Spectrum_data, force.real = TRUE,
                               verbose = verbose)
  Spectrum_data <- begin_info[["Signal_data"]]
  

  
  # NegativeValuesZeroing ----------------------------------------------
  Spectrum_data[Spectrum_data < 0] <- 0
  
  # Data finalisation ----------------------------------------------
  return(endTreatment("NegativeValuesZeroing", begin_info, Spectrum_data, verbose = verbose))
}
