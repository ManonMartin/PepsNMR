#' @export NegativeValuesZeroing
NegativeValuesZeroing <- function(RawSpect_data) {
  # Data initialisation and checks ----------------------------------------------
  begin_info <- beginTreatment("NegativeValuesZeroing", RawSpect_data, force.real = T)
  RawSpect_data <- begin_info[["Signal_data"]]
  
  # NegativeValuesZeroing ----------------------------------------------
  RawSpect_data[RawSpect_data < 0] <- 0
  
  # Data finalisation ----------------------------------------------
  return(endTreatment("NegativeValuesZeroing", begin_info, RawSpect_data))
}
