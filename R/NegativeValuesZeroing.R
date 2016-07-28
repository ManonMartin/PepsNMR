#' @export NegativeValuesZeroing
NegativeValuesZeroing <- function(RawSpect_data) {
  begin_info <- beginTreatment("NegativeValuesZeroing", RawSpect_data, force.real=T)
  RawSpect_data <- begin_info[["Signal_data"]]
  RawSpect_data[RawSpect_data < 0] <- 0
  return(endTreatment("NegativeValuesZeroing", begin_info, RawSpect_data))
}
