#' @export RegionRemoval

RegionRemoval <- function(Spectrum_data, typeofspectra = c("manual", "serum", "urine"), 
                          type.rr = c("zero", "NA"), fromto.rr = list(Water = c(4.5, 5.1))) {
  
  
  # Data initialisation and checks ----------------------------------------------
  
  begin_info <- beginTreatment("RegionRemoval", Spectrum_data, force.real = T)
  Spectrum_data <- begin_info[["Signal_data"]]
  if (!is.list(fromto.rr) & !missing(fromto.rr)) {
    stop(deparse(substitute(fromto.rr)), " is nor a list nor NULL.")
  }
  
  type.rr <- match.arg(type.rr)
  typeofspectra <- match.arg(typeofspectra)
  ppm <- as.numeric(colnames(Spectrum_data))
  
  
  # Region Removal ---------------------------------------------------------------
  # Assign and check fromto.rr values 
  if (typeofspectra == "urine") {
    fromto.rr <- list(Water_Uree_MalAc = c(4.5, 6.1))
  } else if (typeofspectra == "serum") {
    fromto.rr <- list(Water = c(4.5, 5.1))
  } else {
    diff <- diff(unlist(fromto.rr))[1:length(diff(unlist(fromto.rr)))%%2 != 0]
    for (i in 1:length(diff)) {
      if (diff[i] <= 0) {
        stop(paste("Invalid region removal because from > to"))
      }
    }
  }
  
  
  # Defines the intervals bounds
  interval <- vector("list", length(fromto.rr))
  for (i in 1:length(fromto.rr)) {
    interval[[i]] <- indexInterval(ppm, from = fromto.rr[[i]][1], to = fromto.rr[[i]][2], 
      inclusive = TRUE)
  }
  
  # Region removal
  if (type.rr == "zero") {
      Spectrum_data[, unlist(interval)] <- 0
  } else if (type.rr == "NA") {
    Spectrum_data[, unlist(interval)] <- NA
  } else {
    cat("misspeled \"type\" ")
  }
  
  
  # Data finalisation ----------------------------------------------
  
  return(endTreatment("RegionRemoval", begin_info, Spectrum_data))
}
