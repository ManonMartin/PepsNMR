#' @export RegionRemoval
RegionRemoval <- function(Spectrum_data, typeofspectra = c(NULL, "Serum", "Urine"), 
                          type.rr = c( "zero", "NA"), fromto.rr=list(Water =c(4.8, 5.2))) {
  begin_info <- beginTreatment("RegionRemoval", Spectrum_data, force.real=T)
  Spectrum_data <- begin_info[["Signal_data"]]
  if (!is.list(fromto.rr) & !is.null(fromto.rr)) {stop(deparse(substitute(fromto.rr)), " is nor a list nor NULL.")}
  if (!is.null(typeofspectra) & !any(c("Serum", "Urine")  %in% typeofspectra)) { stop(paste("typeofspectra must be NULL,\"Serum\" or \"Urine\" "))}
  type.rr = match.arg(type.rr)
  typeofspectra = match.arg(typeofspectra)
  ppm <- as.numeric(colnames(Spectrum_data))
  
  if (typeofspectra == "Urine")  {
    fromto.rr=list(Water =c(4.8, 5.2), Uree=c(5.5, 6))
    } else if (typeofspectra == "Serum") {
      fromto.rr=list(Water =c(4.8, 5.2))
    } else {
     if (is.null(fromto.rr)) {
       stop(paste("fromto.rr argument is null with no default values"))
     } 
     diff = diff(unlist(fromto.rr))[1:length(diff(unlist(fromto.rr))) %% 2 != 0]
     for (i in 1:length(diff)) {
       if (diff[i]<=0) {
         stop(paste("Invalid region removal because from > to"))
       } 
     }
    }

  
  
  
interval = vector("list", length(fromto.rr))
  for (i in 1:length(fromto.rr)) {
    interval[[i]] <- indexInterval(ppm, from = fromto.rr[[i]][1], to = fromto.rr[[i]][2], inclusive=TRUE)
  }



  if (type.rr== "zero") {
    {Spectrum_data[,unlist(interval)] <- 0}
  } else if (type.rr== "NA")
    {Spectrum_data[,unlist(interval)] <- NA}
  else {cat("misspeled \"type\" ")}
  return(endTreatment("RegionRemoval", begin_info, Spectrum_data))
}
