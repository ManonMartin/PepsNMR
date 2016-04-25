#' @export Scaling
Scaling <- function(Spectrum_data, type.scaling=c("center", "pareto", "uv")) {
  
  begin_info <- beginTreatment("Scaling", Spectrum_data, force.real=T)
  Spectrum_data <- begin_info[["Signal_data"]]
  type.scaling <- match.arg(type.scaling)
  
  stdev<-apply(Spectrum_data,2,sd, na.rm = TRUE)
  
  invalid <- (stdev == 0)
  if (TRUE %in% invalid & type.scaling != "center") {
    invalid_factor <- stdev[invalid]
    warning("The ", type.scaling, "s are 0 for the spectra ", paste(substr(names(invalid_factor),1,4), collapse=", "), " which is null, the scaling will not happen for them.")
    stdev[invalid] <- 1
  }
  
  switch(type.scaling,
         "center" = { 
           Scaled.spectrum <- apply(Spectrum_data,2, (x - mean(x)))
         },
         
         "uv" = { 
           Scaled.spectrum <- apply(Spectrum_data,2, function(x) (x - mean(x)) / sd(x))
         },
         
         "pareto" = { 
           Scaled.spectrum <- apply(Spectrum_data,2, function(x) (x - mean(x)) / sqrt(sd(x)))
         }
  )
return(endTreatment("Scaling", begin_info, Scaled.spectrum))
       
}