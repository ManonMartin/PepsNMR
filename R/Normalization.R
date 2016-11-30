#' @export Normalization
Normalization <- function (Spectrum_data, type.norm=c("mean","pqn", "median","firstquartile", "peak"), fromto.norm = c(3.05, 4.05), ref.norm=1) {
  begin_info <- beginTreatment("Normalization", Spectrum_data, force.real=T)
  Spectrum_data <- begin_info[["Signal_data"]]
  type.norm <- match.arg(type.norm)
  checkArg(fromto.norm, "num")

  if(length(fromto.norm)>2) {
  warning("fromto.norm has a length > 2, only the first two elements are taken into account")
  }
  
  switch(type.norm,
         "mean" = { # mean
           factor <- rowMeans(Spectrum_data, na.rm = TRUE)
         },
         "median" = { # median
#            require("matrixStats")
           factor <- matrixStats::rowMedians(Spectrum_data, na.rm = TRUE)
         },
         "firstquartile" = {
#            require("matrixStats")
           factor <- matrixStats::rowQuantiles(Spectrum_data, 0.25, na.rm=TRUE)[[1]]
         },
         "peak" = {
           ppm <- as.numeric(colnames(Spectrum_data))
           interval <- indexInterval(ppm, fromto.norm[1], fromto.norm[2], inclusive=TRUE)
           Spectrum_dataInZone = Spectrum_data[,interval,drop=F]
           peakInZone <- which.max(colSums(Spectrum_dataInZone))
           factor <- Spectrum_dataInZone[,peakInZone]
         },
          "pqn" = {
            if (missing(ref.norm)) {
              ref.norm <- matrixStats::colMedians(Spectrum_data, na.rm = TRUE)
            } else ref.norm <- Spectrum_data[ref.norm,]
            quotient <- t(Spectrum_data)/ref.norm
            factor <- quotient.median <- matrixStats::colMedians(quotient, na.rm = TRUE)
          }
  )
  invalid <- (factor <= 0)
  if (TRUE %in% invalid) {
    invalid_factor <- factor[invalid]
    warning("The ", type.norm, "s are ", paste(invalid_factor, collapse=", "), " for the spectrums ", paste(names(invalid_factor), collapse=", "), " which is nonpositive, the normalization will not happen for them.")
    factor[invalid] <- 1
  }
  return(endTreatment("Normalization", begin_info, Spectrum_data / factor))
}
