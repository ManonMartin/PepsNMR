#' @export ZoneAggregation
ZoneAggregation <- function (Spectrum_data, fromto.za = list(Citrate =c(2.5, 2.7))) {
  begin_info <- beginTreatment("ZoneAggregation", Spectrum_data, force.real=T)
  Spectrum_data <- begin_info[["Signal_data"]]
  if (!is.list(fromto.za)) {stop(deparse(substitute(fromto.za)), " is not a list.")}
  # Instead of starting the triangle at 0.
  #      /\ ___
  # ----/  \
  # we could start it from the borders values
  #      /\___
  # ----/
  ppm <- as.numeric(colnames(Spectrum_data))
  for (i in 1:length(fromto.za)) {

    interval <- indexInterval(ppm, from = fromto.za[[i]][1], to = fromto.za[[i]][2], inclusive=TRUE)
    p <- length(interval)
    SpectOld <- Spectrum_data[,interval,drop=F]
    S <- rowSums(Re(SpectOld))
    if (p < 3) {
      # Do nothing, the interval is so small that it does not make sense to aggregate it.
      # It does not enter in the general case since we divide by p, p-1 and p-2.
    } else if (p %% 2 == 0) {
      d <- 4*S/((p-2)*p)
      rise <- t(sapply(d*(p/2 - 1),     function (top) seq(0, top, length.out=p/2)))
      triangle <- cbind(rise,            rise[,ncol(rise):1])
    } else {
      d <- 4*S/((p-1)*(p-1))
      rise <- t(sapply(d*((p-1)/2 - 1), function (top) seq(0, top, length.out=(p-1)/2)))
      triangle <- cbind(rise, d*(p-1)/2, rise[,ncol(rise):1,drop=F])
    }
    Spectrum_data[,interval] <- complex(real=triangle, imaginary=Im(SpectOld))
    }
  return(endTreatment("ZoneAggregation", begin_info, Spectrum_data))
}
